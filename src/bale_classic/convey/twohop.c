// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <inttypes.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "convey_impl.h"
#include "private.h"

#define DEBUG2_PRINT(...)
// #define DEBUG2_PRINT(...) mprint(MY_PROC, 0, __VA_ARGS__)


#if MPP_NO_MIMD

convey_t*
convey_new_twohop(size_t capacity, size_t row_procs,
                  const convey_alc8r_t* alloc, uint64_t options)
{
  return convey_new_simple(capacity, alloc, NULL, options);
}

#elif !MPP_USE_MPI && !(MPP_USE_SHMEM && HAVE_SHMEMX_ALLTOALLV && HAVE_SHMEMX_TEAM_ALLTOALLV)

static uint64_t n_calls = 0;

convey_t*
convey_new_twohop(size_t capacity, size_t row_procs,
                  const convey_alc8r_t* alloc, uint64_t options)
{
  if (0 == n_calls++)
    mprint(0, 0, "not supported in this build\n");
  return NULL;
}

#else

#ifdef MPP_USE_MPI
typedef int a2a_off_t;
typedef MPI_Comm pe_team_t;
#else
typedef size_t a2a_off_t;
typedef shmem_team_t pe_team_t;
#endif

typedef struct area {
  uint32_t* next;   // place to write or read next item
  uint32_t* limit;  // upper limit of the area
} area_t;

typedef struct twohop {
  convey_t convey;
  bool overflow;         // has any buffer filled up?
  bool row_quiet;        // nothing more to send from sources?
  bool col_quiet;        // nothing more to pivot?
  bool row_flip;         // which row_sync to use next?
  bool col_flip;         // which col_sync to use next?
  int64_t pull_row;      // draw from which col_recv buffer?
  int64_t pull_col;      // draw from which row_recv buffer?
  void* aligned_item;    // properly aligned staging area for an item
  size_t n_rows;
  size_t n_cols;
  divbymul32_t get_row;  // for division by n_cols
  uint32_t my_col;
  int8_t col_bits;       // number of bits needed to hold column #
  bool odd;              // does n_cols have an odd factor?
  bool dynamic;          // do we deallocate buffers while dormant?
  size_t row_stride;     // stride in 4-byte chunks between row buffers
  size_t col_stride;     // stride in 4-byte chunks between column buffers
  size_t packet_quads;   // how many 4-byte chunks per packet?
  size_t row_quads;      // size of each row buffer (one per column)
  size_t col_quads;      // size of each column buffer (one per row)
  area_t* row_send;
  area_t* row_recv;
  area_t* col_send;
  area_t* col_recv;
  a2a_off_t* send_sizes;    // reused for different all-to-alls
  a2a_off_t* recv_sizes;
  a2a_off_t* row_offsets;   // offsets needed for alltoallv
  a2a_off_t* col_offsets;   // (a.k.a. displacements)
  uint32_t* row_send_bufs;
  uint32_t* row_recv_bufs;
  uint32_t* col_send_bufs;
  uint32_t* col_recv_bufs;
  long* row_sync[2];        // for shmemx_alltoallv
  long* col_sync[2];        // for shmemx_alltoallv
  pe_team_t col_team;    // for communication along columns
  pe_team_t row_team;    // for communication along rows (MPI only)
  convey_alc8r_t alloc;
  mpp_comm_t comm;       // for error checking
  int64_t stats[convey_imp_N_STATS];
} twohop_t;


/*** ALL-TO-ALL COMMUNICATION ***/

static bool
alloc_buffers(twohop_t* twohop)
{
  convey_alc8r_t* alloc = &twohop->alloc;

  size_t row_size = twohop->n_cols * twohop->row_stride * 4;
  twohop->row_send_bufs = alloc->grab(alloc->alc8r, row_size, __FILE__, __LINE__);
  twohop->row_recv_bufs = alloc->grab(alloc->alc8r, row_size, __FILE__, __LINE__);

  size_t col_size = twohop->n_rows * twohop->col_stride * 4;
  twohop->col_send_bufs = alloc->grab(alloc->alc8r, col_size, __FILE__, __LINE__);
  twohop->col_recv_bufs = alloc->grab(alloc->alc8r, col_size, __FILE__, __LINE__);

  return (twohop->row_send_bufs && twohop->row_recv_bufs &&
          twohop->col_send_bufs && twohop->col_recv_bufs);
}

static void
dealloc_buffers(twohop_t* twohop)
{
  convey_alc8r_t* alloc = &twohop->alloc;
  alloc->free(alloc->alc8r, twohop->col_recv_bufs);
  alloc->free(alloc->alc8r, twohop->col_send_bufs);
  alloc->free(alloc->alc8r, twohop->row_recv_bufs);
  alloc->free(alloc->alc8r, twohop->row_send_bufs);
  twohop->row_send_bufs = twohop->row_recv_bufs = NULL;
  twohop->col_send_bufs = twohop->col_recv_bufs = NULL;
}

static void
reset_send_areas(size_t n, area_t areas[n], uint32_t* base,
                 size_t qstride, size_t quads)
{
  for (size_t i = 0; i < n; i++) {
    areas[i].next = base + i * qstride;
    areas[i].limit = areas[i].next + quads;
  }
}

static void
reset_recv_areas(size_t n, area_t areas[n], uint32_t* base, size_t qstride)
{
  for (size_t i = 0; i < n; i++) {
    areas[i].next = base + i * qstride;
    areas[i].limit = areas[i].next;
  }
}

static bool
rows_nonempty(twohop_t* twohop)
{
  const size_t quads = twohop->row_stride;
  for (size_t c = 0; c < twohop->n_cols; c++)
    if (twohop->row_send[c].next > twohop->row_send_bufs + c * quads)
      return true;
  return false;
}

static int
row_advance(twohop_t* twohop)
{
  const size_t n_cols = twohop->n_cols;
  const size_t row_stride = twohop->row_stride;

  // Convert areas to sizes, measured in bytes
  uint64_t n_bytes = 0;
  uint32_t* buffer = twohop->row_send_bufs;
  for (size_t c = 0; c < n_cols; c++) {
    size_t size = (twohop->row_send[c].next - buffer) * 4;
    twohop->send_sizes[c] = size;
    n_bytes += size;
    buffer += row_stride;
  }

  // Call the library function
  {
    CONVEY_PROF_DECL(_sample);
    CONVEY_PROF_START(&_sample);
#ifdef MPP_USE_SHMEM
    const uint32_t my_proc = MY_PROC;
    const int my_start = n_cols * _divbymul32(my_proc, twohop->get_row);
    shmemx_alltoallv(twohop->row_recv_bufs, twohop->row_offsets, twohop->recv_sizes,
                     twohop->row_send_bufs, twohop->row_offsets, twohop->send_sizes,
                     twohop->comm.start + my_start, _trailz(twohop->comm.stride), n_cols,
                     twohop->row_sync[twohop->row_flip]);
    twohop->row_flip ^= 1;
#else
    int error = MPI_Alltoall(twohop->send_sizes, 1, MPI_INT,
                             twohop->recv_sizes, 1, MPI_INT, twohop->row_team);
    if (!error)
      error = MPI_Alltoallv(twohop->row_send_bufs, twohop->send_sizes,
                            twohop->row_offsets, MPI_BYTE,
                            twohop->row_recv_bufs, twohop->recv_sizes,
                            twohop->row_offsets, MPI_BYTE, twohop->row_team);
    if (error)
      return convey_error_COMMS;
#endif
    CONVEY_PROF_STOP(&_sample, PROF_OP_USER_COLLECTIVE, -1, n_bytes);
  }

  // Convert sizes to areas
  buffer = twohop->row_recv_bufs;
  size_t n_quads = 0;
  for (size_t c = 0; c < n_cols; c++) {
    area_t* area = twohop->row_recv + c;
    area->next = buffer;
    area->limit = buffer + (twohop->recv_sizes[c] >> 2);
    n_quads += (twohop->recv_sizes[c] >> 2);
    buffer += row_stride;
  }
  DEBUG2_PRINT("received %zu items\n", n_quads / twohop->packet_quads);

  twohop->stats[convey_COMMS] += n_cols;
  twohop->stats[convey_SYNCS] += 1;
  bool empty = (n_quads == 0);
  twohop->pull_col = empty ? n_cols : 0;
  return empty ? convey_DONE : convey_OK;
}

// Transfer as much data as we can from row buffers to column buffers.
// Returns true if we moved it all (no overflow occurred).
// The header of an input packet is (dest row) << col_bits | (source col).
// The header of an output packet is the source PE number.
static bool
pivot(twohop_t* twohop)
{
  reset_send_areas(twohop->n_rows, twohop->col_send, twohop->col_send_bufs,
                   twohop->col_stride, twohop->col_quads);

  const size_t n_cols = twohop->n_cols;
  const int64_t shift = twohop->col_bits;
  const size_t col_mask = _maskr(shift);
  const size_t row_data = MY_PROC - twohop->my_col;
  const size_t quads = twohop->packet_quads;
  const size_t copy_bytes = 4 * (quads - 1);

  // This process could be optimized!
  bool overflow = false;
  int64_t col = twohop->pull_col;
  while (col < n_cols) {
    area_t* source = twohop->row_recv + col;
    uint32_t* here = source->next;
    if (here < source->limit) {
      uint32_t info = *here;
      area_t* dest = twohop->col_send + (info >> shift);
      if (dest->next >= dest->limit) {
        overflow = true;
        break;
      }
      uint32_t* there = dest->next;
      *there = row_data + (info & col_mask);
      memcpy(there + 1, here + 1, copy_bytes);
      DEBUG2_PRINT("pivot %08x from %u\n", here[1], there[0]);
      source->next = here + quads;
      dest->next = there + quads;
    }
    else
      twohop->pull_col = ++col;
  }

  return !overflow;
}

static int
column_advance(twohop_t* twohop)
{
  const size_t n_rows = twohop->n_rows;

  // Convert areas to sizes, measured in bytes
  uint64_t n_bytes = 0;
  uint32_t* buffer = twohop->col_send_bufs;
  for (size_t r = 0; r < n_rows; r++) {
    size_t size = (twohop->col_send[r].next - buffer) * 4;
    twohop->send_sizes[r] = size;
    n_bytes += size;
    buffer += twohop->col_stride;
  }

  // Call the library function
  {
    CONVEY_PROF_DECL(_sample);
    CONVEY_PROF_START(&_sample);
#ifdef MPP_USE_SHMEM
    if (!twohop->odd) {
      const long my_proc = MY_PROC;
      const int my_start = my_proc & _maskr(twohop->col_bits);
      const int lg_stride = _trailz(twohop->comm.stride) + twohop->col_bits;
      shmemx_alltoallv(twohop->col_recv_bufs, twohop->col_offsets, twohop->recv_sizes,
                       twohop->col_send_bufs, twohop->col_offsets, twohop->send_sizes,
                       twohop->comm.start + my_start, lg_stride, n_rows,
                       twohop->col_sync[twohop->col_flip]);
    }
    else
      shmemx_team_alltoallv(twohop->col_recv_bufs, twohop->col_offsets, twohop->recv_sizes,
                            twohop->col_send_bufs, twohop->col_offsets, twohop->send_sizes,
                            twohop->col_team, twohop->col_sync[twohop->col_flip]);
    twohop->col_flip ^= 1;
#else
    int error = MPI_Alltoall(twohop->send_sizes, 1, MPI_INT,
                             twohop->recv_sizes, 1, MPI_INT, twohop->col_team);
    if (!error)
      error = MPI_Alltoallv(twohop->col_send_bufs, twohop->send_sizes,
                            twohop->col_offsets, MPI_BYTE,
                            twohop->col_recv_bufs, twohop->recv_sizes,
                            twohop->col_offsets, MPI_BYTE, twohop->col_team);
    if (error)
      return convey_error_COMMS;
#endif
    CONVEY_PROF_STOP(&_sample, PROF_OP_A2A, -1, n_bytes);
  }

  // Convert sizes to areas
  buffer = twohop->col_recv_bufs;
  for (size_t r = 0; r < n_rows; r++) {
    area_t* area = twohop->col_recv + r;
    area->next = buffer;
    area->limit = buffer + (twohop->recv_sizes[r] >> 2);
    buffer += twohop->col_stride;
  }

  twohop->pull_row = 0;
  twohop->stats[convey_COMMS] += n_rows;
  twohop->stats[convey_SYNCS] += 1;
  return convey_OK;
}


/*** METHODS (minimal error checking) ***/

static int
twohop_push(convey_t* self, const void* item, int64_t pe)
{
  const size_t item_size = self->item_size;
  uint32_t dest = pe;
  twohop_t* twohop = (twohop_t*) self;
  uint32_t row = _divbymul32(dest, twohop->get_row);
  uint32_t col = dest - (twohop->n_cols * row);
  area_t* area = twohop->row_send + col;

  if (area->next < area->limit) {
    DEBUG2_PRINT("push  %08x to %u\n", *(uint32_t*)item, dest);
    area->next[0] = (row << twohop->col_bits) | twohop->my_col;
    memcpy(area->next + 1, item, item_size);
    area->next += twohop->packet_quads;
    twohop->stats[convey_PUSHES]++;
    return convey_OK;
  }

  twohop->overflow = true;
  return convey_FAIL;
}

// Unaligned pull
static void*
twohop_upull(convey_t* self, int64_t* from)
{
  twohop_t* twohop = (twohop_t*) self;
  const size_t n_rows = twohop->n_rows;
  int64_t row = twohop->pull_row;

  while (row < n_rows) {
    area_t* area = twohop->col_recv + row;
    uint32_t* here = area->next;
    if (here < area->limit) {
      if (from)
        *from = *here;
      area->next = here + twohop->packet_quads;
      DEBUG2_PRINT("pull  %08x from %u\n", here[1], here[0]);
      twohop->stats[convey_PULLS]++;
      return here + 1;
    }
    twohop->pull_row = ++row;
  }

  return NULL;
}

static int
twohop_pull(convey_t* self, void* item, int64_t* from)
{
  void* source = twohop_upull(self, from);
  if (source == NULL)
    return convey_FAIL;
  memcpy(item, source, self->item_size);
  return convey_OK;
}

static void*
twohop_apull(convey_t* self, int64_t* from)
{
  twohop_t* twohop = (twohop_t*) self;
  void* source = twohop_upull(self, from);
  void* temp = twohop->aligned_item;
  if (source == NULL || temp == NULL)
    return source;
  memcpy(temp, source, self->item_size);
  return temp;
}

static int
twohop_unpull(convey_t* self)
{
  twohop_t* twohop = (twohop_t*) self;
  const size_t n_rows = twohop->n_rows;
  int64_t row = twohop->pull_row;
  if (row >= n_rows)
    return convey_FAIL;

  uint32_t* bottom = twohop->col_recv_bufs + row * twohop->col_stride;
  uint32_t* next = twohop->col_recv[row].next;
  if (next <= bottom)
    return convey_FAIL;

  uint32_t* where = next - twohop->packet_quads;
  twohop->col_recv[row].next = where;
  twohop->stats[convey_UNPULLS]++;
  return convey_OK;
}

static int
twohop_advance(convey_t* self, bool done)
{
  twohop_t* twohop = (twohop_t*) self;
  twohop->stats[convey_ADVANCES]++;

  // We are in state WORKING or ENDGAME, and need to decide collectively
  // what to do.  The row_quiet and col_quiet flags are collective state.

  bool need_pull = (twohop->pull_row < twohop->n_rows);
  if (twohop->row_quiet && twohop->col_quiet)
    return need_pull ? convey_NEAR : convey_DONE;

  twohop->stats[convey_SYNCS]++;
  bool cannot_comm = need_pull;
  bool desire_comm = (done | twohop->overflow) ||
    ((self->features & convey_STEADY) && (!twohop->col_quiet || rows_nonempty(twohop)));
  bool further_comm = !done;
  long status = mpp_or_long(cannot_comm << 0 | desire_comm << 1 | further_comm << 2);
  cannot_comm = (status >> 0) & 1;
  desire_comm = (status >> 1) & 1;
  further_comm = (status >> 2) & 1;
  if (cannot_comm || !desire_comm)
    return convey_OK;

  if (twohop->col_quiet) {
    int result = row_advance(twohop);
    reset_send_areas(twohop->n_cols, twohop->row_send, twohop->row_send_bufs,
                     twohop->row_stride, twohop->row_quads);
    if (result < 0)
      return result;
    twohop->overflow = false;
    twohop->col_quiet = false;
    twohop->row_quiet = !further_comm;
    mpp_barrier(1);  // must separate the two shmemx_alltoallv calls
  }
  if (!twohop->col_quiet) {
    bool all_moved = pivot(twohop);
    int result = column_advance(twohop);
    if (result < 0)
      return result;
    twohop->col_quiet = mpp_and_long(all_moved);
  }

  return (twohop->row_quiet && twohop->col_quiet) ? convey_NEAR : convey_OK;
}

static void
reset_areas(twohop_t* twohop)
{
  reset_send_areas(twohop->n_rows, twohop->col_send, twohop->col_send_bufs,
                   twohop->col_stride, twohop->col_quads);
  reset_recv_areas(twohop->n_rows, twohop->col_recv, twohop->col_recv_bufs,
                   twohop->col_stride);
  reset_send_areas(twohop->n_cols, twohop->row_send, twohop->row_send_bufs,
                   twohop->row_stride, twohop->row_quads);
  reset_recv_areas(twohop->n_cols, twohop->row_recv, twohop->row_recv_bufs,
                   twohop->row_stride);
}

static int
twohop_begin(convey_t* self, size_t item_size, size_t align)
{
  twohop_t* twohop = (twohop_t*) self;
  if (!mpp_comm_is_equal(MPP_COMM_CURR, twohop->comm))
    return convey_error_TEAM;
  if (item_size == 0)
    return convey_error_ZERO;
  size_t packet_quads = 1 + ((item_size + 3) >> 2);
  if (packet_quads > twohop->row_stride || packet_quads > twohop->col_stride)
    return convey_error_OFLO;

  twohop->stats[convey_BEGINS]++;
  self->item_size = item_size;
  twohop->packet_quads = packet_quads;
  twohop->row_quads = packet_quads * (twohop->row_stride / packet_quads);
  twohop->col_quads = packet_quads * (twohop->col_stride / packet_quads);

  if (! convey_prep_aligned(&twohop->aligned_item, item_size, 4, align))
    return convey_error_ALLOC;
  if (twohop->dynamic && !alloc_buffers(twohop))
    return convey_error_ALLOC;
  reset_areas(twohop);

  twohop->overflow = false;
  twohop->row_quiet = false;
  twohop->col_quiet = true;
  mpp_barrier(1);
  return convey_OK;
}

static int
twohop_reset(convey_t* self)
{
  twohop_t* twohop = (twohop_t*) self;
  if (!mpp_comm_is_equal(MPP_COMM_CURR, twohop->comm))
    return convey_error_TEAM;
  free(twohop->aligned_item);
  twohop->aligned_item = NULL;
  if (twohop->dynamic)
    dealloc_buffers(twohop);
  convey_imp_update_stats(twohop->stats);
  return convey_OK;
}

static void
free_allocs(twohop_t* twohop)
{
  // Make sure allocations are last-in, first-out
  dealloc_buffers(twohop);
  void* alc8r = twohop->alloc.alc8r;
  void (*dealloc)(void*,void*) = twohop->alloc.free;
  dealloc(alc8r, twohop->col_sync[1]);
  dealloc(alc8r, twohop->row_sync[1]);
  dealloc(alc8r, twohop->col_sync[0]);
  dealloc(alc8r, twohop->row_sync[0]);
  dealloc(alc8r, twohop->col_offsets);
  dealloc(alc8r, twohop->row_offsets);
  dealloc(alc8r, twohop->recv_sizes);
  dealloc(alc8r, twohop->send_sizes);
  dealloc(alc8r, twohop->col_recv);
  dealloc(alc8r, twohop->col_send);
  dealloc(alc8r, twohop->row_recv);
  dealloc(alc8r, twohop->row_send);
  free(twohop->aligned_item);
  free(twohop);
}

static int
twohop_free(convey_t* self)
{
  twohop_t* twohop = (twohop_t*) self;
  if (!mpp_comm_is_equal(MPP_COMM_CURR, twohop->comm))
    return convey_error_TEAM;
#ifdef MPP_USE_MPI
  MPI_Comm_free(&twohop->row_team);
  MPI_Comm_free(&twohop->col_team);
#else
  if (twohop->odd)
    shmemx_team_destroy(&twohop->col_team);
#endif
  free_allocs(twohop);
  return convey_OK;
}

static int64_t
twohop_statistic(convey_t* self, int which)
{
  twohop_t* twohop = (twohop_t*) self;
  return convey_imp_statistic(twohop->stats, which);
}

static const convey_methods_t twohop_fast_methods = {
  .push = &twohop_push,
  .pull = &twohop_pull,
  .unpull = &twohop_unpull,
  .advance = &twohop_advance,
  .begin = &twohop_begin,
  .reset = &twohop_reset,
  .free = &twohop_free,
  .statistic = &twohop_statistic,
  .apull = &twohop_apull,
  .epush = &convey_no_epush,
  .epull = &convey_no_epull,
  .panic = &convey_imp_panic,
};

static const convey_methods_t twohop_methods = {
  ._super_ = &twohop_fast_methods,
  .push = &convey_checked_push,
  .pull = &convey_checked_pull,
  .unpull = &convey_checked_unpull,
  .advance = &twohop_advance,
  .begin = &twohop_begin,
  .reset = &twohop_reset,
  .free = &twohop_free,
  .statistic = &twohop_statistic,
  .apull = &convey_checked_apull,
  .epush = &convey_no_epush,
  .epull = &convey_no_epull,
  .panic = &convey_imp_panic,
};


/*** PUBLIC INTERFACE ***/

convey_t*
convey_new_twohop(size_t capacity, size_t row_procs,
                  const convey_alc8r_t* alloc, uint64_t options)
{
  bool quiet = (options & convey_opt_QUIET);
  if (capacity == 0 || row_procs == 0)
    CONVEY_REJECT(quiet, "invalid arguments");
  if (PROCS % row_procs)
    CONVEY_REJECT(quiet, "row_procs does not divide PROCS");
  if (row_procs == 1)
    return convey_new_simple(capacity, alloc, NULL, options);

  if (alloc == NULL)
    alloc = &convey_imp_alloc_align;
  else if (!alloc->grab || !alloc->free)
    CONVEY_REJECT(quiet, "alloc is missing one or both methods");

  twohop_t* twohop = malloc(sizeof(twohop_t));
  if (twohop == NULL)
    CONVEY_REJECT(quiet, "a small malloc() failed!");

  const long my_proc = MY_PROC;
  size_t n_procs = PROCS;
  size_t n_cols = row_procs;
  size_t n_rows = n_procs / n_cols;
  bool wide = (n_cols >= n_rows);
  size_t row_bytes = wide ? capacity : (capacity * n_rows) / n_cols;
  size_t col_bytes = wide ? (capacity * n_cols) / n_rows : capacity;
  row_bytes = (row_bytes + CONVEY_MAX_ALIGN - 1) & -CONVEY_MAX_ALIGN;
  col_bytes = (col_bytes + CONVEY_MAX_ALIGN - 1) & -CONVEY_MAX_ALIGN;
  int col_bits = 64 - _leadz(n_cols - 1);

  // Initialize some things; set scratch pointers to NULL
  *twohop = (twohop_t) {
    .convey = {
      ._class_ = (options & convey_opt_RECKLESS) ? &twohop_fast_methods : &twohop_methods,
      .features = (options & convey_opt_PROGRESS) ? convey_STEADY : 0,
      .n_procs = n_procs, .state = convey_DORMANT,
      .suppress = UINT64_C(0) - quiet,
    },
    .n_rows = n_rows, .n_cols = n_cols, .col_bits = col_bits,
    .odd = n_cols & (n_cols - 1), .my_col = my_proc % n_cols,
    .row_stride = (row_bytes + 3) >> 2, .col_stride = (col_bytes + 3) >> 2,
    .alloc = *alloc, .dynamic = (options & convey_opt_DYNAMIC),
  };
  twohop->get_row = _divbymul32_prep(n_cols);
  twohop->comm = MPP_COMM_CURR;

  // Make the necessary teams
#ifdef MPP_USE_MPI
  {
    int error;
    error = MPI_Comm_split(mpp_comm_mpi(twohop->comm),
                           my_proc / n_cols, my_proc % n_cols, &twohop->row_team);
    if (!error)
      error = MPI_Comm_split(mpp_comm_mpi(twohop->comm),
                             my_proc % n_cols, my_proc / n_cols, &twohop->col_team);
    if (error) {
      free(twohop);
      CONVEY_REJECT(quiet, "unable to create row and column communicators");
    }
  }
#else
  if (twohop->odd) {
    shmem_team_t full_team;
    shmemx_team_create_strided(twohop->comm.start, twohop->comm.stride, n_procs, &full_team);
    shmemx_team_split(full_team, my_proc % n_cols, my_proc / n_cols, &twohop->col_team);
    shmemx_team_destroy(&full_team);
  }
#endif

  // Try to do all the allocations
#define GRAB(SIZE) alloc->grab(alloc->alc8r, (SIZE), __FILE__, __LINE__)
  twohop->row_send = GRAB(n_cols * sizeof(area_t));
  twohop->row_recv = GRAB(n_cols * sizeof(area_t));
  twohop->col_send = GRAB(n_rows * sizeof(area_t));
  twohop->col_recv = GRAB(n_rows * sizeof(area_t));
  size_t n_max = (n_cols > n_rows) ? n_cols : n_rows;
  twohop->send_sizes = GRAB(n_max * sizeof(a2a_off_t));
  twohop->recv_sizes = GRAB(n_max * sizeof(a2a_off_t));
  twohop->row_offsets = GRAB(n_cols * sizeof(a2a_off_t));
  twohop->col_offsets = GRAB(n_rows * sizeof(a2a_off_t));
  bool ok = (twohop->row_send && twohop->row_recv &&
             twohop->col_send && twohop->col_recv &&
             twohop->send_sizes && twohop->row_offsets &&
             twohop->recv_sizes && twohop->col_offsets);
#ifdef MPP_USE_SHMEM
  for (int i = 0; i < 2; i++) {
    twohop->row_sync[i] = GRAB(_SHMEM_ALLTOALL_SYNC_SIZE * sizeof(long));
    twohop->col_sync[i] = GRAB(_SHMEM_ALLTOALL_SYNC_SIZE * sizeof(long));
    ok &= (twohop->row_sync[i] && twohop->col_sync[i]);
  }
#endif
  if (!twohop->dynamic)
    ok &= alloc_buffers(twohop);
#undef GRAB
  if (!ok) {
    free_allocs(twohop);
    CONVEY_REJECT(quiet, "symmetric allocation failed");
  }

  // Initialize everything else that begin() does not
  for (size_t c = 0; c < n_cols; c++)
    twohop->row_offsets[c] = c * (4 * twohop->row_stride);
  for (size_t r = 0; r < n_rows; r++)
    twohop->col_offsets[r] = r * (4 * twohop->col_stride);
#ifdef MPP_USE_SHMEM
  for (int i = 0; i < 2; i++)
    for (int64_t j = 0; j < _SHMEM_ALLTOALL_SYNC_SIZE; j++) {
      twohop->row_sync[i][j] = _SHMEM_SYNC_VALUE;
      twohop->col_sync[i][j] = _SHMEM_SYNC_VALUE;
    }
#endif

  return &twohop->convey;
}

#endif
