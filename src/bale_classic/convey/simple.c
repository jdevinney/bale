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
#include "simple.h"


/*** UTILITY FUNCTIONS ***/

#if MPP_USE_UPC || (MPP_USE_SHMEM && !HAVE_SHMEMX_ALLTOALLV)
static void
make_perm(int n, int perm[], brand_t* prng)
{
  for (int i = 0; i < n; i++) {
    int j = brand(prng) % (i + 1);
    int k = (j < i) ? perm[j] : i;
    perm[j] = i;
    perm[i] = k;
  }
}

static int*
simple_perm_new(int n_procs, int my_proc)
{
  int* perm = malloc(n_procs * sizeof(int));
  int* pi1 = malloc(n_procs * sizeof(int));
  int* pi2 = malloc(n_procs * sizeof(int));

  if (perm && pi1 && pi2) {
    brand_t _prng;
    brand_init(&_prng, 2017);  // common seed
    make_perm(n_procs, pi1, &_prng);
    make_perm(n_procs, pi2, &_prng);
    for (int i = 0; i < n_procs; i++)
      perm[i] = pi2[(i + pi1[my_proc]) % n_procs];
  }
  else {
    free(perm);
    perm = NULL;
  }

  free(pi2);
  free(pi1);
  return perm;
}
#endif

int
simple_alltoallv(simple_t* simple)
{
  const size_t n_procs = simple->convey.n_procs;
  int error = 0;

  // Convert areas to sizes
  uint64_t n_bytes = 0;
  char* buffer = simple->send_buffers;
  for (size_t p = 0; p < n_procs; p++) {
    size_t size = simple->send[p].next - buffer;
    simple->send_sizes[p] = size;
    n_bytes += size;
    buffer += simple->buffer_bytes;
  }

  // Call the library function
  ACT_START(simple_a2a);
#if MPP_USE_UPC
  error = upcx_alltoallv(simple->all_recv_buffers, simple->all_recv_sizes,
                         simple->all_send_buffers, simple->all_send_sizes,
                         simple->offsets, simple->perm);
#else
# if HAVE_MPP_UTIL
  if (simple->a2a) {
    error = mpp_alltoallv_simple(simple->a2a, simple->recv_buffers, simple->recv_sizes,
                                 simple->buffer_bytes, MPP_A2A_BINNED,
                                 simple->send_buffers, simple->send_sizes);
  }
  else
# endif
  {
    CONVEY_PROF_DECL(_sample);
    CONVEY_PROF_START(&_sample);
# if HAVE_SHMEMX_ALLTOALLV
    shmemx_alltoallv(simple->recv_buffers, simple->offsets, simple->recv_sizes,
                     simple->send_buffers, simple->offsets, simple->send_sizes,
                     simple->comm.start, _trailz(simple->comm.stride), n_procs,
                     simple->sync[simple->flip]);
    simple->flip ^= 1;
# elif !HAVE_MPP_UTIL && MPP_USE_SHMEM
    error = xshmem_alltoallv(simple->recv_buffers, simple->recv_sizes,
                             simple->send_buffers, simple->send_sizes,
                             simple->offsets, simple->perm);
# elif MPP_USE_MPI
    int* send_sizes = malloc(2 * n_procs * sizeof(int));
    for (int64_t p = 0; p < n_procs; p++)
      send_sizes[p] = (int) simple->send_sizes[p];
    int* recv_sizes = send_sizes + n_procs;

    error = MPI_Alltoall(send_sizes, 1, MPI_INT, recv_sizes, 1, MPI_INT,
                         mpp_comm_mpi(simple->comm));
    if (!error) {
      error = MPI_Alltoallv(simple->send_buffers, send_sizes, simple->offsets, MPI_BYTE,
                            simple->recv_buffers, recv_sizes, simple->offsets, MPI_BYTE,
                            mpp_comm_mpi(simple->comm));
      for (int64_t p = 0; p < n_procs; p++)
        simple->recv_sizes[p] = (size_t) recv_sizes[p];
    }
    free(send_sizes);
# elif MPP_NO_MIMD
    memcpy(simple->recv_buffers, simple->send_buffers, simple->send_sizes[0]);
    simple->recv_sizes[0] = simple->send_sizes[0];
# else
    error = 1;
# endif
    CONVEY_PROF_STOP(&_sample, PROF_OP_A2A, -1, n_bytes);
  }
#endif
  ACT_STOP(simple_a2a);
  if (error)
    return convey_error_COMMS;

  // Convert sizes to areas
  buffer = simple->recv_buffers;
  for (size_t p = 0; p < n_procs; p++) {
    area_t* area = simple->recv + p;
    area->next = buffer;
    area->limit = buffer + simple->recv_sizes[p];
    buffer += simple->buffer_bytes;
  }
  simple->pull_from = 0;
  simple->stats[convey_COMMS] += n_procs;
  simple->stats[convey_SYNCS] += 1;
  return convey_OK;
}


/*** SIMPLE METHODS (minimal error checking) ***/

static int
simple_push(convey_t* self, const void* item, int64_t pe)
{
  simple_t* simple = (simple_t*) self;
  simple->nonempty = true;
  if (simple->sorter) {
    // Ensure that we don't push if a buffer might be full
    if (simple->overflow)
      return convey_FAIL;
    if (! sorter_push(simple->sorter, item, pe))
      simple->overflow = true;
    // The sorter_push always succeeds
    return convey_OK;
  }

  // Original implementation:
  const size_t item_size = self->item_size;
  area_t* area = simple->send + pe;
  if (area->next < area->limit) {
    memcpy(area->next, item, item_size);
    area->next += item_size;
    simple->stats[convey_PUSHES]++;
    return convey_OK;
  }
  simple->overflow = true;
  return convey_FAIL;
}

static void*
simple_apull(convey_t* self, int64_t* from)
{
  const int64_t n_procs = self->n_procs;
  simple_t* simple = (simple_t*) self;
  int64_t pe = simple->pull_from;

  while (pe < n_procs) {
    area_t* area = simple->recv + pe;
    char* here = area->next;
    if (here < area->limit) {
      if (from)
        *from = pe;
      area->next = here + self->item_size;
      simple->stats[convey_PULLS]++;
      return here;
    }
    simple->pull_from = ++pe;
  }

  return NULL;
}

static int
simple_pull(convey_t* self, void* item, int64_t* from)
{
  void* source = simple_apull(self, from);
  if (source == NULL)
    return convey_FAIL;
  memcpy(item, source, self->item_size);
  return convey_OK;
}

static int
simple_unpull(convey_t* self)
{
  const int64_t n_procs = self->n_procs;
  simple_t* simple = (simple_t*) self;
  const int64_t pe = simple->pull_from;
  if (pe >= n_procs)
    return convey_FAIL;

  char* bottom = simple->recv_buffers + pe * simple->buffer_bytes;
  char* next = simple->recv[pe].next;
  if (next > bottom) {
    simple->recv[pe].next = next - self->item_size;
    simple->stats[convey_UNPULLS]++;
    return convey_OK;
  }
  return convey_FAIL;
}

void
simple_reset_send_buffers(simple_t* simple)
{
  const size_t n_procs = simple->convey.n_procs;
  const size_t buffer_bytes = simple->buffer_bytes;
  for (size_t p = 0; p < n_procs; p++) {
    simple->send[p].next = simple->send_buffers + p * buffer_bytes;
    simple->send[p].limit = simple->send[p].next + simple->buffer_limit;
  }
  simple->nonempty = false;
  simple->overflow = false;
}

static int
simple_advance(convey_t* self, bool done)
{
  simple_t* simple = (simple_t*) self;
  simple->stats[convey_ADVANCES]++;

  // We are in state WORKING or ENDGAME, and need to decide collectively
  // what to do.  The simple->quiet flag is the same on every process.

  bool need_pull = (simple->pull_from < self->n_procs);
  if (simple->quiet)
    return need_pull ? convey_NEAR : convey_DONE;

  // If using a sorter, we aren't really done until the sorter is flushed.
  bool steady = (self->features & convey_STEADY);
  if (simple->sorter && (done | steady)) {
    if (!simple->overflow) {
      int status = sorter_flush(simple->sorter);
      simple->nonempty |= (status != 0);
      simple->overflow = (status < 0);
      done &= (status >= 0);
    }
    else
      done = false;
  }

  simple->stats[convey_SYNCS]++;
  bool cannot_comm = need_pull;
  bool desire_comm = (done | simple->overflow) || (steady && simple->nonempty);
  bool further_comm = !done;
  long status = mpp_or_long((cannot_comm << 0) + (desire_comm << 1) + (further_comm << 2));
  cannot_comm = (status >> 0) & 1;
  desire_comm = (status >> 1) & 1;
  further_comm = (status >> 2) & 1;
  if (cannot_comm || !desire_comm)
    return convey_OK;

  // This call establishes the recv areas
  int result = simple_alltoallv(simple);
  if (result != convey_OK)
    return result;
  simple_reset_send_buffers(simple);

  simple->quiet = !further_comm;
  return simple->quiet ? convey_NEAR : convey_OK;
}

static void
reset_recv_buffers(simple_t* simple)
{
  const size_t n_procs = simple->convey.n_procs;
  const size_t buffer_bytes = simple->buffer_bytes;
  for (size_t p = 0; p < n_procs; p++) {
    simple->recv[p].next = simple->recv_buffers + p * buffer_bytes;
    simple->recv[p].limit = simple->recv[p].next;
  }
  simple->pull_from = n_procs;
}

static bool
alloc_buffers(simple_t* simple)
{
  size_t n_bytes = simple->convey.n_procs * simple->buffer_bytes;
  convey_alc8r_t* alloc = &simple->alloc;
  PARALLEL_ALLOC(simple, send_buffers, alloc, n_bytes, char);
  PARALLEL_ALLOC(simple, recv_buffers, alloc, n_bytes, char);
  return simple->send_buffers && simple->recv_buffers;
}

static void
dealloc_buffers(simple_t* simple)
{
  convey_alc8r_t* alloc = &simple->alloc;
  PARALLEL_DEALLOC(simple, recv_buffers, alloc);
  PARALLEL_DEALLOC(simple, send_buffers, alloc);
  simple->send_buffers = NULL;
  simple->recv_buffers = NULL;
}

// Alignment can be ignored; as long as the buffer allocations are
// maximally aligned, all items will be fully aligned.
static int
simple_begin(convey_t* self, size_t item_size, size_t align)
{
  simple_t* simple = (simple_t*) self;
  if (!mpp_comm_is_equal(MPP_COMM_CURR, simple->comm))
    return convey_error_TEAM;
  if (item_size > simple->buffer_bytes)
    return convey_error_OFLO;
  size_t old_item_size = self->item_size;

  simple->stats[convey_BEGINS]++;
  self->item_size = item_size;
  simple->capacity = simple->buffer_bytes / item_size;
  simple->buffer_limit = item_size * simple->capacity;

  if (simple->scatter && (!simple->sorter || item_size != old_item_size)) {
    if (simple->sorter)
      sorter_free(simple->sorter);
    const size_t n_procs = simple->convey.n_procs;
#if 0
    // Tree sorter: We want at least 10 items per child, and an odd number overall
    const size_t bucket_depth = 333;
    simple->sorter = sorter_new(n_procs, simple->send, item_size, bucket_depth,
                                alloc, simple->dynamic);
#else
    // Stupid sorter: Use a circular buffer of 16 items; allocator is ignored
    const size_t circle = 16;
    simple->sorter = sorter_new(n_procs, simple->send, item_size, circle, NULL, false);
#endif
    if (!simple->sorter || !sorter_setup(simple->sorter))
      return convey_error_ALLOC;
  }

  if (simple->dynamic && !alloc_buffers(simple))
    return convey_error_ALLOC;
  simple_reset_send_buffers(simple);
  reset_recv_buffers(simple);
  simple->quiet = false;
  mpp_barrier(1);
  return convey_OK;
}

static int
simple_reset(convey_t* self)
{
  simple_t* simple = (simple_t*) self;
  if (!mpp_comm_is_equal(MPP_COMM_CURR, simple->comm))
    return convey_error_TEAM;
  if (simple->sorter)
    sorter_reset(simple->sorter);
  if (simple->dynamic)
    dealloc_buffers(simple);
  convey_imp_update_stats(simple->stats);
  return convey_OK;
}

static void
free_allocs(simple_t* simple)
{
  // Make sure allocations are last-in, first-out
  dealloc_buffers(simple);
  void* alc8r = simple->alloc.alc8r;
  void (*dealloc)(void*,void*) = simple->alloc.free;
  dealloc(alc8r, simple->sync[1]);
  dealloc(alc8r, simple->sync[0]);
  dealloc(alc8r, simple->offsets);
  dealloc(alc8r, simple->recv);
  dealloc(alc8r, simple->send);
  PARALLEL_DEALLOC(simple, recv_sizes, &simple->alloc);
  PARALLEL_DEALLOC(simple, send_sizes, &simple->alloc);
#if HAVE_MPP_UTIL
  if (simple->own_a2a)
    mpp_alltoall_destroy((void*) simple->a2a);
#endif
  free(simple->perm);
  free(simple);
}

static int
simple_free(convey_t* self)
{
  simple_t* simple = (simple_t*) self;
  if (!mpp_comm_is_equal(MPP_COMM_CURR, simple->comm))
    return convey_error_TEAM;
  sorter_free(simple->sorter);
  free_allocs(simple);
  return convey_OK;
}

static int64_t
simple_statistic(convey_t* self, int which)
{
  simple_t* simple = (simple_t*) self;
  return convey_imp_statistic(simple->stats, which);
}

static const convey_methods_t simple_methods = {
  .push = &simple_push,
  .pull = &simple_pull,
  .unpull = &simple_unpull,
  .advance = &simple_advance,
  .begin = &simple_begin,
  .reset = &simple_reset,
  .free = &simple_free,
  .statistic = &simple_statistic,
  .apull = &simple_apull,
  .epush = &convey_no_epush,
  .epull = &convey_no_epull,
  .panic = &convey_imp_panic,
};

static const convey_methods_t debug_methods = {
  ._super_ = &simple_methods,
  .push = &convey_checked_push,
  .pull = &convey_checked_pull,
  .unpull = &convey_checked_unpull,
  .advance = &simple_advance,
  .begin = &simple_begin,
  .reset = &simple_reset,
  .free = &simple_free,
  .statistic = &simple_statistic,
  .apull = &convey_checked_apull,
  .epush = &convey_no_epush,
  .epull = &convey_no_epull,
  .panic = &convey_imp_panic,
};


/*** PUBLIC INTERFACE ***/

convey_t*
convey_new_simple(size_t buffer_bytes,
                  const convey_alc8r_t* alloc, const convey_mpp_a2a_t* a2a,
                  uint64_t options)
{
  bool quiet = (options & convey_opt_QUIET);
  if (buffer_bytes == 0)
    CONVEY_REJECT(quiet, "invalid capacity");
  // Round up to guarantee proper alignment
  buffer_bytes = (buffer_bytes + CONVEY_MAX_ALIGN - 1) & -CONVEY_MAX_ALIGN;
#if MPP_USE_SHMEM
  if (alloc == NULL && !mpp_comm_is_world(MPP_COMM_CURR))
    CONVEY_REJECT(quiet, "SHMEM cannot allocate symmetric memory for a team");
#endif

  convey_mpp_a2a_t* own_a2a = NULL;
#if !HAVE_MPP_UTIL
  a2a = NULL;
#elif MPP_USE_SHMEM && !HAVE_SHMEMX_ALLTOALLV
  if (a2a == NULL) {
    own_a2a = mpp_alltoall_create(5);
    int error = own_a2a ? mpp_alltoall_init(own_a2a, 0, 0, 0, 0, 0) : 1;
    if (error) {
      mpp_alltoall_destroy(own_a2a);
      CONVEY_REJECT(quiet, "unable to create or initialize a default mpp_alltoall");
    }
  }
#elif MPP_USE_MPI
  if (a2a == NULL && buffer_bytes * PROCS > INT_MAX)
    CONVEY_REJECT(quiet, "memory footprint is too large for MPI (and a2a is NULL)");
#endif

  if (alloc == NULL)
    alloc = &convey_imp_alloc_align;
  else if (!alloc->grab || !alloc->free)
    CONVEY_REJECT(quiet, "alloc is missing one or both methods");

  simple_t* simple = malloc(sizeof(simple_t));
  if (simple == NULL)
    CONVEY_REJECT(quiet, "a small malloc() failed!");
  // Initialize some things; set scratch pointers to NULL
  size_t n_procs = PROCS;
  *simple = (simple_t) {
    .convey = {
      ._class_ = (options & convey_opt_RECKLESS) ? &simple_methods : &debug_methods,
      .features = (options & convey_opt_PROGRESS) ? convey_STEADY : 0,
      .n_procs = n_procs, .state = convey_DORMANT, .suppress = UINT64_C(0) - quiet,
    },
    .buffer_bytes = buffer_bytes,
    .alloc = *alloc, .a2a = (own_a2a ? own_a2a : a2a), .own_a2a = (own_a2a != NULL),
    .dynamic = (options & convey_opt_DYNAMIC),
    .scatter = (options & convey_opt_SCATTER),
  };

  // Try to do all the allocations
#define GRAB(SIZE) alloc->grab(alloc->alc8r, (SIZE), __FILE__, __LINE__)
  PARALLEL_ALLOC(simple, send_sizes, alloc, n_procs, size_t);
  PARALLEL_ALLOC(simple, recv_sizes, alloc, n_procs, size_t);
  simple->send = GRAB(n_procs * sizeof(area_t));
  simple->recv = GRAB(n_procs * sizeof(area_t));
  bool ok = (simple->send && simple->send_sizes && simple->recv && simple->recv_sizes);
  if (a2a == NULL) {
    simple->offsets = GRAB(n_procs * sizeof(a2a_off_t));
    ok &= (simple->offsets != NULL);
#if HAVE_SHMEMX_ALLTOALLV
    for (int i = 0; i < 2; i++) {
      simple->sync[i] = GRAB(_SHMEM_ALLTOALL_SYNC_SIZE * sizeof(long));
      ok &= (simple->sync[i] != NULL);
    }
#endif
  }
#undef GRAB
#if MPP_USE_UPC || (MPP_USE_SHMEM && !HAVE_SHMEMX_ALLTOALLV)
  simple->perm = simple_perm_new(n_procs, MY_PROC);
  ok &= (simple->perm != NULL);
#endif

  if (!simple->dynamic)
    ok &= alloc_buffers(simple);
  if (!ok) {
    free_allocs(simple);
    CONVEY_REJECT(quiet, "symmetric allocation failed");
  }

  // Initialize the remaining parts of the conveyor
  if (a2a == NULL) {
    for (size_t p = 0; p < n_procs; p++)
      simple->offsets[p] = p * buffer_bytes;
#if HAVE_SHMEMX_ALLTOALLV
    for (int i = 0; i < 2; i++)
      for (int64_t j = 0; j < _SHMEM_ALLTOALL_SYNC_SIZE; j++)
        simple->sync[i][j] = _SHMEM_SYNC_VALUE;
#endif
  }

  simple->comm = MPP_COMM_CURR;
  return &simple->convey;
}
