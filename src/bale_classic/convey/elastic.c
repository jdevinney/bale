// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#include "convey_impl.h"
#include "private.h"


#if MPP_NO_MIMD

convey_t*
convey_new_etensor(size_t buffer_bytes, size_t monster_bytes,
                   int order, size_t n_local, size_t n_buffers,
                   const convey_alc8r_t* alloc, uint64_t options)
{
  return convey_new_trivial(monster_bytes, alloc, options);
}

#else

#include "porter.h"
#include "tensor.h"

typedef struct ticket {
  size_t size;
  long id;
} ticket_t;

typedef struct elastic {
  tensor_t tensor;
  void* aligned_item;
  size_t max_size;  // maximum item size, in bytes
  size_t big_size;  // above this size, must use these buffers:
  PARALLEL(char*, outgoing);   // holds one large item being pushed
  PARALLEL(char*, incoming);   // holds one large item being pulled
  PARALLEL(long*, bigrcvd);    // how many large items have been received
  long bigsent;     // how many large items have been pushed
  char* backup;     // packet containing most recently pulled item
  bool bigpull;     // was the last pull a monster item?
  bool unpull;      // did we put it back?
  bool dynamic;
#if MPP_USE_MPI
  MPI_Request big_req;
  MPI_Comm comm;
  int convey_id;    // tag for monster messages
#endif
  convey_alc8r_t alloc;
} elastic_t;

static bool
alloc_monster_buffers(elastic_t* elastic)
{
  // We may not need the large buffers at all
  if (elastic->max_size <= elastic->big_size) {
    elastic->incoming = NULL;
    elastic->outgoing = NULL;
    return true;
  }
  size_t n_bytes = elastic->max_size;
  convey_alc8r_t* alloc = &elastic->alloc;
  PARALLEL_ALLOC(elastic, outgoing, alloc, n_bytes, char);
  PARALLEL_ALLOC(elastic, incoming, alloc, n_bytes, char);
  return elastic->outgoing && elastic->incoming;
}

static void
dealloc_monster_buffers(elastic_t* elastic)
{
  PARALLEL_DEALLOC(elastic, incoming, &elastic->alloc);
  PARALLEL_DEALLOC(elastic, outgoing, &elastic->alloc);
  elastic->incoming = NULL;
  elastic->outgoing = NULL;
}


/*** Private Methods ***/

static bool
elastic_pivot_mid(tensor_t* matrix, buffer_t* buffer)
{
  // tag is (y'); we are (x,y'); source is x, tag becomes (x)
  ACT_START(matrix_pivot);
  const uint64_t tag = buffer->source;
  uint32_t* packet = (uint32_t*) &buffer->data[buffer->start];
  uint32_t* limit = (uint32_t*) &buffer->data[buffer->limit];
  while (packet < limit) {
    int dest = packet[0];
    uint32_t descr = packet[1];
    bool ok = porter_epush(matrix->porters[1], tag, descr, &packet[2], dest);
    if (!ok) {
      buffer->start = (char*)packet - buffer->data;
      ACT_STOP(matrix_pivot);
      return false;
    }
    packet += PORTER_QUADS(descr);
  }
  buffer->start = buffer->limit;
  ACT_STOP(matrix_pivot);
  return true;
}

static bool
elastic_pivot_early(tensor_t* tensor, buffer_t* buffer)
{
  // tag is (x',z'); source is z; we are (x,y,y')
  // hop to (x',y,y'), tag becomes (z, z')
  ACT_START(tensor_early);
  const uint32_t source = buffer->source << 8;
  uint32_t* packet = (uint32_t*) &buffer->data[buffer->start];
  uint32_t* limit = (uint32_t*) &buffer->data[buffer->limit];
  while (packet < limit) {
    uint32_t tag = packet[0];
    int dest = tag >> 8;
    tag = (tag & 0xFF) | source;
    uint32_t descr = packet[1];
    bool ok = porter_epush(tensor->porters[1], tag, descr, &packet[2], dest);
    if (!ok) {
      buffer->start = (char*)packet - buffer->data;
      ACT_STOP(tensor_early);
      return false;
    }
    packet += PORTER_QUADS(descr);
  }
  buffer->start = buffer->limit;
  ACT_STOP(tensor_early);
  return true;
}

static bool
elastic_pivot_late(tensor_t* tensor, buffer_t* buffer)
{
  // tag is (z,z'), source is x, we are (x',y',y)
  // hop to (x',y',z'), tag becomes (x,z)
  ACT_START(tensor_late);
  const uint32_t source = buffer->source << 8;
  uint32_t* packet = (uint32_t*) &buffer->data[buffer->start];
  uint32_t* limit = (uint32_t*) &buffer->data[buffer->limit];
  while (packet < limit) {
    uint32_t tag = packet[0];
    int dest = tag & 0xFF;
    tag = (tag >> 8) | source;
    uint32_t descr = packet[1];
    bool ok = porter_epush(tensor->porters[2], tag, descr, &packet[2], dest);
    if (!ok) {
      buffer->start = (char*)packet - buffer->data;
      ACT_STOP(tensor_late);
      return false;
    }
    packet += PORTER_QUADS(descr);
  }
  buffer->start = buffer->limit;
  ACT_STOP(tensor_late);
  return true;
}


/*** MPI Support ***/

#if MPP_USE_MPI

#define CONVEY_BASE_ID 100  // avoid overlap with porter tags
static uint64_t used_convey_ids = 0;

static int
monster_can_send(elastic_t* elastic)
{
  MPI_Request* req = &elastic->big_req;
  if (*req == MPI_REQUEST_NULL)
    return convey_OK;
  MPI_Status _status;
  int ready = 0;
  int rc = MPI_Test(req, &ready, &_status);
  if (rc != MPI_SUCCESS)
    return convey_error_COMMS;
  return ready ? convey_OK : convey_FAIL;
}

#endif


/*** Monster Methods ***/

static int
monster_epush(elastic_t* elastic, size_t bytes, const void* item, int64_t pe)
{
#if MPP_USE_MPI
  int code = monster_can_send(elastic);
  if (code != convey_OK)
    return code;
#else
  if (elastic->bigsent > *elastic->bigrcvd)
    return convey_FAIL;
#endif
  tensor_t* tensor = &elastic->tensor;

  // Try to send a ticket through normal channels, but first, make sure
  // the item is ready to be retrieved by the recipient.
  memcpy(elastic->outgoing, item, bytes);
  route_t _route = tensor->router(tensor, pe);
  ticket_t _ticket = { .size = bytes, .id = elastic->bigsent + 1 };
  bool ok = porter_epush(tensor->porters[0], _route.tag,
                         PORTER_DESCR(sizeof(ticket_t), 1), &_ticket, _route.next);
  if (!ok)
    return convey_FAIL;

#if MPP_USE_MPI
  // mprint(MY_PROC, 0, "send %zu bytes to %ld\n", bytes, pe);
  int rc = MPI_Isend(elastic->outgoing, bytes, MPI_BYTE, pe,
                     elastic->convey_id, elastic->comm, &elastic->big_req);
  if (rc != MPI_SUCCESS)
    return convey_error_COMMS;
#endif

  tensor->stats[convey_PUSHES]++;
  elastic->bigsent++;
  return convey_OK;
}

static int
monster_epull(elastic_t* elastic, ticket_t* ticket, int64_t from, convey_item_t* result)
{
  size_t bytes = ticket->size;

  if (elastic->unpull)
    elastic->unpull = false;
  else {
    // Only do these things once per monster item
#if MPP_USE_SHMEM
    shmem_getmem(elastic->incoming, elastic->outgoing, bytes, from);
    int pe = mpp_rel_to_abs_proc(MPP_COMM_CURR, from);
    shmem_long_p(elastic->bigrcvd, ticket->id, pe);
#elif MPP_USE_UPC
    upc_memget(elastic->incoming, &elastic->all_outgoing[from], bytes);
    elastic->all_bigrcvd[from] = ticket->id;
#else
    MPI_Status _status;
    int rc = MPI_Recv(elastic->incoming, elastic->max_size, MPI_BYTE, from,
                      elastic->convey_id, elastic->comm, &_status);
    // mprint(MY_PROC, 0, "recv %zu bytes from %ld\n", bytes, from);
    if (rc != MPI_SUCCESS)
      return convey_error_COMMS;
    int count;
    MPI_Get_count(&_status, MPI_BYTE, &count);
    if (count != bytes)
      return convey_error_COMMS;
#endif
  }

  result->bytes = bytes;
  result->data = elastic->incoming;
  result->from = from;
  return convey_OK;
}

static int
monster_unpull(elastic_t* elastic)
{
  if (elastic->unpull)
    return convey_FAIL;

  buffer_t* buffer = elastic->tensor.buffer;
  const size_t ticket_quads = PORTER_QUADS(PORTER_DESCR(sizeof(ticket_t), 1));
  if (buffer == NULL || buffer->start < 4 * ticket_quads)
    return convey_FAIL;

  buffer->start -= 4 * ticket_quads;
  elastic->unpull = true;
  elastic->backup = NULL;
  return convey_OK;
}

static int
monster_begin(elastic_t* elastic)
{
  if (elastic->dynamic && !alloc_monster_buffers(elastic))
    return convey_error_ALLOC;
  elastic->bigsent = 0;
  *elastic->bigrcvd = 0;
  elastic->backup = NULL;
  elastic->bigpull = false;
  elastic->unpull = false;
#if MPP_USE_MPI
  elastic->big_req = MPI_REQUEST_NULL;
#endif
  mpp_barrier(1);
  return convey_OK;
}

static int
monster_reset(elastic_t* elastic)
{
  if (elastic->dynamic)
    dealloc_monster_buffers(elastic);
  return convey_OK;
}

static void
monster_free(elastic_t* elastic)
{
#if MPP_USE_MPI
  if (elastic->convey_id >= 0) {
    int bit = elastic->convey_id - CONVEY_BASE_ID;
    used_convey_ids &= ~(UINT64_C(1) << bit);
  }
#endif
  dealloc_monster_buffers(elastic);
  PARALLEL_DEALLOC(elastic, bigrcvd, &elastic->alloc);
}


/*** Public Methods ***/

static int
elastic_epush(convey_t* self, size_t bytes, const void* item, int64_t pe)
{
  elastic_t* elastic = (elastic_t*) self;
  if (bytes > elastic->big_size) {
    if (bytes > elastic->max_size)
      return convey_imp_panic(self, __func__, convey_error_OFLO);
    return monster_epush(elastic, bytes, item, pe);
  }

  tensor_t* tensor = &elastic->tensor;
  route_t _route = tensor->router(tensor, pe);
  bool ok = porter_epush(tensor->porters[0], _route.tag,
                         PORTER_DESCR(bytes, 0), item, _route.next);
  tensor->stats[convey_PUSHES] += ok;
  return ok ? convey_OK : convey_FAIL;
}

static int
elastic_epull(convey_t* self, convey_item_t* result)
{
  elastic_t* elastic = (elastic_t*) self;
  tensor_t* tensor = (tensor_t*) self;
  buffer_t* buffer = tensor->buffer;
  const int order = tensor->order;

  if (buffer && buffer->start == buffer->limit) {
    porter_return(tensor->porters[order - 1]);
    elastic->backup = NULL;
    buffer = NULL;
  }
  if (!buffer) {
    buffer = porter_borrow(tensor->porters[order - 1]);
    tensor->buffer = buffer;
    if (!buffer) {
      result->data = NULL;
      return convey_FAIL;
    }
  }

  uint32_t* packet = (uint32_t*) &buffer->data[buffer->start];
  uint32_t descr = packet[1];
  buffer->start += 4 * PORTER_QUADS(descr);
  int64_t from = origin_from_tag(tensor, order, packet[0], buffer->source);
  bool monster = PORTER_TICKET(descr);
  elastic->bigpull = monster;
  elastic->backup = (char*) packet;

  if (monster) {
    ticket_t _ticket;
    memcpy(&_ticket, &packet[2], sizeof(ticket_t));
    return monster_epull(elastic, &_ticket, from, result);
  }
  uint32_t n_bytes = PORTER_BYTES(descr);
  result->from = from;
  result->bytes = n_bytes;
  result->data = pull_pointer(&packet[2], elastic->aligned_item, n_bytes);
  return convey_OK;
}

static int
elastic_push(convey_t* self, const void* item, int64_t pe)
{
  return elastic_epush(self, self->item_size, item, pe);
}

static int
elastic_pull(convey_t* self, void* item, int64_t* from)
{
  convey_item_t _item;
  int rc = self->_class_->epull(self, &_item);
  if (rc == convey_OK) {
    if (_item.bytes != self->item_size)
      return convey_imp_panic(self, __func__, convey_error_MISFIT);
    if (from)
      *from = _item.from;
    memcpy(item, _item.data, _item.bytes);
  }  
  return rc;
}

static void*
elastic_apull(convey_t* self, int64_t* from)
{
  convey_item_t _item;
  int rc = self->_class_->epull(self, &_item);
  if (rc <= 0)
    return NULL;
  if (_item.bytes != self->item_size) {
    convey_imp_panic(self, __func__, convey_error_MISFIT);
    return NULL;
  }
  if (from)
    *from = _item.from;
  return _item.data;
}

static int
elastic_unpull(convey_t* self)
{
  elastic_t* elastic = (elastic_t*) self;
  if (elastic->bigpull)
    return monster_unpull(elastic);

  tensor_t* tensor = (tensor_t*) self;
  buffer_t* buffer = tensor->buffer;
  if (buffer && buffer->start > 0 && elastic->backup) {
    buffer->start = elastic->backup - buffer->data;
    elastic->backup = NULL;
    tensor->stats[convey_UNPULLS]++;
    return convey_OK;
  }
  return convey_FAIL;
}

static int
elastic_advance(convey_t* self, bool done)
{
  int code = tensor_advance(self, done);

  elastic_t* elastic = (elastic_t*) self;
#if MPP_USE_MPI
  if (code == convey_DONE || code == convey_NEAR) {
    int ready = monster_can_send(elastic);
    if (ready != convey_OK) {
      mprint(MY_PROC, 0, "tensor is %d but ready is %d\n", code, ready);
      code = (ready < 0) ? ready : convey_OK;
    }
  }
#endif
  if (code == convey_DONE && elastic->unpull) 
    code = convey_NEAR;
  return code;
}

static int
elastic_begin(convey_t* self, size_t item_size, size_t align)
{
  elastic_t* elastic = (elastic_t*) self;
  size_t max_size = (elastic->big_size + align - 1) & -align;
  if (! convey_prep_aligned(&elastic->aligned_item, max_size, 4, align))
    return convey_error_ALLOC;

  int a = tensor_begin(self, item_size, align);
  int b = monster_begin((elastic_t*) self);
  return (a < 0) ? a : b;
}

static int
elastic_reset(convey_t* self)
{
  elastic_t* elastic = (elastic_t*) self;
  free(elastic->aligned_item);
  elastic->aligned_item = NULL;
  mpp_barrier(1);

  int a = monster_reset((elastic_t*) self);
  int b = tensor_reset(self);
  return (a < 0) ? a : b;
}

static int
elastic_free(convey_t* self)
{
  // FIXME: check for correct team?
  elastic_t* elastic = (elastic_t*) self;
  free(elastic->aligned_item);
  mpp_barrier(1);
  monster_free((elastic_t*) self);
  return tensor_free(self);
}

// FIXME: Do we want to keep statistics on monster packets?
static int64_t
elastic_statistic(convey_t* self, int which)
{
  return tensor_statistic(self, which);
}

static const convey_methods_t elastic_fast_methods = {
  .push = &elastic_push,
  .pull = &elastic_pull,
  .unpull = &elastic_unpull,
  .advance = &elastic_advance,
  .begin = &elastic_begin,
  .reset = &elastic_reset,
  .free = &elastic_free,
  .statistic = &elastic_statistic,
  .apull = &elastic_apull,
  .epush = &elastic_epush,
  .epull = &elastic_epull,
  .panic = &convey_imp_panic,
};

static const convey_methods_t elastic_debug_methods = {
  ._super_ = &elastic_fast_methods,
  .push = &convey_checked_push,
  .pull = &convey_checked_pull,
  .unpull = &convey_checked_unpull,
  .advance = &elastic_advance,
  .begin = &elastic_begin,
  .reset = &elastic_reset,
  .free = &elastic_free,
  .statistic = &elastic_statistic,
  .apull = &convey_checked_apull,
  .epush = &convey_checked_epush,
  .epull = &convey_checked_epull,
  .panic = &convey_imp_panic,
};


/*** Constructor ***/

convey_t*
convey_new_etensor(size_t buffer_bytes, size_t monster_bytes,
                   int order, size_t n_local, size_t n_buffers,
                   const convey_alc8r_t* alloc, uint64_t options)
{
  const bool quiet = (options & convey_opt_QUIET);
  if (buffer_bytes < 16)
    buffer_bytes = 16;

  // Must ensure that the buffers are at least of size buffer_bytes.
  uint64_t packet_quads = 4;
  uint64_t capacity = 1 + (buffer_bytes - 1) / (4 * packet_quads);
  if ((capacity * packet_quads) >> 29)
    CONVEY_REJECT(quiet, "buffer_bytes is too large");
  size_t big_bytes = (capacity * packet_quads - 2) * sizeof(uint32_t);
  big_bytes = MIN(big_bytes, monster_bytes);

  // Build the underlying tensor conveyor, without compression
  uint64_t tensor_opts = options - (options & convey_opt_COMPRESS);
  tensor_opts |= convey_opt_STANDARD;
  convey_t* self = convey_new_tensor(4 * packet_quads * capacity, order, n_local,
                                     n_buffers, alloc, tensor_opts);
  if (self == NULL)
    return NULL;

  // Extend it to an elastic conveyor
  elastic_t* elastic = realloc(self, sizeof(elastic_t));
  if (elastic == NULL) {
    convey_free(self);
    CONVEY_REJECT(quiet, "a small realloc() failed");
  }
  tensor_t* tensor = &elastic->tensor;
  self = &tensor->convey;
  order = tensor->order;

  // Adjust the fields of the tensor conveyor
  bool reckless = (options & convey_opt_RECKLESS);
  self->_class_ = reckless ? &elastic_fast_methods : &elastic_debug_methods;
  self->features |= convey_ELASTIC;
  if (order == 2)
    tensor->pivots[0] = &elastic_pivot_mid;
  else if (order == 3) {
    tensor->pivots[0] = &elastic_pivot_early;
    tensor->pivots[1] = &elastic_pivot_late;
  }

  // Fill in the unchanging fields of the elastic conveyor
  bool dynamic = (options & convey_opt_DYNAMIC);
  // This allocator is used for monster buffers, not for porters
  if (alloc == NULL)
    alloc = &convey_imp_alloc_align;

  elastic->aligned_item = NULL;
  elastic->max_size = monster_bytes;
  elastic->big_size = big_bytes;
  PARALLEL_NULLIFY(elastic, outgoing);
  PARALLEL_NULLIFY(elastic, incoming);
  elastic->bigrcvd = NULL;
  elastic->dynamic = dynamic;
  elastic->alloc = *alloc;

#if MPP_USE_MPI
  elastic->comm = mpp_comm_mpi(MPP_COMM_CURR);
  uint64_t avail = mpp_and_long(~used_convey_ids);
  if (avail) {
    int bit = _trailz(avail);
    used_convey_ids |= UINT64_C(1) << bit;
    elastic->convey_id = bit + CONVEY_BASE_ID;
  }
  else {
    elastic->convey_id = -1;
    monster_free(elastic);
    tensor_free(self);
    CONVEY_REJECT(quiet, "conveyor IDs used up");
  }
#endif

  // Do the necessary allocations
  PARALLEL_ALLOC(elastic, bigrcvd, alloc, 1, long);
  bool ok = (elastic->bigrcvd != NULL);
  if (!dynamic)
    ok &= alloc_monster_buffers(elastic);
  if (!ok) {
    monster_free(elastic);
    tensor_free(self);
    CONVEY_REJECT(quiet, "symmetric allocation failed");
  }

  // Everything else will be set up by elastic_begin()
  return self;
}

#endif
