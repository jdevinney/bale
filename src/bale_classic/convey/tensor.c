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
convey_new_tensor(size_t capacity, int order, size_t n_local,
                  size_t n_buffers, const convey_alc8r_t* alloc, uint64_t options)
{
  return convey_new_simple(capacity, alloc, NULL, options);
}

#else

/*** INTERNAL METHODS ***/

#include "porter.h"
#include "tensor.h"

#define ROUTER_HINT  // empty: routers not inlined here
#include "router.h"

static inline porter_t*
remote_porter(tensor_t* tensor)
{
  return tensor->porters[tensor->order > (2 - MATRIX_REMOTE_HOP)];
}

static bool
setup_methods(tensor_t* tensor)
{
  convey_t* self = &tensor->convey;
  const int order = tensor->order;
  const size_t item_bytes = self->item_size;
  const uint8_t* tags = tensor->tag_bytes;
  
  push_f* push = tensor_select_push(order, tags[0], item_bytes);
  pull_f* pull = tensor_select_pull(order, tags[order - 1], item_bytes);
  if (push || pull) {
    convey_methods_t* xlr8 = malloc(sizeof(convey_methods_t));
    if (xlr8 == NULL)
      return false;
    *xlr8 = *(self->_class_);
    xlr8->_super_ = self->_class_;
    xlr8->dynamic = true;
    if (push)
      xlr8->push = push;
    if (pull)
      xlr8->pull = pull;
    self->_class_ = xlr8;
  }

  if (order == 2) {
    size_t t = tags[MATRIX_REMOTE_HOP];
    tensor->pivots[0] = tensor_select_pivot_mid(t, item_bytes);
  }
  else if (order == 3) {
    tensor->pivots[0] = tensor_select_pivot_early(tags[1], item_bytes);
    tensor->pivots[1] = tensor_select_pivot_late(tags[1], item_bytes);
  }

  return true;
}

static void
reset_methods(tensor_t* tensor)
{
  convey_t* self = &tensor->convey;
  const convey_methods_t* class = self->_class_;
  if (class->dynamic) {
    const convey_methods_t* standard = class->_super_;
    free((void*) class);
    self->_class_ = standard;
  }
  tensor->pivots[0] = NULL;
  tensor->pivots[1] = NULL;
}


/*** TENSOR METHODS ***/

static int
tensor_push(convey_t* self, const void* item, int64_t pe)
{
  tensor_t* tensor = (tensor_t*) self;
  route_t _route = tensor->router(tensor, pe);
  bool ok = porter_push(tensor->porters[0], _route.tag, item, _route.next);
  tensor->stats[convey_PUSHES] += ok;
  return ok ? convey_OK : convey_FAIL;
}

static void*
tensor_upull(convey_t* self, int64_t* from)
{
  tensor_t* tensor = (tensor_t*) self;
  buffer_t* buffer = tensor->buffer;
  const int order = tensor->order;

  if (buffer && buffer->start == buffer->limit) {
    porter_return(tensor->porters[order - 1]);
    buffer = NULL;
  }
  if (!buffer) {
    buffer = porter_borrow(tensor->porters[order - 1]);
    tensor->buffer = buffer;
    if (!buffer)
      return NULL;
  }

  char* packet = &buffer->data[buffer->start];
  buffer->start += tensor->packet_bytes;
  if (from) {
    uint32_t source = buffer->source;
    if (order == 1)
      *from = source;
    else {
      uint32_t tag;
#if MATRIX_REMOTE_HOP == 1
      if (order == 2 && tensor->item_offset == 1)
        tag = *(uint8_t*)packet;
      else
#endif
        tag = *(uint32_t*)packet;
      *from = origin_from_tag(tensor, order, tag, source);
    }
  }
  tensor->stats[convey_PULLS]++;
  return packet + tensor->item_offset;
}

static int
tensor_pull(convey_t* self, void* item, int64_t* from)
{
  void* source = tensor_upull(self, from);
  if (source == NULL)
    return convey_FAIL;
  memcpy(item, source, self->item_size);
  return convey_OK;
}

static void*
tensor_apull(convey_t* self, int64_t* from)
{
  tensor_t* tensor = (tensor_t*) self;
  void* source = tensor_upull(self, from);
  return pull_pointer(source, tensor->aligned_item, self->item_size);
}

static int
tensor_unpull(convey_t* self)
{
  tensor_t* tensor = (tensor_t*) self;
  buffer_t* buffer = tensor->buffer;
  if (buffer && buffer->start > 0) {
    size_t prev = buffer->start - tensor->packet_bytes;
    buffer->start = prev;
    tensor->stats[convey_UNPULLS]++;
    return convey_OK;
  }
  return convey_FAIL;
}

int
tensor_advance(convey_t* self, bool done)
{
  tensor_t* tensor = (tensor_t*) self;
  tensor->stats[convey_ADVANCES]++;
  // We won't be called if state is already COMPLETE

  const int order = tensor->order;
  // No shortcuts: must advance every porter to ensure progress
  for (int k = tensor->n_complete; k < order; k++) {
    bool go = porter_advance(tensor->porters[k], done);
    if (!go) {
      tensor->n_complete++;
      continue;  // done is true
    }
    done = false;
    if (k == order - 1)
      break;

    // Try to interleave work in a reasonable way
    for (int loop = 0; loop < 5; loop++) {
      buffer_t* buffer = porter_borrow(tensor->porters[k]);
      if (!buffer)
        break;
      go = (tensor->pivots[k])(tensor, buffer);
      porter_return(tensor->porters[k]);
      if (!go)
        break;
    }
  }

  if (done) {
    porter_t* porter = remote_porter(tensor);
    tensor->stats[convey_COMMS] = porter_get_stats(porter, 0);
    tensor->stats[convey_SYNCS] = porter_get_stats(porter, 1);
    tensor->stats[convey_BYTES] = porter_get_stats(porter, 2);
  }
  return done ? convey_DONE : convey_OK;
}

int
tensor_begin(convey_t* self, size_t item_size, size_t align)
{
  tensor_t* tensor = (tensor_t*) self;
  if (!mpp_comm_is_equal(MPP_COMM_CURR, tensor->comm))
    return convey_error_TEAM;
  if (item_size == 0)
    return convey_error_ZERO;
  size_t packet_bytes = porter_packet_size(tensor->item_offset, item_size);
  if (packet_bytes > tensor->max_bytes)
    return convey_error_OFLO;

  self->item_size = item_size;
  size_t header_bytes = tensor->item_offset;
  if (header_bytes == 0)
    header_bytes = sizeof(buffer_t);
  bool ok = true;

  if (tensor->accelerate)
    ok &= setup_methods(tensor);
  ok &= convey_prep_aligned(&tensor->aligned_item, item_size, header_bytes, align);
  for (int i = 0; i < tensor->order; i++)
    ok &= porter_setup(tensor->porters[i], item_size);
  tensor->n_complete = 0;
  tensor->stats[convey_BEGINS]++;
  tensor->packet_bytes = packet_bytes;
  return ok ? convey_OK : convey_error_ALLOC;
}

int
tensor_reset(convey_t* self)
{
  tensor_t* tensor = (tensor_t*) self;
  if (!mpp_comm_is_equal(MPP_COMM_CURR, tensor->comm))
    return convey_error_TEAM;
  for (int i = tensor->order - 1; i >= 0; i--)
    porter_reset(tensor->porters[i]);
  free(tensor->aligned_item);
  tensor->aligned_item = NULL;
  if (tensor->accelerate)
    reset_methods(tensor);
  convey_imp_update_stats(tensor->stats);
  return convey_OK;
}

int
tensor_set_codec(convey_t* self, const convey_codec_t* codec, void* arg)
{
  tensor_t* tensor = (tensor_t*) self;
  bool ok = porter_set_codec(remote_porter(tensor), codec, arg);
  return ok ? convey_OK : convey_FAIL;
}

int
tensor_free(convey_t* self)
{
  tensor_t* tensor = (tensor_t*) self;
  if (!mpp_comm_is_equal(MPP_COMM_CURR, tensor->comm))
    return convey_error_TEAM;
  for (int i = tensor->order - 1; i >= 0; i--)
    porter_free(tensor->porters[i]);
  free(tensor->aligned_item);
  free(self);
  return convey_OK;
}

int64_t
tensor_statistic(convey_t* self, int which)
{
  tensor_t* tensor = (tensor_t*) self;
  return convey_imp_statistic(tensor->stats, which);
}


/*** METHOD TABLES ***/

static const convey_methods_t tensor_fast_methods = {
  .push = &tensor_push,
  .pull = &tensor_pull,
  .unpull = &tensor_unpull,
  .advance = &tensor_advance,
  .begin = &tensor_begin,
  .reset = &tensor_reset,
  .free = &tensor_free,
  .statistic = &tensor_statistic,
  .apull = &tensor_apull,
  .epush = &convey_no_epush,
  .epull = &convey_no_epull,
  .set_codec = &tensor_set_codec,
  .panic = &convey_imp_panic,
};

static const convey_methods_t tensor_debug_methods = {
  ._super_ = &tensor_fast_methods,
  .push = &convey_checked_push,
  .pull = &convey_checked_pull,
  .unpull = &convey_checked_unpull,
  .advance = &tensor_advance,
  .begin = &tensor_begin,
  .reset = &tensor_reset,
  .free = &tensor_free,
  .statistic = &tensor_statistic,
  .apull = &convey_checked_apull,
  .epush = &convey_no_epush,
  .epull = &convey_no_epull,
  .set_codec = &tensor_set_codec,
  .panic = &convey_imp_panic,
};


/*** PRIVATE CONSTRUCTORS ***/

// options are reduced to the ones that porters understand
typedef tensor_t* (make_f)(convey_t*, size_t, size_t, size_t, size_t,
                           const convey_alc8r_t*, uint64_t, bool);

static tensor_t*
vector_new(convey_t* base, size_t capacity, size_t n_procs,
           size_t n_local, size_t n_buffers, const convey_alc8r_t* alloc,
           uint64_t options, bool quiet)
{
  const bool shrink = !(options & convey_opt_STANDARD);
  tensor_t* vector = malloc(sizeof(tensor_t));
  int32_t* friends = malloc(n_procs * sizeof(uint32_t));
  size_t my_proc = MY_PROC;
  bool ok = vector && friends;

  if (ok) {
    const size_t t = (shrink ? 0 : 4);  // omit routing tags if permitted
    *vector = (tensor_t) {
      .convey = *base, .order = 1, .tag_bytes = { t },
      .router = &vector_route,
    };
    for (int i = 0; i < n_procs; i++)
      friends[i] = i;
    if (n_procs <= n_local)
      options |= porter_opt_LOCAL;
    vector->porters[0] = porter_new(n_procs, friends, my_proc, t, capacity, n_buffers,
                                    alloc, options, CONVEY_SEND_1);
  }

  if (!ok || !vector->porters[0]) {
    free(friends);
    free(vector);
    CONVEY_REJECT(quiet, ok ? "construction failed" : "a small malloc() failed");
  }
  return vector;
}

static tensor_t*
matrix_new(convey_t* base, size_t capacity, size_t n_procs,
           size_t n_local, size_t n_buffers, const convey_alc8r_t* alloc,
           uint64_t options, bool quiet)
{
  const size_t n_rows = n_procs / n_local;
  const bool shrink = ((MATRIX_REMOTE_HOP ? n_rows : n_local) <= 256) &&
    !(options & (convey_opt_COMPRESS | convey_opt_STANDARD));

  tensor_t* matrix = malloc(sizeof(tensor_t));
  int32_t* friends[2];
  friends[0] = malloc(n_local * sizeof(uint32_t));
  friends[1] = malloc(n_rows * sizeof(uint32_t));
  if (! (matrix && friends[0] && friends[1])) {
    free(friends[1]);
    free(friends[0]);
    free(matrix);
    CONVEY_REJECT(quiet, "a small malloc() failed");
  }

  size_t my_proc = MY_PROC;
  uint32_t me[2] = { my_proc % n_local, my_proc / n_local };
  const size_t t = (shrink ? 1 : 4);
  *matrix = (tensor_t) {
    .convey = *base, .order = 2,
    .n_local = n_local, .router = &matrix_route,
  };
  matrix->tag_bytes[MATRIX_REMOTE_HOP] = t;
  matrix->tag_bytes[MATRIX_REMOTE_HOP ^ 1] = 4;
  matrix->div_local = _divbymul32_prep(n_local);
  matrix->pivots[0] = tensor_select_pivot_mid(t, 0);

  for (int i = 0; i < n_local; i++)
    friends[0][i] = my_proc + (i - me[0]);
  for (int j = 0; j < n_rows; j++)
    friends[1][j] = my_proc + (j - me[1]) * n_local;
  // must build porter[0] first because it gets freed last!
  for (int i = 0; i < 2; i++) {
    matrix->porters[i] = (i == MATRIX_REMOTE_HOP)
      ? porter_new(n_rows, friends[1], me[1], t, capacity, n_buffers, alloc,
               options, CONVEY_SEND_1)
      : porter_new(n_local, friends[0], me[0], 4, capacity, n_buffers, alloc,
               options | porter_opt_LOCAL, CONVEY_SEND_0);
  }    

  if (! (matrix->porters[0] && matrix->porters[1])) {
    for (int k = 0; k < 2; k++)
      if (!matrix->porters[k])
        free(friends[k]);
    tensor_free(&matrix->convey);
    CONVEY_REJECT(quiet, "construction failed");
  }

  return matrix;
}

static tensor_t*
tensor_new(convey_t* base, size_t capacity, size_t n_procs,
           size_t n_local, size_t n_buffers, const convey_alc8r_t* alloc,
           uint64_t options, bool quiet)
{
  const size_t squared = n_local * n_local;
  const size_t n_middle = 1 + (n_procs - 1) / squared;
  if (n_local > 256)
    CONVEY_REJECT(quiet, "number of local PEs is too large");
  const bool shrink = !(options & (convey_opt_COMPRESS | convey_opt_STANDARD));

  tensor_t* tensor = malloc(sizeof(tensor_t));
  int32_t* friends[3];
  friends[0] = malloc(n_local * sizeof(uint32_t));
  friends[1] = malloc(n_middle * sizeof(uint32_t));
  friends[2] = malloc(n_local * sizeof(uint32_t));
  if (! (tensor && friends[0] && friends[1] && friends[2])) {
    for (int i = 0; i < 3; i++)
      free(friends[i]);
    free(tensor);
    CONVEY_REJECT(quiet, "a small malloc() failed");
  }

  // Compute my coordinates
  size_t my_proc = MY_PROC;
  uint32_t me[3];
  me[0] = my_proc % n_local;
  me[1] = (my_proc / n_local) % n_local;
  me[2] = my_proc / squared;

  const size_t t = (shrink ? 1 + (n_local > 16) : 4);
  *tensor = (tensor_t) {
    .convey = *base, .order = 3, .tag_bytes = { 4, t, 4 },
    .n_local = n_local, .squared = squared, .router = &tensor_route,
  };
  tensor->div_local = _divbymul32_prep(n_local);
  tensor->div_square = _divbymul32_prep(squared);
  tensor->pivots[0] = tensor_select_pivot_early(t, 0);
  tensor->pivots[1] = tensor_select_pivot_late(t, 0);

  for (int i = 0; i < n_local; i++)
    friends[0][i] = friends[2][i] = my_proc - me[0] + i;
  for (int j = 0; j < n_middle; j++)
    friends[1][j] = j * squared + n_local * me[0] + me[1];
  tensor->porters[0] = porter_new(n_local, friends[0], me[0], 4, capacity, n_buffers,
                                  alloc, options | porter_opt_LOCAL, CONVEY_SEND_0);
  tensor->porters[1] = porter_new(n_middle, friends[1], me[2], t, capacity, n_buffers,
                                  alloc, options, CONVEY_SEND_1);
  tensor->porters[2] = porter_new(n_local, friends[2], me[0], 4, capacity, n_buffers,
                                  alloc, options | porter_opt_LOCAL, CONVEY_SEND_2);

  if (! (tensor->porters[0] && tensor->porters[1] && tensor->porters[2])) {
    for (int k = 0; k < 3; k++)
      if (!tensor->porters[k])
        free(friends[k]);
    tensor_free(&tensor->convey);
    CONVEY_REJECT(quiet, "construction failed");
  }

  return tensor;
}


/*** PUBLIC CONSTRUCTOR ***/

convey_t*
convey_new_tensor(size_t capacity, int order, size_t n_local, size_t n_buffers,
                  const convey_alc8r_t* alloc, uint64_t options)
{
  const bool quiet = (options & convey_opt_QUIET);
  const bool reckless = (options & convey_opt_RECKLESS);
  const bool steady = (options & convey_opt_PROGRESS);
  const bool compress = (options & convey_opt_COMPRESS);
  const bool standard = (options & convey_opt_STANDARD);
  const uint64_t porter_options = options &
    (convey_opt_DYNAMIC | convey_opt_PROGRESS | convey_opt_ALERT |
     convey_opt_COMPRESS | convey_opt_STANDARD | convey_opt_BLOCKING);

  if (capacity == 0 || order < 1 || order > 3 ||
      n_buffers == 0 || (n_buffers & n_buffers-1))
    CONVEY_REJECT(quiet, "invalid arguments");
  if (alloc && (!alloc->grab || !alloc->free))
    CONVEY_REJECT(quiet, "alloc is missing one or both methods");
  if (options &~ (convey_opt_RECKLESS | convey_opt_DYNAMIC | convey_opt_QUIET |
                  convey_opt_PROGRESS | convey_opt_ALERT | convey_opt_COMPRESS |
                  convey_opt_STANDARD | convey_opt_BLOCKING))
    CONVEY_REJECT(quiet, "unrecognized option(s)");

  size_t n_procs = PROCS;
  if (order > 1) {
    if (n_local == 0 || n_procs % n_local)
      CONVEY_REJECT(quiet, "n_local must divide PROCS if order > 1");
    if (order == 3 && n_procs <= n_local * n_local)
      order = 2;
    if (order == 2 && n_procs <= n_local)
      order = 1;
    if (n_local == 1)
      order = 1;
  }

  convey_t _base = {
    ._class_ = reckless ? &tensor_fast_methods : &tensor_debug_methods,
    .features = steady * convey_STEADY + compress * convey_THRIFTY,
    .suppress = UINT64_C(0) - quiet,
    .n_procs = n_procs, .state = convey_DORMANT,
  };
  make_f* maker = ((make_f* []) { &vector_new, &matrix_new, &tensor_new }) [order - 1];
  tensor_t* tensor = (*maker)(&_base, capacity, n_procs, n_local, n_buffers,
                              alloc, porter_options, quiet);
  if (tensor == NULL)
    return NULL;

  tensor->accelerate = reckless && !standard;
  tensor->item_offset = tensor->tag_bytes[order - 1];
  tensor->max_bytes = capacity;
  tensor->comm = MPP_COMM_CURR;
  if (compress)
    convey_set_codec(&tensor->convey, &convey_standard_codec, NULL);

  return &tensor->convey;
}

#endif
