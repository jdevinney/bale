// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include "convey_impl.h"
#if HAVE_CONFIG_H
#include "config.h"
#endif
#define PANIC(ERR) convey_panic(self, __func__, ERR)

const char convey_version[] = PACKAGE_VERSION;


/*** Fast dispatch functions that don't necessarily check for errors ***/

int
convey_push(convey_t* self, const void* item, int64_t pe)
{
  return self->_class_->push(self, item, pe);
}

int
convey_pull(convey_t* self, void* item, int64_t* from)
{
  return self->_class_->pull(self, item, from);
}

int
convey_unpull(convey_t* self)
{
  if (self->state == convey_COMPLETE)
    return convey_FAIL;
  return self->_class_->unpull(self);
}

void*
convey_apull(convey_t* self, int64_t* from)
{
  return self->_class_->apull(self, from);
}

int
convey_epush(convey_t* self, size_t bytes, const void* item, int64_t pe)
{
  return self->_class_->epush(self, bytes, item, pe);
}

int
convey_epull(convey_t* self, convey_item_t* result)
{
  return self->_class_->epull(self, result);
}


/*** Checking methods ***/

int
convey_checked_push(convey_t* self, const void* item, int64_t pe)
{
  if (self->state != convey_WORKING)
    return PANIC(convey_error_STATE);
  if (pe < 0 || pe >= self->n_procs)
    return PANIC(convey_error_RANGE);
  return self->_class_->_super_->push(self, item, pe);
}

int
convey_checked_pull(convey_t* self, void* item, int64_t* from)
{
  if (self->state == convey_DORMANT)
    return PANIC(convey_error_STATE);
  return self->_class_->_super_->pull(self, item, from);
}

int
convey_checked_unpull(convey_t* self)
{
  if (self->state == convey_DORMANT)
    return PANIC(convey_error_STATE);
  return self->_class_->_super_->unpull(self);
}

void*
convey_checked_apull(convey_t* self, int64_t* from)
{
  if (self->state == convey_DORMANT)
    return NULL;
  return self->_class_->_super_->apull(self, from);
}

int
convey_checked_epush(convey_t* self, size_t bytes, const void* item, int64_t pe)
{
  if (self->state != convey_WORKING)
    return PANIC(convey_error_STATE);
  if (pe < 0 || pe >= self->n_procs)
    return PANIC(convey_error_RANGE);
  return self->_class_->_super_->epush(self, bytes, item, pe);
}

int
convey_checked_epull(convey_t* self, convey_item_t* result)
{
  if (self->state == convey_DORMANT) {
    result->data = NULL;
    return convey_FAIL;
  }
  return self->_class_->_super_->epull(self, result);
}


/*** Default null methods ***/

int
convey_no_epush(convey_t* self, size_t bytes, const void* item, int64_t pe)
{
  return PANIC(convey_error_RIGID);
}

int
convey_no_epull(convey_t* self, convey_item_t* result)
{
  return PANIC(convey_error_RIGID);
}


/*** Slow dispatch functions that check stuff and update the state ***/

int
convey_begin(convey_t* self, size_t item_bytes, size_t align)
{
  if (self == NULL)
    return convey_error_NULL;
  if (self->state != convey_DORMANT)
    return PANIC(convey_error_STATE);
  if (item_bytes == 0)
    return PANIC(convey_error_ZERO);

  if (align == 0 || (align & (align - 1)) || (item_bytes & (align - 1)))
    align = 1 + ((item_bytes - 1) & ~item_bytes);  // 1 << _trailz(item_bytes)
  if (align > CONVEY_MAX_ALIGN)
    align = CONVEY_MAX_ALIGN;

  int result = self->_class_->begin(self, item_bytes, align);
  self->state = convey_WORKING;
  return (result < 0) ? PANIC(result) : result;
}

int
convey_advance(convey_t* self, bool done)
{
  if (self == NULL)
    return convey_error_NULL;
  if (self->state == convey_DORMANT)
    return PANIC(convey_error_STATE);
  // Rejecting a wrong 'done' value might cause deadlock
  done |= (self->state != convey_WORKING);

  int result = self->_class_->advance(self, done);

  if (result == convey_NEAR)
    self->state = convey_CLEANUP;
  else if (result == convey_DONE)
    self->state = convey_COMPLETE;
  else if (done && self->state == convey_WORKING)
    self->state = convey_ENDGAME;
  return (result < 0) ? PANIC(result) : result;
}

int
convey_reset(convey_t* self)
{
  if (self == NULL)
    return convey_error_NULL;
  if (self->state == convey_DORMANT)
    return convey_OK;
  if (self->state != convey_COMPLETE)
    return PANIC(convey_error_STATE);
  int result = self->_class_->reset(self);
  self->state = convey_DORMANT;
  return (result < 0) ? PANIC(result) : result;
}

int
convey_free(convey_t* self)
{
  if (self == NULL)
    return convey_OK;
  if (self->state != convey_COMPLETE && self->state != convey_DORMANT)
    return PANIC(convey_error_STATE);
  int result = convey_OK;
  if (self->features & convey_THRIFTY)
    result = self->_class_->set_codec(self, NULL, NULL);
  if (result >= 0)
    result = (self->_class_->free)(self);
  return (result < 0) ? PANIC(result) : result;
}

int
convey_set_codec(convey_t* self, const convey_codec_t* methods, void* arg)
{
  if (self == NULL)
    return convey_error_NULL;
  if (self->state != convey_DORMANT)
    return convey_error_STATE;
  if (!(self->features & convey_THRIFTY))
    return convey_FAIL;
  if (methods != NULL && methods->plan == NULL)
    return convey_FAIL;
  int result = self->_class_->set_codec(self, methods, arg);
  return (result < 0) ? PANIC(result) : result;
}

size_t
convey_item_size(convey_t* self)
{
  if (self && self->state != convey_DORMANT)
    return self->item_size;
  return (size_t) 0;
}

uint64_t
convey_features(convey_t* self)
{
  return self ? self->features : UINT64_C(0);
}

int64_t
convey_statistic(convey_t* self, int which)
{
  return self ? self->_class_->statistic(self, which) : INT64_C(0);
}


/*** Error Messages ***/

static const char * const error_strings[] = {
  [-convey_error_STATE] = "operation is illegal in current state",
  [-convey_error_USAGE] = "misuse of conveyor API detected",
  [-convey_error_NULL] = "conveyor is NULL",
  [-convey_error_RANGE] = "PE number is out of range",
  [-convey_error_TEAM] = "current communicator is wrong",
  [-convey_error_COMMS] = "internal communication error",
  [-convey_error_ALLOC] = "memory allocation failed",
  [-convey_error_RIGID] = "conveyor is not elastic",
  [-convey_error_MISFIT] = "item has wrong size for pull",
  [-convey_error_OFLO] = "item is too large for conveyor",
  [-convey_error_ZERO] = "item size is zero",
  [-convey_error_TOOBIG] = "number of buffered items is too large",
  [-convey_error_NOFUNC] = "function pointer is NULL",
};

const char*
convey_error_string(convey_t* self, int error)
{
  if (error == 0)
    return "ordinary failure";
  else if (error > 0)
    return "value indicates success";
  const char* string = NULL;
  if (-error < sizeof(error_strings) / sizeof(const char*))
    string = error_strings[-error];
  else if (self != NULL && self->_class_->error_string)
    string = self->_class_->error_string(self, error);
  return string ? string : "unknown error";
}

int
convey_panic(convey_t* self, const char* where, int error)
{
  if (-error > 0 && -error < 64 && (self->suppress >> -error & 1))
    return error;
  if (self->_class_->panic)
    return self->_class_->panic(self, where, error);
  return error;
}


/*** Generic Constructors ***/

#include <stdlib.h>
#include "private.h"

static size_t
porter_usage(size_t capacity, size_t n_procs, size_t n_buffers)
{
  return n_procs * (2 * n_buffers * (2 * sizeof(size_t) + capacity) +
                    sizeof(uint64_t) + sizeof(long long) * n_buffers);
}

// Compute the amount of "symmetric" memory needed by a conveyor
// with the given parameters.  (Can be a slight overestimate.)
size_t
convey_memory_usage(size_t capacity, bool sync, int order,
		    size_t n_procs, size_t n_local, size_t n_buffers)
{
  // Simple conveyors
  if (sync && order == 1)
    return n_procs * 2 * (capacity + sizeof(size_t));

  // Twohop conveyors
  if (sync) {
    size_t n_rows = n_procs / n_local;
    size_t big_cap = (n_local > n_rows) ?
      (capacity * n_local) / n_rows : (capacity * n_rows) / n_local;
    return 2 * (capacity + big_cap) + 56 * (n_local + n_rows);
  }

  // Tensor conveyors
  if (order == 1)
    return porter_usage(capacity, n_procs, n_buffers);
  else if (order == 2)
    return porter_usage(capacity, n_local, n_buffers) +
      porter_usage(capacity, n_procs / n_local, n_buffers);
  else {
    size_t n_middle = 1 + (n_procs - 1) / (n_local * n_local);
    return 2 * porter_usage(capacity, n_local, n_buffers) +
      porter_usage(capacity, n_middle, n_buffers);
  }
}

void
convey_parameters(size_t max_bytes, size_t n_local,
                  size_t* capacity_, size_t* n_buffers_, int* sync_, int* order_)
{
  // Set up some default values
  bool sync = false;
  bool twohop = true;
  int order = 1;
  size_t n_buffers = 2;
#if !MPP_USE_MPI && !(MPP_USE_SHMEM && HAVE_SHMEMX_ALLTOALLV && HAVE_SHMEMX_TEAM_ALLTOALLV)
  twohop = false;
#endif

  // Find the preferred capacity
  size_t capacity = 0;
  char* number = getenv("CONVEY_BUFFER_SIZE");
  if (number) {
    char* endptr;
    capacity = strtoul(number, &endptr, 0);
    if (*endptr == 'k' || *endptr == 'K')
      capacity *= 1000;
  }
  if (capacity == 0)
    capacity = CONVEY_BUFFER_SIZE;

  // Try to choose the fastest conveyor, ignoring max_bytes
  const size_t n_procs = PROCS;
  if (n_local > 1 && n_local < n_procs) {
    if (!sync || twohop)
      order = 2;
    size_t middle = n_procs / (n_local * n_local);
    if (!sync && middle * middle >= n_local)
      order = 3;
  }

  if (max_bytes < SIZE_MAX)
    for (int step = 1; true; step++) {
      size_t usage = convey_memory_usage(capacity, sync, order, n_procs, n_local, n_buffers);
      if (usage <= max_bytes)
        break;

      // Step 1: Try to reduce the number of buffers per link
      if (step == 1)
        n_buffers = 1;

      // Step 2: Try to balance the hops
      else if (step == 2) {
        if (!sync || twohop) {
          size_t new_local = n_local;
          while (new_local % 2 == 0 && new_local * new_local > 4 * n_procs)
            new_local >>= 1;
          if (new_local != n_local && n_procs > new_local) {
            size_t try = convey_memory_usage(capacity, true, 2, n_procs, new_local, 1);
            if (try < usage) {
              order = 2;
              n_local = new_local;
              usage = try;
            }
          }
        }
        if (!sync) {
          size_t new_local = n_local;
          while (new_local % 2 == 0 && new_local * new_local * new_local > 8 * n_procs)
            new_local >>= 1;
          if (new_local != n_local && n_procs > new_local * new_local) {
            size_t try = convey_memory_usage(capacity, false, 3, n_procs, new_local, 1);
            if (try < usage) {
              order = 3;
              n_local = new_local;
            }
          }
        }
      }

      // Step 3: Set the buffer capacity to an appropriate value
      else if (step == 3) {
        capacity = (capacity * usage) / max_bytes;
        capacity = (capacity > 1) ? capacity - 1 : 1;
      }

      // Steps 4+: Try to reduce the buffer capacity further
      else {
        if (capacity == 1)
          break;
        capacity = (capacity + 1) >> 1;
      }
    }

  *capacity_ = capacity;
  *n_buffers_ = n_buffers;
  *sync_ = sync;
  *order_ = order;
}

static const uint64_t common_options = convey_opt_RECKLESS | convey_opt_DYNAMIC |
  convey_opt_QUIET | convey_opt_PROGRESS | convey_opt_ALERT;

convey_t*
convey_new(size_t max_bytes, size_t n_local,
           const convey_alc8r_t* alloc, uint64_t options)
{
  if (n_local == 0)
    n_local = convey_procs_per_node();

  size_t capacity, n_buffers;
  int sync, order;
  convey_parameters(max_bytes, n_local, &capacity, &n_buffers, &sync, &order);

  // Build the chosen conveyor
#if 0
  if (sync && order > 1)
    return convey_new_twohop(capacity, n_local, alloc,
                             options & common_options);
#endif
  if (sync)
    return convey_new_simple(capacity, alloc, NULL,
                             options & (common_options | convey_opt_SCATTER));
  else
    return convey_new_tensor(capacity, order, n_local, n_buffers, alloc,
                             options & (common_options | convey_opt_STANDARD |
                                        convey_opt_COMPRESS | convey_opt_BLOCKING));
}

convey_t*
convey_new_elastic(size_t item_bound, size_t max_bytes, size_t n_local,
                   const convey_alc8r_t* alloc, uint64_t options)
{
#if MPP_NO_MIMD
  return convey_new_trivial(item_bound, alloc, options);
#else
  // Deduct space for the monster buffers, even if they aren't needed
  if (max_bytes < SIZE_MAX) {
    size_t big = 2 * item_bound;
    max_bytes = (big > max_bytes) ? 0 : max_bytes - big;
  }

  size_t buffer_bytes, n_buffers;
  int sync, order;
  convey_parameters(max_bytes, n_local, &buffer_bytes, &n_buffers, &sync, &order);
  if (sync)
    return NULL;

  return convey_new_etensor(buffer_bytes, item_bound, order, n_local, n_buffers,
                            alloc, options & (common_options | convey_opt_BLOCKING));

#endif
}

#undef PANIC
