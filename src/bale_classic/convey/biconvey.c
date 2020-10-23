// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include "biconvey_impl.h"
#include "private.h"
#define PANIC(ERR) biconvey_panic(self, __func__, ERR)

static int
biconvey_panic(biconvey_t* self, const char* where, int error)
{
  bool in_range = (-error > 0 && -error < 64);
  if (in_range && (self->suppress >> -error & 1))
    return error;
  mprint(MY_PROC, 0, "%s: %s\n", where, convey_error_string(NULL, error));
  if (in_range)
    self->suppress |= UINT64_C(1) << -error;
  return error;
}

int
biconvey_push(biconvey_t* self, const void* query, int64_t to)
{
  return self->_class_->push(self, query, to);
}

int
biconvey_pull(biconvey_t* self, void* reply)
{
  return self->_class_->pull(self, reply);
}

int
biconvey_advance(biconvey_t* self, bool done)
{
  if (self == NULL)
    return convey_error_NULL;
  // FIXME: do more error checking...
  int result = self->_class_->advance(self, done);
  return (result < 0) ? PANIC(result) : result;
}

int
biconvey_begin(biconvey_t* self, size_t query_bytes, size_t reply_bytes,
	       void (*answer)(const void* query, void* reply, void* context),
	       void* context)
{
  if (self == NULL)
    return convey_error_NULL;
  if (answer == NULL)
    return PANIC(convey_error_NOFUNC);
  if (query_bytes == 0 || reply_bytes == 0)
    return PANIC(convey_error_ZERO);

  self->query_bytes = query_bytes;
  self->reply_bytes = reply_bytes;
  self->answer = answer;
  self->context = context;
  int result = self->_class_->begin(self, query_bytes, reply_bytes);
  return (result < 0) ? PANIC(result) : result;
}

int
biconvey_reset(biconvey_t* self)
{
  if (self == NULL)
    return convey_error_NULL;
  int err = self->_class_->reset(self);
  return (err < 0) ? PANIC(err) : err;
}

int
biconvey_free(biconvey_t* self)
{
  if (self == NULL)
    return convey_OK;
  self->answer = NULL;
  int err = self->_class_->free(self);
  return (err < 0) ? PANIC(err) : err;
}


/*** Generic Constructors ***/

biconvey_t*
biconvey_new(size_t max_bytes, size_t n_local,
	     const convey_alc8r_t* alloc, uint64_t options)
{
  if (max_bytes < SIZE_MAX)
    max_bytes /= 3;
  if (n_local == 0)
    n_local = convey_procs_per_node();

  size_t capacity, n_buffers;
  int sync, order;
  convey_parameters(max_bytes, n_local, &capacity, &n_buffers, &sync, &order);

  if (sync)
    return biconvey_new_simple(capacity, alloc, NULL, options);
  else
    return biconvey_new_tensor(capacity, order, n_local, n_buffers, alloc, options);
}
