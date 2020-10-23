// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#ifndef CONVEY_IMPL_H
#define CONVEY_IMPL_H

#include "convey.h"

// Negative error codes:
enum convey_error {
  convey_error_NOFUNC = -13,  // function pointer is NULL
  convey_error_TOOBIG,     // total buffer capacity is too large
  convey_error_ZERO,       // item size is zero
  convey_error_OFLO,       // item is too large for conveyor
  convey_error_MISFIT,     // next item is wrong size for pull
  convey_error_RIGID,      // conveyor is not elastic
  convey_error_ALLOC,      // a memory allocation failed
  convey_error_COMMS,      // an internal communication operation failed
  convey_error_TEAM,       // collective method called with wrong communicator
  convey_error_RANGE,      // PE number out of range
  convey_error_NULL,       // method called on a null conveyor
  convey_error_USAGE,      // caller has violated the contract
  convey_error_STATE,      // method called in the wrong state
};

enum convey_state {
  convey_DORMANT = 0,
  convey_WORKING = 1,
  convey_ENDGAME = 2,
  convey_CLEANUP = 3,
  convey_COMPLETE = 4,
};

typedef struct conveyor_methods convey_methods_t;

struct conveyor_methods {
  const convey_methods_t* _super_;
  int (*push)(convey_t* self, const void* item, int64_t pe);
  int (*pull)(convey_t* self, void* item, int64_t* from);
  int (*unpull)(convey_t* self);
  int (*advance)(convey_t* self, bool done);
  int (*begin)(convey_t* self, size_t item_bytes, size_t align);
  int (*reset)(convey_t* self);
  int (*free)(convey_t* self);
  int64_t (*statistic)(convey_t* self, int which);
  const char* (*error_string)(convey_t* self, int error);
  void* (*apull)(convey_t* self, int64_t* from);
  int (*epush)(convey_t* self, size_t bytes, const void* item, int64_t pe);
  int (*epull)(convey_t* self, convey_item_t* result);
  int (*set_codec)(convey_t* self, const convey_codec_t* methods, void* arg);
  // Private method
  int (*panic)(convey_t* self, const char* where, int error);
  // Should this structure be freed when the conveyor is freed?
  bool dynamic;
};

struct conveyor {
  const convey_methods_t* _class_;
  uint64_t features;
  size_t item_size;
  size_t n_procs;
  uint64_t suppress;  // bitmask: errors to pass through
  int64_t state;
};

int convey_checked_push(convey_t* self, const void* item, int64_t pe);
int convey_checked_pull(convey_t* self, void* item, int64_t* from);
int convey_checked_unpull(convey_t* self);
void* convey_checked_apull(convey_t* self, int64_t* from);
int convey_checked_epush(convey_t* c, size_t bytes, const void* item, int64_t pe);
int convey_checked_epull(convey_t* c, convey_item_t* result);

int convey_no_epush(convey_t* c, size_t bytes, const void* item, int64_t pe);
int convey_no_epull(convey_t* c, convey_item_t* result);

int convey_panic(convey_t* self, const char* where, int error);

// choose good parameters using heuristics and environment variables
void convey_parameters(size_t max_bytes, size_t n_local,
		       size_t* capacity_, size_t* n_buffers_, int* sync_, int* order_);
size_t convey_memory_usage(size_t capacity, bool sync, int order,
			   size_t n_procs, size_t n_local, size_t n_buffers);

#endif
