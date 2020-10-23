// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#ifndef CONVEY_TENSOR_H
#define CONVEY_TENSOR_H

#include <stddef.h>
#include <stdint.h>
#include <string.h>

typedef struct tensor tensor_t;
typedef bool (pivot_f)(tensor_t*, buffer_t*);
typedef struct { uint32_t tag; int64_t next; } route_t;
typedef route_t (route_f)(tensor_t*, int64_t);

struct tensor {
  convey_t convey;
  int n_complete;
  int order;  // 1=vector, 2=matrix, 3=tensor
  // Fixed data for dealing with routing tags
  uint8_t tag_bytes[3];
  bool accelerate;
  uint32_t n_local;
  uint32_t squared;
  divbymul32_t div_local;
  divbymul32_t div_square;
  // Other stuff
  size_t item_offset;   // tag_bytes[order - 1]
  size_t packet_bytes;  // size of pulled packets
  void* aligned_item;
  buffer_t* buffer;
  porter_t* porters[3];
  route_f* router;
  pivot_f* pivots[2];
  size_t max_bytes;
  mpp_comm_t comm;
  int64_t stats[convey_imp_N_STATS];
};

static inline int64_t
origin_from_tag(tensor_t* tensor, int order, uint32_t tag, uint32_t source)
{
  if (order == 1)
    return source;
  else if (order == 2)
#ifndef MATRIX_REMOTE_HOP
#error "MATRIX_REMOTE_HOP needs to be defined"
#elif MATRIX_REMOTE_HOP == 0
    return tag * tensor->n_local + source;
#else
    return source * tensor->n_local + tag;
#endif
  else
    return (tag >> 8) * tensor->squared
      + source * tensor->n_local + (tag & 0xFF);
}

static inline void*
pull_pointer(uint32_t* item, void* temp, size_t item_size)
{
  return (temp && item) ? memcpy(temp, item, item_size) : (void*)item;
}

/*** Functions that etensor conveyors need to call ***/

int tensor_advance(convey_t* self, bool done);
int tensor_begin(convey_t* self, size_t item_bytes, size_t align);
int tensor_reset(convey_t* self);
int tensor_free(convey_t* self);
int64_t tensor_statistic(convey_t* self, int which);

/*** Functions for retrieving codelets (accel.c) ***/

typedef int (push_f)(convey_t*, const void*, int64_t);
typedef int (pull_f)(convey_t*, void*, int64_t*);
// FIXME: accelerate tensor_apull() also, for completeness

push_f* tensor_select_push(int order, size_t tag_bytes, size_t item_bytes);
pull_f* tensor_select_pull(int order, size_t tag_bytes, size_t item_bytes);
// mid: either tag_bytes => 4 or (if MATRIX_REMOTE_HOP) 4 => tag bytes
pivot_f* tensor_select_pivot_mid(size_t tag_bytes, size_t item_bytes);
// early: 4 => tag_bytes; late: tag_bytes => 4
pivot_f* tensor_select_pivot_early(size_t tag_bytes, size_t item_bytes);
pivot_f* tensor_select_pivot_late(size_t tag_bytes, size_t item_bytes);

#endif
