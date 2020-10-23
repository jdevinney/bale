// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: Accelerated push, pull, and pivot functions for tensor conveyors
// that exploit knowledge of tag size, item size, and porter internals.


#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "convey_impl.h"
#include "porter_impl.h"
#include "private.h"
#include "tensor.h"
#define ROUTER_HINT inline
#include "router.h"


#defcases Item = 4, 8, 12, 16, 20, 24, 28, 32

#defcases (Tag,Type) = (0,void), (1,uint8_t), (2,uint16_t), (4,uint32_t)
#calc Packet = (Tag == 0) ? Item : Tag * (1 + (Item + Tag - 1) / Tag)

static inline bool
porter_push_##Tag##_##Item(porter_t* self, uint64_t tag, const void* item, int dest)
{
  area_t* area = &self->send_areas[dest];
  bool room = (area->next < area->limit);
  if (room) {
    _prefetch_x(area->next + 96);
#if Tag > 0
    *(Type*)(area->next) = (Type) tag;
#endif
    memcpy(area->next + Tag, item, Item);
    area->next += Packet;
    if (area->next >= area->limit) {
      porter_close_buffer(self, dest, area);
      porter_try_send(self, dest);
    }
  }
  else
    porter_try_send(self, dest);
  return room;
}

#endcalc
#endcases


// Conveyor push methods

#defcases Tag = 0, 4
static int
vector_push_##Tag##_##Item(convey_t* self, const void* item, int64_t pe)
{
  tensor_t* vector = (tensor_t*) self;
  route_t _route = vector_route(vector, pe);
  bool ok = porter_push_##Tag##_##Item
    (vector->porters[0], _route.tag, item, _route.next);
  vector->stats[convey_PUSHES] += ok;
  return ok;
}
#endcases

#defcases Tag = 1, 4
static int
matrix_push_##Tag##_##Item(convey_t* self, const void* item, int64_t pe)
{
  tensor_t* matrix = (tensor_t*) self;
  route_t _route = matrix_route(matrix, pe);
  bool ok = porter_push_##Tag##_##Item
    (matrix->porters[0], _route.tag, item, _route.next);
  matrix->stats[convey_PUSHES] += ok;
  return ok;
}
#endcases

static int
tensor_push_4_##Item(convey_t* self, const void* item, int64_t pe)
{
  tensor_t* tensor = (tensor_t*) self;
  route_t _route = tensor_route(tensor, pe);
  bool ok = porter_push_4_##Item
    (tensor->porters[0], _route.tag, item, _route.next);
  tensor->stats[convey_PUSHES] += ok;
  return ok;
}


// Conveyor pull methods

#defcases (Order,Flavor,Tag) = (1,vector,0), (2,matrix,1), (2,matrix,4), (3,tensor,4)
#calc Packet = (Tag == 0) ? Item : Tag * (1 + (Item + Tag - 1) / Tag)

static int
Flavor##_pull_##Tag##_##Item(convey_t* self, void* item, int64_t* from)
{
  tensor_t* tensor = (tensor_t*) self;
  buffer_t* buffer = tensor->buffer;

  if (buffer && buffer->start == buffer->limit) {
    porter_return(tensor->porters[Order - 1]);
    buffer = NULL;
  }
  if (!buffer) {
    buffer = porter_borrow(tensor->porters[Order - 1]);
    tensor->buffer = buffer;
    if (!buffer)
      return convey_FAIL;
  }

  char* packet = &buffer->data[buffer->start];
  buffer->start += Packet;
  if (from) {
    uint32_t source = buffer->source;
#if Order == 1
    *from = source;
#else
    uint32_t tag = (Tag == 1) ? *(uint8_t*)packet : *(uint32_t*)packet;
    *from = origin_from_tag(tensor, Order, tag, source);
#endif
  }

  tensor->stats[convey_PULLS]++;
  memcpy(item, packet + Tag, Item);
  return convey_OK;
}

#endcalc
#endcases


// Specialized pivot functions
#encase "pivot.h"

#endcases


// Standard pivot functions
#defcases Item = 0
#encase "pivot.h"
#endcases


/*** Selector Functions ***/

static push_f* const push_functions[][5] = {
#defcases Item = 4, 8, 12, 16, 20, 24, 28, 32
  { &vector_push_0_##Item, &vector_push_4_##Item,
      &matrix_push_1_##Item, &matrix_push_4_##Item, 
      &tensor_push_4_##Item },
#endcases
};

static pull_f* const pull_functions[][4] = {
#defcases Item = 4, 8, 12, 16, 20, 24, 28, 32
  { &vector_pull_0_##Item, &matrix_pull_4_##Item, &tensor_pull_4_##Item,
    &matrix_pull_1_##Item },
#endcases
};

static pivot_f* const pivot_functions[][8] = {
#defcases Item = 4, 8, 12, 16, 20, 24, 28, 32
  { &pivot_mid_1_##Item, &pivot_mid_4_##Item,
    &pivot_early_1_##Item, &pivot_early_2_##Item, &pivot_early_4_##Item,
    &pivot_late_1_##Item, &pivot_late_2_##Item, &pivot_late_4_##Item },
#endcases
};

static int
item_index(size_t item_bytes)
{
  if (item_bytes > 32 || (item_bytes & 3))
    return -1;
  return (item_bytes / 4) - 1;
}

push_f*
tensor_select_push(int order, size_t tag_bytes, size_t item_bytes)
{
  int index = item_index(item_bytes);
  if (index < 0)
    return NULL;
  if (order == 1 && tag_bytes == 0)
    return push_functions[index][0];
  if (order == 1 && tag_bytes == 4)
    return push_functions[index][1];
  if (order == 2 && tag_bytes == 1)
    return push_functions[index][2];
  if (order == 2 && tag_bytes == 4)
    return push_functions[index][3];
  if (order == 3 && tag_bytes == 4)
    return push_functions[index][4];
  return NULL;
}

pull_f*
tensor_select_pull(int order, size_t tag_bytes, size_t item_bytes)
{
  int index = item_index(item_bytes);
  if (index < 0)
    return NULL;
#if MATRIX_REMOTE_HOP == 1
  if (order == 2 && tag_bytes == 1)
    return pull_functions[index][3];
#endif
  if (tag_bytes != ((order == 1) ? 0 : 4))
    return NULL;
  return pull_functions[index][order - 1];
}

pivot_f*
tensor_select_pivot_mid(size_t tag_bytes, size_t item_bytes)
{
  int index = item_index(item_bytes);
  if (tag_bytes != 1 && tag_bytes != 4)
    return NULL;
  if (index < 0)
    return (tag_bytes == 4) ? &pivot_mid_4_0 : &pivot_mid_1_0;
  return pivot_functions[index][tag_bytes == 4];
}

pivot_f*
tensor_select_pivot_early(size_t tag_bytes, size_t item_bytes)
{
  int index = item_index(item_bytes);
  int logtag = _trailz(tag_bytes);
  if (_popcnt(tag_bytes) != 1 || logtag > 2)
    return NULL;
  if (index < 0)
    return ((pivot_f* [])
      { &pivot_early_1_0, &pivot_early_2_0, &pivot_early_4_0 }) [logtag];
  return pivot_functions[index][2 + logtag];
}

pivot_f*
tensor_select_pivot_late(size_t tag_bytes, size_t item_bytes)
{
  int index = item_index(item_bytes);
  int logtag = _trailz(tag_bytes);
  if (_popcnt(tag_bytes) != 1 || logtag > 2)
    return NULL;
  if (index < 0)
    return ((pivot_f* [])
      { &pivot_late_1_0, &pivot_late_2_0, &pivot_late_4_0 }) [logtag];
  return pivot_functions[index][5 + logtag];
}
