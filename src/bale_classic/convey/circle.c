// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: A simple sorter implementation that keeps a small circular
// buffer of recently pushed items in order to prefetch the buckets they
// belong in.


#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "sorter.h"
#include "bolite.h"


struct sorter {
  area_t* areas;
  char* circular;     // items in transit
  int* targets;       // destination PEs
  size_t item_bytes;
  uint64_t mask;
  // Mutable fields:
  uint64_t head;      // number of items pushed to us
  uint64_t tail;      // number of items we have written out
};

sorter_t*
sorter_new(int n, area_t areas[n], size_t item_bytes, size_t capacity,
           const convey_alc8r_t* alloc, bool dynamic)
{
  if (n == 0 || item_bytes == 0 || capacity == 0 || (capacity & (capacity - 1)))
    return NULL;

  sorter_t* self = malloc(sizeof(sorter_t));
  if (self == 0)
    return NULL;
  *self = (sorter_t) { .areas = areas, .item_bytes = item_bytes, .mask = capacity - 1 };
  self->circular = malloc(capacity * item_bytes);
  self->targets = malloc(capacity * sizeof(int));
  if (!self->circular || !self->targets) {
    sorter_free(self);
    self = NULL;
  }

  return self;
}

bool
sorter_setup(sorter_t* self)
{
  self->head = 0;
  self->tail = 0;
  return true;
}

bool
sorter_push(sorter_t* self, const void* item, int dest)
{
  _prefetch_x(self->areas[dest].next);

  const uint64_t mask = self->mask;
  const size_t size = self->item_bytes;
  uint64_t head = self->head;
  uint64_t tail = self->tail;
  uint64_t index = head & mask;
  char* position = self->circular + index * size;
  bool space = true;

  if (head > tail + mask) {
    int target = self->targets[index];
    area_t* area = &self->areas[target];
    memcpy(area->next, position, size);
    area->next += size;
    space = (area->next < area->limit);
    self->tail = tail + 1;
  }

  self->targets[index] = dest;
  memcpy(position, item, size);
  self->head = head + 1;
  return space;
}

int
sorter_flush(sorter_t* self)
{
  if (self->tail == self->head)
    return 0;

  const size_t size = self->item_bytes;
  while (self->tail < self->head) {
    uint64_t index = self->tail & self->mask;
    int target = self->targets[index];
    area_t* area = &self->areas[target];
    memcpy(area->next, self->circular + index * size, size);
    area->next += size;
    self->tail++;
    if (area->next >= area->limit)
      return -1;
  }

  return 1;
}

void
sorter_reset(sorter_t* self)
{
  // nothing to do
}

void
sorter_free(sorter_t* self)
{
  if (self == NULL)
    return;
  free(self->targets);
  free(self->circular);
  free(self);
}
