// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <stdlib.h>
#include <string.h>

#include "bolite.h"
#include "sorter.h"



typedef bool (sorter_push_f)(sorter_t*, const void*, int);

struct sorter {
  sorter_push_f* push;
  size_t item_bytes;
  area_t* areas;
  int degree;    // number of areas this node sorts into; positive
  // The following fields are relevant only for internal nodes:
  int shift;     // the top bits of a tag determine the child
  sorter_t** children;
  size_t packet_bytes;
  size_t bucket_bytes;
  char* buckets;     // points to the array of buckets
  size_t n_buckets;  // # used by this node and its descendants
  // The following fields are only meaningful at the root:
  bool dynamic;
  const convey_alc8r_t* alloc;
};


/*** Internal Methods for Sorting ***/

static bool
sorter_push_leaf(sorter_t* self, const void* item, int dest)
{
  area_t* area = &self->areas[dest];
  _prefetch_x(area->next + 64);
  memcpy(area->next, item, self->item_bytes);
  area->next += self->item_bytes;
  return (area->next < area->limit);
}

// Distribute items into sub-buckets.  Return true if no buffer filled up.
static bool
scatter(sorter_t* self, int index)
{
  area_t* area = &self->areas[index];
  char* const limit = area->next;
  char* const start = self->buckets + index * self->bucket_bytes;
  if (start == limit)
    return true;

  sorter_t* child = self->children[index];
  sorter_push_f* push = child->push;
  const size_t step = self->packet_bytes;

  int n = child->degree;
  for (int i = 0; i < n; i++)
    _prefetch_x(child->areas[i].next);
  for (char* p = start; p < limit; ) {
    uint32_t* packet = (uint32_t*) p;
    bool ok = (*push)(child, packet + 1, packet[0]);
    p += step;
    if (!ok) {
      memmove(start, p, limit - p);
      area->next = start + (limit - p);
      return false;
    }
  }

  area->next = start;
  return true;
}

static bool
sorter_push_node(sorter_t* self, const void* item, int dest)
{
  int shift = self->shift;
  int index = dest >> shift;
  area_t* area = &self->areas[index];

  uint32_t* packet = (uint32_t*) area->next;
  _prefetch_x(packet + 16);
  packet[0] = dest - (index << shift);
  memcpy(packet + 1, item, self->item_bytes);
  area->next += self->packet_bytes;
  if (area->next < area->limit)
    return true;

  return scatter(self, index);
}


/*** Construction and Destruction ***/

static void
sorter_setup_areas(sorter_t* self)
{
  int n = self->degree;
  if (self->children) {
    for (int i = 0; i < n; i++)
      sorter_setup_areas(self->children[i]);
    for (int i = 0; i < n; i++)
      self->areas[i].next = self->buckets + i * self->bucket_bytes;
  }
}

static char*
sorter_build_buckets(sorter_t* self, char* memory)
{
  if (self->children) {
    int n = self->degree;
    self->buckets = memory;
    for (int i = 0; i < n; i++) {
      memory += self->bucket_bytes;
      self->areas[i].limit = memory;
    }
    for (int i = 0; i < n; i++)
      memory = sorter_build_buckets(self->children[i], memory);
  }
  return memory;
}

static bool
sorter_setup_buckets(sorter_t* self)
{
  if (!self->children)
    return true;

  size_t memory_bytes = self->n_buckets * self->bucket_bytes;
  char* memory = self->alloc->grab(self->alloc->alc8r, memory_bytes, __FILE__, __LINE__);
  if (memory)
    sorter_build_buckets(self, memory);
  return (memory != NULL);
}

static void
sorter_free_buckets(sorter_t* self)
{
  if (self->buckets) {
    self->alloc->free(self->alloc->alc8r, self->buckets);
    self->buckets = NULL;
  }
}

static sorter_t*
sorter_tree(int n, area_t areas[n], size_t item_bytes, size_t capacity)
{
  sorter_t* node = malloc(sizeof(sorter_t));
  if (node == NULL)
    return NULL;

  if (n <= 128)
    *node = (sorter_t) {
      .push = &sorter_push_leaf, .degree = n,
      .item_bytes = item_bytes, .areas = areas,
    };
  else {
    // Decide how many bits to bite off
    int span = 64 - _leadz(n - 1);   // at least 8
    int levels = (span + 6) / 7;     // is this a good idea?
    int bits = span / levels;
    int shift = span - bits;
    int degree = 1 + (n-1 >> shift);

    size_t packet_bytes = 4 * (1 + (item_bytes + 3 >> 2));
    size_t bucket_bytes = packet_bytes * capacity;

    *node = (sorter_t) {
      .push = &sorter_push_node, .degree = degree, .shift = shift,
      .item_bytes = item_bytes, .packet_bytes = packet_bytes, .bucket_bytes = bucket_bytes,
    };
    node->areas = malloc(degree * sizeof(area_t));
    node->children = malloc(degree * sizeof(sorter_t*));
    if (node->children)
      for (int i = 0; i < degree; i++)
        node->children[i] = NULL;

    size_t n_buckets = degree;
    bool ok = node->areas && node->children;
    for (int i = 0; ok && i < degree; i++) {
      int m = MIN(n, i+1 << shift) - (i << shift);
      sorter_t* child = sorter_tree(m, &areas[i << shift], item_bytes, capacity);
      if (child) {
        node->children[i] = child;
        n_buckets += child->n_buckets;
      }
      else
        ok = false;
    }
    node->n_buckets = n_buckets;

    if (!ok) {
      sorter_free(node);
      node = NULL;
    }
  }

  return node;
}


/*** External Functions ***/

bool
sorter_push(sorter_t* self, const void* item, int dest)
{
  return self->push(self, item, dest);
}

bool
sorter_flush(sorter_t* self)
{
  bool done = true;
  if (self->children) {
    int n = self->degree;
    for (int i = 0; i < n; i++)
      done &= scatter(self, i) && sorter_flush(self->children[i]);
  }
  return done;
}

bool
sorter_setup(sorter_t* self)
{
  if (self->dynamic && !sorter_setup_buckets(self))
    return false;
  sorter_setup_areas(self);
  return true;
}

void
sorter_reset(sorter_t* self)
{
  if (self->dynamic)
    sorter_free_buckets(self);
}

void
sorter_free(sorter_t* self)
{
  if (self == NULL)
    return;
  if (self->alloc)
    sorter_free_buckets(self);
  if (self->children) {
    for (int i = self->degree - 1; i >= 0; i--)
      sorter_free(self->children[i]);
    free(self->children);
    free(self->areas);
  }
  free(self);
}

sorter_t*
sorter_new(int n, area_t areas[n], size_t item_bytes, size_t capacity,
           const convey_alc8r_t* alloc, bool dynamic)
{
  if (n <= 0 || areas == NULL || item_bytes == 0 || capacity == 0 || alloc == NULL)
    return NULL;

  sorter_t* self = sorter_tree(n, areas, item_bytes, capacity);
  if (self) {
    self->alloc = alloc;
    self->dynamic = dynamic;
    if (!dynamic && !sorter_setup_buckets(self)) {
      sorter_free(self);
      self = NULL;
    }
  }

  return self;
}
