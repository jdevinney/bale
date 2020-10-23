// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: An elastic conveyor to handle the non-parallel case (MPP_NO_MIMD).

#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#include "convey_impl.h"
#include "private.h"


typedef struct trivial {
  convey_t convey;
  convey_item_t item;
  bool full;
  bool endgame;
  bool dynamic;
  size_t capacity;
  convey_alc8r_t alloc;
  int64_t stats[convey_imp_N_STATS];
} trivial_t;


/*** Memory ***/

static void
alloc_buffer(trivial_t* trivial)
{
  convey_alc8r_t* alloc = &trivial->alloc;
  void* buffer = alloc->grab(alloc->alc8r, trivial->capacity, __FILE__, __LINE__);
  trivial->item.data = buffer;
}

static void
dealloc_buffer(trivial_t* trivial)
{
  convey_alc8r_t* alloc = &trivial->alloc;
  (*alloc->free)(alloc->alc8r, trivial->item.data);
  trivial->item.data = NULL;
}

static void
free_allocs(trivial_t* trivial)
{
  dealloc_buffer(trivial);
  free(trivial);
}


/*** Methods ***/

static int
unchecked_epush(trivial_t* trivial, size_t bytes, const void* item)
{
  if (trivial->full)
    return convey_FAIL;
  trivial->item.bytes = bytes;
  memcpy(trivial->item.data, item, bytes);
  trivial->full = true;
  trivial->stats[convey_PUSHES]++;
  return convey_OK;
}

static int
trivial_push(convey_t* self, const void* item, int64_t pe)
{
  return unchecked_epush((trivial_t*) self, self->item_size, item);
}

static int
trivial_epush(convey_t* self, size_t bytes, const void* item, int64_t pe)
{
  trivial_t* trivial = (trivial_t*) self;
  if (bytes > trivial->capacity)
    return convey_imp_panic(self, __func__, convey_error_OFLO);
  return unchecked_epush(trivial, bytes, item);
}

static int
trivial_pull(convey_t* self, void* item, int64_t* from)
{
  trivial_t* trivial = (trivial_t*) self;
  if (!trivial->full)
    return convey_FAIL;
  if (trivial->item.bytes != self->item_size)
    return convey_imp_panic(self, __func__, convey_error_MISFIT);
  if (from)
    *from = 0;
  trivial->full = false;
  trivial->stats[convey_PULLS]++;
  memcpy(item, trivial->item.data, self->item_size);
  return convey_OK;
}

static void*
trivial_apull(convey_t* self, int64_t* from)
{
  trivial_t* trivial = (trivial_t*) self;
  if (!trivial->full)
    return NULL;
  if (trivial->item.bytes != self->item_size) {
    convey_imp_panic(self, __func__, convey_error_MISFIT);
    return NULL;
  }
  if (from)
    *from = 0;
  trivial->full = false;
  trivial->stats[convey_PULLS]++;
  return trivial->item.data;
}

static int
trivial_epull(convey_t* self, convey_item_t* result)
{
  trivial_t* trivial = (trivial_t*) self;
  if (!trivial->full) {
    *result = (convey_item_t) { .data = NULL };
    return convey_FAIL;
  }
  *result = trivial->item;
  trivial->full = false;
  trivial->stats[convey_PULLS]++;
  return convey_OK;
}

static int
trivial_unpull(convey_t* self)
{
  trivial_t* trivial = (trivial_t*) self;
  if (trivial->full)
    return convey_FAIL;
  trivial->full = true;
  trivial->stats[convey_UNPULLS]++;
  return convey_OK;
}

static int
trivial_advance(convey_t* self, bool done)
{
  trivial_t* trivial = (trivial_t*) self;
  trivial->stats[convey_ADVANCES]++;
  trivial->endgame = done;
  if (done)
    return trivial->full ? convey_NEAR : convey_DONE;
  return convey_OK;
}

static int
trivial_begin(convey_t* self, size_t item_size, size_t align)
{
  trivial_t* trivial = (trivial_t*) self;
  if (item_size > trivial->capacity)
    return convey_error_OFLO;
  if (trivial->dynamic)
    alloc_buffer(trivial);
  self->item_size = item_size;
  trivial->full = false;
  trivial->endgame = false;
  trivial->stats[convey_BEGINS]++;
  return convey_OK;
}

static int
trivial_reset(convey_t* self)
{
  trivial_t* trivial = (trivial_t*) self;
  if (trivial->dynamic)
    dealloc_buffer(trivial);
  convey_imp_update_stats(trivial->stats);
  return convey_OK;
}

static int
trivial_free(convey_t* self)
{
  trivial_t* trivial = (trivial_t*) self;
  free_allocs(trivial);
  return convey_OK;
}

static int64_t
trivial_statistic(convey_t* self, int which)
{
  trivial_t* trivial = (trivial_t*) self;
  return convey_imp_statistic(trivial->stats, which);
}

static const convey_methods_t trivial_methods = {
  .push = &trivial_push,
  .pull = &trivial_pull,
  .unpull = &trivial_unpull,
  .advance = &trivial_advance,
  .begin = &trivial_begin,
  .reset = &trivial_reset,
  .free = &trivial_free,
  .statistic = &trivial_statistic,
  .apull = &trivial_apull,
  .epush = &trivial_epush,
  .epull = &trivial_epull,
  .panic = &convey_imp_panic,
};

static const convey_methods_t debug_methods = {
  ._super_ = &trivial_methods,
  .push = &convey_checked_push,
  .pull = &convey_checked_pull,
  .unpull = &convey_checked_unpull,
  .advance = &trivial_advance,
  .begin = &trivial_begin,
  .reset = &trivial_reset,
  .free = &trivial_free,
  .statistic = &trivial_statistic,
  .apull = &convey_checked_apull,
  .epush = &convey_checked_epush,
  .epull = &convey_checked_epull,
  .panic = &convey_imp_panic,
};


/*** Constructor ***/

convey_t*
convey_new_trivial(size_t monster_size, const convey_alc8r_t* alloc,
                   uint64_t options)
{
  bool reckless = (options & convey_opt_RECKLESS);
  bool quiet = (options & convey_opt_QUIET);
  if (monster_size == 0)
    CONVEY_REJECT(quiet, "invalid arguments");

  if (alloc == NULL)
    alloc = &convey_imp_alloc_align;
  else if (!alloc->grab || !alloc->free)
    CONVEY_REJECT(quiet, "alloc is missing one or both methods");

  trivial_t* trivial = malloc(sizeof(trivial_t));
  if (trivial == NULL)
    CONVEY_REJECT(quiet, "a small malloc() failed!");
  *trivial = (trivial_t) {
    .convey = { ._class_ = (reckless ? &trivial_methods : &debug_methods),
                .features = convey_ELASTIC | convey_STEADY,
                .n_procs = 1, .state = convey_DORMANT,
                .suppress = UINT64_C(0) - quiet, },
    .item = { .data = NULL, .from = 0 },
    .dynamic = (options & convey_opt_DYNAMIC),
    .capacity = monster_size, .alloc = *alloc,
  };

  if (!trivial->dynamic) {
    alloc_buffer(trivial);
    if (trivial->item.data == NULL)
      CONVEY_REJECT(quiet, "memory allocation failed");
  }

  return &trivial->convey;
}
