// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "biconvey_impl.h"
#include "private.h"
#include "simple.h"

typedef struct bisimple bisimple_t;

// FIXME: replay can share memory with reorder!
struct bisimple {
  biconvey_t biconvey;
  convey_t* convey;
  PARALLEL(uint32_t*, replay);
  PARALLEL(char*, reorder);
  size_t reorder_bytes;
  size_t limit;      // max replies in reorder buffer
  bool dynamic;
  size_t n_pushed;   // indexes into replay array
  size_t eldest;     // where to pull next
  size_t n_replies;  // how many replies are waiting
};


/*** Working Methods ***/

static inline size_t
easy_mod(size_t x, size_t m)
{
  return (x >= m) ? x - m : x;
}

static int
bisimple_push(biconvey_t* self, const void* query, int64_t to)
{
  bisimple_t* b = (bisimple_t*) self;
  if (b->n_replies + b->n_pushed == b->limit)
    return convey_FAIL;
  int result = convey_push(b->convey, query, to);
  if (result == convey_OK) {
    b->replay[b->n_pushed] = (uint32_t) to;
    b->n_pushed += 1;
  }
  return result;
}

static int
bisimple_pull(biconvey_t* self, void* reply)
{
  bisimple_t* b = (bisimple_t*) self;
  if (b->n_replies == 0)
    return convey_FAIL;
  size_t item_size = b->biconvey.reply_bytes;
  // FIXME: could use eldest pointer rather than index
  memcpy(reply, b->reorder + item_size * b->eldest, item_size);
  b->eldest = easy_mod(b->eldest + 1, b->limit);
  b->n_replies -= 1;
  return convey_OK;
}

static int
bisimple_advance(biconvey_t* self, bool done)
{
  bisimple_t* b = (bisimple_t*) self;
  simple_t* simple = (simple_t*) b->convey;

  // We need to know whether an exchange happened.  We do this in a very
  // hacky way, by watching the communication count.
  int64_t old_stat = simple->stats[convey_COMMS];
  int result = convey_advance(&simple->convey, done);
  int64_t new_stat = simple->stats[convey_COMMS];
  if (result < 0 || new_stat == old_stat)
    return result;

  // Now we handle all the queries.  The send buffers have been reset.
  const size_t n_procs = simple->convey.n_procs;
  const size_t query_bytes = b->biconvey.query_bytes;
  const size_t reply_bytes = b->biconvey.reply_bytes;
  for (size_t i = 0; i < n_procs; i++) {
    area_t* area = simple->recv + i;
    char* r = simple->send[i].next;
    for (char* q = area->next; q < area->limit; q += query_bytes) {
      b->biconvey.answer(q, r, b->biconvey.context);
      r += reply_bytes;
    }
    simple->send[i].next = r;
    // send[i].limit doesn't matter!
  }

  // Now exchange back, which sets up the recv areas
  int err = simple_alltoallv(simple);
  if (err != convey_OK)
    return err;
  simple_reset_send_buffers(simple);

  // Gather from recv buffers into reorder buffer
  size_t start = easy_mod(b->eldest + b->n_replies, b->limit);
  char* reply = b->reorder + start * reply_bytes;
  char* limit = b->reorder + b->limit * reply_bytes;
  for (size_t i = 0; i < b->n_pushed; i++) {
    area_t* area = simple->recv + b->replay[i];
    assert(area->next + reply_bytes <= area->limit);
    memcpy(reply, area->next, reply_bytes);
    area->next += reply_bytes;
    reply += reply_bytes;
    if (reply == limit)
      reply = b->reorder;
  }

  // Update our state: nothing more to pull; replay array is empty
  simple->pull_from = simple->convey.n_procs;
  b->n_replies += b->n_pushed;
  b->n_pushed = 0;

  // The result of the first exchange is the correct return value
  return result;
}


/*** Setup and Teardown ***/

static bool
alloc_reorder(bisimple_t* b)
{
  simple_t* simple = (simple_t*) b->convey;
  convey_alc8r_t* alloc = &simple->alloc;
  PARALLEL_ALLOC(b, reorder, alloc, b->reorder_bytes, char);
  return (b->reorder != NULL);
}

static void
dealloc_reorder(bisimple_t* b)
{
  simple_t* simple = (simple_t*) b->convey;
  convey_alc8r_t* alloc = &simple->alloc;
  PARALLEL_DEALLOC(b, reorder, alloc);
  b->reorder = NULL;
}

static int
bisimple_begin(biconvey_t* self, size_t query_bytes, size_t reply_bytes)
{
  bisimple_t* b = (bisimple_t*) self;
  simple_t* simple = (simple_t*) b->convey;
  int result = convey_begin(b->convey, query_bytes, 0);
  if (result < 0)
    return result;

  if (reply_bytes > query_bytes) {
    simple->capacity = simple->buffer_bytes / reply_bytes;
    if (simple->capacity == 0)
      return convey_error_OFLO;
    simple->buffer_limit = query_bytes * simple->capacity;
    simple_reset_send_buffers(simple);
  }

  if (b->dynamic && !alloc_reorder(b))
    return convey_error_ALLOC;
  b->limit = b->reorder_bytes / reply_bytes;

  size_t max_pushes = simple->convey.n_procs * simple->capacity;
  if (max_pushes > UINT32_MAX)
    return convey_error_TOOBIG;
  convey_alc8r_t* alloc = &simple->alloc;
  PARALLEL_ALLOC(b, replay, alloc, max_pushes, uint32_t);
  if (b->replay == NULL)
    return convey_error_ALLOC;
  return convey_OK;
}

static int
bisimple_reset(biconvey_t* self)
{
  bisimple_t* b = (bisimple_t*) self;
  simple_t* simple = (simple_t*) b->convey;
  convey_alc8r_t* alloc = &simple->alloc;
  PARALLEL_DEALLOC(b, replay, alloc);
  b->replay = NULL;
  if (b->dynamic)
    dealloc_reorder(b);
  return convey_reset(b->convey);
}

static int
bisimple_free(biconvey_t* self)
{
  bisimple_t* b = (bisimple_t*) self;
  if (b->convey->state != convey_DORMANT) {
    int result = bisimple_reset(self);
    if (result < 0)
      return result;
  }

  if (!b->dynamic)
    dealloc_reorder(b);
  int result = convey_free(b->convey);
  b->convey = NULL;
  free(b);
  return result;
}


/*** Constructor ***/

static const biconvey_methods_t bisimple_methods = {
  .push = &bisimple_push,
  .pull = &bisimple_pull,
  .advance = &bisimple_advance,
  .begin = &bisimple_begin,
  .reset = &bisimple_reset,
  .free = &bisimple_free,
};

biconvey_t*
biconvey_new_simple(size_t capacity, const convey_alc8r_t* alloc,
		    const convey_mpp_a2a_t* a2a, uint64_t options)
{
  bool quiet = (options & convey_opt_QUIET);
  size_t n_procs = PROCS;

  bisimple_t* b = malloc(sizeof(bisimple_t));
  if (b == NULL)
    CONVEY_REJECT(quiet, "a small malloc() failed!");
  *b = (bisimple_t) {
     .biconvey = { ._class_ = &bisimple_methods },
     .reorder_bytes = capacity * n_procs,
     .dynamic = (options & convey_opt_DYNAMIC),
  };
  options |= convey_opt_PROGRESS;
  b->convey = convey_new_simple(capacity, alloc, a2a, options);
  if (b->convey == NULL) {
    free(b);
    CONVEY_REJECT(quiet, "construction of simple conveyor failed");
  }
  if (!b->dynamic && !alloc_reorder(b)) {
    convey_free(b->convey);
    free(b);
    CONVEY_REJECT(quiet, "symmetric allocation failed");
  }

  return &b->biconvey;
}
