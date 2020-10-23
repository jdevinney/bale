// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "biconvey_impl.h"
#include "private.h"

typedef struct bitensor bitensor_t;
typedef struct packet packet_t;

struct packet {
  uint32_t token;
  char item[];
};

struct bitensor {
  biconvey_t biconvey;
  convey_t* queries;
  convey_t* replies;
  size_t reorder_bytes;
  PARALLEL(char*, reorder);
  uint64_t* present;
  packet_t* query;
  packet_t* reply;  // reply->item should be strongly aligned
  uint32_t await;
  uint32_t inflight;
  uint32_t limit;
  bool dynamic;
  convey_alc8r_t alloc;
};


/*** Default Allocator for Reorder Buffer ***/

static void*
local_alloc(void* alc8r, size_t size, const char* tag, uint64_t value)
{
  return malloc(size);
}

static void
local_free(void* alc8r, void* ptr)
{
  free(ptr);
}

static const convey_alc8r_t local_alc8r = {
  .alc8r = NULL, .grab = &local_alloc, .free = &local_free,
};


/*** Working Methods ***/

static inline uint32_t
easy_mod(uint32_t x, uint32_t m)
{
  return (x >= m) ? x - m : x;
}

static int
bitensor_push(biconvey_t* self, const void* query, int64_t to)
{
  bitensor_t* b = (bitensor_t*) self;
  uint32_t inflight = b->inflight;
  if (inflight + 1 == b->limit)
    return convey_FAIL;

  convey_t* q = b->queries;
  packet_t* packet = b->query;
  packet->token = easy_mod(b->await + inflight, b->limit);
  memcpy(packet->item, query, b->biconvey.query_bytes);

  int result = convey_push(q, packet, to);
  if (result == convey_OK)
    b->inflight = inflight + 1;
  return result;
}

static int
bitensor_pull(biconvey_t* self, void* item)
{
  bitensor_t* b = (bitensor_t*) self;
  uint32_t await = b->await;
  uint64_t word = b->present[await >> 6];
  uint64_t bit = UINT64_C(1) << (await & 63);
  if (! (word & bit))
    return convey_FAIL;

  b->present[await >> 6] = word & ~bit;
  b->await = easy_mod(await + 1, b->limit);
  b->inflight -= 1;
  size_t item_size = b->biconvey.reply_bytes;
  memcpy(item, b->reorder + await * item_size, item_size);
  return convey_OK;
}

static int
bitensor_advance(biconvey_t* self, bool done)
{
  bitensor_t* b = (bitensor_t*) self;
  convey_t* q = b->queries;
  convey_t* r = b->replies;

  int result = convey_advance(q, done);
  if (result < 0)
    return result;
  done = (result == convey_DONE);
  result = convey_advance(r, done);
  if (result < 0)
    return result;

  if (!done) {
    int64_t from;
    packet_t* query;
    packet_t* reply = b->reply;
    while ((query = convey_apull(q, &from)) != NULL) {
      memcpy(&reply->token, &query->token, sizeof(query->token));
      assert(b->biconvey.answer != NULL);
      b->biconvey.answer(query->item, reply->item, b->biconvey.context);
      // FIXME: can we ensure this push will succeed?
      if (convey_push(r, reply, from) != convey_OK) {
	convey_unpull(q);
	break;
      }
    }
  }

  packet_t* reply;
  size_t item_size = b->biconvey.reply_bytes;
  while ((reply = convey_apull(r, NULL)) != NULL) {
    uint32_t token;
    memcpy(&token, &reply->token, sizeof(token));
    b->present[token >> 6] |= UINT64_C(1) << (token & 63);
    memcpy(b->reorder + token * item_size, reply->item, item_size);
  }

  if (result == convey_DONE && b->inflight)
    result = convey_NEAR;
  return result;
}


/*** Setup and Teardown ***/

static bool
alloc_reorder(bitensor_t* b)
{
  convey_alc8r_t* alloc = &b->alloc;
  PARALLEL_ALLOC(b, reorder, alloc, b->reorder_bytes, char);
  return (b->reorder != NULL);
}

static void
dealloc_reorder(bitensor_t* b)
{
  convey_alc8r_t* alloc = &b->alloc;
  PARALLEL_DEALLOC(b, reorder, alloc);
  b->reorder = NULL;
}

static int
bitensor_begin(biconvey_t* self, size_t query_bytes, size_t reply_bytes)
{
  bitensor_t* b = (bitensor_t*) self;
  if (b->queries->state != convey_DORMANT ||
      b->replies->state != convey_DORMANT)
    return convey_error_STATE;

  // the caller has checked that query_bytes and reply_bytes are positive
  size_t capacity = b->reorder_bytes / reply_bytes;
  if (capacity < 2)
    return convey_error_OFLO;
  if (b->dynamic && !alloc_reorder(b))
    return convey_error_ALLOC;

  b->present = calloc((capacity + 63) / 64, sizeof(uint64_t));
  b->query = malloc(sizeof(packet_t) + query_bytes);
  // FIXME: align this pointer more strongly?
  b->reply = malloc(sizeof(packet_t) + reply_bytes);
  if (! (b->present && b->query && b->reply))
    return convey_error_ALLOC;

  b->limit = capacity;
  b->inflight = 0;
  b->await = 0;

  const int token_bytes = sizeof(uint32_t);
  int result = convey_begin(b->queries, token_bytes + query_bytes, 1);
  if (result < 0)
    return result;
  return convey_begin(b->replies, token_bytes + reply_bytes, 1);
}

static int
bitensor_reset(biconvey_t* self)
{
  bitensor_t* b = (bitensor_t*) self;
  int err = convey_reset(b->queries);
  if (err != convey_OK)
    return err;
  err = convey_reset(b->replies);
  if (err != convey_OK)
    return err;

  free(b->reply);
  free(b->query);
  free(b->present);
  b->reply = NULL;
  b->query = NULL;
  b->present = NULL;
  if (b->dynamic)
    dealloc_reorder(b);
  return convey_OK;
}

static int
bitensor_free(biconvey_t* self)
{
  bitensor_t* b = (bitensor_t*) self;
  int err = bitensor_reset(&b->biconvey);
  if (err < 0)
    return err;

  err = convey_free(b->replies);
  if (err != convey_OK)
    return err;
  b->replies = NULL;
  err = convey_free(b->queries);
  if (err != convey_OK)
    return err;
  b->queries = NULL;

  if (!b->dynamic)
    dealloc_reorder(b);
  free(b);
  return convey_OK;
}


/*** Constructor ***/

static const biconvey_methods_t bitensor_methods = {
  .push = &bitensor_push,
  .pull = &bitensor_pull,
  .advance = &bitensor_advance,
  .begin = &bitensor_begin,
  .reset = &bitensor_reset,
  .free = &bitensor_free,
};

biconvey_t*
biconvey_new_tensor(size_t capacity, int order, size_t n_local, size_t n_buffers,
		    const convey_alc8r_t* alloc, uint64_t options)
{
  bool quiet = (options & convey_opt_QUIET);
  size_t reorder_bytes =
    convey_memory_usage(capacity, false, order, PROCS, n_local, n_buffers);
  mprint(0, 0, "biconveyor will use %zu bytes for reordering\n", reorder_bytes);

  bitensor_t* b = malloc(sizeof(bitensor_t));
  if (b == NULL)
    CONVEY_REJECT(quiet, "a small malloc() failed!");
  *b = (bitensor_t) {
    .biconvey = { ._class_ = &bitensor_methods },
    .reorder_bytes = reorder_bytes,
    .dynamic = (options & convey_opt_DYNAMIC),
    .alloc = (alloc ? *alloc : local_alc8r),
  };

  options |= convey_opt_PROGRESS;
  b->queries = convey_new_tensor(capacity, order, n_local, n_buffers, alloc, options);
  b->replies = convey_new_tensor(capacity, order, n_local, n_buffers, alloc, options);
  bool ok = (b->queries && b->replies);
  if (ok && !b->dynamic)
    ok = alloc_reorder(b);
  if (!ok) {
    convey_free(b->replies);
    convey_free(b->queries);
    free(b);
    CONVEY_REJECT(quiet, "symmetric allocation failed");
  }

  return &b->biconvey;
}
