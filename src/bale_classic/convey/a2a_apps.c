// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "alltoallv.h"


/*** Two Conveyors, Gather ***/

// Idea: table[pos] = (c * pos) % p
// We sum up (i+1) * table[pos[i]] modulo p, and should get
// c * (sum of (i+1) * pos[i]) modulo p.

#define PRIME UINT64_C(3690948119)
#define MULTIPLIER UINT64_C(2646891621)

// Can't cast from local pointer to shared, so save the shared pointer :(
#if MPP_USE_UPC
static shared uint64_t* all_table;
#endif

uint64_t*
global_table_init(size_t echo_size, size_t entries, brand_t* prng)
{
  const size_t w = (echo_size - 1) >> 3;
  const size_t n = entries + w - 1;
#if MPP_USE_UPC
  all_table = upc_all_alloc(n * THREADS, sizeof(uint64_t));
  uint64_t* table = (uint64_t*) &all_table[MYTHREAD];
#else
  uint64_t* table = mpp_alloc(n * sizeof(uint64_t));
#endif

  // Fill the top 32 bits of each table entry with incompressible junk
  for (size_t j = 0; j < n; j++) {
    size_t i = (MY_PROC + j * PROCS) % PRIME;
    table[j] = (i * MULTIPLIER) % PRIME | (brand(prng) << 32);
  }

  mpp_barrier(1);
  return table;
}

void
global_table_free(uint64_t* table)
{
#if MPP_USE_UPC
  upcx_all_free(all_table);
#else
  mpp_free(table);
#endif
}

checksum_t
indexgather(convey_t* request, convey_t* reply, size_t reply_size, brand_t* prng,
            double load, size_t entries, const uint64_t source[entries])
{
  typedef struct { int32_t slot; uint32_t local; } query_t;
  typedef struct { int32_t slot; uint32_t data[]; } packet_t;
  assert(sizeof(query_t) == 8);
  assert(reply_size >= 12);
  packet_t* packet = malloc(reply_size);
  const size_t reply_extra = reply_size - 12;
  int rv;

  rv = convey_begin(request, sizeof(query_t), 1);
  if (rv != convey_OK)
    conveyor_bug(request, "convey_begin", rv);
  rv = convey_begin(reply, reply_size, 1);
  if (rv != convey_OK)
    conveyor_bug(reply, "convey_begin", rv);

  const size_t n_procs = PROCS;
  checksum_t sums = { 0, 0 };
  bool pending = false;
  uint64_t index = 0;

  size_t sent = 0;
  size_t bulk = ceil(load * n_procs);
  bool more;
  while (1) {
    rv = convey_advance(request, sent == bulk);
    if (rv < 0)
      conveyor_bug(request, "convey_advance", rv);
    more = (rv != convey_DONE);
    rv = convey_advance(reply, !more);
    if (rv < 0)
      conveyor_bug(reply, "convey_advance", rv);
    if (!more && rv == convey_DONE)
      break;

    for (; sent < bulk; sent++) {
      if (!pending)
        index = dbrand(prng) * (entries * n_procs);
      uint64_t local = index / n_procs;
      query_t query = { .slot = sent, .local = local };
      int64_t pe = index - local * n_procs;

      rv = convey_push(request, &query, pe);
      pending = (rv == convey_FAIL);
      if (pending)
        break;
      if (rv != convey_OK)
        conveyor_bug(request, "convey_push", rv);

      sums.sent += (sent + 1) * index;
      if (sums.sent >> 63)
        sums.sent %= PRIME;
    }

    query_t query;
    int64_t from;
    while ((rv = convey_pull(request, &query, &from)) == convey_OK) {
      int64_t local = query.local;
      packet->slot = query.slot;
      memcpy(packet->data, &source[local], 8);
      if (reply_extra > 0)
        memcpy(&packet->data[2], &source[local + 1], reply_extra);
      rv = convey_push(reply, packet, from);
      if (rv == convey_FAIL) {
        int rv1 = convey_unpull(request);
        if (rv1 != convey_OK)
          conveyor_bug(request, "convey_unpull", rv1);
        break;
      }
      if (rv != convey_OK)
        conveyor_bug(reply, "convey_push", rv);
    }
    if (rv != convey_FAIL)
      conveyor_bug(request, "convey_pull", rv);

    while ((rv = convey_pull(reply, packet, NULL)) == convey_OK) {
      uint64_t value;
      memcpy(&value, packet->data, 8);
      value &= UINT64_C(0xFFFFFFFF);
      sums.rcvd += (uint64_t)(packet->slot + 1) * value;
      if (sums.rcvd >> 63)
        sums.rcvd %= PRIME;
    }
    if (rv != convey_FAIL)
      conveyor_bug(reply, "convey_pull", rv);
  }

  rv = convey_reset(reply);
  if (rv != convey_OK)
    conveyor_bug(reply, "convey_reset", rv);
  rv = convey_reset(request);
  if (rv != convey_OK)
    conveyor_bug(request, "convey_reset", rv);
  free(packet);

  // Convert the accumulated indices into the correct value
  sums.sent %= PRIME;
  sums.sent = (sums.sent * MULTIPLIER) % PRIME;
  sums.rcvd %= PRIME;
  return sums;
}


/*** One Conveyor, Scatter ***/

checksum_t
histogram(convey_t* conveyor, size_t size, brand_t* prng, double load,
          size_t entries, uint64_t target[entries])
{
  typedef struct { int32_t slot; uint32_t data[]; } packet_t;
  assert(size >= 8);
  const size_t n_words = (size + 7) >> 3;
  packet_t* packet = malloc(sizeof(packet_t) + n_words * 2 * sizeof(uint32_t));
  packet_t* update = malloc(sizeof(packet_t) + n_words * 2 * sizeof(uint32_t));
  checksum_t sums = { 0, 0 };

  const size_t n_procs = PROCS;
  const size_t my_proc = MY_PROC;
  const size_t bulk = ceil(load * n_procs);
  size_t i, sent = 0;
  uint64_t index = 0;
  bool pending = false;
  int rv;

  rv = convey_begin(conveyor, sizeof(packet_t) + size, 1);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_begin", rv);

  while (1) {
    rv = convey_advance(conveyor, sent == bulk);
    if (rv == convey_DONE)
      break;
    if (rv < 0)
      conveyor_bug(conveyor, "convey_advance", rv);

    for (; sent < bulk; sent++) {
      if (!pending) {
        index = dbrand(prng) * (entries * n_procs);
        for (i = 0; i < n_words; i++) {
          uint64_t send = brand(prng);
          memcpy(&packet->data[2*i], &send, 8);
        }
      }
      uint64_t slot = index / n_procs;
      int64_t pe = index % n_procs;
      packet->slot = slot;

      rv = convey_push(conveyor, packet, pe);
      pending = (rv == convey_FAIL);
      if (pending)
        break;
      if (rv != convey_OK)
        conveyor_bug(conveyor, "convey_push", rv);

      uint64_t sent;
      memcpy(&sent, packet->data, 8);
      sums.sent += sent * (index + 1);
    }

    while ((rv = convey_pull(conveyor, update, NULL)) == convey_OK) {
      uint64_t rcvd;
      memcpy(&rcvd, update->data, 8);
      target[update->slot] += rcvd;
      uint64_t globl = my_proc + update->slot * n_procs;
      sums.rcvd += rcvd * (globl + 1);
    }
    if (rv != convey_FAIL)
      conveyor_bug(conveyor, "convey_pull", rv);
  }

  rv = convey_reset(conveyor);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_reset", rv);
  free(update);
  free(packet);

  return sums;
}
