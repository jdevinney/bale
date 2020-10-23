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


static checksum_t
rigidalign(convey_t* conveyor, brand_t* prng, double load, size_t size)
{
  const size_t n_procs = PROCS;
  const size_t n_words = (size + 7) >> 3;
  uint64_t send[n_words], valid[n_words];
  memset(send, 0, 8 * n_words);
  memset(valid, 0, 8 * n_words);
  memset(valid, 0xFF, size);

  const size_t mask = _maskr(64 - _leadz(n_procs - 1));
  const size_t bulk = ceil(load * n_procs);
  uintptr_t ptrbits = 0;
  size_t i, sent = 0;
  int64_t pe = 0;
  bool pending = false;
  int rv;

  rv = convey_begin(conveyor, size, 0);
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
	for (i = 0; i < n_words; i++)
	  send[i] = brand(prng);
	do pe = brand(prng) & mask;
	while (pe >= n_procs);
      }

      rv = convey_push(conveyor, send, pe);
      pending = (rv == convey_FAIL);
      if (pending)
	break;
      if (rv != convey_OK)
	conveyor_bug(conveyor, "convey_push", rv);
    }

    void* packet;
    while ((packet = convey_apull(conveyor, NULL)) != NULL)
      ptrbits |= (uintptr_t) packet;
  }

  rv = convey_reset(conveyor);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_reset", rv);

  if (sent != bulk) {
    mprint(MY_PROC, 0, "convey_advance stopped early (%zu of %zu sent)\n",
	   sent, bulk);
    mpp_exit(1);
  }

  ptrbits &= (~size & (size - 1) & (CONVEY_MAX_ALIGN - 1));
  return (checksum_t) { 0, ptrbits };
}

static checksum_t
flexalign(convey_t* conveyor, brand_t* prng, double load, size_t max_size)
{
  const size_t max_words = (max_size + 7) >> 3;
  uint64_t* send = malloc(max_words * sizeof(uint64_t));

  // all item sizes will be a multiple of 'align'
  const size_t align = 1 + ((max_size - 1) & ~max_size);
  const size_t n_procs = PROCS;
  const size_t mask = _maskr(64 - _leadz(n_procs - 1));
  const size_t bulk = ceil(load * n_procs);

  uint64_t ptrbits = 0;
  size_t i, n_bytes = 0, sent = 0;
  int64_t pe = 0;
  bool pending = false;
  int rv;

  rv = convey_begin(conveyor, align, align);
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
	// For now, use a uniform distribution
	n_bytes = align + (size_t)(max_size * dbrand(prng));
	n_bytes &= -align;
	assert(n_bytes <= max_size);
	size_t n_words = (n_bytes + 7) >> 3;
	for (i = 0; i < n_words; i++)
	  send[i] = brand(prng);
	do pe = brand(prng) & mask;
	while (pe >= n_procs);
      }

      rv = convey_epush(conveyor, n_bytes, send, pe);
      pending = (rv == convey_FAIL);
      if (pending)
	break;
      if (rv != convey_OK)
	conveyor_bug(conveyor, "convey_epush", rv);
    }

    convey_item_t _item;
    while ((rv = convey_epull(conveyor, &_item)) == convey_OK)
      ptrbits |= (uintptr_t) _item.data;
    if (rv != convey_FAIL)
      conveyor_bug(conveyor, "convey_epull", rv);
  }

  rv = convey_reset(conveyor);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_reset", rv);

  if (sent != bulk) {
    mprint(MY_PROC, 0, "convey_advance stopped early (%zu of %zu sent)\n",
	   sent, bulk);
    mpp_exit(1);
  }

  free(send);
  ptrbits &= (~align & (align - 1) & (CONVEY_MAX_ALIGN - 1));
  return (checksum_t) { 0, ptrbits };
}


/*** Wrapper Function ***/

checksum_t
aligntest(convey_t* conveyor, brand_t* prng, double load,
	  size_t size, bool elastic)
{
  return (elastic ? &flexalign : &rigidalign)(conveyor, prng, load, size);
}
