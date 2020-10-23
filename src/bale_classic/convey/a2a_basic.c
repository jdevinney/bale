// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <assert.h>
#include <math.h>
#include <stdalign.h>
#include <stdlib.h>
#include <string.h>

#include "alltoallv.h"


/*** Helper Functions ***/

static int*
create_perm(int n, brand_t* prng)
{
  int* perm = malloc(n * sizeof(int));
  for (int i = 0; i < n; i++) {
    int j = brand(prng) % (i + 1);
    if (j < i)
      perm[i] = perm[j];
    perm[j] = i;
  }
  return perm;
}

static void
randomize_mask(brand_t* prng, int entropy, size_t n, void* mask)
{
  uint8_t* bytes = mask;
  memset(bytes, 0, n);
  // randomize either bits or bytes
  if (entropy % 8 == 0) {
    for (int j = 0; j < n; j++)
      if (dbrand(prng) * (n - j) < 0.125 * entropy) {
        bytes[j] = 0xFF;
        entropy -= 8;
      }
  } else {
    int m = 8 * n;
    for (int i = 0; i < m; i++)
      if (dbrand(prng) * (m - i) < 1.0 * entropy) {
        bytes[i >> 3] |= 1 << (i & 7);
        entropy--;
      }
  }
}

// Compare point-to-point checksums and report all errors
static void
compare_checksums(size_t n_procs, size_t n_words,
		  uint64_t sendsums[n_procs][n_words],
		  uint64_t recvsums[n_procs][n_words],
		  brand_t* prng)
{
  uint64_t (*source)[n_words];
  uint64_t (*target)[n_words];
  const size_t item_bytes = n_words * sizeof(uint64_t);
  const size_t local_bytes = n_procs * item_bytes;
#if MPP_USE_UPC
  shared char* all_source = upc_all_alloc(THREADS * local_bytes, sizeof(char));
  shared char* all_target = upc_all_alloc(THREADS * local_bytes, sizeof(char));
  source = (uint64_t (*)[n_words]) &all_source[MYTHREAD];
  target = (uint64_t (*)[n_words]) &all_target[MYTHREAD];
#else
  source = mpp_alloc(local_bytes);
  target = mpp_alloc(local_bytes);
#endif

  memcpy(source, &sendsums[0][0], local_bytes);

  // Now we need a basic all-to-all (not alltoallv)
  int error = 0;
#if MPP_USE_UPC
  int* perm = create_perm(n_procs, prng);
  error = upcx_alltoall(all_target, all_source, item_bytes, perm);
  free(perm);
#elif HAVE_MPP_UTIL
  mpp_alltoall_t* a2a = mpp_alltoall_create(1);
  mpp_alltoall_init(a2a, 0, 0, 0, 0, 0);
  mpp_alltoall(a2a, target, source, item_bytes);
  mpp_alltoall_destroy(a2a);
#elif MPP_USE_SHMEM
  int* perm = create_perm(n_procs, prng);
  error = xshmem_alltoall((char*) target, (char*) source, item_bytes, perm);
  free(perm);
#elif MPP_USE_MPI
  error = MPI_Alltoall(source, item_bytes, MPI_BYTE,
		       target, item_bytes, MPI_BYTE, MPI_COMM_WORLD);
#else
  memcpy(target, source, local_bytes);
#endif
  if (error) {
    mprint(MY_PROC, 0, "failed transpose of checksums: error %d\n", error);
    mpp_exit(6);  // -convey_error_COMMS
  }

  for (size_t pe = 0; pe < n_procs; pe++)
    for (size_t i = 0; i < n_words; i++)
      if (target[pe][i] != recvsums[pe][i])
	mprint(MY_PROC, 0, "data from PE %zu damaged: [%zu] expected %016lx got %016lx\n",
	       pe, i, target[pe][i], recvsums[pe][i]);

#if MPP_USE_UPC
  upcx_all_free(all_target);
  upcx_all_free(all_source);
#else
  mpp_free(target);
  mpp_free(source);
#endif
}


/*** One Conveyor, Fixed Size ***/

static checksum_t
communicate(convey_t* conveyor, brand_t* prng, double load, size_t size,
            int entropy, convey_t* tally, double reject, bool p2p_sums)
{
  const size_t n_procs = PROCS;
  const size_t n_words = (size + 7) >> 3;
  uint64_t (*sendsums)[n_words] = NULL;
  uint64_t (*recvsums)[n_words] = NULL;
  if (p2p_sums) {
    sendsums = calloc(n_procs, n_words * sizeof(uint64_t));
    recvsums = calloc(n_procs, n_words * sizeof(uint64_t));
  }

  checksum_t sums = { 0, 0 };
  uint64_t send[n_words], recv[n_words], valid[n_words], varying[n_words];
  memset(send, 0, 8 * n_words);
  memset(recv, 0, 8 * n_words);
  memset(valid, 0, 8 * n_words);
  memset(valid, 0xFF, size);
  varying[n_words - 1] = 0;
  if (convey_features(conveyor) & convey_THRIFTY)
    randomize_mask(prng, entropy, size, varying);
  else
    memcpy(varying, valid, 8 * n_words);

  const size_t bulk = ceil(load * n_procs);
  const uint64_t self = MY_PROC;
  const uint64_t mask = _maskr(64 - _leadz(n_procs - 1));
  // make sure data generation is deterministic
  size_t i, sent = 0;
  int64_t pe = 0;
  bool pending = false;
  int rv;

  rv = convey_begin(conveyor, size, 1);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_begin", rv);

  if (tally == NULL)
    while (1) {
      rv = convey_advance(conveyor, sent == bulk);
      if (rv == convey_DONE)
        break;
      if (rv < 0)
        conveyor_bug(conveyor, "convey_advance", rv);

      for (; sent < bulk; sent++) {
        if (!pending) {
          for (i = 0; i < n_words; i++)
            send[i] = brand(prng) & varying[i];
          do pe = brand(prng) & mask;
          while (pe >= n_procs);
        }

        rv = convey_push(conveyor, send, pe);
        pending = (rv == convey_FAIL);
        if (pending)
          break;
        if (rv != convey_OK)
          conveyor_bug(conveyor, "convey_push", rv);

        for (i = 0; i < n_words; i++)
          sums.sent += (2*i + 2*(self ^ pe) + 1) * (valid[i] & send[i]);
	if (sendsums != NULL)
	  for (i = 0; i < n_words; i++)
	    sendsums[pe][i] += (valid[i] & send[i]);
      }

      int64_t from;
      while ((rv = convey_pull(conveyor, recv, &from)) == convey_OK) {
        for (i = 0; i < n_words; i++)
          sums.rcvd += (2*i + 2*(from ^ self) + 1) * (valid[i] & recv[i]);
	if (recvsums != NULL)
	  for (i = 0; i < n_words; i++)
	    recvsums[from][i] += (valid[i] & recv[i]);
      }
      if (rv != convey_FAIL)
        conveyor_bug(conveyor, "convey_pull", rv);
    }

  else {
    // To test a steady conveyor, we can yoke it with a "tally" conveyor to
    // tell each PE how many incoming messages to expect.  We wait for all
    // messages to arrive before starting the main conveyor's endgame.
    size_t rcvd = 0, goal = 0;
    bool tally_done = false;
    bool stuck = false;  // tally conveyor has a pending push?
    uint32_t one = 1;

    while (1) {
      // Advance the main conveyor
      rv = convey_advance(conveyor, tally_done && (rcvd == goal));
      if (rv == convey_DONE)
        break;
      if (rv < 0)
        conveyor_bug(conveyor, "convey_advance", rv);

      // Advance the tally conveyor unless it is already done
      if (!tally_done) {
        rv = convey_advance(tally, (sent == bulk) && !stuck);
        tally_done = (rv == convey_DONE);
        if (rv < 0)
          conveyor_bug(tally, "convey_advance", rv);
      }

      // Push loop.  Stop if we fail to push into either conveyor.
      for (; sent < bulk; sent++) {
        if (!stuck) {
          if (!pending) {
            for (i = 0; i < n_words; i++)
              send[i] = brand(prng) & varying[i];
            do pe = brand(prng) & mask;
            while (pe >= n_procs);
          }

          rv = convey_push(conveyor, send, pe);
          pending = (rv == convey_FAIL);
          if (pending)
            break;
          if (rv != convey_OK)
            conveyor_bug(conveyor, "convey_push", rv);

          for (i = 0; i < n_words; i++)
            sums.sent += (2*i + 2*(self ^ pe) + 1) * (valid[i] & send[i]);
	  if (sendsums != NULL)
	    for (i = 0; i < n_words; i++)
	      sendsums[pe][i] += (valid[i] & send[i]);
        }

        rv = convey_push(tally, &one, pe);
        if (rv < 0)
          conveyor_bug(tally, "convey_push", rv);
        stuck = (rv == convey_FAIL);
        if (stuck)
          break;
      }

      // Pull from the main conveyor and count what we get
      int64_t from;
      while ((rv = convey_pull(conveyor, recv, &from)) == convey_OK) {
        for (i = 0; i < n_words; i++)
          sums.rcvd += (2*i + 2*(self ^ from) + 1) * (valid[i] & recv[i]);
	if (recvsums != NULL)
	  for (i = 0; i < n_words; i++)
	    recvsums[from][i] += (valid[i] & recv[i]);
        rcvd++;
      }
      if (rv != convey_FAIL)
        conveyor_bug(conveyor, "convey_pull", rv);

      // Pull the tallies also
      if (!tally_done)
        while (convey_apull(tally, NULL))
          goal++;
    }
  }

  rv = convey_reset(conveyor);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_reset", rv);

  if (sent != bulk) {
    mprint(MY_PROC, 0, "convey_advance stopped early (%zu of %zu sent)\n",
           sent, bulk);
    mpp_exit(1);
  }

  if (p2p_sums) {
    compare_checksums(n_procs, n_words, sendsums, recvsums, prng);
    free(recvsums);
    free(sendsums);
  }

  return sums;
}


/*** One Conveyor, Varying Sizes ***/

static checksum_t
elasticomm(convey_t* conveyor, brand_t* prng, double load, size_t max_size,
           int entropy, convey_t* tally, double reject, bool p2p_sums)
{
  checksum_t sums = { 0, 0 };
  size_t max_words = (max_size + 7) >> 3;
  uint64_t* send = malloc(max_words * sizeof(uint64_t));
  uint64_t* recv = malloc(max_words * sizeof(uint64_t));

  uint64_t valid[8];
  memset(valid, 0, 8 * sizeof(uint64_t));
  for (int i = 0; i < 8; i++)
    memset(&valid[i], 0xFF, (i ? i : 8));

  const size_t n_procs = PROCS;
  const size_t mask = _maskr(64 - _leadz(n_procs - 1));
  const size_t bulk = ceil(load * n_procs);
  const uint64_t self = MY_PROC;

  size_t i, n_bytes = 0, sent = 0;
  int64_t pe = 0;
  bool pending = false;
  int rv;

  rv = convey_begin(conveyor, 1, 1);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_begin", rv);

  if (tally == NULL)
    while (1) {
      rv = convey_advance(conveyor, sent == bulk);
      if (rv == convey_DONE)
        break;
      if (rv < 0)
        conveyor_bug(conveyor, "convey_advance", rv);

      for (; sent < bulk; sent++) {
        if (!pending) {
          // For now, use a uniform distribution
          n_bytes = 1 + (size_t)(max_size * dbrand(prng));
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

        size_t n_words = (n_bytes + 7) >> 3;
        for (i = 0; i < n_words - 1; i++)
          sums.sent += (2*i + 2*self + 1) * send[i];
        sums.sent += (2*i + 2*self + 1) * (send[i] & valid[n_bytes & 7]);
      }

      convey_item_t _item;
      while ((rv = convey_epull(conveyor, &_item)) == convey_OK) {
        if (reject > 0.0 && dbrand(prng) < reject) {
          rv = convey_unpull(conveyor);
          if (rv != convey_OK)
            conveyor_bug(conveyor, "convey_unpull", rv);
          rv = convey_FAIL;
          break;
        }
        size_t n_words = (_item.bytes + 7) >> 3;
        recv[n_words - 1] = 0;
        memcpy(recv, _item.data, _item.bytes);
        uint64_t from = _item.from;
        for (i = 0; i < n_words; i++)
          sums.rcvd += (2*i + 2*from + 1) * recv[i];
      }
      if (rv != convey_FAIL)
        conveyor_bug(conveyor, "convey_epull", rv);
    }

  else {
    // See communicate() for an explanation of this code.
    size_t rcvd = 0, goal = 0;
    bool tally_done = false;
    bool stuck = false;
    uint32_t one = 1;

    while (1) {
      rv = convey_advance(conveyor, tally_done && (rcvd == goal));
      if (rv == convey_DONE)
        break;
      if (rv < 0)
        conveyor_bug(conveyor, "convey_advance", rv);

      if (!tally_done) {
        rv = convey_advance(tally, (sent == bulk) && !stuck);
        tally_done = (rv == convey_DONE);
        if (rv < 0)
          conveyor_bug(tally, "convey_advance", rv);
      }

      for (; sent < bulk; sent++) {
        if (!stuck) {
          if (!pending) {
            n_bytes = 1 + (size_t)(max_size * dbrand(prng));
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

          size_t n_words = (n_bytes + 7) >> 3;
          for (i = 0; i < n_words - 1; i++)
            sums.sent += (2*i + 2*self + 1) * send[i];
          sums.sent += (2*i + 2*self + 1) * (send[i] & valid[n_bytes & 7]);
        }

        rv = convey_push(tally, &one, pe);
        if (rv < 0)
          conveyor_bug(tally, "convey_push", rv);
        stuck = (rv == convey_FAIL);
        if (stuck)
          break;
      }

      convey_item_t _item;
      while ((rv = convey_epull(conveyor, &_item)) == convey_OK) {
        if (reject > 0.0 && dbrand(prng) < reject) {
          rv = convey_unpull(conveyor);
          if (rv != convey_OK)
            conveyor_bug(conveyor, "convey_unpull", rv);
          rv = convey_FAIL;
          break;
        }
        size_t n_words = (_item.bytes + 7) >> 3;
        recv[n_words - 1] = 0;
        memcpy(recv, _item.data, _item.bytes);
        uint64_t from = _item.from;
        for (i = 0; i < n_words; i++)
          sums.rcvd += (2*i + 2*from + 1) * recv[i];
        rcvd++;
      }
      if (rv != convey_FAIL)
        conveyor_bug(conveyor, "convey_epull", rv);

      if (!tally_done)
        while (convey_apull(tally, NULL))
          goal++;
    }
  }

  rv = convey_reset(conveyor);
  if (rv != convey_OK)
    conveyor_bug(conveyor, "convey_reset", rv);

  if (sent != bulk) {
    mprint(MY_PROC, 0, "convey_advance stopped early (%zu of %zu sent)\n",
           sent, bulk);
    mpp_exit(1);
  }

  free(recv);
  free(send);
  return sums;
}


/*** Wrapper Function ***/

checksum_t
basictest(convey_t* conveyor, brand_t* prng, double load, size_t size,
          int entropy, convey_t* tally, bool elastic, double reject,
	  bool p2p_sums)
{
  if (tally) {
    int rv = convey_begin(tally, sizeof(uint32_t), alignof(uint32_t));
    if (rv != convey_OK)
      conveyor_bug(tally, "convey_begin", rv);
  }

  checksum_t sums = (elastic ? &elasticomm : &communicate)
    (conveyor, prng, load, size, entropy, tally, reject, p2p_sums);

  if (tally) {
    int rv = convey_reset(tally);
    if (rv != convey_OK)
      conveyor_bug(tally, "convey_reset", rv);
  }

  return sums;
}

