/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
// 
//
//  All rights reserved.
//  
//   This file is a part of Bale.  For license information see the
//   LICENSE file in the top level directory of the distribution.
//  
// 
 *****************************************************************/ 
#include <assert.h>
#include <shmem.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "bolite.h"

#if HAVE_CONFIG_H
#include "config.h"
#endif
#if !HAVE_SHMEM_FREE
#define shmem_malloc(S) shmalloc(S)
#define shmem_free(P) if (P) shfree(P)
#endif
#if !HAVE_SHMEM_PUTMEM_NBI
#define shmem_put64_nbi shmem_put64
#define shmem_get64_nbi shmem_get64
#endif


static double now(void)
{
  struct timeval _t;
  gettimeofday(&_t, NULL);
  return _t.tv_sec + 1.0e-6 * _t.tv_usec;
}

static void usage(const char* prog)
{
  if (shmem_my_pe() == 0)
    fprintf(stderr, "usage: %s min_words max_words reps\n", prog);
  exit(2);
}

int
main(int argc, char* argv[])
{
  shmem_init();
  int npes = shmem_n_pes();
  int mype = shmem_my_pe();

  const size_t lg_words = 24;
  const size_t lg_iters = 20;
  brand_t _prng;
  if (argc < 4)
    usage(argv[0]);
  const size_t min_words = strtoul(argv[1], NULL, 10);
  const size_t max_words = strtoul(argv[2], NULL, 10);
  const size_t n_reps = strtoul(argv[3], NULL, 10);
  if (min_words == 0 || min_words > max_words || n_reps == 0)
    usage(argv[0]);
  brand_init(&_prng, 1);
  const divbymul32_t by_npes = _divbymul32_prep(npes);

  uint64_t* arena = shmem_malloc((sizeof(uint64_t) << lg_words) + 8 * (max_words - 1));
  assert(arena != NULL);
  uint64_t* buffer = malloc(sizeof(uint64_t) * max_words << lg_iters);
  assert(buffer != NULL);

  for (size_t w = min_words; w <= max_words; w++) {
    size_t bulk = sizeof(uint64_t) * w << lg_iters;
    for (size_t i = 0; i < ((size_t)1 << lg_words); i++)
      arena[i] = i;
    for (size_t i = 0; i < (w << lg_iters); i++)
      buffer[i] = i;

    // Benchmark puts
    double start = 0.0, elapsed;
    for (int loop = 0; loop < n_reps + 1; loop++) {
      if (loop == 1)
        start = now();
      for (size_t i = 0; i < ((size_t)1 << lg_iters); i++) {
        uint64_t r = brand(&_prng);
        uint32_t pe = r >> 32;
        pe -= npes * _divbymul32(pe, by_npes);
        uint64_t* addr = arena + (r & _maskr(lg_words));
        shmem_put64_nbi(addr, buffer + i * w, w, pe);
      }
      shmem_barrier_all();
    }
    elapsed = now() - start;
    if (mype == 0) {
      printf("put(%zu) bandwidth: %g bytes/sec\n", 8 * w, bulk * n_reps / elapsed);
      fflush(stdout);
    }

    // Benchmark gets
    for (int loop = 0; loop < n_reps + 1; loop++) {
      if (loop == 1)
        start = now();
      for (size_t i = 0; i < ((size_t)1 << lg_iters); i++) {
        uint64_t r = brand(&_prng);
        uint32_t pe = r >> 32;
        pe -= npes * _divbymul32(pe, by_npes);
        uint64_t* addr = arena + (r & _maskr(lg_words));
        shmem_get64_nbi(buffer + i * w, addr, w, pe);
      }
      shmem_barrier_all();
    }
    elapsed = now() - start;
    if (mype == 0) {
      printf("get(%zu) bandwidth: %g bytes/sec\n", 8 * w, bulk * n_reps / elapsed);
      fflush(stdout);
    }
  }

  shmem_finalize();
  return 0;
}
