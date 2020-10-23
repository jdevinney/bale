// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include "example.h"

int
main(int argc, char* argv[])
{
  example_start();

  // Parse command line and environment
  long bins = 10000;
  long load = 100000L;
  uint64_t seed = 1;
  if (argc > 1)
    bins = strtol(argv[1], NULL, 0);
  if (argc > 2)
    load = strtol(argv[2], NULL, 0);
  if (argc > 3)
    seed = strtoull(argv[3], NULL, 0);
  if (MY_PROC == 0)
    printf("command: %s %ld %ld %" PRIu64 "\n"
           "(parameters are: bins, load, seed)\n",
           argv[0], bins, load, seed);

  // Initialize local data
  brand_t _prng;
  brand_init(&_prng, (seed << 32) + MY_PROC);
  long* counts = calloc(bins, sizeof(long));
  long area = PROCS * bins;

  // Build, use, and release the conveyor
  int status = EXIT_FAILURE;
  convey_t* conveyor = convey_new(SIZE_MAX, 0, NULL,
                                  convey_opt_SCATTER | convey_opt_ALERT);
  if (conveyor && counts) {
    convey_begin(conveyor, sizeof(long), alignof(long));

    /*** START OF CONVEYOR LOOP ***/
    long n = 0;
    long index = brand(&_prng) % area;
    while (convey_advance(conveyor, n == load)) {
      for (; n < load; n++) {
        long payload = index / PROCS;
        long pe = index % PROCS;
        if (! convey_push(conveyor, &payload, pe))
          break;
        index = brand(&_prng) % area;
      }

      long* local;
      while ((local = convey_apull(conveyor, NULL)) != NULL)
        counts[*local] += 1;
    }
    /*** END OF CONVEYOR LOOP ***/

    convey_reset(conveyor);
    status = EXIT_SUCCESS;

    // Produce a modest amount of output without further communication
    long peak = 0, where = 0;;
    for (long i = 0; i < bins; i++)
      if (counts[i] > peak) {
        peak = counts[i];
        where = i;
      }
    double lambda = load * 1.0 / bins;
    double tail = 0.0;
    for (long i = 0; i < 10; i++)
      tail += exp(-lambda + (peak + i) * log(lambda) - lgamma(peak + i + 1));
    if (tail * area < 4.0) {
      printf("RESULT: %ld[%ld] = %ld\n", (long)(MY_PROC), where, peak);
      fflush(stdout);
    }
  }
  convey_free(conveyor);
  free(counts);

  example_end();
  exit(status);
}
