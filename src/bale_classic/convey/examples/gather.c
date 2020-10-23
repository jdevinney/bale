// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include "example.h"

typedef struct {
  long slot;
  long value;
} packet_t;

int
main(int argc, char* argv[])
{
  example_start();

  // Parse command line and environment
  long dim = 10000;
  if (argc > 1)
    dim = strtol(argv[1], NULL, 0);
  if (MY_PROC == 0)
    printf("command: %s %ld\n(parameter is: local_array_length)\n",
           argv[0], dim);

  // Prepare local data
  int status = EXIT_FAILURE;
  long* index = malloc(dim * sizeof(long));
  long* xedni = malloc(dim * sizeof(long));
  long* source = malloc(dim * sizeof(long));
  long* target = malloc(dim * sizeof(long));
  long* original = malloc(dim * sizeof(long));
  long n_procs = PROCS;
  long my_proc = MY_PROC;
  if (index && source && target) {
    brand_t _prng;
    brand_init(&_prng, 1 + my_proc);
    for (long i = 0; i < dim; i++) {
      original[i] = brand(&_prng);
      source[i] = original[i];
      index[i] = dim * my_proc + i;
      long k = n_procs * i + my_proc;
      xedni[i] = n_procs * (k % dim) + (k / dim);
    }
  }

  convey_t* request = convey_new(SIZE_MAX, 0, NULL, convey_opt_ALERT | convey_opt_SCATTER);
  convey_t* reply = convey_new(SIZE_MAX, 0, NULL, convey_opt_ALERT);

  if (request && reply && index && source && target) {
    // Transpose, then perform inverse transpose
    for (int loop = 0; loop < 2; loop++) {
      convey_begin(request, sizeof(packet_t), alignof(packet_t));
      convey_begin(reply, sizeof(packet_t), alignof(packet_t));

      /*** START OF CONVEYOR LOOP ***/
      long n = 0;
      bool more;
      while (more = convey_advance(request, n == dim),
             more | convey_advance(reply, !more)) {\
        for (; n < dim; n++) {
          packet_t packet = { .slot = n, .value = index[n] / n_procs };
          long pe = index[n] % n_procs;
          if (! convey_push(request, &packet, pe))
            break;
        }

        packet_t* p;
        int64_t from;
        while ((p = convey_apull(request, &from)) != NULL) {
          packet_t packet = { .slot = p->slot, .value = source[p->value] };
          if (! convey_push(reply, &packet, from)) {
            convey_unpull(request);
            break;
          }
        }

        while ((p = convey_apull(reply, NULL)) != NULL)
          target[p->slot] = p->value;
      }
      /*** END OF CONVEYOR LOOP ***/

      convey_reset(reply);
      convey_reset(request);

      long* temp;
      temp = source, source = target, target = temp;
      temp = index, index = xedni, xedni = temp;
    }

    // Check that the values are back where they started
    bool ok = true;
    for (long i = 0; ok && i < dim; i++)
      if (source[i] != original[i]) {
        printf("ERROR: %ld[%ld] is wrong\n", my_proc, i);
        ok = false;
      }
    if (ok && my_proc == 0)
      printf("no errors on PE 0\n");
    fflush(stdout);

    if (ok)
      status = EXIT_SUCCESS;
  }

  convey_free(reply);
  convey_free(request);
  free(original);
  free(target);
  free(source);
  free(xedni);
  free(index);

  example_end();
  exit(status);
}
