// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#ifndef CONVEY_PRIVATE_H
#define CONVEY_PRIVATE_H

#include <stddef.h>
#include <stdint.h>
#include "common.h"
#include "convey.h"
#include "bolite.h"

#ifndef MPP_USE_SHMEM
# undef HAVE_SHMEMX_ALLTOALLV
# undef HAVE_SHMEMX_PUTMEM_SIGNAL
# undef HAVE_SHMEMX_TEAM_ALLTOALLV
#else
# if (HAVE_SHMEMX_ALLTOALLV || HAVE_SHMEMX_TEAM_ALLTOALLV || HAVE_SHMEMX_PUTMEM_SIGNAL)
#  include <shmemx.h>
# endif
# if !HAVE_SHMEM_FREE
// Fall back to old-style names
#  define shmem_free shfree
#  define shmem_malloc shmalloc
#  define shmem_align shmemalign
# endif
#endif

enum convey_imp_timers {
  simple_a2a = 0,
  matrix_pivot,
  tensor_early,
  tensor_late,
};


/*** Memory Allocation ***/

extern const convey_alc8r_t convey_imp_alloc, convey_imp_alloc_align;

// Compute the required alignment from the item size and the alignment
// request (coming from convey_begin, and validated).  If it exceeds the
// alignment automatically guaranteed by the item size and the header size,
// then allocate an aligned buffer of item_size bytes and store its address
// in *handle.  On entry, *handle must be NULL.  The return value is false
// iff allocation failed.
bool
convey_prep_aligned(void** handle, size_t item_size, size_t header_size,
                    size_t align_wanted);


/*** Processes Per Node ***/

size_t
convey_procs_per_node(void);


/*** Error Handling ***/

#define CONVEY_REJECT(QUIET,MSG) \
  do { if (!(QUIET)) mprint(0, 0, "%s\n", (MSG)); return NULL; } while(0)

int
convey_imp_panic(convey_t* self, const char* where, int error);


/*** Statistics ***/

enum convey_imp_statistic {
  convey_imp_N_STATS = 2 * convey_CUMULATIVE
};

// Retrieve a statistic if it exists
int64_t
convey_imp_statistic(int64_t stats[], int which);

// Add transient stats into persistent stats, then reset them
void
convey_imp_update_stats(int64_t stats[]);


/*** Fallback Conveyor ***/

convey_t*
convey_new_trivial(size_t monster_size, const convey_alc8r_t* alloc,
                   uint64_t options);


#endif
