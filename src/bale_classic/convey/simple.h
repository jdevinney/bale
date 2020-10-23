// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#ifndef CONVEY_SIMPLE_H
#define CONVEY_SIMPLE_H

#include "convey_impl.h"
#include "sorter.h"

#ifdef MPP_USE_MPI
typedef int a2a_off_t;
#else
typedef size_t a2a_off_t;
#endif

typedef struct simple {
  convey_t convey;
  sorter_t* sorter;      // for distributing items into buffers
  bool nonempty;         // is there anything in any buffer?
  bool overflow;         // has any buffer filled up?
  bool quiet;            // are communications finished?
  bool own_a2a;          // did we create the a2a?
  bool dynamic;          // deallocate buffers when dormant?
  bool scatter;          // do we want a sorter?
  bool flip;             // which sync array should we use next?
  int64_t pull_from;     // draw from which recv buffer?
  size_t capacity;
  size_t buffer_bytes;
  size_t buffer_limit;   // item_size * capacity
  area_t* send;          // current state of send buffers
  area_t* recv;          // current state of recv buffers
  PARALLEL(size_t*, send_sizes);
  PARALLEL(size_t*, recv_sizes);
  PARALLEL(char*, send_buffers);
  PARALLEL(char*, recv_buffers);
  // Symmetric (in SHMEM case) but not shared (in UPC case):
  a2a_off_t* offsets;   // offsets (displacements) for alltoallv
  // Miscellaneous fields
  convey_alc8r_t alloc;
  const mpp_alltoall_t* a2a;
  int* perm;             // for xshmem_alltoallv or upcx_alltoallv
  long* sync[2];         // for shmemx_alltoallv
  mpp_comm_t comm;       // for error checking and MPI_Alltoallv
  int64_t stats[convey_imp_N_STATS];
} simple_t;


/*** Functions that bisimple conveyors need to call ***/

int simple_alltoallv(simple_t* simple);

void simple_reset_send_buffers(simple_t* simple);

#endif
