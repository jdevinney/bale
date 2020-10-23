// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: Include configuration settings and enable actimer.

#ifndef CONVEY_COMMON_H
#define CONVEY_COMMON_H

#if HAVE_CONFIG_H
#include "config.h"
#endif
#if HAVE_STDATOMIC_H
# include <stdatomic.h>
#endif

#define MATRIX_REMOTE_HOP 1
#if (MPP_USE_UPC || HAVE_SHMEM_PTR) && HAVE_STDATOMIC_H && \
    HAVE__ATOMIC_UINT64_T && (ATOMIC_LLONG_LOCK_FREE >= 2)
# define CONVEY_INTRANODE 1
#endif

#if MPP_USE_UPC
# include "mpp2upc.h"
#else
# if HAVE_MPP_UTIL
#  include "mpp_utilV4.h"
# elif MPP_RAW_MPI
#  include "mpp2mpi.h"
#  define MPP_USE_MPI 1
# elif MPP_RAW_SHMEM
#  include "mpp2shmem.h"
#  define MPP_USE_SHMEM 1
# else
#  include "mpp2nil.h"
#  define MPP_NO_MIMD 1
# endif
# define PARALLEL(TYPE,FIELD) TYPE FIELD
# define PARALLEL_NULLIFY(OBJECT,FIELD) \
  (OBJECT)->FIELD = NULL
# define PARALLEL_ALLOC(OBJECT,FIELD,ALLOC,SIZE,TYPE) \
  (OBJECT)->FIELD = (TYPE*) (ALLOC)->grab((ALLOC)->alc8r, \
    (SIZE) * sizeof(TYPE), __FILE__, __LINE__)
# define PARALLEL_DEALLOC(OBJECT,FIELD,ALLOC) \
  (ALLOC)->free((ALLOC)->alc8r, (void*) (OBJECT)->FIELD)
#endif

#if ENABLE_PROFILING
// this is redundant because mpp_utilV4.h includes it:
# include "mpp_utilV4_profile.h"
# define CONVEY_PROF_DECL(x) mpp_profile_t x
# define CONVEY_PROF_START mpp_profile_start
# define CONVEY_PROF_STOP mpp_profile_stop
# define CONVEY_SEND_0 PROF_OP_USER(0)
# define CONVEY_SEND_1 PROF_OP_USER_P2P
# define CONVEY_SEND_2 PROF_OP_USER(16)
#else
# define CONVEY_PROF_DECL(x)
# define CONVEY_PROF_START(sample)
# define CONVEY_PROF_STOP(sample,opcode,other_pe,data)
# define CONVEY_SEND_0 (0)
# define CONVEY_SEND_1 (1)
# define CONVEY_SEND_2 (2)
#endif

#if HAVE_MPP_UTIL
# define ACTIMER_MODULE_NAME conveyors
# define ACTIMER_SHARED
# include "actimer.h"
// Use detail level 3 for now, to tidy up mpp_util reporting
# define ACT_START(timer) actimer_start(timer, 3)
# define ACT_STOP(timer) actimer_stop(timer, 3, 1)
#else
# define ACT_START(timer)
# define ACT_STOP(timer)
#endif

#if MPP_USE_MPI
# if HAVE_MPP_UTIL
#  define mpp_comm_mpi(C) ((C).internal->mpi_comm)
# else
#  define mpp_comm_mpi(C) C
# endif
#endif

#endif
