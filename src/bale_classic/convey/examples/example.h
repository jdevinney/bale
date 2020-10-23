// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <inttypes.h>
#include <math.h>
#include <stdalign.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "convey.h"
#if HAVE_CONFIG_H
# include "config.h"
#endif
#include "bolite.h"

#if MPP_USE_UPC
# define PROCS THREADS
# define MY_PROC MYTHREAD
# define example_start()
# define example_end()
#elif HAVE_MPP_UTIL
# include "mpp_utilV4.h"
# define example_start() argc = mpp_util_init(argc, argv, NULL)
# define example_end() mpp_util_fin()
#elif MPP_RAW_MPI
# include <mpi.h>
extern long xmpi_n_procs, xmpi_my_proc;
extern int xmpi_init(int argc, char* argv[]);
# define PROCS xmpi_n_procs
# define MY_PROC xmpi_my_proc
# define example_start() xmpi_init(argc,argv)
# define example_end() MPI_Finalize()
#elif MPP_RAW_SHMEM
# include <shmem.h>
# define PROCS shmem_n_pes()
# define MY_PROC shmem_my_pe()
# define example_start() shmem_init()
# define example_end() shmem_finalize()
#else
# define PROCS (1L)
# define MY_PROC (0L)
# define example_start()
# define example_end()
#endif
