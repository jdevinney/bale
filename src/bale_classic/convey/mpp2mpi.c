// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpp2mpi.h"

long xmpi_my_proc = 0;
long xmpi_n_procs = 0;

int
xmpi_init(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  xmpi_my_proc = rank;
  xmpi_n_procs = size;
  MPI_Barrier(MPI_COMM_WORLD);

  return argc;
}

void
xmpi_exit(void)
{
  MPI_Finalize();
  exit(0);
}

void*
xmpi_alloc_align(size_t align, size_t size)
{
  void* ptr = NULL;
  posix_memalign(&ptr, align, size);
  return ptr;
}

#define REDUCER(Name,Op) \
long \
xmpi_ ## Name ## _long(long myval) \
{ \
  long result; \
  MPI_Allreduce(&myval, &result, 1, MPI_LONG, Op, MPI_COMM_WORLD); \
  return result; \
}

REDUCER(accum, MPI_SUM)
REDUCER(and, MPI_BAND)
REDUCER(or, MPI_BOR)
REDUCER(max, MPI_MAX)
REDUCER(min, MPI_MIN)
#undef REDUCER

void
xmpi_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...)
{
  va_list args;
  va_start(args, format);
  if (prefix)
    fprintf(stream, "PE %d: %s> ", (int) xmpi_my_proc, func);
  vfprintf(stream, format, args);
  va_end(args);
  fflush(stream);
}

void
xmpi_broadcast64(void* data, size_t count, int root)
{
  MPI_Bcast(data, count, MPI_UINT64_T, root, MPI_COMM_WORLD);
}
