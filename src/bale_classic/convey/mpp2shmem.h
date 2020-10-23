// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: Translate a limited set of mpp_utilV4 functions directly into
// SHMEM. Requires that config.h has already been included (if HAVE_CONFIG_H).

#ifndef CONVEY_MPP2SHMEM_H
#define CONVEY_MPP2SHMEM_H

#include <stdio.h>
#include <shmem.h>

#define mpp_accum_long(L) xshmem_accum_long(L)
#define mpp_and_long(L) xshmem_and_long(L)
#define mpp_or_long(L) xshmem_or_long(L)
#define mpp_max_long(L) xshmem_max_long(L)
#define mpp_min_long(L) xshmem_min_long(L)

#if HAVE_SHMEM_FREE
# define mpp_alloc(S) shmem_malloc(S)
# define mpp_free(P) if (P) shmem_free(P)
#else
# define mpp_alloc(S) shmalloc(S)
# define mpp_free(P) if (P) shfree(P)
#endif

#if HAVE_SHMEM_ALIGN
# define mpp_alloc_align(S,A,B) shmem_align(A,S)
#else
# define mpp_alloc_align(S,A,B) mpp_alloc(S)
#endif

#if HAVE_SHMEM_GLOBAL_EXIT
# define mpp_exit(N) shmem_global_exit(N)
#else
# define mpp_exit(N) xshmem_exit(N)
#endif

#define mpp_broadcast64(P,N,R) xshmem_broadcast64(P,N,R)
#define mpp_comm_is_equal(X,Y) (1)
#define mpp_comm_is_world(X) (1)
#define mpp_barrier(X) shmem_barrier_all()
#define mpp_fence() shmem_fence()
#define mpp_quiet() shmem_quiet()
#define mpp_rel_to_abs_proc(C,I) (I)
#define mpp_util_end() xshmem_final()
#define mpp_util_init(ARGC, ARGV, X) xshmem_init(ARGC)
#define mprint(WHO, LEVEL, ...) mfprint(stdout, WHO, 1, __VA_ARGS__)
#define mprint_np(WHO, LEVEL, ...) mfprint(stdout, WHO, 0, __VA_ARGS__)
#define mfprint(FILE, WHO, PREFIX, ...) \
  do { if ((WHO) == MY_PROC) xshmem_mfprint(FILE, PREFIX, __func__, __VA_ARGS__); } while (0)

typedef void mpp_alltoall_t;
typedef struct {
  long start;
  long stride;
} mpp_comm_t;

#define MY_PROC shmem_my_pe()
#define PROCS shmem_n_pes()
#define MPP_COMM_CURR xshmem_comm_world
extern const mpp_comm_t xshmem_comm_world;

// Function prototypes
int xshmem_init(int argc);
long xshmem_accum_long(long myval);
long xshmem_and_long(long myval);
long xshmem_or_long(long myval);
long xshmem_max_long(long myval);
long xshmem_min_long(long myval);
void xshmem_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...);
void xshmem_broadcast64(void* data, size_t count, int root);
// The arrays recv, recv_bytes, send, send_bytes must be symmetric
int xshmem_alltoallv(char* recv, size_t recv_bytes[],
                     char* send, size_t send_bytes[],
                     const size_t offsets[], const int perm[]);
// Likewise recv and send must be symmetric
int xshmem_alltoall(char* recv, char* send, size_t n_bytes, const int perm[]);
void xshmem_final(void);
void xshmem_exit(int status);

#endif
