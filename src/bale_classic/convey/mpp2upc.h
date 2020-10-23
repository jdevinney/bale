// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: Emulate a limited set of mpp_utilV4 and shmem functions within
// UPC.  Requires that config.h already be included (if HAVE_CONFIG_H).

#ifndef CONVEY_MPP2UPC_H
#define CONVEY_MPP2UPC_H

#include <stdio.h>
#include <stdlib.h>
#include <upc.h>
#if HAVE_UPC_CASTABLE_H
# include <upc_castable.h>
#endif

// A structure member that is a PARALLEL pointer provides access to a
// shared array.  The member itself is a ordinary pointer to the local
// slice of the shared array, but the corresponding shared pointer is also
// present as a member whose name is prefixed by all_.
#define PARALLEL(TYPE,FIELD) TYPE FIELD; shared TYPE all_ ## FIELD
#define PARALLEL_NULLIFY(OBJECT,FIELD) \
  (OBJECT)->all_ ## FIELD = NULL; \
  (OBJECT)->FIELD = NULL
#define PARALLEL_ALLOC(OBJECT,FIELD,ALLOC,SIZE,TYPE) \
  (OBJECT)->all_ ## FIELD = upc_all_alloc(THREADS * (SIZE), sizeof(TYPE)); \
  (OBJECT)->FIELD = ((OBJECT)->all_ ## FIELD == NULL) ? NULL : \
    (TYPE*) & (OBJECT)->all_ ## FIELD [MYTHREAD]
#define PARALLEL_DEALLOC(OBJECT,FIELD,ALLOC) \
  upcx_all_free((OBJECT)->all_ ## FIELD); \
  (OBJECT)->all_ ## FIELD = NULL

#define mpp_accum_long(L) upcx_accum_long(L)
#define mpp_and_long(L) upcx_and_long(L)
#define mpp_or_long(L) upcx_or_long(L)
#define mpp_max_long(L) upcx_max_long(L)
#define mpp_min_long(L) upcx_min_long(L)
#define mpp_alloc(S) malloc(S)
#define mpp_alloc_align(S,A,B) malloc(S)
#define mpp_free(P) free(P)
#define mpp_broadcast64(P,N,R) upcx_broadcast64(P,N,R)
#define mpp_comm_is_equal(X,Y) (1)
#define mpp_comm_is_world(X) (1)
#define mpp_barrier(X) upc_barrier
#define mpp_exit(N) upc_global_exit(N)
#define mpp_fence() upc_fence
#define mpp_rel_to_abs_proc(C,I) (I)
#define mpp_util_end() exit(0)
#define mpp_util_init(ARGC, ARGV, X) (ARGC)
#define mprint(WHO, LEVEL, ...) mfprint(stdout, WHO, 1, __VA_ARGS__)
#define mprint_np(WHO, LEVEL, ...) mfprint(stdout, WHO, 0, __VA_ARGS__)
#define mfprint(FILE, WHO, PREFIX, ...) \
  do { if ((WHO) == MYTHREAD) upcx_mfprint(FILE, PREFIX, __func__, __VA_ARGS__); } while (0)

#define MY_PROC MYTHREAD
#define PROCS THREADS
#define MPP_COMM_CURR (0)

#define shmem_fence() upc_fence
#define shmem_malloc(S) malloc(S)
#define shmem_free(P) free(P)
#define shmem_ptr(P) upc_cast(P)

typedef void mpp_alltoall_t;
typedef int mpp_comm_t;

// These reduce-to-all functions usually have names like upc_all_rdc_add_l
long upcx_accum_long(long myval);
long upcx_and_long(long myval);
long upcx_or_long(long myval);
long upcx_max_long(long myval);
long upcx_min_long(long myval);

// This collective deallocator is usually called upc_all_free()
void upcx_all_free(shared void* ptr);

// Need minimal support for printing errors
void upcx_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...);

// Broadcast function...
void upcx_broadcast64(void* data, size_t count, int root);

// upcx_alltoallv: For each pair (i,j), thread i sends a block of data from
// its local slice of send[], starting at byte offset offsets[i], to the
// part of recv[] local to thread j starting at byte offset offsets[j].
// (The offsets[] array is local but should be identical across threads.)
// The number of bytes sent is given by send_bytes[THREADS*j + i], and
// recv_bytes[THREADS*i + j] is set to this value.

int upcx_alltoallv(shared char* recv, shared size_t* recv_bytes,
                   shared char* send, shared size_t* send_bytes,
                   const size_t offsets[], const int perm[]);

int upcx_alltoall(shared char* recv, shared char* send,
		  size_t n_bytes, const int perm[]);

#endif
