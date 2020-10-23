/******************************************************************
//
//
//  Copyright(C) 2019-2020, Institute for Defense Analyses
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
/*! \file libgetput.upc
 * \brief some standard parallel programming support functions
 */
#include "libgetput.h"

#if __UPC__ && __UPC_ATOMIC__ && !( __cray__ || _CRAYC )
// this is relevant for BUPC or GUPC
#include <upc_atomic.h>
upc_atomicdomain_t * lgp_atomic_domain;
#endif


/*!
 * \brief Wrapper for atomic add to help with ifdef noise
 * \ingroup libgetputgrp
 */
void lgp_atomic_add(SHARED int64_t * ptr, int64_t index, int64_t value) {
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //shmem_int64_atomic_add(&ptr[lindex], value, pe);
  shmem_atomic_add(&ptr[lindex], value, (int)pe);
#elif __cray__ || _CRAYC
  _amo_aadd(&ptr[index], value);
#elif __BERKELEY_UPC_RUNTIME__
  bupc_atomicI64_fetchadd_relaxed(&ptr[index], value);
#elif __UPC__
  upc_atomic_relaxed(lgp_atomic_domain, NULL, UPC_ADD, &ptr[index], &value, NULL);
#endif
}

/*!
 * \brief Wrapper for non-blocking atomic add to help with ifdef noise
 * \ingroup libgetputgrp
 */
void lgp_atomic_add_async(SHARED int64_t * ptr, int64_t index, int64_t value){
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //shmem_int64_atomic_add(&ptr[lindex], value, pe);
  shmem_atomic_add(&ptr[lindex], value, (int)pe);
#elif __cray__ || _CRAYC
#pragma pgas defer_sync
  _amo_aadd_upc(&ptr[index], value);
#elif __BERKELEY_UPC_RUNTIME__
  bupc_atomicI64_fetchadd_relaxed(&ptr[index], value);
#endif
}

/*!
 * \brief Wrapper for atomic fetch and inc to help with ifdef noise
 * \ingroup libgetputgrp
 */
int64_t lgp_fetch_and_inc(SHARED int64_t * ptr, int64_t index) {
  int64_t ret;
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //ret = shmem_int64_atomic_fetch_inc(&ptr[lindex], pe);
  ret = shmem_atomic_fetch_inc(&ptr[lindex], (int)pe);
#elif __cray__ || _CRAYC
  ret = _amo_afadd(&ptr[index], 1L);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_fetchadd_relaxed(&ptr[index], 1L);
#elif __UPC__
  upc_atomic_relaxed(lgp_atomic_domain, &ret, UPC_INC, &ptr[index], NULL, NULL);
#endif
  return(ret);
}

/*!
 * \brief Wrapper for atomic fetch and inc to help with ifdef noise
 * \ingroup libgetputgrp
 */
int64_t lgp_fetch_and_add(SHARED int64_t * ptr, int64_t index, int64_t value) {
  int64_t ret;
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //ret = shmem_int64_atomic_fetch_add(&ptr[lindex], value, pe);
  ret = shmem_atomic_fetch_add(&ptr[lindex], value, (int)pe);
#elif __cray__ || _CRAYC
  ret = _amo_afadd(&ptr[index], value);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_fetchadd_relaxed(&ptr[index], value);
#elif __UPC__
  upc_atomic_relaxed(lgp_atomic_domain, &ret, UPC_ADD, &ptr[index], &value, NULL);
#endif
  return(ret);
}

/*!
 * \brief Wrapper for atomic compare and swap to help with ifdef noise
 * \return the old value
 * \ingroup libgetputgrp
 */
int64_t lgp_cmp_and_swap(SHARED int64_t * ptr, int64_t index, int64_t cmp_val, int64_t swap_val) {
  int64_t ret;
#if USE_SHMEM
  int64_t lindex = index/shmem_n_pes();
  int64_t pe = index % shmem_n_pes();
  //ret = shmem_int64_atomic_compare_swap(&ptr[lindex], cmp_val, swap_val, pe);
  ret = shmem_atomic_compare_swap(&ptr[lindex], cmp_val, swap_val, (int)pe);
#elif __cray__ || _CRAYC
  ret = _amo_acswap_upc(&ptr[index], cmp_val, swap_val);
#elif __BERKELEY_UPC_RUNTIME__
  ret = bupc_atomicI64_cswap_relaxed(&ptr[index], cmp_val, swap_val);
#elif __UPC__
  upc_atomic_relaxed(lgp_atomic_domain, &ret, UPC_CSWAP, &ptr[index], &cmp_val, &swap_val);
#endif
  return(ret);
}


/******************************************************************************************/
/* COLLECTIVES */
/******************************************************************************************/
#if __UPC__

/*!
 * \brief function to print a banner and do initialization if needed.
 * \param argc  from main
 * \param argv  from main
 * \return the old value
 * \ingroup libgetputgrp
 */
void lgp_init(int argc, char *argv[]) {
  
  setlocale(LC_NUMERIC,"");

#if __UPC_ATOMIC__ && !( __cray__ || _CRAYC )
  lgp_atomic_domain = upc_all_atomicdomain_alloc(UPC_INT64, UPC_ADD | UPC_INC | UPC_MAX | UPC_MIN | UPC_CSWAP, 0);
#endif
}


/*!
 * \brief function to shutdown a model if needed
 * \ingroup libgetputgrp
 */
void lgp_finalize(){
  return;
}

/*! \brief macro magic to generate all the reduction operation of all the types */
#define Define_Reducer( NAME, XTYPE, STYPE, RED_FUNC, UPC_FUNC)                    \
  XTYPE NAME (XTYPE myval) {                                                       \
    static shared STYPE *dst=NULL, * src;                                          \
    if (dst==NULL) {                                                               \
      dst = upc_all_alloc(THREADS, sizeof(STYPE));                                 \
      src = upc_all_alloc(THREADS, sizeof(STYPE));                                 \
    }                                                                              \
    src[MYTHREAD] = myval;                                                         \
    upc_barrier;                                                                   \
    RED_FUNC(dst,src,UPC_FUNC, THREADS, 1, NULL, UPC_IN_NOSYNC || UPC_OUT_NOSYNC); \
    upc_barrier;                                                                   \
    return dst[0]; }

Define_Reducer(lgp_reduce_add_l, int64_t, int64_t, upc_all_reduceL, UPC_ADD)
Define_Reducer(lgp_reduce_min_l, int64_t, int64_t, upc_all_reduceL, UPC_MIN)
Define_Reducer(lgp_reduce_max_l, int64_t, int64_t, upc_all_reduceL, UPC_MAX)

Define_Reducer(lgp_reduce_or_int, int, int, upc_all_reduceL, UPC_OR)

Define_Reducer(lgp_reduce_add_d, double, double, upc_all_reduceD, UPC_ADD)
Define_Reducer(lgp_reduce_min_d, double, double, upc_all_reduceD, UPC_MIN)
Define_Reducer(lgp_reduce_max_d, double, double, upc_all_reduceD, UPC_MAX)




/******************************************************************************************/
/******************************************************************************************/
#elif USE_SHMEM

void lgp_init(int argc, char *argv[]) {
  shmem_init();
  setlocale(LC_NUMERIC,"");
}

void lgp_finalize(){
  shmem_finalize();
}

static void *setup_shmem_reduce_workdata(long **psync, size_t xsize) {
  int *work;
  int i;
  
  work=shmem_malloc(_SHMEM_REDUCE_MIN_WRKDATA_SIZE*xsize);
  *psync=shmem_malloc(_SHMEM_REDUCE_SYNC_SIZE*sizeof(long));
  for(i=0;i<_SHMEM_REDUCE_SYNC_SIZE;++i) {
    (*psync)[i]=_SHMEM_SYNC_VALUE;
  }
  shmem_barrier_all();
  return work;
}

/* Macro to define wrapper around shmem reduce functions */
#define Define_Reducer( NAME, XTYPE, STYPE, CTYPE, SHMEM_FUNC)      \
  XTYPE NAME (XTYPE myval) {                                        \
    static STYPE *buff=NULL, *work; static long *sync;              \
    assert(sizeof(STYPE) == sizeof(CTYPE));                         \
    if (buff==NULL) {                                               \
      buff=shmem_malloc(2*sizeof(STYPE));                           \
      work=setup_shmem_reduce_workdata(&sync,sizeof(STYPE));        \
    }                                                               \
    buff[0]=myval;                                                  \
    SHMEM_FUNC(&buff[1],buff,1,0,0,shmem_n_pes(),work,sync);        \
    shmem_barrier_all();                                            \
    return buff[1]; }

Define_Reducer(lgp_reduce_add_l, int64_t, int64_t, long, shmem_long_sum_to_all)
Define_Reducer(lgp_reduce_min_l, int64_t, int64_t, long, shmem_long_min_to_all)
Define_Reducer(lgp_reduce_max_l, int64_t, int64_t, long, shmem_long_max_to_all)

Define_Reducer(lgp_reduce_add_d, double, double, double, shmem_double_sum_to_all)
Define_Reducer(lgp_reduce_min_d, double, double, double, shmem_double_min_to_all)
Define_Reducer(lgp_reduce_max_d, double, double, double, shmem_double_max_to_all)

Define_Reducer(lgp_reduce_or_int, int, int, int, shmem_int_or_to_all)

/*!
* \ingroup libgetputgrp
*/
void lgp_shmem_write_upc_array_int64(SHARED int64_t *addr, size_t index, size_t blocksize, int64_t val) {
  int pe;
  size_t local_index;
  int64_t *local_ptr;


  pe = index % shmem_n_pes();
  local_index = (index / shmem_n_pes())*blocksize;

  local_ptr =(int64_t*)(( (char*)addr ) + local_index);

  shmem_int64_p ( local_ptr, val, pe );
}

void lgp_shmem_write_upc_array_double(SHARED double *addr, size_t index, size_t blocksize, double val) {
  int pe;
  size_t local_index;
  double *local_ptr;


  pe = index % shmem_n_pes();
  local_index = (index / shmem_n_pes())*blocksize;

  local_ptr =(double*)(( (char*)addr ) + local_index);

  shmem_double_p ( local_ptr, val, pe );
}

/*!
* \ingroup libgetputgrp
*/
int64_t lgp_shmem_read_upc_array_int64(const SHARED int64_t *addr, size_t index, size_t blocksize) {
  int pe;
  size_t local_index;
  int64_t *local_ptr;


  pe = index % shmem_n_pes();
  local_index = (index / shmem_n_pes())*blocksize;

  local_ptr =(int64_t*)(( (char*)addr ) + local_index);

  return shmem_int64_g ( local_ptr, pe );
}

/*!
* \ingroup libgetputgrp
*/
double lgp_shmem_read_upc_array_double(const SHARED double *addr, size_t index, size_t blocksize) {
  int pe;
  size_t local_index;
  double *local_ptr;

  /* asupc_init tests that (long long) == (int64_t) */

  pe = index % shmem_n_pes();
  local_index = (index / shmem_n_pes())*blocksize;

  local_ptr =(double*)(( (char*)addr ) + local_index);

  return shmem_double_g ( local_ptr, pe );
}

#endif

/******************************************************************************************/
/******************************************************************************************/

/*!
\brief Compute partial sums across threads.  
\param x input value \f$x_m\f$
\return \f$\sum_{i<=m} x_i\f$  where \f$m\f$ is MYTHREAD

NB: This function must be called on all threads.
* \ingroup libgetputgrp
*/
int64_t lgp_partial_add_l(int64_t x) {

  SHARED int64_t * tmp = lgp_all_alloc(THREADS, sizeof(int64_t));
  int64_t out = 0;
  
  lgp_put_int64(tmp, MYTHREAD, x);

  lgp_barrier();

  for (int i = 0; i <= MYTHREAD; i++) {
    out += lgp_get_int64(tmp, i);
  }

  lgp_barrier();

  return out;
}

/*! 
\brief Compute prior partial sums (not including this value) across threads.  
\param x the value on <tt>MYTHREAD</tt>
\return \f$\sum_{i<m} x_i\f$, where \f$x_m\f$ is <tt>MYTHREAD</tt>
NB: This function must be called on all threads.
\ingroup libgetputgrp
*/
int64_t lgp_prior_add_l(int64_t x) {
  return lgp_partial_add_l(x) - x;
}


/*!
 * \brief This routine finds the min average and max of a collection 
 *  of myval's (int64_t's) across all threads
 *
 * This routine is collective and fills in a struct that holds the min average max
 * of the values given by each thread.
 *
 * \param s struct hold the computed reductions 
 * \param myval The value to contributed by MYTHREAD
 * \param dem   the denominator for the average.
 * \return to the minavgmax struct the computed reductions
 * \ingroup libgetputgrp
 *
 */
int64_t lgp_min_avg_max_l(minavgmaxL_t *s, int64_t myval, int64_t dem) {
    long retval=0;

    s->min = lgp_reduce_min_l(myval);
    s->max = lgp_reduce_max_l(myval);
    s->avg = lgp_reduce_add_l(myval) / dem;
    return( retval );
}

/*!
 * \brief This routine finds the min average and max of a collection 
 *  of myval's (int64_t's) across all threads
 *
 * This routine is collective and fills in a struct that holds the min average max
 * of the values given by each thread.
 * \param s struct hold the computed reductions 
 * \param myval The value to contributed by MYTHREAD
 * \param dem   the denominator for the average.
 * \return to the minavgmax struct the computed reductions
 * \ingroup libgetputgrp
 *
 */
int64_t lgp_min_avg_max_d(minavgmaxD_t * s, double myval, int64_t dem){
  long retval = 0;
  s->min = lgp_reduce_min_d(myval);
  s->max = lgp_reduce_max_d(myval);
  s->avg = lgp_reduce_add_d(myval)/dem;  
  return( retval );
}

/*! 
 * \brief This routine uses gettimeofday routine to give access 
 *  to the wall clock timer on most UNIX-like systems.
 * \ingroup libgetputgrp
*/
double wall_seconds() {
  struct timeval tp;
  int retVal = gettimeofday(&tp,NULL);
  if (retVal == -1) { perror("gettimeofday:"); fflush(stderr); }
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}



#define USE_KNUTH   /*!< Default define to set whether we use the Knuth random number generator or rand48 */
#ifdef USE_KNUTH
#define LGP_RAND_MAX 2251799813685248  /*!< max random number depends on which rng we use */
#include "knuth_rng_double_2019.h"
#else
#define LGP_RAND_MAX 281474976710656
#endif

/*! 
 * \brief seed for the random number generator
 * \param seed the seed
 * Note: if all thread call this with the same seed they actually get different seeds.
 */
void lgp_rand_seed(int64_t seed){
#ifdef USE_KNUTH
  ranf_start(seed + 1 + MYTHREAD);
#else
  srand48(seed + 1 + MYTHREAD);
#endif
}

/*! 
 * \brief return a random integer mod N.
 * \param N the modulus
 */
int64_t lgp_rand_int64(int64_t N){
  assert(N < LGP_RAND_MAX);
#ifdef USE_KNUTH
  return((int64_t)(ranf_arr_next()*N));
#else
  return((int64_t)(drand48()*N));
#endif
}

/*! 
 * \brief return a random double in the interval (0,1]
 */
double lgp_rand_double(){
#ifdef USE_KNUTH
  return(ranf_arr_next());
#else
  return(drand48());
#endif
}



/*!
 * \brief This routine finds the pe and local index on that pe for given global index.
 * We support a BLOCK and a CYCLIC distribution of indices to pe's.
 * 
 * \param pe address into which we write the pe
 * \param lidx address into which we write the local index
 * \param gidx the global index
 * \param n the global number of indices (needed for the block distribution calculation)
 * \param layout the enum for BLOCK or CYCLIC
 * \ingroup libgetputgrp
 */
void global_index_to_pe_and_offset(int64_t *pe, int64_t *lidx, int64_t gidx, int64_t n, layout layout)
{
  if(layout == CYCLIC){
    *pe = gidx % THREADS;
    *lidx = gidx / THREADS;
    return;
  }
  int64_t idx;
  int64_t upper_points_per_pe = n / THREADS + ((n % THREADS > 0) ? 1 : 0);
  int64_t rem = n % THREADS;

  if( (rem == 0) || (gidx / upper_points_per_pe < rem) ){
    *pe = gidx / upper_points_per_pe;
    *lidx = gidx % upper_points_per_pe;
  }else{
    *pe = (gidx - rem) / (upper_points_per_pe - 1);
    *lidx= (gidx - rem) % (upper_points_per_pe - 1);
  }
  return;
}

/*!
 * \brief This routine finds the global index for a given pe and local index.
 * We support a BLOCK and a CYCLIC distribution of indices to pe's.
 * 
 * \param pe the given pe
 * \param lidx the local index on that pe
 * \param n the global number of indices (needed for the block distribution calculation)
 * \param layout the enum for BLOCK or CYCLIC
 * \return the global index
 * \ingroup libgetputgrp
 */
int64_t pe_and_offset_to_global_index(int64_t pe, int64_t lidx, int64_t n, layout layout)
{
  if(layout == CYCLIC)
    return(lidx*THREADS + pe);
  else{
    int64_t i, gidx = 0;
    for(i = 0; i < pe; i++)
      gidx += (n + THREADS - i - 1)/THREADS;
    gidx += lidx;
    return(gidx);
  }
}

