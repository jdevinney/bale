/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
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

/*! \file triangle_exstack2.upc
 * \brief Implementation of triangle counting algorithm using exstack2.
 */
#include "triangle.h"

/*! \brief make it the same */
typedef struct pkg_tri_t {
  int64_t w;     //!< w
  int64_t vj;    //!< vj
}pkg_tri_t;

/*!
 * \brief routine to handle the exstack2 push of local rows to remote rows
 * \param *c a place to return the number of hits.
 * \param *ex2 the extack2 buffers
 * \param *mat the input sparse matrix 
 * \param done the signal to exstack2_proceed that this thread is done
 * \return the return value from exstack2_proceed
 * NB. The matrix must be tidy.
 */
static int64_t 
tri_exstack2_push_process(int64_t *c, exstack2_t *ex2, sparsemat_t *mat, int64_t done) {
  int64_t k, cnt = 0;
  struct pkg_tri_t pkg;
    
  while( exstack2_pop(ex2, &pkg, NULL) ){
    for(k = mat->loffset[pkg.vj]; k < mat->loffset[pkg.vj + 1]; k++){
       if( pkg.w == mat->lnonzero[k]){
         cnt ++;
         break;
       }
       if( pkg.w < mat->lnonzero[k]) // requires that nonzeros are increasing
         break;
    }
  }
  *c += cnt;

  return( exstack2_proceed(ex2, done) );
}


/*!
 * \brief This routine implements the exstack2 variant of triangle counting.
 * where one pushes the appropriate part of the local row to the remote row. 
 * \param *count a place to return the number hits
 * \param *sr a place to return the number of pushes 
 * \param *L lower triangle of the input matrix
 * \param *U upper triangle of the input matrix
 * \param alg 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
 * \param buf_cnt the size of the exstack2 buffers
 * \return average run time
 * NB. The matrix must be tidy.
 */
double triangle_exstack2_push(int64_t * count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg, int64_t buf_cnt) {

  exstack2_t * ex2 = exstack2_init(buf_cnt, sizeof(pkg_tri_t));
  if( ex2 == NULL ) return(-1.0);

  int64_t cnt = 0;
  int64_t numpushed = 0;
  
  int64_t k, kk, pe;
  int64_t l_i, L_i, L_j;

  pkg_tri_t pkg;

  double t1 = wall_seconds();
  if(alg == 0){
    // foreach nonzero in L
    for(l_i=0; l_i < L->lnumrows; l_i++) { 
      for(k=L->loffset[l_i]; k< L->loffset[l_i + 1]; k++) {
        L_i = l_i * THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_i > L_j );
        
        pe = L_j % THREADS;
        pkg.vj = L_j / THREADS;
        for(kk = L->loffset[l_i]; kk < L->loffset[l_i + 1]; kk++) {
          pkg.w = L->lnonzero[kk]; 
          if( pkg.w > L_j) 
            break;
          numpushed++;
          if(!exstack2_push(ex2, &pkg, pe)){
            tri_exstack2_push_process(&cnt, ex2, L, 0); 
            numpushed--;
            kk--;
          }
        }
      }
    }
    while ( tri_exstack2_push_process(&cnt, ex2, L, 1) )
      ;
  }else{
    // foreach nonzero in L
    for(l_i=0; l_i < L->lnumrows; l_i++) { 
      for(k=L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
        L_i = l_i * THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_i > L_j );
        
        pe = L_j % THREADS;
        pkg.vj = L_j / THREADS;
        for(kk = U->loffset[l_i]; kk < U->loffset[l_i + 1]; kk++) {
          pkg.w = U->lnonzero[kk]; 
          numpushed++;
          if(!exstack2_push(ex2, &pkg, pe)){
            tri_exstack2_push_process(&cnt, ex2, U, 0); 
            numpushed--;
            kk--;
          }
        }
      }
    }
    while ( tri_exstack2_push_process(&cnt, ex2, U, 1) )
      ;
  }

  lgp_barrier();
  *sr = numpushed;
  *count = cnt;
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  return(stat->avg);
}


