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
/*! \file triangle_agp.upc
 * \brief The intuitive implementation of triangle counting 
 * that uses generic global references
 */

#include "triangle.h"

/*!
 * \brief This routine implements the AGP variant of triangle counting. We are performing
 *                             L & L * L^T (if alg == 0)
 *                             L & U * U^T (if alg == 1)
 * These two computations are very similar and could have been merged into one
 * block of code but for clarity we have separated them into two code blocks.
 *
 * \param *count a place to return the counts from this thread
 * \param *sr a place to return the number of shared references
 * \param *L lower triangle of the input matrix
 * \param *U upper triangle of the input matrix
 * \param alg 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
 * \return average run time
 */

double triangle_agp(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg) {
  int64_t cnt=0;
  int64_t numpulled=0;
  int64_t l_i, ii, k, kk, w, L_i, L_j;
   
  double t1 = wall_seconds();

  if(!L){ T0_printf("ERROR: triangle_agp: NULL L!\n"); return(-1); }
  
  if(alg == 0){    
    //foreach nonzero (i, j) in L
    for(l_i = 0; l_i < L->lnumrows; l_i++){ 
      for(k = L->loffset[l_i] + 1; k < L->loffset[l_i + 1]; k++) {
        L_i = l_i*THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_j < L_i );

        // NOW: get L[L_j,:] and count intersections with L[L_i,:]
        int64_t start = L->loffset[l_i];
        int64_t kbegin = lgp_get_int64(L->offset,L_j); 
        int64_t kend   = lgp_get_int64(L->offset, L_j + THREADS);
        numpulled+=2;
        for( kk = kbegin; kk < kend; kk++){
          w = lgp_get_int64(L->nonzero, L_j%THREADS + kk*THREADS);
          numpulled++;
          assert( w < L_j );
          for(ii = start; ii < L->loffset[l_i + 1]; ii++) {
            if( w ==  L->lnonzero[ii] ){ 
              cnt++;
              start = ii + 1;
              break;
            }
            if( w < L->lnonzero[ii] ){ // the rest are all bigger because L is tidy
              start = ii;
              break;
            }
          }
        }
      }
    }
  }else{
    if(!U){ T0_printf("ERROR: triangle_agp: NULL U!\n"); return(-1); }
    
    //foreach nonzero (i, j) in L
    for(l_i = 0; l_i < L->lnumrows; l_i++){ 
      for(k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
        L_i = l_i*THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_j < L_i );
        
        // NOW: get U[L_j,:] and count intersections with U[L_i,:]
        int64_t start = U->loffset[l_i + 1] - 1;
        int64_t smallest_col_in_row_i = U->lnonzero[U->loffset[l_i]];
        numpulled+=2;
        int64_t kbegin = lgp_get_int64(U->offset, L_j + THREADS) - 1;
        int64_t kend   = lgp_get_int64(U->offset,L_j);
        for( kk = kbegin; kk >= kend; kk--){
          w = lgp_get_int64(U->nonzero, L_j%THREADS + kk*THREADS);
          numpulled++;
          if (w < smallest_col_in_row_i) break; // there can't be any more intersections in these rows
          assert( w >= L_j );
          for(ii = start; ii >= U->loffset[l_i]; ii--) {
            if( w ==  U->lnonzero[ii] ){ 
              cnt++;
              start = ii-1;
              break;
            }
            if( w > U->lnonzero[ii] ){ // the rest are all smaller because U is tidy
              start = ii;
              break;
            }
          }
        }
      }
    }
    
  }

  lgp_barrier();
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );

  *sr = numpulled; 
  *count = cnt;
  return(stat->avg);
}

