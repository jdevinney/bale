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
/*! \file triangle_agp_opt1.upc
 * \brief Another intuitive implementation of triangle counting 
 * that uses generic global references. This implementation differs from
 * triangle_agp only in that the shared references to the row_start and row_end
 * are pulled out of the for loop.
 */

#include "triangle.h"

/*!
 * \brief This routine implements another AGP variant of triangle counting
 * \param *count a place to return the counts from this thread
 * \param *sr a place to return the number of shared references
 * \param *L the lower triangle of the input sparse matrix 
 * \param *U the upper triangle of the input sparse matrix 
 * \param alg 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
 * \return average run time
 */
double triangle_agp_opt1(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg) {
  int64_t cnt=0;
  int64_t numpulled=0;
  int64_t l_i, ii, k, kk, w, L_i, L_j;
   
  double t1 = wall_seconds();
  //foreach nonzero (i, j) in L
  if(alg == 0){
    for(l_i = 0; l_i < L->lnumrows; l_i++){ 
      for(k = L->loffset[l_i] + 1; k < L->loffset[l_i + 1]; k++) {
        L_i = l_i*THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_j < L_i );
        
        numpulled+=2;  // to get the offsets in their shared form
        int64_t row_start = lgp_get_int64(L->offset, L_j);
        int64_t row_end   = lgp_get_int64(L->offset, L_j + THREADS);
        numpulled += row_end - row_start;
        int64_t start = L->loffset[l_i];

        for( kk = row_start; kk < row_end; kk++){
          numpulled++;
          w = lgp_get_int64(L->nonzero, L_j%THREADS + kk*THREADS);
          assert( w < L_j );          

          for(ii = start; ii < L->loffset[l_i + 1]; ii++) {
            if( w ==  L->lnonzero[ii] ){ 
              cnt++;
              start = ii + 1;
              break;
            }
            if( w < L->lnonzero[ii] ){
              start = ii;
              break;
            }
          }
        }
      }
    }
  }else{
    for(l_i = 0; l_i < L->lnumrows; l_i++){ 
      for(k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
        L_i = l_i*THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        assert( L_j < L_i );

        numpulled+=2;  // to get the offsets in their shared form
        int64_t row_start = lgp_get_int64(U->offset, L_j + THREADS) - 1;
        int64_t row_end   = lgp_get_int64(U->offset, L_j);
        int64_t start = U->loffset[l_i + 1] - 1;

        for( kk = row_start; kk >= row_end; kk--){
          numpulled++;
          w = lgp_get_int64(U->nonzero, L_j%THREADS + kk*THREADS);
          if (w < L_i) break; // there can't be any more intersections in these rows
          assert( w >= L_j );
          for(ii = start; ii >= U->loffset[l_i]; ii--) {
            if( w == U->lnonzero[ii] ){ 
              cnt++;
              start = ii - 1;
              break;
            }
            if( w > U->lnonzero[ii] ){
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
