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
/*! \file triangle_agp_opt2.upc
 * \brief Another intuitive implementation of triangle counting 
 * that uses generic global references. This implementation differs from
 * triangle_agp_opt1 only in that we do some buffering when reading remote rows.
 */

#include "triangle.h"

// this depends on the row nonzeros being sorted
// optimizations:
// 1. fetches row nonzeros in buffers (rather than one word at a time)
// 2. search for intersections in two rows is done more intelligently (since nonzeros are sorted)
// I think one thing that helps this model over conveyors/exstack is that it is more cache friendly
// since it is crawling over the l_i row many times all in one go, before moving on to row l_i + 1.
/*!
 * \brief This routine implements another AGP variant of triangle counting
 * \param *count a place to return the counts from this thread
 * \param *sr a place to return the number of shared references
 * \param *L the lower triangle of the input sparse matrix 
 * \param *U the upper triangle of the input sparse matrix 
 * \param alg 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
 * \return average run time
 */
#define PULL_BUF_DEPTH 64        //!< the number of things one pulls at a time
double triangle_agp_opt2(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg) {
  int64_t cnt=0;
  int64_t numpulled=0;
  int64_t l_i, ii, k, kk, kt, num2pull, L_i, L_j;

/*! \brief pull remote gets in a buffer with this many nonzeros */
  int64_t w[PULL_BUF_DEPTH];

  double t1 = wall_seconds();
  //foreach nonzero (i, j) in L

  for(l_i = 0; l_i < L->lnumrows; l_i++){ 
    for(k = L->loffset[l_i]; k < L->loffset[l_i + 1]; k++) {
      L_i = l_i*THREADS + MYTHREAD;
      L_j = L->lnonzero[k];
      assert( L_j < L_i );

      numpulled+=2;  // to get the offsets in their shared form
      int64_t row_start, row_end, inc, num2pull, num2pull_now, start;
      if(alg == 0){
        row_start = lgp_get_int64(L->offset, L_j);
        row_end   = lgp_get_int64(L->offset, L_j + THREADS);
        inc = PULL_BUF_DEPTH;
        num2pull = row_end - row_start;
        start = L->loffset[l_i];
      }else{
        row_start = lgp_get_int64(U->offset, L_j + THREADS) - 1;
        row_end   = lgp_get_int64(U->offset, L_j);
        inc = -PULL_BUF_DEPTH;
        num2pull = row_start - row_end;
        start = U->loffset[l_i + 1] - 1;
      }
      numpulled += num2pull;
      
      for( kk = row_start; kk >= row_end; kk += inc ){
        //num2pull = ((row_end - kk) <= PULL_BUF_DEPTH) ? row_end - kk : PULL_BUF_DEPTH;
        num2pull_now = (PULL_BUF_DEPTH > num2pull ? num2pull : PULL_BUF_DEPTH);
        num2pull -= num2pull_now;
        if(alg == 0){
          lgp_memget(w,
                     L->nonzero,
                     num2pull_now*sizeof(int64_t),
                     kk*THREADS + L_j % THREADS);
          
          for( kt = 0; kt < num2pull_now; kt++){ 
            
            for(ii = start; ii < L->loffset[l_i + 1]; ii++) {
              if( w[kt] ==  L->lnonzero[ii] ){ 
                cnt++;
                start = ii + 1;
                break;
              }
              if( w[kt] < L->lnonzero[ii] ){ // the rest are all bigger too, cause L is tidy 
                start = ii; // since pulled row is sorted, we can change start
                break;
              }
            }
          }
        
        }else{
          lgp_memget(w,
                     U->nonzero,
                     num2pull_now*sizeof(int64_t),
                     L_j%THREADS + (kk - num2pull_now)*THREADS);
          for( kt = num2pull_now - 1; kt >= 0; kt--){
            
            //if( w[kt] < L_i){kk = row_end; break;}
            
            for(ii = start; ii >= U->loffset[l_i]; ii--) {
              if( w[kt] == U->lnonzero[ii] ){ 
                cnt++;
                start = ii - 1;
                break;
              }
              if( w[kt] > U->lnonzero[ii] ){ // the rest are all bigger too, cause U is tidy 
                start = ii; // since pulled row is sorted, we can change start
                break;
              }
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
