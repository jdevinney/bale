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
/*! \file triangle_agp_oo.upc
 * \brief The intuitive implementation of triangle counting 
 * that uses generic global references
 */

#include "triangle.h"

/*!
 * \brief This routine implements another AGP variant of triangle counting
 * \param *count a place to return the counts from this thread
 * \param *sr a place to return the number of shared references
 * \param *mat the input sparse matrix 
 *        NB: This must be the tidy lower triangular matrix from the adjacency matrix
 * \return average run time
 */
double triangle_agp_oo(int64_t *count, int64_t *sr, sparsemat_t * mat) {
  int64_t cnt=0;
  int64_t numpulled=0;
  int64_t l_i, L_i, L_k, L_j;
  nxnz_t *nxnz_I = init_nxnz(mat); 
  nxnz_t *nxnz_J = init_nxnz(mat); 
  nxnz_t *nxnz_K = init_nxnz(mat); 

  double t1 = wall_seconds();
  for(l_i = 0; l_i < mat->lnumrows; l_i++){ 
    L_i = l_i*THREADS + MYTHREAD;
    for(first_l_nxnz(nxnz_I, l_i); has_l_nxnz(nxnz_I, l_i); incr_l_nxnz(nxnz_I, l_i) ){
      L_j = nxnz_I->col;

      numpulled += 2;
      for(first_S_nxnz(nxnz_J, L_j); has_S_nxnz(nxnz_J, L_j); incr_S_nxnz(nxnz_J, L_j) ){
        L_k = nxnz_J->col;
        numpulled++;

        for(first_l_nxnz(nxnz_K, l_i); has_l_nxnz(nxnz_K, l_i); incr_l_nxnz(nxnz_K, l_i) ){
          if( L_k == nxnz_K->col ){ 
            cnt++;
            break;
          }
          if( L_k < nxnz_K->col ){ // the rest are all bigger too, cause mat is tidy 
            break;
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
  free(nxnz_I);
  free(nxnz_J);
  free(nxnz_K);
  return(stat->avg);
}
