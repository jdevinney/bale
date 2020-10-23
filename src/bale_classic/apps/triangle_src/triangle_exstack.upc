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

/*! \file triangle_exstack.upc
 * \brief Implementation of triangle counting algorithm using exstack.
 */

#include "triangle.h"

/*! \brief same as the others */
typedef struct pkg_tri_t {
  //int64_t vi;    
  int32_t w; //!< w
  int32_t vj; //!< vj
}pkg_tri_t;

/*!
\brief routine to handle the exstack push of local rows to remote rows
\param *c a place to return the number of hits.
\param *ex the extack buffers
\param *mat the input sparse matrix 
\param done the signal to exstack_proceed that this thread is done
\return the return value from exstack_proceed

NB. The matrix must be tidy.
 */
static int64_t tri_exstack_push_process(int64_t *c, exstack_t *ex, sparsemat_t * mat, int64_t done) {
  int64_t k, fromth, cnt = 0;
  struct pkg_tri_t pkg;

  exstack_exchange(ex);
  
  while(exstack_pop(ex, &pkg, &fromth)){
    //Is pkg.w on row pkg.vj
    for(k = mat->loffset[pkg.vj]; k < mat->loffset[pkg.vj + 1]; k++){
      if( pkg.w == mat->lnonzero[k] ) {
        cnt ++;
        break;
      }
      if( pkg.w < mat->lnonzero[k] ) // requires that nonzeros are increasing
        break;
    }
  }
  
  *c += cnt;

  return( exstack_proceed(ex, done) );
}

/*!
 * \brief This routine implements the exstack variant of triangle counting,
 * where one pushes the appropriate part of the local row to the remote row. 
 * \param *count a place to return the number hits
 * \param *sr a place to return the number of pushes 
 * \param *L lower triangle of the input matrix
 * \param *U upper triangle of the input matrix
 * \param alg 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
 * \param buf_cnt the number of packets in the exstack buffers
 * \return average run time
 * NB. The matrix must be tidy.
 */
double triangle_exstack_push(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg, int64_t buf_cnt) {
  exstack_t * ex = exstack_init(buf_cnt, sizeof(pkg_tri_t));
  if( ex == NULL ){return(-1.0);}
  
  int64_t cnt = 0;
  int64_t numpushed = 0;

  int64_t k, kk, pe;
  int64_t l_i, L_i, L_j;

  pkg_tri_t pkg;

  double t1 = wall_seconds();
  if(alg == 0){
    // foreach nonzero in L
    for(l_i=0; l_i < L->lnumrows; l_i++){
      for(k = L->loffset[l_i]; k < L->loffset[l_i+1]; k++){
        L_i = l_i * THREADS + MYTHREAD; // shared name for row i
        L_j = L->lnonzero[k];         // shared name for col j
        assert( L_i > L_j );
        
        //pkg.vi = L_i;
        pe = L_j % THREADS;             // thread with affinity to L_j
        pkg.vj = L_j / THREADS;         // local name for col j on pe
        for(kk = L->loffset[l_i]; kk < L->loffset[l_i + 1]; kk++){ 
          pkg.w = L->lnonzero[kk];  // possible wedge vertex
          if( pkg.w == L_j )           // all cols on L_j are less than L_j
            break;
          numpushed++;
          if( exstack_push(ex, &pkg, pe) == 0 ) {
            tri_exstack_push_process(&cnt, ex, L, 0);
            numpushed--;
            kk--;
          } 
        }
      }
    }
    while( tri_exstack_push_process(&cnt, ex, L, 1))
      ;
  }else{
    // foreach nonzero in L
    for(l_i=0; l_i < L->lnumrows; l_i++){
      for(k = L->loffset[l_i]; k < L->loffset[l_i+1]; k++){
        L_i = l_i * THREADS + MYTHREAD; // shared name for row i
        L_j = L->lnonzero[k];         // shared name for col j
        assert( L_i > L_j );

        //send U[L_i,:] to the PE who owns row U[L_j,:] and ask it to do the intersection
        pe = L_j % THREADS;             // thread with affinity to L_j
        pkg.vj = L_j / THREADS;         // local name for col j on pe
        for(kk = U->loffset[l_i]; kk < U->loffset[l_i + 1]; kk++){ 
          pkg.w = U->lnonzero[kk];  // possible wedge vertex
          numpushed++;
          if( exstack_push(ex, &pkg, pe) == 0 ) {
            tri_exstack_push_process(&cnt, ex, U, 0);
            numpushed--;
            kk--;
          } 
        }
      }
    }
    while( tri_exstack_push_process(&cnt, ex, U, 1))
      ;
  }

  lgp_barrier();
  *sr = numpushed;
  t1 = wall_seconds() - t1;
  minavgmaxD_t stat[1];
  lgp_min_avg_max_d( stat, t1, THREADS );

  *count = cnt;
  return(stat->avg);
}



/*!
\brief this routine implements the exstack pull variant of triangle counting,
       where one pulls the remote row to the local row.
\param count a place to return the number hits
\param sr a place to return the number of pushes 
\param mat the input matrix
\param alg pull from remote rows or push to remote rows
\param buf_cnt the number of packets in the exstack buffers
\return average run time

NB. The matrix must be tidy
*/
double triangle_exstack_pull(int64_t *count, int64_t *sr, sparsemat_t * mat, int64_t alg, int64_t buf_cnt) {

  exstack_t * ex_req = exstack_init(buf_cnt, sizeof(pkg_tri_t));
  exstack_t * ex_resp = exstack_init(buf_cnt, sizeof(pkg_tri_t));
  if( !ex_req || !ex_resp ){return(-1.0);}

  int64_t exstack_pushes = 0;
  int64_t l_i, L_i, L_j, k, kk, ii, pe, fromth;
  pkg_tri_t pkg;
  int64_t cnt = 0;
  double t1 = wall_seconds();
  int ret;
  
  // foreach nonzero in L
  l_i = 0;
  L_i = MYTHREAD;
  k = 0;
  while(k == mat->loffset[l_i + 1]){
    l_i++;
    L_i += THREADS;
  }
  
  while(exstack_proceed(ex_req, (k == mat->lnnz))){

    while(k < mat->lnnz){
      L_j = mat->lnonzero[k];       // shared name for col j
      assert( L_i > L_j );
      pkg.w = l_i;
      pkg.vj = L_j/THREADS;      
      pe = L_j % THREADS;
      ret = exstack_push(ex_req, &pkg, pe);
      exstack_pushes++;
      
      k++;
      while(k == mat->loffset[l_i + 1]){
        l_i++;
        L_i += THREADS;
      }

      if(ret == 1)
        break;      
    }

    exstack_exchange(ex_req);

    int mid_request = 0;
    int done_popping = 0;
    int64_t row, label, requstr;
    
    exstack_reset(ex_resp);
    
    while(exstack_proceed(ex_resp, done_popping)){

      while(!done_popping){
        if(!mid_request){
          ret = exstack_pop(ex_req, &pkg, &fromth);
          if(ret == 0){
            done_popping = 1;
          }
          else{
            //PE 'fromth' wants me to send him all columns in row 'pkg.vj' each with label 'pkg.w'
            row = pkg.vj;
            label = pkg.w;
            requstr = fromth;
            kk = mat->loffset[row];
          }
        }
      
        if(!done_popping && kk < mat->loffset[row+1]){
          pkg.w = label;
          pkg.vj = mat->lnonzero[kk++];
          mid_request = (kk == mat->loffset[row + 1]) ? 0 : 1;
          exstack_pushes++;
          if(exstack_push(ex_resp, &pkg, requstr) == 1)
            break;
        }
      }
      
      exstack_exchange(ex_resp);
      
      while(exstack_pop(ex_resp, &pkg, &fromth)){
        for(ii = mat->loffset[pkg.w]; ii < mat->loffset[pkg.w+1]; ii++){
          if(mat->lnonzero[ii] == pkg.vj)
            cnt++;
          else if(mat->lnonzero[ii] > pkg.vj)
            break;
        }
      }
    }
  }

  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  *sr = exstack_pushes;
  *count = cnt;
  return(stat->avg);
  
}

