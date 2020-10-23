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

/*! \file triangle_conveyor.upc
 * \brief Implementations of triangle counting algorithm using conveyors.
 */

#include "triangle.h"
/*! \brief same as the others */
typedef struct pkg_tri_t{
  int64_t w;    //!< w
  int64_t vj;   //!< vj
}pkg_tri_t;

/*!
 * \brief routine to handle the conveyor pushes to the remote thread
 * \param *c a place to return the number of hits.
 * \param *conv the conveyor
 * \param *mat the input sparse matrix 
 * \param done the signal to convey_advance that this thread is done
 * \return the return value from convey_advance
 * NB. The matrix must be tidy.
 */
static int64_t tri_convey_push_process(int64_t *c, convey_t *conv, sparsemat_t *mat, int64_t done) {
  int64_t k, cnt = 0;
  struct pkg_tri_t pkg;
    
  while(convey_pull(conv, &pkg, NULL) == convey_OK){
    //Is pkg.w on row pkg.vj
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

  return( convey_advance(conv, done) );
}

/*!
 * \brief This routine implements the conveyor variant of triangle counting,
 * where one pushes the appropriate part of the local row to the remote row. 
 * \param *count  a place to return the counts from each thread
 * \param *sr a place to return the number of shared references
 * \param *L lower triangle of the input matrix
 * \param *U upper triangle of the input matrix
 * \param alg 0,1: 0 to compute (L & L * U), 1 to compute (L & U * L).
 * \return average run time
 * NB. The matrix must be tidy.
 */
double triangle_convey_push(int64_t *count, int64_t *sr, sparsemat_t * L, sparsemat_t * U, int64_t alg) {
  
  convey_t * conv = convey_new(SIZE_MAX, 0, NULL, 0);
  if(conv == NULL){return(-1);}
  if(convey_begin( conv, sizeof(pkg_tri_t), 0 ) != convey_OK){return(-1);}

  int64_t cnt = 0;
  int64_t numpushed = 0;
  double t1 = wall_seconds();

  pkg_tri_t pkg;
  int64_t k,kk, pe;
  int64_t l_i, L_i, L_j;

  if(alg == 0){
    // foreach nonzero (i,j) in L
    for(l_i=0; l_i < L->lnumrows; l_i++) { 
      for(k=L->loffset[l_i]; k< L->loffset[l_i + 1]; k++) {
        L_i = l_i * THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        
        pe = L_j % THREADS;
        pkg.vj = L_j / THREADS;
        for(kk = L->loffset[l_i]; kk < L->loffset[l_i + 1]; kk++) {
          pkg.w = L->lnonzero[kk]; 
          if( pkg.w > L_j) 
            break;
          numpushed++;
          if(convey_push(conv, &pkg, pe) != convey_OK){
            tri_convey_push_process(&cnt, conv, L, 0); 
            kk--;
            numpushed--;
          }
        }
      }
    }
    while ( tri_convey_push_process(&cnt, conv, L, 1) ) // keep popping til all threads are done
      ;
  }else{
    for(l_i=0; l_i < L->lnumrows; l_i++) { 
      for(k=L->loffset[l_i]; k< L->loffset[l_i + 1]; k++) {
        L_i = l_i * THREADS + MYTHREAD;
        L_j = L->lnonzero[k];
        
        pe = L_j % THREADS;
        pkg.vj = L_j / THREADS;
        for(kk = U->loffset[l_i]; kk < U->loffset[l_i + 1]; kk++) {
          pkg.w = U->lnonzero[kk]; 
          numpushed++;
          if(convey_push(conv, &pkg, pe) != convey_OK){
            tri_convey_push_process(&cnt, conv, U, 0); 
            kk--;
            numpushed--;
          }
        }
      }
    }
    while ( tri_convey_push_process(&cnt, conv, U, 1) ) // keep popping til all threads are done
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

/*!
\brief This routine implements the conveyor variant of triangle counting,
       where one pulls the remote row to the local row.
\param count  a place to return the counts from each thread
\param sr a place to return the number of shared references
\param mat lower triangle of the input matrix
\return average run time

NB. The matrix must be tidy.
*/
double triangle_convey_pull(int64_t *count, int64_t *sr, sparsemat_t *mat) {

  int64_t cnt = 0, row, col, i, ii, j, k, rowcnt, pos, pos2;
  int64_t L_i, l_i, L_j, pe, frompe;

  // determine the largest nnz in a row
  int64_t max_row = 0;
  for(i = 0; i < mat->lnumrows; i++)
    max_row = ((max_row < mat->loffset[i+1] - mat->loffset[i]) ? mat->loffset[i+1] - mat->loffset[i]: max_row);
  max_row = lgp_reduce_max_l(max_row);

  struct pkg_tri_t pkg;
  int64_t * buf = calloc(max_row + 2, sizeof(int64_t));
  
  convey_t * conv_req = convey_new(SIZE_MAX, 0, NULL, 0);
  convey_t * conv_resp = convey_new_elastic((max_row + 2)*sizeof(int64_t), SIZE_MAX, 0, NULL, 0);
  
  if(!conv_req || !conv_resp){return(-1);}
  if(convey_begin( conv_req, sizeof(pkg_tri_t), 0 ) != convey_OK){return(-1);}
  if(convey_begin( conv_resp, sizeof(int64_t), 0 ) != convey_OK){return(-1);}
  
  int64_t numpushed = 0;
  bool more, new_pull = 1;
  double t1 = wall_seconds();
  convey_item_t * convitem = calloc(1, sizeof(convey_item_t));
  int64_t lastrow = -1, start;
  
  lgp_barrier();
  
  k = 0;
  l_i = 0;
  L_i = MYTHREAD;
  
  while(more = convey_advance(conv_req, k == mat->lnnz),
        more | convey_advance(conv_resp, !more)){

    /* push requests: each request is a (i,j) tuple where i is the row 
       this PE will need to look at once the response is sent back. 
       j is the row this PE is asking for. The response will be the 
       nonzeros in local row j on the PE the request is sent to */
    
    while(k < mat->loffset[l_i + 1]){
      L_j = mat->lnonzero[k]; // shared name for col j
      assert( L_i > L_j );
      pkg.w = l_i;
      pkg.vj = L_j/THREADS;
      pe = L_j % THREADS;
      if(convey_push(conv_req, &pkg, pe) != convey_OK)
        break;      
      numpushed++;
      k++;
    }

    /* we just hit the end of a row */
    if(k == mat->loffset[l_i + 1]){
      l_i++;
      L_i += THREADS;
    }

    /* pop requests and send responses */
    while(1){
      if(new_pull){
        if(convey_pull(conv_req, &pkg, &frompe) != convey_OK)
          break;
        pos = 0;
        row = pkg.vj;
        buf[pos++] = pkg.w;
        buf[pos++] = mat->loffset[row + 1] - mat->loffset[row];
        for(j = mat->loffset[row]; j < mat->loffset[row+1]; j++){
          buf[pos++] = mat->lnonzero[j];
        }
      }

      // this push can still fail, so we have to save some state, or unpop.
      // how much code messiness is the elastic conveyor really saving us?
      if(convey_epush(conv_resp, pos*sizeof(int64_t), buf, frompe) != convey_OK){
        new_pull = 0;
        break;
      }
      new_pull = 1;
      numpushed+=pos;
    }

    /* pop responses and search for collisions */
    // this depends on the row nonzeros being sorted
    while(convey_epull(conv_resp, convitem) == convey_OK){
      int64_t * data = (int64_t *)convitem->data;
      pos2 = 0;
      row = data[pos2++];
      rowcnt = data[pos2++];
      start = mat->loffset[row];
      for(i = 0; i < rowcnt; i++, pos2++){
        for(ii = start; ii < mat->loffset[row+1]; ii++){
          if(mat->lnonzero[ii] == data[pos2]){
            cnt++;
            start = ii + 1;
            break;
          }else if(mat->lnonzero[ii] > data[pos2]){
            start = ii; // since pulled row is sorted, we can change start
            break;
          }
        }
      }
    }
    
  }
  
  
  *sr = numpushed;
  *count = cnt;
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );

  free(buf);
  free(convitem);
  convey_free(conv_req);
  convey_free(conv_resp);
  
  return(stat->avg);
}
