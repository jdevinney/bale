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

/*! \file toposort_conveyor.upc
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */

#include "toposort.h"
/*!
 * \brief This routine implements the exstack2 variant of toposort
 * \param *rperm returns the row permutation that is found
 * \param *cperm returns the column permutation that is found
 * \param *mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
 * \param *tmat the transpose of mat
 * \return average run time
 */

double toposort_matrix_convey(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat) {
  
  typedef struct pkg_topo_t{
    int64_t row;    
    int64_t col;
    int64_t level;
  }pkg_topo_t;
  typedef struct pkg_cperm_t{
    int64_t pos;
    int64_t col;
  }pkg_cperm_t;


  //T0_printf("Running Toposort with conveyors ...");
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
  int64_t lnr = (nr + THREADS - MYTHREAD - 1)/THREADS;
  int64_t lnc = (nc + THREADS - MYTHREAD - 1)/THREADS;
  int64_t i,j,row,col,curr_col,pe,fromth,ret, pos;

  int64_t * lrperm = lgp_local_part(int64_t, rperm);
  int64_t * lcperm = lgp_local_part(int64_t, cperm);

  uint64_t type_mask = 0x8000000000000000;
  int64_t * lrowqueue  = calloc(lnr, sizeof(int64_t));
  int64_t * lcolqueue  = calloc(lnc, sizeof(int64_t));
  int64_t * lcolqueue_level  = calloc(lnc, sizeof(int64_t));
  int64_t * lrowsum    = calloc(lnr, sizeof(int64_t));
  int64_t * lrowcnt    = calloc(lnr, sizeof(int64_t));
  int64_t * level      = calloc(lnr, sizeof(int64_t));
  int64_t * matched_col= calloc(lnr, sizeof(int64_t));

  convey_t * conv = convey_new(SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
  if(conv == NULL){return(-1);}

  if(convey_begin( conv, sizeof(pkg_topo_t), 0 ) != convey_OK){return(-1);}
  
  /* initialize rowsum, rowcnt, and queue (queue holds degree one rows) */
  int64_t rownext, rowlast;
  int64_t colnext, collast;
  int64_t colstart,colend;
  rownext = rowlast = colnext = collast = colstart = colend = 0;

  for(i = 0; i < mat->lnumrows; i++){
    lrowsum[i] = 0L;
    lrowcnt[i] = mat->loffset[i+1] - mat->loffset[i];
    if(lrowcnt[i] == 1){
      lrowqueue[rowlast++] = i;
      level[i] = 0;
    }
    for(j = mat->loffset[i]; j < mat->loffset[i+1]; j++)
      lrowsum[i] += mat->lnonzero[j];
  }
  
  lgp_barrier();    
  
  double t1 = wall_seconds();

  int64_t r_and_c_done = 0;
  int64_t max_level = 0, col_level;
  pkg_topo_t pkg, pkg_ptr;
  while(convey_advance(conv, (r_and_c_done == (lnr + lnc)))){
    
    int64_t full = 0;
    
    /* push the column for each degree one row */
    while(rownext < rowlast){
      row = pkg.row = lrowqueue[rownext];
      pkg.row |= type_mask;
      pkg.col = lrowsum[row];
      pkg.level = level[row];
      matched_col[row] = pkg.col;
      pe = pkg.col % THREADS;

      if(convey_push(conv, &pkg, pe) != convey_OK){
        full = 1;
        break;
      }
      r_and_c_done++;
      rownext++;
    }

    if(!full){
      while(colnext <= collast){
        if(colstart == colend){
          if(colnext == collast){break;}
          curr_col = lcolqueue[colnext];
          col_level = lcolqueue_level[colnext++];
          colstart = tmat->loffset[curr_col];
          colend = tmat->loffset[curr_col + 1];
        }
        
        row = tmat->lnonzero[colstart];
        pkg.row = row/THREADS;
        pkg.col = curr_col*THREADS + MYTHREAD;
        pkg.level = col_level;
        pe = row % THREADS;
        if(convey_push(conv, &pkg, pe) != convey_OK)
          break;
        colstart++;
        if(colstart == colend)
          r_and_c_done++;
      }
    }
    
    while(convey_pull(conv, &pkg_ptr, NULL) == convey_OK){
      if(pkg_ptr.row & type_mask){
        lcolqueue[collast] = (pkg_ptr.col)/THREADS;
        lcolqueue_level[collast++] = pkg_ptr.level;
      }else{
        lrowsum[pkg_ptr.row] -= pkg_ptr.col;
        lrowcnt[pkg_ptr.row]--;
        /* update the level for this row 
         If this ancestor was at level i, 
         it must be at least level i+1. */
        if(pkg_ptr.level + 1 > level[pkg_ptr.row]){
          level[pkg_ptr.row] = pkg_ptr.level + 1;
          if((pkg_ptr.level+1) > max_level)
            max_level = pkg_ptr.level + 1;
        }
        if(lrowcnt[pkg_ptr.row] == 1){
          lrowqueue[rowlast++] = pkg_ptr.row;
        }
      }
    }
  }
  
  //convey_reset(conv);
  convey_free(conv);
  
  /* at this point we know for each row its level and the column it was matched with.
     we need to create cperm and rperm from this information */
  max_level = lgp_reduce_max_l(max_level);

  int64_t num_levels = max_level + 1;
  int64_t * level_sizes = calloc(num_levels, sizeof(int64_t));
  int64_t * level_start = calloc(num_levels, sizeof(int64_t));
  
  int64_t total = 0;
  for(i = 0; i < lnr; i++){
    level_sizes[level[i]]++;
  }

  for(i = 0; i < num_levels; i++){
    level_start[i] = total + lgp_prior_add_l(level_sizes[i]);
    level_sizes[i] = lgp_reduce_add_l(level_sizes[i]);
    total += level_sizes[i];
  }

  lgp_barrier();

  for(i = 0; i < lnr; i++){
    lrperm[i] = (nr - 1) - level_start[level[i]]++;
  }
  
  convey_t * conv2 = convey_new(SIZE_MAX, 0, NULL, 0);
  convey_begin( conv2, sizeof(pkg_cperm_t), 0 );
  if(conv == NULL){return(-1);}
  pkg_cperm_t pkg2, pkg2_ptr;

  /* create cperm */
  i = 0;
  while(convey_advance( conv2, (i==lnr) )) {
    for( ; i < lnr; i++){
      pkg2.pos = lrperm[i];
      pkg2.col = matched_col[i];
      pe  = pkg2.col % THREADS;
      if(convey_push(conv2, &pkg2, pe) != convey_OK)
        break;
    }

    while(convey_pull(conv2, &pkg2_ptr, NULL) == convey_OK){
      lcperm[pkg2_ptr.col/THREADS] = pkg2_ptr.pos;
    }
  }
  convey_free(conv2);

  lgp_barrier();
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  free(lrowcnt);
  free(lrowsum);
  free(lrowqueue);
  free(lcolqueue);
  T0_fprintf(stderr, "num levels = %ld ", num_levels);
  return(stat->avg);
}

