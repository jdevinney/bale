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

/*! \file toposort_exstack.upc
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */

#include "toposort.h"

/*!
\brief Implements the exstack variant of toposort. It does not use
       any atomic ops, but rather sends extra data around to construct the permutations.
\param rperm returns the row permutation that is found
\param cperm returns the column permutation that is found
\param mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
\param tmat the transpose of mat
\param buf_cnt the number of packages in the exstack buffer
\return average run time
*/
double toposort_matrix_exstack(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat, int64_t buf_cnt) {
  typedef struct pkg_topo_t{
    int64_t row;    
    int64_t col;
    int64_t level;
  }pkg_topo_t;

  typedef struct pkg_cperm_t{
    int64_t pos;
    int64_t col;
  }pkg_cperm_t;

  //T0_printf("Running Toposort with exstack ...");
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
  int64_t lnr = (nr + THREADS - MYTHREAD - 1)/THREADS;
  int64_t lnc = (nc + THREADS - MYTHREAD - 1)/THREADS;
  int64_t i,j,row,col,curr_col,pe,fromth,ret, pos;

  int64_t * lrperm = lgp_local_part(int64_t, rperm);
  int64_t * lcperm = lgp_local_part(int64_t, cperm);

  int64_t items_processed = 0;
  uint64_t type_mask = 0x8000000000000000;
  int64_t start, end;
  int64_t * lrowqueue  = calloc(lnr, sizeof(int64_t));
  int64_t * lcolqueue  = calloc(lnc, sizeof(int64_t));
  int64_t * lcolqueue_level  = calloc(lnc, sizeof(int64_t));
  int64_t * lrowsum    = calloc(lnr, sizeof(int64_t));
  int64_t * lrowcnt    = calloc(lnr, sizeof(int64_t));
  int64_t * level      = calloc(lnr, sizeof(int64_t));
  int64_t * matched_col= calloc(lnr, sizeof(int64_t));
  pkg_topo_t pkg;

  exstack_t * ex = exstack_init(buf_cnt, sizeof(pkg_topo_t));
  if( ex == NULL ){return(-1.0);}
  
  /* initialize rowsum, rowcnt, and queue (queue holds degree one rows) */
  int64_t rownext, rowlast;
  int64_t colnext, collast;
  int64_t colstart,colend, col_level;
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
  
  SHARED int64_t * pivots = lgp_all_alloc(THREADS, sizeof(int64_t));
  int64_t * lpivots = lgp_local_part(int64_t, pivots); 
  lpivots[0] = 0L;
  //lgp_put_int64(pivots, 0, 0);

  lgp_barrier();

  double t1 = wall_seconds();
  
  int64_t num_levels = 0;
  while(exstack_proceed(ex, (items_processed == (lnr+lnc)))){
    
    int full = 0;
    
    /* push the column for each degree one row */
    while(rownext < rowlast){
      row = pkg.row = lrowqueue[rownext++];
      pkg.row |= type_mask;
      pkg.col = lrowsum[row];
      pkg.level = level[row];
      matched_col[row] = pkg.col;
      pe = pkg.col % THREADS;

      items_processed++;
      if(exstack_push(ex, &pkg, pe) == 1){
        full = 1;
        break;
      }
    }
    
    if(!full){
      while(colnext <= collast){
        if(colstart == colend){
          if(colnext == collast) break;
          items_processed++;
          curr_col = lcolqueue[colnext];
          col_level = lcolqueue_level[colnext];
          colnext++;
          colstart = tmat->loffset[curr_col];
          colend = tmat->loffset[curr_col + 1];          
        }
        row = tmat->lnonzero[colstart];
        pkg.row = row/THREADS;
        pkg.col = curr_col*THREADS + MYTHREAD;
        pkg.level = col_level;
        pe = row % THREADS;
        if(exstack_push(ex, &pkg, pe) == 0)
          break;
        colstart++;
      
      }
    }
    
    exstack_exchange(ex);
    
    while(exstack_pop(ex, &pkg, &fromth)){
      if(pkg.row & type_mask){
        lcolqueue[collast] = pkg.col/THREADS;
        lcolqueue_level[collast] = pkg.level;
        collast++;
      }else{
        lrowsum[pkg.row] -= pkg.col;
        lrowcnt[pkg.row]--;
        /* update the level for this row */
        if(pkg.level >= level[pkg.row]){
          level[pkg.row] = pkg.level + 1;
          if((pkg.level+1) > num_levels)
            num_levels = pkg.level + 1;
        }
        if(lrowcnt[pkg.row] == 1){
          //printf("found a d1 row with level %ld\n", level[pkg.row]);
          lrowqueue[rowlast++] = pkg.row;
        }
      }
    }
  }

  num_levels++;
  /* at this point we know for each row its level and the column it was matched with.
     we need to create cperm and rperm from this information */
  num_levels = lgp_reduce_max_l(num_levels);

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

  lgp_barrier();

  exstack_clear(ex);
  ex = exstack_init(buf_cnt, sizeof(pkg_cperm_t));
  pkg_cperm_t pkg2;
  
  /* create cperm */
  i = 0;
  while(exstack_proceed( ex, (i==lnr) )) {
    for( ; i < lnr; i++){
      pkg2.pos = lrperm[i];
      pkg2.col = matched_col[i];
      pe  = pkg2.col % THREADS;
      if(exstack_push(ex, &pkg2, pe) == 0)
        break;
    }

    exstack_exchange(ex);

    while(exstack_pop(ex, &pkg2, &fromth)){
      lcperm[pkg2.col/THREADS] = pkg2.pos;
    }
  }
  
  lgp_barrier();
  
  t1 = wall_seconds() - t1;
  minavgmaxD_t stat[1];
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  exstack_clear(ex);

  free(lrowcnt);
  free(lrowsum);
  free(lrowqueue);
  free(lcolqueue);
  free(level_start);
  free(level_sizes);
  T0_fprintf(stderr, "num levels = %ld ", num_levels);
  return(stat->avg);
  
}

