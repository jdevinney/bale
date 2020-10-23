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

/*! \file toposort_exstack2.upc
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 * using the exstack2 buffering model.
 */
#include "toposort.h"

/*!
 * \brief This routine implements the exstack2 variant of toposort
 * \param *rperm place to return the row permutation that is found
 * \param *cperm place to return the column permutation that is found
 * \param *mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
 * \param *tmat the transpose of mat
 * \param buf_cnt number of package in the exstack2 buffers
 * \return average run time
 */
double toposort_matrix_exstack2(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat, int64_t buf_cnt) {
  typedef struct pkg_rowcol_t{
    int64_t row;    
    int64_t col;
  }pkg_rowcol_t;
  //T0_printf("Running Toposort with exstack2 ...");
  int64_t col2;
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
  int64_t lnr = (nr + THREADS - MYTHREAD - 1)/THREADS;
  int64_t lnc = (nc + THREADS - MYTHREAD - 1)/THREADS;
  int64_t * lrowqueue  = calloc(lnr, sizeof(int64_t));
  int64_t * lrowsum    = calloc(lnr, sizeof(int64_t));
  int64_t * lrowcnt    = calloc(lnr, sizeof(int64_t));
  int64_t * lpivots     = calloc(lnr, sizeof(int64_t));
  int64_t * lrperm = lgp_local_part(int64_t, rperm);
  int64_t * lcperm = lgp_local_part(int64_t, cperm);
  
  int64_t i,j, start, end;

//X#if __UPC_ATOMIC__
//X  upc_atomicdomain_t * domain = upc_all_atomicdomain_alloc(UPC_INT64, UPC_INC, 0);
//X#else
//X  void * domain;
//X#endif

  /* initialize lrowsum, lrowcnt, and lrowqueue (lrowqueue holds degree one rows) */
  start = end = 0;
  for(i = 0; i < mat->lnumrows; i++){
    lrowsum[i] = 0L;
    lrowcnt[i] = mat->loffset[i+1] - mat->loffset[i];
    if(lrowcnt[i] == 1)
      lrowqueue[end++] = i;
    for(j = mat->loffset[i]; j < mat->loffset[i+1]; j++)
      lrowsum[i] += mat->lnonzero[j];
  }
  
  //SHARED volatile int64_t * pivots = lgp_all_alloc(THREADS, sizeof(int64_t));
  //pivots[MYTHREAD] = 0L;

  lgp_barrier();
  
  // we a pick a row with a single nonzero = col.
  // setting tri_rperm[pos] = row and tri_cprem[pos] = col
  // moves that nonzero to the diagonal.
  // Now, cross out that row and col by decrementing 
  //  the lrowcnt for any row that contains that column
  // repeat

  int64_t pe, fromth, pos, row, col, l_col, thislevel, total_rows = 0;
  int64_t lrows_this_level = end - start;
  int64_t rows_this_level = lgp_reduce_add_l(lrows_this_level);
  int64_t perm_pos = nr - 1 - lgp_prior_add_l(lrows_this_level);
  int64_t rows_per_level = 0L;

  exstack2_t * ex3 =  exstack2_init(buf_cnt, sizeof(int64_t));
  pkg_rowcol_t pkg_nz;
  exstack2_t * ex4 = exstack2_init(buf_cnt, sizeof(pkg_rowcol_t));

  // Do it in levels:
  // Per level -- everyone picks their degree one rows
  //              and claims their spots in the permutations,
  //              then they all go back and update the 
  //              appropriate lrowcnts for that level
  //                 
  // cool trick: lrowsum[i] is the sum of all the column indices in row i,
  //             so lrowsum[i] is the column indice we want when lrowcnt = 1;
  int64_t curr_k=0, curr_end=0, curr_col;
  int64_t level = 0;
  double t1 = wall_seconds();
  
  while(rows_this_level){

    thislevel = end;
    rows_per_level += thislevel-start;
    //T0_printf("level -- %ld =%ld\n", level, thislevel-start);
    
    while(exstack2_proceed(ex3, (start==thislevel)) || exstack2_proceed(ex4, ex3->all_done)){

      /* try to send out the column from a degree one row in your lrowqueue */
      while(start < thislevel){
        row = lrowqueue[start];
        col = lrowsum[row];  // see cool trick
        pe = col % THREADS;
        // push this col out to the owner of this row in the transpose matrix
        // looking at the columns in the transpose matrix will give us the 
        // other rows that contain this column
        if( !exstack2_push(ex3, &col, pe) )
          break;
        start++;
        //pos = my_fetch_and_inc(&pivots[0], domain);
        //rperm[row*THREADS + MYTHREAD] = pos;
        //cperm[col] = pos;
        lrperm[row] = perm_pos--;
        lpivots[row] = col;
      }
      
      /* see if you have any columns to process */
      do {
        if(curr_k == curr_end ){
          if(exstack2_pop(ex3, &col2, &fromth)) {
             curr_col = col2/THREADS;
             curr_k   =  tmat->loffset[curr_col];
             curr_end =  tmat->loffset[curr_col+1];
          } else 
            break;
        }
        while( curr_k < curr_end ) { 
          pkg_nz.row = tmat->lnonzero[curr_k] / THREADS;
          pe = tmat->lnonzero[curr_k] % THREADS;
          pkg_nz.col = curr_col * THREADS + MYTHREAD;
          if( !exstack2_push(ex4, &pkg_nz, pe) )
            break;
          curr_k++;
        }
      }while(curr_k == curr_end);
      
      // pop the rows and decrement lrowcnt
      while(exstack2_pop(ex4, &pkg_nz, &fromth)){
        lrowcnt[pkg_nz.row]--;
        lrowsum[pkg_nz.row] -= pkg_nz.col;

        if(lrowcnt[pkg_nz.row] == 1){  // became 1 in level, will used in the next level
          lrowqueue[end++] = pkg_nz.row;
        }
      }
    }
    
    total_rows += rows_this_level;
    lrows_this_level = end - start;
    rows_this_level = lgp_reduce_add_l(lrows_this_level);
    perm_pos = (nr - 1 - total_rows) - lgp_prior_add_l(lrows_this_level);

    exstack2_reset(ex3);
    exstack2_reset(ex4);
    level++;
  }
  
  level = lgp_reduce_max_l(level);
  lgp_barrier();
  exstack2_reset(ex4);

  /* create cperm */
  i = 0;
  while(exstack2_proceed( ex4, (i==lnr) )) {
    for( ; i < lnr; i++){
      pkg_nz.row = lrperm[i];
      pkg_nz.col = lpivots[i];
      pe  = pkg_nz.col % THREADS;
      if(exstack2_push(ex4, &pkg_nz, pe) == 0)
        break;
    }
    while(exstack2_pop(ex4, &pkg_nz, &fromth)){
      lcperm[pkg_nz.col/THREADS] = pkg_nz.row;
      //printf("cperm[%ld] = %ld\n", pkg_nz.col, pkg_nz.row);
    }
  }
  
  lgp_barrier();
  
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );

  //rows_per_level /= level;
  //lgp_min_avg_max_d( stat, rows_per_level, THREADS );
  //T0_fprintf(stderr,"rows_per_level : %'16ld %'16ld %'16ld\n", stat->min, stat->avg, stat->max);

  exstack2_clear(ex3);
  exstack2_clear(ex4);
  free(ex3);
  free(ex4);
  
  free(lpivots);
  free(lrowqueue);
  free(lrowsum);
  free(lrowcnt);
  T0_fprintf(stderr, "num levels = %ld ", level+1);
  //T0_printf( "done\n");
  return(stat->avg);
}

