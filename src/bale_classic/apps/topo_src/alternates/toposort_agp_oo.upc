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
/*! \file toposort_agp_oo.upc
 * \brief The intuitive implementation of toposort that does uses 
 * atomics and generic global references
 */

#include "toposort.h"

/* To play with some higher level abstractions for queues and 
 * moving across the rows of a sparse matrix, we have written 
   toposort_matrix_upc_oo which uses the following support functions.
 * the toposort_matrix_upc_orig  uses more the explicit code.
 * the main toposort functions just calls toposort_matrix_upc,
   hence the macros above.
 */

/*!
\brief A structure to hold a queue that can be filled by other threads and emptied only by MYTHREAD

We used shared arrays for the items and both shared and a local variable to hold the 
indices for the heads of the queues. Since one dequeues items locally, we only need
a local variable for the tails.
The race condition to add items to a queue is handled with a fetch_and_add to the shared head.
The read after write race condition to dequeue items up to head is side stepped by using a 
snap shot of the queues and only updating the snap shots between barriers. 
*/  
typedef struct swlrQ {
   SHARED int64_t * S_queue;  //!< space for all the queues
   int64_t * l_queue;         //!< localized pointers to the individual queues
   SHARED int64_t *S_head;    //!< shared version of the heads of each queue (inherently suffers races conditions)
   int64_t l_head;            //!< snap shot version of the head of the local queue (this is only set by qrab_swlrQ)
   int64_t l_tail;            //!< tail of the local queue
} swlrQ_t;

/*!
\brief initialize the shared write local read queue
\param nitems the number of items on each local queue
*/
swlrQ_t * init_swlrQ(int64_t nitems) {
  swlrQ_t * sq = calloc(1, sizeof(swlrQ_t));
  sq->S_queue  = lgp_all_alloc(((nitems/THREADS+1)*THREADS), sizeof(int64_t));
  sq->l_queue  = lgp_local_part(int64_t, sq->S_queue);

  sq->S_head = lgp_all_alloc(THREADS, sizeof(int64_t));
  lgp_put_int64(sq->S_head, MYTHREAD, 0);
  sq->l_head = 0;
  sq->l_tail = 0;

  lgp_barrier();
  return(sq);
}

/*!
\brief clears the shared write local read queue
\param sq pointer to the shared queue
*/
void * clear_swlrQ( swlrQ_t * sq ) {
  if( !sq ) return(NULL);
  lgp_all_free(sq->S_queue);
  lgp_all_free(sq->S_head);
  return(NULL);
}

/*!
\brief grabs a local snap shot version of the queue by setting the "local head" of the queue
\param sq pointer to the shared queue
\param nitems place to hold the number of items on the local queue
\return the total number of items in all the threads queues
This allows one to work with local queue while other threads are use S_head to add more items to the queue.
*/
int64_t grab_swlrQ(swlrQ_t * sq, int64_t *nitems) {
   lgp_barrier();   // finish inflight pushes, updates to rowcnt and rowsum
   sq->l_head = lgp_get_int64(sq->S_head, MYTHREAD);
   *nitems = sq->l_head - sq->l_tail;
   return( lgp_reduce_add_l(*nitems) ); // implicit barrier here is necessary
}

/*!
\brief pushs an item (int64_t) on to another threads queue
\param sq pointer to the shared queue
\param owner the target thread for the item
\param item the int64_t we are pushing
\return true (might return false if we did some error checking)
*/
bool en_swlrQ(swlrQ_t * sq, int64_t owner, int64_t item) {
   int64_t l_pos;
   
   // the race for the spot on the another threads queue
   // is handled with an atomic fetch_and_add to its head
   l_pos = lgp_fetch_and_inc(sq->S_head, owner);
   lgp_put_int64(sq->S_queue, l_pos*THREADS + owner, item);
   //printf("%d: >> %d to %d into %d\n",MYTHREAD, item, owner, l_pos);
   return true;
}

/*!
\brief pull the tail element from a thread's queue.
\param sq pointer to the shared write local read queue
\param ret_item a place to put the int64_t item from the queue
\return false if the queue was empty, else true
*/
bool de_swlrQ(swlrQ_t * sq, int64_t *ret_item ) {
  if( sq->l_head > sq->l_tail ) {
     *ret_item = sq->l_queue[sq->l_tail];
     sq->l_tail++;
     //printf("%d: << %d\n",MYTHREAD, *ret_item);
     return true;
  } else {
     //printf("%d: << empty\n",MYTHREAD);
     return false;
  }
}

/*!
\brief This routine implements an AGP variant of toposort with code to encapsulate working 
  with the queue of degree one rows and moving across the rows a transpose matrix.
\param rperm place to return the row permutation that is found
\param cperm place to return the column permutation that is found
\param mat the input sparse matrix NB. it must be a permuted upper triangular matrix 
\param tmat the transpose of mat
\return average run time
*/
double toposort_matrix_oo(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat) {
  //T0_fprintf(stderr,"Running Toposort with UPC ....");
  int64_t l_row, S_col, S_row;
  int64_t old_row_sum, old_row_cnt;
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;

  swlrQ_t * sq = init_swlrQ(nr);

  int64_t *l_rperm = lgp_local_part(int64_t, rperm);
  SHARED int64_t * rowsum = lgp_all_alloc(nr+THREADS, sizeof(int64_t));
  int64_t *lrowsum = lgp_local_part(int64_t, rowsum);
  SHARED int64_t * rowcnt = lgp_all_alloc(nr+THREADS, sizeof(int64_t));
  int64_t *lrowcnt = lgp_local_part(int64_t, rowcnt);
  nxnz_t *lnxz  = init_nxnz( mat );
  nxnz_t *Snxzt = init_nxnz( tmat );

  double t1 = wall_seconds();
   
  /* initialize rowsum, rowcnt, and queue holding the degree one rows */
  for(int64_t i = 0; i < mat->lnumrows; i++){
    lrowcnt[i] = rowcount_l(mat, i);
    if(lrowcnt[i] == 1)
      en_swlrQ(sq, MYTHREAD, i);
    lrowsum[i] = 0L;
    for(first_l_nxnz(lnxz,i); has_l_nxnz(lnxz,i); incr_l_nxnz(lnxz,i))
      lrowsum[i] += lnxz->col;
  }

  int64_t mypivs, n_pivs; // number of pivots on my queue, all queues (done when n_pivs == 0)
  int64_t pivs_found = 0; // total pivots found so far
  int64_t pos;
  while( (n_pivs = grab_swlrQ(sq, &mypivs)) != 0 ) {  // getting a local copy of the queue defines a level

    pos = pivs_found + lgp_prior_add_l(mypivs); // pivots in this level and thread start here
    while( de_swlrQ(sq, &l_row) ) {
      S_col = lrowsum[l_row];  // see cool trick

      l_rperm[l_row]  = nr - 1 - pos;
      lgp_put_int64(cperm, S_col, nc - 1 - pos);
      pos++;

      // for every S_row in the S_col'th row of tmat
      for(first_S_nxnz(Snxzt,S_col); has_S_nxnz(Snxzt,S_col); incr_S_nxnz(Snxzt,S_col)){
        S_row = Snxzt->col;
        old_row_cnt = lgp_fetch_and_add(rowcnt, S_row, -1L);       // race condition between these two 
        old_row_sum = lgp_fetch_and_add(rowsum, S_row, (-1L)*S_col); // will be resolved before the next grab
        if( old_row_cnt == 2L ) {
          en_swlrQ(sq, S_row%THREADS, S_row/THREADS);
        }
      }
    }
    pivs_found += n_pivs;
  }
  
  minavgmaxD_t stat[1];
  t1 = wall_seconds() - t1;
  lgp_min_avg_max_d( stat, t1, THREADS );
  
  if(pivs_found != nr){
    printf("ERROR! toposort_matrix_upc_oo: found %"PRId64" pivots but expected %"PRId64"!\n", pivs_found, nr);
    exit(1);
  }

  clear_swlrQ(sq); free(sq);
  free(lnxz);
  free(Snxzt);
  lgp_all_free(rowsum);
  lgp_all_free(rowcnt);
  return(stat->avg);
}

