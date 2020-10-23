/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

#include "spmat_utils.h"

/*! \file sssp_dijsktra.c
\brief These are the function that implement and support two implementations of Dijsktra's Alg.
*/

/*!
\brief Implementation of naive version of Dijkstra's algorithm
\param weight array report the weight of the lightest path to the node
\param mat sparsemat_t that holds the graph. 
\param v0 is the starting vertex
\return runtime

This implementation uses a simple array to hold the
tentative weights and uses tricks with the weights to flag
the transition of vertices from unvisited to tentative to resolved.
It also does a linear search to find the smallest tentative weight.
 */
double sssp_dijsktra_linear(d_array_t *weight, sparsemat_t * mat, int64_t v0)
{
  double tm = wall_seconds();
  int64_t i, k, rn;
  int64_t numrows = mat->numrows;
  double minwt, curwt;
  int64_t minidx;

  set_d_array(weight, INFINITY);
  double *tent = weight->entry;
  for(k = mat->offset[v0];  k < mat->offset[v0+1]; k++)
    tent[ mat->nonzero[k] ] = mat->value[k];
  tent[v0] = -0.0;                          // trick: -0.0 is different than 0.0
  //printf("tent[v0] = %lg\n", tent[v0]);

  while( 1 ) { // find the smallest tentative distance
    minwt = INFINITY;
    minidx = numrows;
    for(i=0; i<numrows; i++){
      if(tent[i] > 0 && tent[i] < minwt){
        minwt = tent[i];
        minidx = i;
      }
    }
    // done if all connected vertices have been resolved
    if( (minidx == numrows) || (minwt == INFINITY) ) 
      break;

    // update tentative distance from the current vertex
    curwt = tent[minidx];
    for(k = mat->offset[minidx];  k < mat->offset[minidx+1]; k++){
      rn = mat->nonzero[k];
      if(tent[rn] > curwt + mat->value[k])
        tent[rn] = curwt + mat->value[k];
    }
    //printf("tent[%"PRId64"] = %lg\n", minidx, curwt);
    tent[minidx] = -tent[minidx];
  }

  for(i=0; i<numrows; i++){
    if(tent[i] != INFINITY)
      tent[i] = -tent[i];
  }
  return(wall_seconds() - tm);
}

// TODO move to README   rename node to heap or hidx 
/* The min-heap implementation of Dijkstra's alg.
 * 
 * In the heap implementation the tentative weights are stored in
 * a min-heap data structure that is coupled with the rows of the matrix.
 * We will refer to the entries in the heap as nodes and the vertices 
 * in the graph as rows (of the matrix).  
 * The coupling is maintained by the arrays row[] and node[]. 
 * row[] and node[] are inverses of each other. 
 * row[n] gives a way look up the row when the heap node n reaches the root of the heap
 * node[r] gives a way look up the heap node that currently handles row r's tentative weight
 */

/*!
\brief holds the simple binary heap that we use as our priority queue to find the lightest tent[v]
*/
typedef struct PQ_t {
  int64_t numrows; //!< the number of vertices in the graph, rows in the matrix
  int64_t tail;    //!< the first node not actively in the queue;
  double  *wt;     //!< the corresponding weight that is being prioritized
  int64_t *row;    //!< the index of the row being represented by the given node in the heap
  int64_t *node;   //!< the index into the heap where the weight for the given row is stored
}PQ_t;

PQ_t * init_pqueue(int64_t numrows); //!< allocate and initialize the heap
void bubble_up(PQ_t *pq, int64_t k); //!< starting a position k in the heap swap children for a parent to heapify the element before position k
void delete_root(PQ_t * pq);         //!< pop the root element from the heap

void  heapify_pqueue(PQ_t * pq);     //!< starting with the heap filled with arbitrary values, repeatedly bubble_up element until the heap property
void print_queue(PQ_t * pq);         //!< print the elements in the active part of the heap
int64_t check_pqueue(PQ_t * pq);     //!< check that the heap property holds

/*!
\brief initialize the heap
\param numrows the number in the matrix (nodes in the graph)

Note the indexing into the queue will be "1-up" so the parent and children are easy to calculate.
A given parent, p, has kids 2*p and 2*p+1, and the parent of kid, k, is \f$\lfloor k/2 \rfloor\f$.
As a convenience, we will use the zero node in the heap for the initial row (starting vertex).
We allocate an extra node so that numrows is a safe index.
This allows us to use certain positions in the heap as flags to the rows: 
  node[row] == numrows means the row is unvisited
  node[row] == 0 means the row is resolved
  otherwise node[row] is the index of the entry in the heap that holds the tentative weight of row.
*/  
PQ_t * init_pqueue(int64_t numrows)
{
  PQ_t * pq = calloc(1, sizeof(PQ_t));
  pq->numrows = numrows;
  pq->tail = 1;
  pq->wt  = (double  *)calloc(numrows+1, sizeof(double));
  pq->row  = (int64_t *)calloc(numrows+1, sizeof(int64_t));
  pq->node = (int64_t *)calloc(numrows+1, sizeof(int64_t));
  return(pq);
}

/*!
\brief bubble up a given node to return the heap to a legal state
\param pq the priority queue
\param nd the index of the node that has changed

Since a tentative weight only changes if it gets smaller
and we change only one node at a time, it is enough 
to bubble up the changed node until it is not smaller than its parent.
*/  
void bubble_up(PQ_t *pq, int64_t nd)
{
   double w;
   int64_t kid_row, par_row; 
   while( nd > 1 ){
     if( pq->wt[ nd/2 ] <= pq->wt[nd] )
       return; 
     //printf("swap kid %"PRId64" and parent %"PRId64"\n", nd/2, nd);
     w = pq->wt[nd];
     pq->wt[nd] = pq->wt[nd/2];
     pq->wt[nd/2] = w;
     
     kid_row = pq->row[nd];
     par_row = pq->row[nd/2];
     assert( kid_row != pq->numrows); 
     assert( par_row != pq->numrows); 
     pq->row[nd] = par_row;
     pq->row[nd/2] = kid_row;
     
     pq->node[kid_row] = nd/2;
     pq->node[par_row] = nd;
     
     nd = nd/2;
  }    
  return;
}

/*!
\brief remove the root node of the heap 
\param pq the priority queue
\return 0 if the heap is empty, 1 otherwise

Replace it with the last node in the heap (making the heap one node shorter).
Then restore the heap property by bubbling the root node down until
it is not bigger than either of its children.
*/  
void delete_root(PQ_t * pq)
{
  double w;
  int64_t kid_row, par_row; 
  int64_t par_nd, kid_nd;

  pq->tail--;                 // now the index of the last active node in the heap
  if( pq->tail == 1)          // the heap is now empty
    return;
  pq->node[pq->row[1]] = 0;           // mark this row as resolved
  pq->row[1] = pq->row[ pq->tail];
  pq->wt[1] = pq->wt[ pq->tail];
  pq->wt[pq->tail] = 99.0;            // pq->tail is now available
  pq->row[pq->tail] = pq->numrows;
  //printf("delete root:\n");
  //print_queue(pq);

  // recursively, if the parent node is not less than both it's children, then swap parent with smaller child
  // note: if right child is not in the heap (ie. pq->tail is the right child), 
  // then that nodes exists (even though it is not active) and its wt is INFINITY.  
  // We don't have to determine and make a special case if there is only a left child.
  par_nd = 1;
  //printf("bubble down: %lg : %lg, %lg\n", pq->wt[par_nd], pq->wt[2*par_nd], pq->wt[2*par_nd + 1]);
  while( (2*par_nd + 1 <= pq->tail) && 
         ((pq->wt[ 2*par_nd ] <  pq->wt[par_nd]) || (pq->wt[ 2*par_nd + 1 ] <  pq->wt[par_nd]) ) ) {
    if( pq->wt[2*par_nd] < pq->wt[2*par_nd + 1] ){
      kid_nd = 2*par_nd;
    } else {
      kid_nd = 2*par_nd + 1;
    }
    //printf("bubble down: swap parent  %"PRId64" with kid %"PRId64"\n", par_nd, kid_nd);
    w = pq->wt[par_nd];
    pq->wt[par_nd] = pq->wt[kid_nd];
    pq->wt[kid_nd] = w;
    
    par_row = pq->row[par_nd];
    kid_row = pq->row[kid_nd];
    pq->row[par_nd] = kid_row;
    pq->row[kid_nd] = par_row;

    pq->node[par_row] = kid_nd;
    pq->node[kid_row] = par_nd;
   
    par_nd = kid_nd;
  }
  //print_queue(pq);
}


/*! \brief DRPT control the printing of debugging printfs */
#define DPRT 0
/*!
\brief The implementation of Dijkstra's algorithm that uses a heap to prioritize the tentative vertices.
\param weight array report the weight of the lightest path to the node
\param mat sparsemat_t that holds the graph. 
\param r0 is the starting row (vertex)
\return runtime
*/
double sssp_dijsktra_heap(d_array_t *weight, sparsemat_t * mat, int64_t r0)
{
  double tm = wall_seconds();
  int64_t i, k;
  int64_t numrows = mat->numrows;

  double *tent = weight->entry;

  // initialize the heap 
  PQ_t * pq = init_pqueue(numrows);
  for(i=0; i<numrows+1; i++){
    pq->wt[i] = INFINITY;
    pq->row[i] = numrows;
    pq->node[i] = numrows;
  }
  pq->wt[0] = 0.0;
  pq->row[0] = r0;
  pq->node[r0] = 0;

  tent[r0] = 0.0;
  if(DPRT){printf(">>>>>>>>>>>>> tent[%"PRId64"] = %lg\n", r0, tent[r0]);}

  int64_t rn, row, nd;
  double e_wt, vn_wt;
  // start the heap with weights to vertices adjacent to r0
  pq->tail = 1;
  for(k = mat->offset[r0];  k < mat->offset[r0+1]; k++){
    row = mat->nonzero[k];
    e_wt = mat->value[k];
    nd = pq->node[row] = pq->tail;
    pq->wt[nd] = e_wt;
    pq->row[nd] = row;
    pq->tail++;
    if(DPRT){printf("explore (%2"PRId64",%2"PRId64"): new row %"PRId64" with wt %lg at node %"PRId64"\n", r0, row, row, pq->wt[nd], nd);}
    bubble_up(pq, nd);
    if(DPRT){print_queue(pq);}
  }
  if(DPRT){print_queue(pq);printf("\n");}

  while(pq->tail > 1){       // pq-tail == 1 means the queue is empty
    rn     = pq->row[1];
    vn_wt  = pq->wt[1];
    tent[rn] = vn_wt;
    if(DPRT){printf(">>>>>>>>>>>>> tent[%"PRId64"] = %lg\n", rn, vn_wt);}
    for(k = mat->offset[rn];  k < mat->offset[rn+1]; k++){
      row = mat->nonzero[k];
      e_wt  = mat->value[k];
      nd   = pq->node[row];
      if(nd == 0){                  // row is done
        if(DPRT){printf("explore (%2"PRId64",%2"PRId64"): done with row (%"PRId64")\n", rn, row, row);}
        continue;
      }
      if( nd == numrows ) {                      // row is new
        nd = pq->tail;
        pq->node[row] = nd;
        pq->wt[nd] = vn_wt + e_wt;
        pq->row[nd] = row;
        pq->tail++;
        if(DPRT){printf("explore (%2"PRId64",%2"PRId64"): new row %"PRId64" with wt %lg at node %"PRId64"\n", rn, row, row, pq->wt[nd], nd);}
        bubble_up(pq, nd);
        if(DPRT){print_queue(pq);}
      } else if(pq->wt[nd] > (vn_wt + e_wt)) {          // improved the weight
        pq->wt[nd] = vn_wt + e_wt;
        if(DPRT){printf("explore (%2"PRId64",%2"PRId64"): updated row %"PRId64" with wt %lg at node %"PRId64"\n", rn, row, row, pq->wt[nd], nd);}
        bubble_up(pq, nd);
        if(DPRT){print_queue(pq);}
      } else {
        if(DPRT){printf("explore (%2"PRId64",%2"PRId64"): no change\n", rn, row);}
      }
    }
    delete_root(pq);    // remove the root and decrement pq->tail
  };

  return(wall_seconds() - tm);
}


/*
 * These functions were useful during the development of the heap based implementation.
 * We left them here for debugging and exploring the code.
 */


#if 0 // test the heap stuff
    int64_t i;
    PQ_t * pq = init_pqueue(numrows);
    for(i=1; i<numrows; i++) {
      pq->val[i] = (double)(numrows-i);
      pq->row[i] = i;
      pq->node[i] = i;
    }
    pq->tail = numrows;
    print_queue(pq);

    heapify_pqueue(pq);
    exit(1);
#endif


// start at the top and bubble down
// going to assume that a bunch of the node are  out of order
void heapify_pqueue(PQ_t * pq)
{
  int64_t k;

  //printf("heapify:\n");
  for(k=2;  k < pq->tail; k++) {
    bubble_up(pq, k);
  }
  //print_queue(pq);
}


// check that the wt of a parent is never greater than either child
int64_t check_pqueue(PQ_t * pq)
{
  int64_t p, ret_ok=1;
  for(p=1; p <= (pq->tail)/2; p++){
    if(  ((2*p) < pq->tail && pq->wt[p] > pq->wt[(2*p)]) 
       ||((2*p+1) < pq->tail &&  pq->wt[p] > pq->wt[(2*p+1)] ) ) {
      ret_ok = 0;
      break;
    }
  }
  return(ret_ok);
}

// just dump the entries in the heap it three consecutive lines.
// numrows needs to be small
void print_queue(PQ_t * pq)
{
  int i;
  printf("wt:   "); 
  for(i=0; i<pq->numrows+1; i++){
    printf("%3lg ", pq->wt[i]);
  }
  printf("\nrow:  "); 
  for(i=0; i<pq->numrows+1; i++){
    printf("%3"PRId64" ", pq->row[i]);
  }
  printf("\nnode: ");
  for(i=0; i<pq->numrows+1; i++){
    printf("%3"PRId64" ", pq->node[i]);
  }
  printf("\n\n");
}
  
