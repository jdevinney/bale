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
/*! \file sssp_delta.c
 * \brief Serial implementation of delta stepping for SSSP
 */
#include "spmat_utils.h"

/*! \file sssp_delta.c
\brief An implementation of the delta-stepping algorithm
*/

/*!  \brief a struct to hold (pointers to) all the state needed for the Delta-Stepping algorithm

The storage required 6 times the number of nodes plus a little for space for B.
Note we need a doubly link list because we delete nodes from any where on the list.
*/
typedef struct ds_t {
  int64_t *next;       //!< array of next "pointers" for the linked list storing a bucket
  int64_t *prev;       //!< array of prev "pointers" for the linked list storing a bucket
  int64_t *in_bucket;  //!< record which bucket a node is in, -1 if it is in no bucket
  int64_t *B;          //!< array of Buckets: B[i] is the index of the first node on the list, or -1 if empty
  int64_t num_buckets; //!< number of buckets needed to cover all possible live buckets
  double  *tent;       //!< the tentative weight for the vertex
  int64_t *deleted;    //!< deleted means resolved?
  int64_t *R;          //!< queue to hold tail vertices that need to relax their heavy edges
  double  delta;       //!< the algorithmic parameter delta
} ds_t;

/*! \brief debugging function that dumps a bucket
\param ds points to main struct
\param i_m the given bucket
*/
static void dump_bucket(ds_t *ds, int64_t i_m)
{
  int64_t v, w ;
  printf("Bucket[%ld] =", i_m);
  if( ds->B[i_m] == -1 ){
    printf("empty \n");
    return;
  }
  v = w = ds->B[i_m];
  do {
    printf("%ld ", w);
  } while( (w = ds->next[w]) != v );
  printf("\n");
  return;
}

/*! \brief Prepend node v into a bucket i_m
\param ds the main struct
\param v the given node
\param i_m the desired bucket
*/
void insert_node_in_bucket(ds_t *ds, int64_t v, int64_t i_m)
{
  int64_t w;                   // node on list, given by ds->B[i_m]
  Dprintf("Adding %"PRId64" to bucket %"PRId64" of %"PRId64"\n", v, i_m, ds->num_buckets);
  
  assert(i_m >= -1 && i_m < ds->num_buckets);
  
  if(ds->in_bucket[v] == i_m){      // it is ok if this node is already in this bucket
    return; 
  }
  assert(ds->in_bucket[v] == -1);   // better not be in a different bucket

  ds->next[v] = v;                  // buckets are circular lists, so a single element points to itself
  ds->prev[v] = v;
  if(ds->B[i_m] != -1){    
    w = ds->B[i_m];                 // w is "first" on the list, insert v before w              
    Dprintf("non-empty: w=%ld, prev=%ld\n", w, ds->prev[w]);
    ds->prev[v] = ds->prev[w];            
    ds->next[ds->prev[w]] = v;            
    ds->prev[w] = v;
    ds->next[v] = w;
    Dprintf("   v=%ld, prev=%ld, next %ld, (%ld,%ld)\n", v, ds->prev[v], ds->next[v], ds->prev[ds->next[v]], ds->next[ds->prev[v]] );
  }
  ds->B[i_m] = v;                       // set v to be the "first" on the list
  ds->in_bucket[v] = i_m;
}


/*! \brief Remove node v from its bucket (need not be in any bucket)
\param ds the main struct
\param v the given node

It eases the logic and causes no harm to "remove" a node 
from its bucket and find out here that it is actually not in a bucket.
*/
static void remove_node_from_bucket(ds_t *ds, int64_t v)
{
  int64_t i_m, w;

  i_m = ds->in_bucket[v];
  assert(i_m >= -1 && i_m < ds->num_buckets);
  Dprintf("Removing %"PRId64" from bucket %"PRId64"\n", v, i_m);

  if(i_m == -1){
    Dprintf("   %"PRId64" was not in bucket %"PRId64"\n", v, i_m);
    return;
  }

  ds->in_bucket[v] = -1;
  if((ds->next[v] == v) && (ds->prev[v] == v)){  // the only thing on the list
    ds->B[i_m] = -1;
    return;
  }

  // move the pointers to remove v and make the next thing in the bucket the new "first" thing
  w = ds->next[v];  
  ds->prev[w] = ds->prev[v];
  ds->next[ds->prev[v]] = w;
  ds->B[i_m] = w;
  return;
}

/*!  \brief relax the edge (possible change the tentative weight of the head of a edge)
\param ds the main struct
\param w the head of the given edge
\param cand_wt the new weight = tent[v] + c(v,w) where v is the tail of the edge

We don't actually care who is the tail of the edge.
*/
void relax_edge(ds_t *ds, int64_t w, double cand_wt)
{
  int64_t iold, inew;
  Dprintf("relax head %"PRId64" cand_wt = %lf < %lf?\n", w, cand_wt, ds->tent[w]);
  if ( cand_wt < ds->tent[w] ){
    iold = ds->in_bucket[w];
    inew = ((int64_t)floor(cand_wt/ds->delta)) % (ds->num_buckets);

    assert(iold >= -1 && iold < ds->num_buckets);
    assert(inew >= -1 && inew < ds->num_buckets);
    Dprintf("winner: %"PRId64"  move from bucket %"PRId64" to %"PRId64"\n", w, iold, inew);
    if( iold != inew ){
      if(iold >= 0)
        remove_node_from_bucket(ds, w);
      insert_node_in_bucket(ds, w, inew);
    }
    ds->tent[w] = cand_wt;
  }
}

/*!  \brief the array based implementation of the delta stepping algorithm
\param weight pointer to the array that holds the results
\param mat the graph (matrix)
\param r0 the given single source node
\param del the algorithm parameter delta
*/
double sssp_delta_stepping(d_array_t *weight, sparsemat_t * mat, int64_t r0, double del)
{
  int64_t i, i_m, k;
  int64_t v;

  double tm = wall_seconds();

  assert((r0 >= 0) && (r0<mat->numrows));
  
  
  /* calculate delta and set tentative distances to infinity */  
  double delta = 0.0;
  int64_t max_degree = 0;
  for(i = 0; i < mat->numrows; i++){
    if(max_degree < (mat->offset[i+1] - mat->offset[i]))
      max_degree = (mat->offset[i+1] - mat->offset[i]);
  }
  assert(max_degree > 0);
  delta = 1.0/max_degree;
  
  double max_edge_weight = 0.0;
  for(i = 0; i < mat->nnz; i++)
    if(max_edge_weight < mat->value[i])
      max_edge_weight = mat->value[i];
  Dprintf("max edge weight = %lf\n", max_edge_weight);
  int64_t num_buckets = (int64_t)ceil(max_edge_weight/delta) + 1;
  
  // allocate and initialize the all the arrays
  ds_t * ds = (ds_t *)calloc(1,sizeof(ds_t)); assert(ds != NULL);
  ds->next = (int64_t *)malloc(mat->numrows * sizeof(int64_t)); assert(ds->next != NULL);
  ds->prev = (int64_t *)malloc(mat->numrows * sizeof(int64_t)); assert(ds->prev != NULL);
  ds->in_bucket = (int64_t *)malloc(mat->numrows * sizeof(int64_t)); assert(ds->in_bucket != NULL);
  ds->R = malloc(mat->numrows * sizeof(int64_t)); assert(ds->R != NULL);
  ds->deleted = (int64_t *)malloc(mat->numrows * sizeof(int64_t)); assert(ds->deleted != NULL);
  ds->tent = (double *)malloc(mat->numrows * sizeof(double)); assert(ds->tent != NULL);
  for(i = 0; i < mat->numrows; i++){
    ds->next[i]      = i;
    ds->prev[i]      = i;
    ds->in_bucket[i] = -1;
    ds->deleted[i]   =  0;
    ds->R[i]         =  0;
    ds->tent[i]      = INFINITY;
  }
  Dprintf("Allocate buckets\n");
  ds->num_buckets = num_buckets;
  ds->B = (int64_t *)calloc(num_buckets, sizeof(int64_t)); assert(ds->B != NULL);
  // bucket indices are mod num_buckets.
  for(i = 0; i < num_buckets; i++){
    ds->B[i] = -1;
  }
  ds->delta = delta;

  int64_t light_nnz=0;
  int64_t heavy_nnz=0;
  for(i = 0; i < mat->nnz; i++){
    if(mat->value[i] <= delta)
      light_nnz++;
  }
  heavy_nnz = mat->nnz - light_nnz;
  sparsemat_t *light = init_matrix(mat->numrows, mat->numcols, light_nnz, 1);
  sparsemat_t *heavy = init_matrix(mat->numrows, mat->numcols, heavy_nnz, 1);
  light->offset[0] = 0;
  heavy->offset[0] = 0;
  int64_t lpos=0, hpos=0;
  light->nnz = light_nnz;
  heavy->nnz = heavy_nnz;
  Dprintf("nnz %ld =  %ld + %ld (light + heavy)\n",mat->nnz, light->nnz, heavy->nnz);
  for(i=0; i<mat->numrows; i++){
    for(k=mat->offset[i]; k<mat->offset[i+1]; k++){
      if(mat->value[k] <= delta){
        light->nonzero[lpos] = mat->nonzero[k];
        light->value[lpos] = mat->value[k];
        lpos++;
      } else {
        heavy->nonzero[hpos] = mat->nonzero[k];
        heavy->value[hpos] = mat->value[k];
        hpos++;
      }
    }
    light->offset[i+1] = lpos;
    heavy->offset[i+1] = hpos;
  }
  printf("nnz %ld =  %ld + %ld (light + heavy)\n", k, lpos, hpos);

  // set the distance to r0 equal to 0.0
  Dprintf("Set source node %ld\n", r0);
  ds->tent[r0] = 0.0;
  insert_node_in_bucket(ds, r0, 0);
  Dprintf("putting source r0=%ld in bucket %ld\n", r0, 0L);
  
  int64_t rmb = 0;      // the real min bucket (smallest non-empty bucket)
  int64_t start = 0;
  int64_t end = 0;
  /* main loop */
  // MAIN LOOP: exists when all buckets are empty
  while(1){
    Dprintf("Outer loop rmb = %"PRId64"\n", rmb);
    // find the minimum indexed non-empty bucket, given that we use them mod num_buckets 
    for(i = 0; i < ds->num_buckets; i++){
      if(ds->B[(rmb + i) % ds->num_buckets] != -1)
        break;
    }
    if(i == num_buckets)                //All buckets are empty, we are done
      break;
    rmb = rmb + i;
    i_m = rmb % ds->num_buckets;
    Dprintf("Starting inner loop: working on bucket %"PRId64" in %"PRId64"\n", rmb, i_m);
    if(DEBUG) dump_bucket(ds, i_m);

    start = 0;
    end = 0;
    while( ds->B[i_m] >= 0 ) {     // removes vertices from this bucket possibly adding more vertices to the bucket

      v = ds->B[i_m];
      Dprintf("Processing Node %"PRId64" in Bucket %"PRId64"\n", v, i_m);
      remove_node_from_bucket(ds, v);
      for(k = light->offset[v]; k < light->offset[v + 1]; k++){   // relax light edges from v 
          relax_edge(ds, light->nonzero[k], ds->tent[v] + light->value[k]);
      } 
      
      if(ds->deleted[v] == 0){   // insert v into R if it is not already there/
        ds->deleted[v] = 1;
        ds->R[end++] = v;
        Dprintf("deleted %"PRId64"\n", v);
      }
    }

    for(start=0; start<end; start++){ // relax heavy requests edges for everything in R
      v = ds->R[start];
      for(k = heavy->offset[v]; k < heavy->offset[v + 1]; k++){
          relax_edge(ds, heavy->nonzero[k], ds->tent[v] + heavy->value[k]);
      }
    }
  }

  for(v = 0; v < mat->numrows; v++){          // Copy the answers to dist array
    weight->entry[v] = ds->tent[v];
  }

  free(ds->next);
  free(ds->prev);
  free(ds->in_bucket);
  free(ds->deleted);
  free(ds->tent);
  free(ds->R);
  free(ds);
  
  return(wall_seconds() - tm);
}
