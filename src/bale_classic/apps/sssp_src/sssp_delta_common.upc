/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//
//  All rights reserved.
//  
//  This file is a part of Bale.  For license information see the
//  LICENSE file in the top level directory of the distribution.
// 
 *****************************************************************/

/*! \file sssp_delta_common.upc
 * \brief Routines common to all the delta stepping implementations.
 */

#include "sssp.h"
#include "sssp_delta_common.h"


/*!
\brief A debugging routine that dumps the contents of a bucket
\param ds the data structure that hold or points to everything
\param i_m the ith bucket (mod num_buckets)
*/
void dump_bucket(ds_t *ds, int64_t i_m)
{
  int64_t v, w ;
  char lineout[2048], wordout[32];

  lineout[0] = '\0';
  sprintf(lineout, "%02d: Bucket[%ld] =", MYTHREAD, i_m);
  if( ds->B[i_m] == -1 ){
    strncat(lineout, " empty", 2047);
    puts(lineout);
    return;
  }
  v = w = ds->B[i_m];
  do {
    sprintf(wordout, "%3ld ", w);
    strncat(lineout, wordout, 2047);
  } while( (w = ds->next[w]) != v );
  puts(lineout);
  return;
}


/*!
\brief Adds a vertex to the ith bucket
\param ds the data structure that holds or points to everything
\param v the node
\param i_m the ith bucket (mod num_buckets)

The prepends v to the list of nodes in the bucket
*/
void insert_node_in_bucket(ds_t *ds, int64_t v, int64_t i_m)
{
  int64_t w;                   // node on list, given by ds->B[i_m]
  //if(DPRT){printf("%02d: Adding %"PRId64" to bucket %"PRId64" of %"PRId64"\n", MYTHREAD, v, i_m, ds->num_buckets);}
  
  assert(i_m >= -1 && i_m < ds->num_buckets);
  
  if(ds->in_bucket[v] == i_m){      // it is ok if this node is already in this bucket
    return; 
  }
  assert(ds->in_bucket[v] == -1);   // better not be in a different bucket

  ds->next[v] = v;
  ds->prev[v] = v;
  if(ds->B[i_m] != -1){    
    w = ds->B[i_m];             // w is "first" on the list, insert v before w              
    //if(DPRT){printf("%02d: non-empty: w=%ld, prev=%ld\n", MYTHREAD, w, ds->prev[w]);}
    ds->prev[v] = ds->prev[w];            
    ds->next[ds->prev[w]] = v;            
    ds->prev[w] = v;
    ds->next[v] = w;
    //if(DPRT){printf("%02d:    v=%ld, prev=%ld, next %ld, (%ld,%ld)\n", MYTHREAD, v, ds->prev[v], ds->next[v], ds->prev[ds->next[v]], ds->next[ds->prev[v]] );}
  }
  ds->B[i_m] = v;                       // set v to be the "first" on the list
  ds->in_bucket[v] = i_m;
}

/*! \brief Removes a vertex from the ith bucket
 * \param ds the data structure that holds or points to everything
 * \param v the node
 * Finds v in the list and removes it (v need not be in a bucket, this avoids checking in the main loop).
 */
void remove_node_from_bucket(ds_t *ds, int64_t v)
{
  int64_t i_m, w;

  i_m = ds->in_bucket[v];
  assert(i_m >= -1 && i_m < ds->num_buckets);
  //if(DPRT){printf("%02d: Removing %"PRId64" from bucket %"PRId64"\n", MYTHREAD, v, i_m);}
  if(00 && DPRT && !( ((ds->next[v] == v) && (ds->prev[v] == v)) || (ds->next[v] != ds->prev[v]) )){
    printf("%02d: ERROR:  v, next, prev = %ld %ld %ld\n", MYTHREAD, v, ds->next[v], ds->prev[v] );
  }

  if(i_m == -1)     // v wasn't in a bucket
    return;

  ds->in_bucket[v] = -1;
  if((ds->next[v] == v) && (ds->prev[v] == v)){  // the only thing on the list
    ds->B[i_m] = -1;
    return;
  }

  w = ds->next[v];
  ds->prev[w] = ds->prev[v];
  ds->next[ds->prev[v]] = w;
  ds->B[i_m] = w;         //move the B[i] pointer to w, 
  return;
}

/*! \brief Relax an edge, specifically given the head, w, and cand_dist, update tent[w].
 * \param ds the data structure that holds or points to everything
 * \param w the head vertex of the given edge
 * \param cand_dist = tent[v]+c(v,w) the candidate new weight
 */
void local_relax(ds_t *ds, int64_t w, double cand_dist)
{
  int64_t iold, inew;
  //if(DPRT){printf("%02d: relax head %"PRId64" cand_dist = %lf < %lf?\n", MYTHREAD, w, cand_dist, ds->tent[w]);}
  if ( cand_dist < ds->tent[w] ){
    iold = ds->in_bucket[w];
    inew = ((int64_t)floor(cand_dist/ds->delta)) % (ds->num_buckets);

    assert(iold >= -1 && iold < ds->num_buckets);
    assert(inew >= -1 && inew < ds->num_buckets);
    //if(DPRT){printf("%02d: winner: %"PRId64"  move from bucket %"PRId64" to %"PRId64"\n", MYTHREAD, w, iold, inew);}
    if( iold != inew ){
      if(iold >= 0)
        remove_node_from_bucket(ds, w);
      insert_node_in_bucket(ds, w, inew);
    }
    ds->tent[w] = cand_dist;
  }
}

/*!
\brief Allocate all the arrays need to hold pointer and values
\param ds the data structure that holds or points to everything
\param lnumrows the local number of rows on this thread
\param num_buckets the number of buckets 
\param delta delta
*/
void allocate_and_initialize_delta_stepping_struct(ds_t *ds, int64_t lnumrows, int64_t num_buckets, double delta)
{
  int64_t i;
  ds->next = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->next != NULL);
  ds->prev = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->prev != NULL);
  ds->in_bucket = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->in_bucket != NULL);
  ds->deleted = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->deleted != NULL);
  ds->R = (int64_t *)malloc(lnumrows * sizeof(int64_t)); assert(ds->R != NULL);
  ds->tent = (double *)malloc(lnumrows * sizeof(double)); assert(ds->tent != NULL);
  for(i = 0; i < lnumrows; i++){
    ds->next[i]      = i;
    ds->prev[i]      = i;
    ds->in_bucket[i] = -1;
    ds->deleted[i]   =  0;
    ds->tent[i]      = INFINITY;
  }
  if(D0PRT){printf("Allocate buckets\n");}
  ds->num_buckets = num_buckets;
  ds->B = (int64_t *)calloc(num_buckets, sizeof(int64_t)); assert(ds->B != NULL);

  for(i = 0; i < num_buckets; i++){
    ds->B[i] = -1;
  }
  ds->delta = delta;
}

/*! \brief Clear the memory used by ds
\param ds the data structure that holds or points to everything
*/
void clear_ds_struct(ds_t *ds)
{
  free(ds->next);
  free(ds->prev);
  free(ds->in_bucket);
  free(ds->deleted);
  free(ds->R);
  free(ds->tent);
}

/*!
\brief Determines delta and number of buckets, based on commandline options and the graph
\param delta the place to write the computed delta
\param num_buckets the place to write the number of buckets 
\param mat the sparse matrix that holds the graph
\param opt_delta the delta given on the command line, delta=0.0 means compute it here
*/
void calculate_delta_and_num_buckets(double *delta, int64_t *num_buckets, sparsemat_t *mat, double opt_delta)
{ 
  double del;
  int64_t nb;
  int64_t i;

  if( opt_delta > 0.0 && opt_delta <= 1.0) {
    del = opt_delta;
  } else {
    int64_t max_degree = 0;
    for(i = 0; i < mat->lnumrows; i++){
      if(max_degree < (mat->loffset[i+1] - mat->loffset[i]))
        max_degree = (mat->loffset[i+1] - mat->loffset[i]);
    }
    max_degree = lgp_reduce_max_l(max_degree);
    assert(max_degree > 0);
    del = 1.0/max_degree;
  }
  if(D0PRT){printf("%02d:del = %lf\n", MYTHREAD, del);}
  
  double max_edge_weight = 0.0;
  for(i = 0; i < mat->lnnz; i++)
    if(max_edge_weight < mat->lvalue[i])
      max_edge_weight = mat->lvalue[i];
  max_edge_weight = lgp_reduce_max_d(max_edge_weight);
  if(D0PRT){printf("%02d:max edge weight = %lf\n", MYTHREAD, max_edge_weight);}

  if(D0PRT){printf("%02d:Init Buckets\n", MYTHREAD);}
  nb = (int64_t)ceil(max_edge_weight/del) + 1;
  *delta = del;
  *num_buckets = nb;
}


/*! \brief Produces a sparse matrix with only the light edges (weights are less < delta)
 * \param mat the sparse matrix that holds the whole graph
 * \param delta  delta
 * \return the light sparse matrix
 */
sparsemat_t * get_light_edges(sparsemat_t *mat, double delta)
{
  int64_t i, k;
  int64_t light_lnnz=0;

  for(i = 0; i < mat->lnnz; i++){
    if(mat->lvalue[i] <= delta)
      light_lnnz++;
  }

  sparsemat_t *retmat = init_matrix(mat->numrows, mat->numcols, light_lnnz, 1);
  retmat->loffset[0] = 0;
  int64_t lpos=0;
  for(i=0; i<mat->lnumrows; i++){
    for(k=mat->loffset[i]; k<mat->loffset[i+1]; k++){
      if(mat->lvalue[k] <= delta){
        retmat->lnonzero[lpos] = mat->lnonzero[k];
        retmat->lvalue[lpos] = mat->lvalue[k];
        lpos++;
      }
    }
    retmat->loffset[i+1] = lpos;
  }
  if(D0PRT){printf("%02d: nnz %ld =  light  %ld\n", MYTHREAD, k, lpos);}
  return(retmat);
}

/*! \brief Produces a sparse matrix with only the light edges (weights are less >= delta)
 * \param mat the sparse matrix that holds the whole graph
 * \param delta  delta
 * \return the light sparse matrix
 */
sparsemat_t * get_heavy_edges(sparsemat_t *mat, double delta)
{
  int64_t i, k;
  int64_t heavy_lnnz=0;

  for(i = 0; i < mat->lnnz; i++){
    if(mat->lvalue[i] > delta)
      heavy_lnnz++;
  }

  sparsemat_t *retmat = init_matrix(mat->numrows, mat->numcols, heavy_lnnz, 1);
  retmat->loffset[0] = 0;
  int64_t hpos=0;
  for(i=0; i<mat->lnumrows; i++){
    for(k=mat->loffset[i]; k<mat->loffset[i+1]; k++){
      if(mat->lvalue[k] > delta){
        retmat->lnonzero[hpos] = mat->lnonzero[k];
        retmat->lvalue[hpos] = mat->lvalue[k];
        hpos++;
      }
    }
    retmat->loffset[i+1] = hpos;
  }
  if(D0PRT){printf("%02d: nnz %ld = heavy %ld\n", MYTHREAD, k, hpos);}
  return(retmat);
}


