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

/*! \file sssp_delta_common.h
 * \brief Header for the data structure and common functions common to all the delta stepping implementations.
 */

#ifndef sssp_delta_common_INCLUDED
#define sssp_delta_common_INCLUDED

/*! \brief set to 1 for debug print statements on all threads */
#define DPRT 0
/*! \brief set to 1 for debug print statements on thread 0 */
#define D0PRT (!MYTHREAD && DPRT)

/*! 
\brief a struct to hold (pointers to) all the state needed for the Delta-Stepping algorithm

The storage required is 6 times the number of nodes plus a little for B.
Need a doubly link list because we delete nodes from any where on the list.
Note: It is interesting the all of this state is held in local arrays.
*/
typedef struct ds_t{
  int64_t *next;        /*!< next in linked list */
  int64_t *prev;        /*!< prev in linked list */
  int64_t *in_bucket;   /*!< which bucket it is in, else -1 */
  int64_t *B;           /*!< array of Buckets: B[i] is the index of the first node on the list, or -1 if empty */
  int64_t num_buckets;  /*!< num_buckets needed for "recycling" */
  double  *tent;        /*!< the tentative weight for the vertex */
  int64_t *deleted;     /*!< deleted means resolved? */
  int64_t *R;           /*!< queue to hold tail vertices that need to relax their heavy edges */
  double  delta;        /*!< delta parameter */
}ds_t;

void dump_bucket(ds_t *ds, int64_t i_m);                           /*!< debugging aid that dumps the bucket */
void insert_node_in_bucket(ds_t *ds, int64_t v, int64_t i_m);      /*!< insert a node into given bucket */
void remove_node_from_bucket(ds_t *ds, int64_t v);                 /*!< remove a given node from a bucket (if it is in bucket) */
void local_relax(ds_t *ds, int64_t w, double cand_dist);           /*!< locally relax the heads of edges (exactly the serial code) */
sparsemat_t *get_light_edges(sparsemat_t *mat, double delta);      /*!< sparse matrix to hold the "light" edges */
sparsemat_t *get_heavy_edges(sparsemat_t *mat, double delta);      /*!< sparse matrix to hold the "heavy" edges */


void calculate_delta_and_num_buckets(double *delta, int64_t *num_buckets, sparsemat_t *mat, double opt_delta); /* pick or use the given delta and find the number of buckets needed for that delta */
void allocate_and_initialize_delta_stepping_struct(ds_t *ds, int64_t lnumrows, int64_t num_buckets, double delta); /* all the data needed for the algorithm is in this one struct */
void clear_ds_struct(ds_t *ds); /*!< free the arrays within the struct */

#endif  // sssp_delta_common_INCLUDED
