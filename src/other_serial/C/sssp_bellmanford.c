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

/*! \file sssp_bellmanford.c
\brief Implementations of several versions of Bellman-Ford Alg
*/

/*!  \brief This routine implements "relax edge" for all the implementation of Bellman-Ford.
\param tent_head pointer to the tent weight of the head vertex
\param tent_tail pointer to the tent weight of the tail vertex
\param edge_wt the edge weight
\return 1 if the tentative weight of the head vertex has improved, 0 otherwise
*/
static int relax_edge( double *tent_head, double *tent_tail, double edge_wt )
{
   if( *tent_head > *tent_tail + edge_wt ){
     if(0)printf("relaxing %lg to %lg\n", *tent_head, (*tent_tail + edge_wt)); 
     *tent_head = *tent_tail + edge_wt;
     return(1);
   }
   return(0);
}


/*! \brief This routine implements the most naive version of the Bellman-Ford algorithm
\param weight array report the weight of the lightest path to the node
\param mat sparsemat_t that holds the graph. 
\param v0 is the starting vertex
\return runtime

NB: One should only run this on very small cases because it takes so long
*/
double sssp_bellmanford_simple(d_array_t *weight, sparsemat_t * mat, int64_t v0)
{
  double tm = wall_seconds();
  int64_t i, k, loop;
  int64_t numrows = mat->numrows;
  
  for(i=0; i<numrows; i++)
    weight->entry[i] = INFINITY;
  weight->entry[v0] = 0.0;

  for(loop=0; loop<numrows; loop++){
    for(i=0; i<numrows; i++){ 
      for(k = mat->offset[i]; k < mat->offset[i+1]; k++){
        relax_edge( &(weight->entry[ mat->nonzero[k] ]), &(weight->entry[i]), mat->value[k]);
      }
    }
  }

  return(wall_seconds() - tm);
}


/*! \brief This routine implements the textbook version of the 
 Bellman-Ford algorithm as a Dynamic Programming algorithm.
\param weight array report the weight of the lightest path to the node
\param mat sparsemat_t that holds the graph. 
\param v0 is the starting vertex

NB: Only run this on small problems because storage is also bad.
*/
double sssp_bellmanford_dynprog(d_array_t *weight, sparsemat_t * mat, int64_t v0)
{
  double tm = wall_seconds();
  int64_t i,j,k, loop;
  int64_t numrows = mat->numrows;
  int64_t changed;

  double ** tent =  (double**) malloc( numrows * sizeof(double*) );
  for(i=0; i<numrows; i++)
    tent[i] =  (double*) malloc( numrows * sizeof(double) );

  loop = 0;
	for(i=0; i<numrows; i++){
		tent[loop][i] = INFINITY;
  }
  tent[loop][v0] = 0.0;
  if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent[loop][i]); printf("\n");}

  loop = 1;
	for(i=0; i<numrows; i++){
		tent[loop][i] = tent[loop-1][i];
  }
  for(k = mat->offset[v0]; k < mat->offset[v0+1]; k++){
    j = mat->nonzero[k];
    tent[loop][j] = mat->value[k];
  }
  if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent[loop][i]); printf("\n");}

  for(loop=2; loop<numrows; loop++){
    changed = 0;
    for(i=0; i<numrows; i++)
      tent[loop][i] = tent[loop-1][i];
    for(i=0; i<numrows; i++){
      if( tent[loop-2][i] == tent[loop-1][i] )
        continue;
      for(k = mat->offset[i]; k < mat->offset[i+1]; k++){
        changed |= relax_edge( &(tent[loop][mat->nonzero[k]]), &(tent[loop-1][i]), mat->value[k] );
      }
    }
    if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent[loop][i]); printf("\n");}
    if(changed == 0){
      break;
    }
  }
  assert( loop <= numrows );
	for(i=0; i<numrows; i++)
	  weight->entry[i] = tent[loop][i];

  for(i=0; i<numrows; i++)
    free(tent[i]);
  free(tent);

  return(wall_seconds() - tm);
}

/*! \brief This routine implements a version of Bellman-Ford as a Dynamic Programming algorithm 
that "recycles" storage and tries to do less redundant relaxations.
\param weight array report the weight of the lightest path to the node
\param mat sparsemat_t that holds the graph. 
\param v0 is the starting vertex
*/
double sssp_bellmanford(d_array_t *weight, sparsemat_t * mat, int64_t v0)
{
  double tm = wall_seconds();
  int64_t i,k, loop;
  int64_t numrows = mat->numrows;
  int64_t changed;
  double *tent_old, *tent_cur, *tent_new, *tent_temp;

  double *tent0 =  (double*) malloc( numrows * sizeof(double) );
  double *tent1 =  (double*) malloc( numrows * sizeof(double) );
  double *tent2 =  (double*) malloc( numrows * sizeof(double) );
  if( tent0 == NULL || tent1 == NULL || tent2 == NULL) return(-1.0);

  loop = 0;
	for(i=0; i<numrows; i++)
		tent0[i] = INFINITY;
  tent0[v0] = 0.0;
  if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent0[i]); printf("\n");}

  loop = 1;
	for(i=0; i<numrows; i++)
		tent1[i] = tent0[i];
  for(k = mat->offset[v0]; k < mat->offset[v0+1]; k++){
    tent1[mat->nonzero[k]] = mat->value[k];
  }
  if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent1[i]); printf("\n");}

  tent_old = tent0;
  tent_cur = tent1;
  tent_new = tent2;
  for(loop=2; loop<numrows; loop++){
    changed = 0;
    for(i=0; i<numrows; i++)
      tent_new[i] = tent_cur[i];
    for(i=0; i<numrows; i++){
      if( tent_old[i] == tent_cur[i] )
        continue;
      for(k = mat->offset[i]; k < mat->offset[i+1]; k++){
        changed |= relax_edge( &(tent_new[mat->nonzero[k]]), &(tent_cur[i]), mat->value[k] );
      }
    }
    if(0){printf("Bell %02ld : ",loop); for(i=0; i<numrows; i++) printf("%lg ",tent_new[i]); printf("\n");}
    if(changed == 0){
      break;
    }
    tent_temp = tent_old;
    tent_old = tent_cur;
    tent_cur = tent_new;
    tent_new = tent_temp;
  }
  assert( loop <= numrows );
	for(i=0; i<numrows; i++)
	  weight->entry[i] = tent_new[i];

  free(tent0);
  free(tent1);
  free(tent2);

  return(wall_seconds() - tm);
}

