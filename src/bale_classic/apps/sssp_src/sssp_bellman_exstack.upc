/******************************************************************
//
//
//  Copyright(C) 2018-2020, Institute for Defense Analyses
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

/*! \file sssp_bellman_exstack.upc
 * \brief Implementation of sssp using exstack
 */

#include "sssp.h"

/*!
 * \brief Relax the head of the edges delivered by an exstack buffer
 * \param tent pointer to the tentative distances array
 * \param *ex the extack buffers
 * \param done the signal to exstack_proceed that this thread is done
 * \return the return value from exstack_proceed
 */
static int64_t bellman_exstack_relax_process(d_array_t *tent, exstack_t *ex, int64_t done) 
{
  int64_t fromth;
  sssp_pkg_t pkg;

  exstack_exchange(ex);
  
  while(exstack_pop(ex, &pkg, &fromth)){
    if( tent->lentry[pkg.lj] > pkg.tw ){
      if(0){printf("Ex: replace %ld weight %lg with %lg\n", pkg.lj*THREADS + MYTHREAD, tent->lentry[pkg.lj], pkg.tw);}
      tent->lentry[pkg.lj] = pkg.tw;
    }
  }
  return( exstack_proceed(ex, done) );
}

/*!
 * \brief Push the potentially improved weight to the thread handling the head of the edge
 * \param ex the extack buffers
 * \param tent pointer to the tentative distances array to be passed thru to bellman_exstack_relax_process
 * \param J the head of the edge, given by it global name
 * \param tw the new weight
 * \return the value from the push
 */
static int64_t bellman_exstack_push(exstack_t *ex, d_array_t *tent, int64_t J, double tw)
{
  int64_t ret, pe;
  sssp_pkg_t pkg;
  pe     = J % THREADS;
  pkg.lj = J / THREADS;
  pkg.tw = tw;
  if((ret = exstack_push(ex, &pkg, pe)) == 0){
    bellman_exstack_relax_process(tent, ex, 0);
  }
  return(ret);
}


/*!
\brief This routine implements the Bellman-Ford algorithm using exstack
\param dist pointer to the array for return the distances (weights) to each vertex
\param mat the input matrix that holds the graph
\param buf_cnt the number of packages in an exstack buffer
\param v0 is the the staring row (vertex)
\return average run time
*/
double sssp_bellman_exstack(d_array_t *dist, sparsemat_t * mat, int64_t buf_cnt, int64_t v0)
{
  double t1 = wall_seconds();

  if(!mat){ T0_printf("ERROR: sssp_bellman_exstack: NULL L!\n"); return(-1); }
  
  int64_t k, li, J, pe;
  int64_t pe_v0, li_v0;
  int64_t loop;
  int64_t changed;
  sssp_pkg_t pkg;
  d_array_t *tent0, *tent1, *tent2;
  d_array_t *tent_old, *tent_cur, *tent_new, *tent_temp;

  exstack_t * ex = exstack_init(buf_cnt, sizeof(sssp_pkg_t));
  if( ex == NULL) return(-1.0);


  pe_v0 = v0 % THREADS;
  li_v0 = v0 / THREADS;

  tent0 = init_d_array(dist->num);
  set_d_array(tent0, INFINITY);
  lgp_barrier();
  if(pe_v0 == MYTHREAD){
    tent0->lentry[li_v0] = 0.0;
  }
  lgp_barrier();
  tent1 = copy_d_array(tent0);
  if(pe_v0 == MYTHREAD ){
    for(k = mat->loffset[li_v0]; k < mat->loffset[li_v0+1]; k++){
      lgp_put_double(tent1->entry, mat->lnonzero[k], mat->lvalue[k]);
    }
  }

  lgp_barrier();

  tent2 = init_d_array(tent1->num);

  tent_old = tent0;
  tent_cur = tent1;
  tent_new = tent2;
  for(loop=2; loop<mat->numrows; loop++){
    changed = 0;
    for(li=0; li < mat->lnumrows; li++)
      tent_new->lentry[li] = tent_cur->lentry[li];

    for(li=0; li < mat->lnumrows; li++){
      if(tent_old->lentry[li] == tent_cur->lentry[li])
        continue;
      changed = 1;
      for(k=mat->loffset[li]; k< mat->loffset[li+1]; k++){
        if(0){printf("%ld %d: relaxing (%ld,%ld)   %lg %lg\n", loop, MYTHREAD, li*THREADS + MYTHREAD, J, tent_cur->lentry[li],  pkg.tw);}
        if( bellman_exstack_push(ex, tent_new, mat->lnonzero[k],  tent_cur->lentry[li] + mat->lvalue[k]) == 0)
            k--;
      }
    }
    while(bellman_exstack_relax_process(tent_new, ex, 1))  // keep popping til all threads are done
      ;

    lgp_barrier();
    if( lgp_reduce_add_l(changed) == 0 ){                  // keep looping until nothing changes
      replace_d_array(dist, tent_new);
      break;
    }

    //dump_tent("EXSTACK:  ", tent_new);

    tent_temp = tent_old;
    tent_old = tent_cur;
    tent_cur = tent_new;
    tent_new = tent_temp;

    lgp_barrier();
    exstack_reset(ex);
  }

  //dump_tent("Bellman EXSTACK: ", tent);
  lgp_barrier();
  clear_d_array(tent0); free(tent0);
  clear_d_array(tent1); free(tent1);
  clear_d_array(tent2); free(tent2);
  t1 = wall_seconds() - t1;
  minavgmaxD_t stat[1];
  lgp_min_avg_max_d( stat, t1, THREADS );

  return(stat->avg);
}

