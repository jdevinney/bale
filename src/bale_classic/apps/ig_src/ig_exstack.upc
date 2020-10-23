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

/*! \file ig_exstack.upc
 * \brief The classic exstack implementation of indexgather.
 */
#include "ig.h"

/*!
 * \brief This routine implements the exstack classic variant of indexgather.
 * \param *tgt array of target locations for the gathered values
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param l_num_req the length of the pcindx array
 * \param *ltable localized pointer to the count array.
 * \param buf_cnt the exstack buffer size in packets
 * \return average run time
 *
 */
double ig_exstack(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *ltable, int64_t buf_cnt) {
  double tm;
  int imdone;
  int64_t ret;
  int64_t room;
  int64_t l_indx, idx, i, i0, j;
  int64_t pe, fromth;
  minavgmaxD_t stat[1];

  exstack_t * ex = exstack_init(buf_cnt, sizeof(int64_t));
  if( ex == NULL ) return(-1.0);
  
  lgp_barrier();
  tm = wall_seconds();
  i=0;
  while( exstack_proceed(ex, (i==l_num_req)) ) {
    i0 = i;  
    while(i < l_num_req) {
      l_indx = pckindx[i] >> 16;
      pe  = pckindx[i] & 0xffff;
      if(!exstack_push(ex, &l_indx, pe)) 
        break; 
      i++;
    }   

    exstack_exchange(ex);
    
    while(exstack_pop(ex, &idx , &fromth)) {
      idx  = ltable[idx];
      exstack_push(ex, &idx, fromth);   // don't need check for room 
    }   
    lgp_barrier();
    exstack_exchange(ex);
  
    for(j=i0; j<i; j++) {  // retrace the requests
      fromth = pckindx[j] & 0xffff;
      exstack_pop_thread(ex, &idx, (uint64_t)fromth);
      tgt[j] = idx;
    }
    lgp_barrier();
  }
  tm = wall_seconds() - tm;
  lgp_barrier();
  lgp_min_avg_max_d( stat, tm, THREADS );

  exstack_clear(ex);
  free(ex);
  return( stat->avg );
}

