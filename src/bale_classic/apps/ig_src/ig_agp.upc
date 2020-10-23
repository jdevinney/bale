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

/*! \file ig_agp.upc
 * \brief The intuitive implementation of indexgather that uses single word gets to shared addresses.
 */
#include "ig.h"

/*!
 * \brief This routine implements the single word get version indexgather
 * \param *tgt array of target locations for the gathered values
 * \param *index array of indices into the global array of counts
 * \param l_num_req the length of the index array
 * \param *table shared pointer to the shared table array.
 * \return average run time
 *
 */
double ig_agp(int64_t *tgt, int64_t *index, int64_t l_num_req,  SHARED int64_t *table) {
  int64_t i;
  double tm;
  minavgmaxD_t stat[1];

  lgp_barrier();
  tm = wall_seconds();

  for(i = 0; i < l_num_req; i++){
    #if __cray__ || _CRAYC
    #pragma pgas defer_sync
    #endif
    tgt[i] = lgp_get_int64(table, index[i]);
  }

  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  return( stat->avg );
}

