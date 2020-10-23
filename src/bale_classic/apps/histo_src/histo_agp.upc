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
/*! \file histo_agp.upc
 * \brief The intuitive implementation of histogram that uses global atomics.
 */
#include "histo.h"

/*!
 * \brief This routine implements straight forward, single word atomic updates to implement histogram.
 * \param data the histo_t struct that carries all the parameters for the implementations
 * \return average run time
 */
double histo_agp(histo_t * data){
  double tm;
  int64_t i;
  minavgmaxD_t stat[1];

  lgp_barrier();
  tm = wall_seconds();

  for(i = 0; i < data->l_num_ups; i++) {
    #if __cray__ || _CRAYC
      #pragma pgas defer_sync
    #endif
    //_amo_aadd(&counts[index[i]], 1);
    //counts[index[i]] += 1;
    assert(data->index[i] < data->num_counts);
    lgp_atomic_add(data->counts, data->index[i], 1L);
  }
  
  lgp_barrier();
  tm = wall_seconds() - tm;
  lgp_min_avg_max_d( stat, tm, THREADS);
  
  return( stat->avg );
}

