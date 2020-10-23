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

/*! \file histo_exstack2_cyclic.upc
 * \brief The exstack2 implementation of histogram looping over single pushes and pops.
 */
#include "histo.h"

/*!
 * \brief This routine implements the exstack2 variant of histogram where one 
 *  pushes and pops, pushes and pops, .... instead of 
 *  pushing until you can't, then popping until you can't.
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param l_num_ups the length of the pcindx array
 * \param *lcounts localized pointer to the count array.
 * \param buf_cnt the size of the exstack buffers in packages
 * \return average run time
 */
double histo_exstack2_cyclic(int64_t *pckindx, int64_t l_num_ups,  int64_t *lcounts, int64_t buf_cnt) {
  int ret;
  double tm;
  int64_t i;
  int64_t pe, col, idx;
  minavgmaxD_t stat[1];
  exstack2_t * ex2 = exstack2_init(buf_cnt, sizeof(int64_t));
  if( ex2 == NULL ) return(-1.0);

  lgp_barrier();
  tm = wall_seconds();
  for( i=0;  exstack2_proceed( ex2, (i==l_num_ups) );  ) {
    if( i < l_num_ups ) {
      col = pckindx[i] >> 16;
      pe  = pckindx[i] & 0xffff;
      if( exstack2_push(ex2, &col, pe) )
        i++;
    }
    if( exstack2_pop(ex2, &idx, NULL) )
      lcounts[idx]++;
  }
  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  lgp_barrier();
  exstack2_clear(ex2);
  free(ex2);
  return( stat->avg );
}
