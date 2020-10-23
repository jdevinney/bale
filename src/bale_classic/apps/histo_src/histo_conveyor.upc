/******************************************************************
//
//
//  Copyright(C) 2019-2020, Institute for Defense Analyses
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

/*! \file histo_conveyor.upc
 * \brief A conveyor implementation of histogram.
 */
#include "histo.h"

/*!
 * \brief This routine implements histogram using conveyors
 * \param data the histo_t struct that carries all the parameters for the implementations
 * \return average run time
 */
double histo_conveyor(histo_t * data){
  int ret;
  int64_t i;
  double tm;
  int64_t pe, col;
  int64_t pop_col;
  
  minavgmaxD_t stat[1];

  int status = EXIT_FAILURE;
  convey_t* conveyor = convey_new(SIZE_MAX, 0, NULL, convey_opt_SCATTER);
  if(!conveyor){printf("ERROR: histo_conveyor: convey_new failed!\n"); return(-1.0);}
  
  ret = convey_begin(conveyor, sizeof(int64_t), 0);
  if(ret < 0){printf("ERROR: histo_conveyor: begin failed!\n"); return(-1.0);}

  lgp_barrier();  
  tm = wall_seconds();
  i = 0UL;
  while(convey_advance(conveyor, i == data->l_num_ups)) {
    for(; i< data->l_num_ups; i++){
      col = data->pckindx[i] >> 20;
      pe  = data->pckindx[i] & 0xfffff;
      assert(pe < THREADS);
      if( !convey_push(conveyor, &col, pe))
	    break;
    }
    while( convey_pull(conveyor, &pop_col, NULL) == convey_OK){
      assert(pop_col < data->lnum_counts);
      data->lcounts[pop_col] += 1;
    }
  }

  lgp_barrier();
  tm = wall_seconds() - tm;

  lgp_min_avg_max_d( stat, tm, THREADS );

  lgp_barrier();
  convey_free(conveyor);
  return( stat->avg );
}
