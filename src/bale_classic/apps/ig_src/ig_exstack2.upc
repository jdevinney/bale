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

/*! \file ig_exstack2.upc
 * \brief An exstack2 implementation of indexgather.
 */
#include "ig.h"

/*!
 * \brief This routine implements the exstack2 variant of indexgather.
 * \param *tgt array of target locations for the gathered values
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param l_num_req the length of the pcindx array
 * \param *ltable localized pointer to the count array.
 * \param buf_cnt the exstack buffer size in packets
 * \return average run time
 *
 */
double ig_exstack2(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *ltable, int64_t buf_cnt) {
  double tm;
  int64_t more_req = 1;
  int64_t ret;
  int64_t i;
  int64_t pe, fromth, fromth2;
  minavgmaxD_t stat[1];


  typedef struct index_val_t{
    int64_t i;    
    int64_t idx;
  }index_val_t;

  index_val_t pkg, pkge, pkgr;

  exstack2_t * ex2r = exstack2_init(buf_cnt , sizeof(index_val_t));
  exstack2_t * ex2e = exstack2_init(buf_cnt , sizeof(index_val_t));
  if( (ex2r == NULL) || (ex2e == NULL) ) return(-1.0);

  lgp_barrier();
  
  tm = wall_seconds();
  i=0;
  while( more_req || exstack2_proceed(ex2e, !more_req) ){
    while(i < l_num_req){
      pkgr.i = i;
      pkgr.idx = pckindx[i] >> 16;
      pe = pckindx[i] & 0xffff;

      if( !exstack2_push(ex2r, &pkgr, pe) )
        break;
      i++; 
    }

    while(exstack2_pop(ex2r, &pkg, &fromth2)){
      pkge.i = pkg.i;
      pkge.idx = ltable[pkg.idx];
      if( !exstack2_push(ex2e, &pkge, fromth2) ){ // if we can't push, then we must
        exstack2_unpop(ex2r);                     // unpop the last request
        break;
      }
    }

    while(exstack2_pop(ex2e, &pkg, NULL)){      // process responses to your trg_len
      tgt[pkg.i] = pkg.idx;
    }

    more_req = exstack2_proceed( ex2r, (i==l_num_req) );
  }

  tm = wall_seconds() - tm;
  lgp_barrier();

  lgp_min_avg_max_d( stat, tm, THREADS );
  exstack2_clear(ex2r);
  exstack2_clear(ex2e);
  free(ex2r);
  free(ex2e);
  return( stat->avg );
}

