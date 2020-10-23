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

/*! \file ig_exstack_function.upc
 * \brief The classic exstack implementation of indexgather using functions to hide the loops.
 */
#include "ig.h"
  
/*! \brief the package struct for elements of the exstack buffers */
typedef struct index_val_t {
  int64_t i;          /*!< this local index */
  int64_t idx;        /*!< the local index on remote thread OR the return value value */
} index_val_t;

/*!
 * \brief  routine to service the pops
 * \param *ex the exstack struct
 * \param *tgt the target array
 * \param *ltable localized pointer to the count array
 * \param imdone flag for this thread
 * \return average run time
 */
void ig_exstack_service(exstack_t *ex, int64_t *tgt, int64_t *ltable, int64_t imdone ) {
  int64_t fromth;
  index_val_t pkg;

  do {
    exstack_exchange(ex);

    while(exstack_pop(ex, &pkg , &fromth)) {
      pkg.idx  = ltable[pkg.idx];
      exstack_push(ex, &pkg, fromth); 
    }   
    lgp_barrier();
    exstack_exchange(ex);

    while(exstack_pop(ex, &pkg, NULL)){ 
      tgt[pkg.i] = pkg.idx;
    }
  } while (imdone && exstack_proceed( ex , imdone ));
}


/*!
 * \brief  routine to generate the requests
 * \param *ex the exstack struct
 * \param *pckindx array of "packed" indices
 * \param i the ith request
 * \param *tgt the target array
 * \param *ltable localized pointer to the count array
 * \return average run time
 *
 */
void ig_exstack_request(exstack_t *ex, int64_t *pckindx, int64_t i,  int64_t *tgt, int64_t *ltable) {
  int64_t l_indx, pe;
  index_val_t pkg;

  pkg.i = i;
  pkg.idx = pckindx[i] >> 16;
  pe = pckindx[i] & 0xffff;
  if( exstack_push(ex, &pkg, pe) > 1 ) 
    return;

  ig_exstack_service( ex, tgt, ltable, 0 );
}


/*!
 * \brief This routine implements the exstack classic variant of indexgather.
 * \param *tgt array of target locations for the gathered values
 * \param *pckindx array of packed indices for the distributed version of the global array of counts.
 * \param l_num_req the length of the pcindx array
 * \param *ltable localized pointer to the count array
 * \param buf_cnt the exstack buffer size in packages
 * \return average run time
 *
 */
double ig_exstack_function(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *ltable, int64_t buf_cnt) {
  double tm;
  int imdone;
  int64_t ret;
  int64_t room;
  int64_t l_indx, idx, i, i0, j;
  int64_t pe, fromth;
  minavgmaxD_t stat[1];

  index_val_t pkg;

  exstack_t * ex = exstack_init(buf_cnt, sizeof(index_val_t));
  if( ex == NULL ) return(-1.0);
  
  lgp_barrier();
  tm = wall_seconds();

  for(i=0; i<l_num_req; i++) {
    ig_exstack_request( ex, pckindx, i, tgt, ltable ); 
  }
  ig_exstack_service(ex, tgt, ltable, 1);

  lgp_barrier();
  
  tm = wall_seconds() - tm;
  lgp_barrier();
  lgp_min_avg_max_d( stat, tm, THREADS );

  exstack_clear(ex);
  free(ex);
  return( stat->avg );
}
