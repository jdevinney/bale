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

/*! \file histo.h
 * \brief header file for the histogram app.
 */
#ifndef HISTO_H
#define HISTO_H
#include <libgetput.h>
#include <exstack.h>
#include <convey.h>
#include <locale.h>

/*!
\brief A structure to carry all the histogram arrays, counts to different implementations,
and aids in error checking
*/ 
typedef struct histo_t {
  SHARED int64_t * counts;  /*!< the shared array that holds the histogram counts */
  int64_t * lcounts;        /*!< the local pointer to the per thread parts of counts */
  int64_t num_counts;       /*!< the global size of the counts array */
  int64_t lnum_counts;      /*!< the local size of the counts array */
  int64_t * index;          /*!< the local index array */
  int64_t * pckindx;        /*!< the packed index with the divmod calculation already done */
  int64_t l_num_ups;        /*!< the local number of update to do */
} histo_t;

double histo_agp(histo_t * data);                       /*!< The AGP implementation */
double histo_exstack(histo_t * data, int64_t buf_cnt);  /*!< The EXSTACK implementation */
double histo_exstack2(histo_t * data, int64_t buf_cnt); /*!< The EXSTACK2 implementation */
double histo_conveyor(histo_t * data);                  /*!< The CONVEYOR implementation */

#include "alternates/histo_alternates.h"

#endif
