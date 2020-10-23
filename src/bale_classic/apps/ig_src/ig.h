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

/*! \file ig.h
 * \brief Demo program that computes an indexed gather of elements from a 
 *  shared array into local arrays.
 *  The size of the source array should be large enough that the elements 
 *  need to be spread across the whole machine
 */

#ifndef IG_H
#define IG_H
#include <exstack.h>
#include <convey.h>

double ig_agp(int64_t *tgt, int64_t *index, int64_t l_num_req,  SHARED int64_t *table);
double ig_exstack(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *ltable, int64_t buf_cnt); 
double ig_exstack2(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *ltable, int64_t buf_cnt); 
double ig_conveyor(int64_t *tgt, int64_t *pckindx, int64_t l_num_req,  int64_t *ltable);

#include "alternates/ig_alternates.h"

#endif
