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

/*! \file histo_alternates.h
 * \brief header file for the alternate models for histo
 */

double histo_exstack2_goto(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t buf_cnt);
double histo_exstack2_cyclic(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t buf_cnt);
double histo_exstack_function(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t buf_cnt);
double histo_exstack2_function(int64_t *pckindx, int64_t T,  int64_t *lcounts, int64_t buf_cnt);
