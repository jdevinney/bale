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

/*! \file toposort_alternates.h
 * \brief header file for the alternate models for toposort
 */

double toposort_matrix_oo(SHARED int64_t *tri_rperm, SHARED int64_t *tri_cperm, sparsemat_t *mat, sparsemat_t *tmat);
double toposort_matrix_cooler(SHARED int64_t *tri_rperm, SHARED int64_t *tri_cperm, sparsemat_t *mat, sparsemat_t *tmat);
double toposort_matrix_exstack_orig(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat);


