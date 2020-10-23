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

/*! \file toposort.h
 * \brief Demo application that does a toposort on a permuted upper triangular matrix
 */
#ifndef TOPOSORT_H
#define TOPOSORT_H
#include <libgetput.h>
#include <exstack.h>
#include <convey.h>
#include <spmat.h>
#include <locale.h>

double toposort_matrix_agp(SHARED int64_t *tri_rperm, SHARED int64_t *tri_cperm, sparsemat_t *mat, sparsemat_t *tmat);
double toposort_matrix_exstack(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat, int64_t buf_cnt);
double toposort_matrix_exstack2(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat, int64_t buf_cnt);
double toposort_matrix_convey(SHARED int64_t *rperm, SHARED int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat);

#include "alternates/toposort_alternates.h"

#endif
