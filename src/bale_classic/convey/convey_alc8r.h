// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


/** \file convey_alc8r.h
 * The auxiliary API for symmetric memory allocation.
 */


#ifndef CONVEY_ALC8R_H
#define CONVEY_ALC8R_H

#include <stddef.h>
#include <stdint.h>

/** Memory allocation objects that can be passed to conveyor constructors.
 * 
 * Many conveyor constructors can be told how to allocate and deallocate
 * symmetric memory by passing them a pointer to a structure of this type.
 * The structure does not have to stay in scope; the constructor copies it.
 * If the pointer is \c NULL, the conveyor uses standard memory management
 * functions instead.
 */
typedef struct convey_alc8r {
  /// Passed as the first argument of the following functions.
  void* alc8r;
  /// The function for allocating memory: the \a tag and \a value will be the file name and line number of the call.
  void* (*grab)(void* alc8r, size_t size, const char* tag, uint64_t value);
  /// The function for releasing previously allocated memory.
  void  (*free)(void* alc8r, void* ptr);
} convey_alc8r_t;

#endif
