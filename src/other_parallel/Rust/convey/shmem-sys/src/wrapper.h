//
// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of Convey, a conveyor library for rust.  For
// licence information see the LICENSE file in the top level dirctory
// of the distribution.

// This define is needed so that Cray's shmem.h can be parsed by bindgen
//  The symbol referenced is a #define defined after use, which seems to work
//  in cc but not bindgen
#define _SHMEM_MAX_RADIX 64
#include <shmem.h>
