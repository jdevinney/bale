// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpp2nil.h"


void*
xnil_alloc_align(size_t align, size_t size)
{
  void* ptr = NULL;
  posix_memalign(&ptr, align, size);
  return ptr;
}

void
xnil_mfprint(FILE* stream, int prefix, const char* func, const char* format, ...)
{
  va_list args;
  va_start(args, format);
  if (prefix)
    fprintf(stream, "%s> ", func);
  vfprintf(stream, format, args);
  va_end(args);
  fflush(stream);
}

