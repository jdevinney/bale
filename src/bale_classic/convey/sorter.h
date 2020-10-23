// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#ifndef CONVEY_SORTER_H
#define CONVEY_SORTER_H

#include <stdbool.h>
#include <stddef.h>
#include "convey_alc8r.h"


// A sorter is an object that can efficiently distribute a stream of
// fixed-sized items into an array of buffers.  ...

typedef struct sorter sorter_t;

typedef struct area {
  char* next;   // place to write or read next item
  char* limit;  // upper limit of the area
} area_t;


sorter_t* sorter_new(int n, area_t areas[n], size_t item_bytes, size_t capacity,
                     const convey_alc8r_t* alloc, bool dynamic);

// Returns false if something goes wrong (a memory allocation error).
bool sorter_setup(sorter_t* self);

// Returns false if one of the outgoing buffers is full.  The item is
// always successfully pushed, provided that the caller responds to a
// return value of 'false' by emptying the full buffer(s).
bool sorter_push(sorter_t* self, const void* item, int dest);

// Returns 0 if there is nothing to flush, -1 if an outgoing buffer
// filled up, and 1 otherwise.
int sorter_flush(sorter_t* self);

void sorter_reset(sorter_t* self);
void sorter_free(sorter_t* self);


#endif
