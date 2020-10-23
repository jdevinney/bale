// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#ifndef CONVEY_PORTER_H
#define CONVEY_PORTER_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "convey_alc8r.h"
#include "convey_codec.h"


typedef struct porter porter_t;
// The fields of a buffer are used for different purposes externally
// -- in buffers obtained by porter_borrow() -- and internally.
typedef struct buffer {
  uint32_t start;    // first valid item starts at data + start
  // (internally: may hold a termination flag)
  uint32_t limit;    // data + limit is just beyond last valid item
  uint32_t n_items;  // positive <==> buffer is compressed
  uint32_t source;   // transmitter of most recent hop
  char data[];       // each item is preceded by a tag and may be padded
} buffer_t;

// Compute the packet size in the non-elastic case.
size_t porter_packet_size(size_t tag_bytes, size_t item_bytes);

// A variable-length item starts with a 4-byte descriptor which equals
// 2 * (item length in bytes) + ticket_flag.  A 'ticket' ends its buffer
// and must be transmitted ASAP.

#define PORTER_DESCR(bytes, flag) (2*(bytes) | (flag))
#define PORTER_TICKET(descr) ((descr) & 1)
#define PORTER_BYTES(descr) ((descr) >> 1)
#define PORTER_QUADS(descr) (((descr) + 22) >> 3)
// that's an optimized version of (2 + ((descr)/2 + 3)/4)

// Creation of a porter is a collective operation.  Locally, the porter is
// able to communicate with the n given PEs (whose indices are relative to
// the current communicator), and my index among their friends is my_rank.
// (For 0 <= i < n, if PE friends[i] sent me a copy of its 'friends' array,
// then copy[my_rank] would be MY_PROC.)  The 'friends' array is absorbed
// by the porter.
//
// The 'multiplicity' is the number of buffers per friend; it must be a
// power of two.  The 'opcode' is used for mpp profiling of buffer sends.
// The recognized options are: convey_opt_DYNAMIC, convey_opt_PROGRESS,
// convey_opt_ALERT, and also porter_opt_LOCAL, which indicates that
// communication is thought to be within a node.  The tag_bytes value,
// which must be 0, 1, 2, 4, or 8, determines the format of buffer data.
// Each tag is stored in a field of size tag_bytes, and the space used to
// hold each item is rounded up (if tag_bytes > 0) to a multiple of
// tag_bytes.

enum porter_option {
  porter_opt_LOCAL = 0x10000,
};

porter_t* porter_new(int n, int32_t friends[n], int my_rank,
                     size_t tag_bytes, size_t capacity, size_t multiplicity,
                     const convey_alc8r_t* alloc, uint64_t options, int opcode);

// Tear down any existing codec data, and prepare to set up the given codec.
// Should only fail (return false) if self->codata is NULL.
bool      porter_set_codec(porter_t* self, const convey_codec_t* codec, void* arg);

bool      porter_setup(porter_t* self, size_t item_size);

// Push a fixed-length item of item_size bytes.  Calls to porter_push and
// porter_epush must not be mixed!  A porter may be used for one kind of
// push or the other, but not both.
bool      porter_push(porter_t* self, uint64_t tag, const void* item, int dest);

// Push a variable-length item, preceded by its tag and descriptor.
bool      porter_epush(porter_t* self, uint32_t tag, uint32_t descr,
                       const void* item, int dest);

// Attempt to retrieve a nonempty buffer of items.
buffer_t* porter_borrow(porter_t* self);

// Give back the most recently taken buffer, having adjusted buffer->start.
void      porter_return(porter_t* self);

bool      porter_advance(porter_t* self, bool done);

// Return porter->packet_bytes
size_t    porter_stride(porter_t* self);

// Obtain some statistics (since the last setup).  Currently which=0 means
// buffer sends, which=1 means synchronization operations, and which=2 means
// bytes actually transferred.
int64_t   porter_get_stats(porter_t* self, int which);

void      porter_reset(porter_t* self);

// Tear down a non-NULL porter without waiting for other PEs
// and without freeing the 'relative' array.
void      porter_demolish(porter_t* self);

void      porter_free(porter_t* self);


#endif
