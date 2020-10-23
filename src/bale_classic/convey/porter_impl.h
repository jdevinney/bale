// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#ifndef PORTER_IMPL_H
#define PORTER_IMPL_H

#define PORTER_DEBUG 0

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "common.h"
#include "porter.h"

#if PORTER_DEBUG
# define DEBUG_PRINT(...) mprint(MY_PROC, 0, __VA_ARGS__)
#else
# define DEBUG_PRINT(...)
#endif


// FIXME: Describe buffer states here...

// For a "steady" conveyor, which cycles through the ranks to make sure
// all data eventually gets delivered, this is how many cycles to wait
// before moving things along.  It must fit in a uint8_t.
#define PATIENCE 2

typedef struct area {
  char* next;       // place to write next packet
  char* limit;      // upper limit of the area
} area_t;

// Data kept locally about each outgoing communication channel.
// By definition, produced >= emitted >= delivered.
typedef struct channel {
  uint64_t produced;    // counts outgoing buffers filled and closed
  uint64_t emitted;     // counts outgoing buffers for which sending has begun
  uint64_t delivered;   // counts outgoing buffers actually delivered
  uint64_t urgent;      // if positive, 1 + maximum index of an urgent buffer
  // Must achieve delivered >= urgent in order to send all tickets.
} channel_t;

typedef struct porter_methods porter_methods_t;
typedef struct porter_codata  porter_codata_t;

struct porter_methods {
  // const porter_methods_t* _super_;
  // Public methods.  Note that setup(), reset(), and demolish() are
  // responsible for allocation and deallocation of send_buffers.
  bool (*setup)(porter_t* self);
  buffer_t* (*borrow)(porter_t* self);
  void (*turnin)(porter_t* self);
  void (*reset)(porter_t* self);
  void (*demolish)(porter_t* self);

  // Private methods...
  // 0. Check whether the buffer with sequence number n for this
  //    destination (the next one to be sent) is ready to be received.
  bool (*ready)(porter_t* self, int dest, uint64_t n);

  // 1. Send a buffer; return true if the buffer has been delivered.
  bool (*send)(porter_t*, int, uint64_t, size_t, buffer_t*, uint64_t);

  // 2. Ensure progress; send necessary signals; update channels and n_urgent.
  //    If arg >= 0, make progress sending to that destination (at least).
  //    Else try to deliver of all outstanding sends to all destinations.
  //    Return true if delivery was successful.
  bool (*progress)(porter_t*, int);

  // 3. Free any extra data.
  void (*release)(void*);
};

// Supplementary data structure needed to support compression
struct porter_codata {
  convey_cargo_t cargo;  // item_size==0 means not set up
  convey_layout_t layout;
  size_t work_bytes;     // size of scratch buffer
  void* work;            // scratch buffer for de/compression
  void** compressors;    // array of length n_ranks
  void** decompressors;  // array of length n_ranks
  // Spread n_items packets from zone into self->work at self->layout.stride
  void (*unpack)(porter_codata_t* self, size_t n_items, const void* zone);
  // Repack n_items packets from self->work into zone
  void (*repack)(porter_codata_t* self, size_t n_items, void* zone);
  // FIXME: consider adding a divbymul32_t for item_size + tag_size
};

struct porter {
  const porter_methods_t* _class_;
  size_t tag_bytes;
  size_t item_bytes;
  size_t packet_bytes;  // tag_bytes + item_bytes + padding
  size_t buffer_bytes;  // bytes per buffer, including sizeof(buffer_t)
  size_t buffer_stride; // stride between buffers (rounded up for alignment)
  int n_ranks;
  int my_rank;          // we write to this position in remote arrays
  int abundance;        // log2(# of send buffers per communicant)
  int opcode;           // used for profiling
  // Local memory
  int32_t* relative;    // relative PE numbers
  area_t* send_areas;
  bool* all_sent;       // FIXME: perhaps combine this array with waiting[]
  channel_t* channels;
  uint8_t* waiting;     // tracks how long each link has waited with undelivered data
  // Symmetric memory
  PARALLEL(char*, send_buffers);
  PARALLEL(char*, recv_buffers);
  // State and statistics
  bool endgame;         // no more pushes?
  bool flushed;         // all buffers sent?
  bool drained;         // all buffers received and taken?
  bool dynamic;         // do we alloc on setup and dealloc on reset?
  bool compress;        // are we using a codec?
  bool alert;           // warn about performance issues?
  bool inmult;          // are receive buffers also abundant?
  int n_urgent;         // number of links with undelivered urgent buffers
  int phase;            // (# of advances + my_rank) modulo n_ranks
  int64_t send_count;   // counts sends of buffers
  int64_t sync_count;   // counts calls to quiet and fence
  int64_t byte_count;   // counts bytes transferred
  size_t buffer_align;  // actual alignment of send and receive buffers
  // Compression state
  const convey_codec_t* codec;  // non-NULL if compression is active
  porter_codata_t* codata;      // non-NULL if this porter supports compression
  // Embedded allocator
  convey_alc8r_t alloc;
};

// Functions in porter.c for use by subclasses and optimized functions
bool porter_grab_buffers(porter_t* self);
void porter_free_buffers(porter_t* self);
void porter_record_delivery(porter_t* self, int dest, uint64_t delivered);
void porter_close_buffer(porter_t* self, int dest, area_t* area);
bool porter_try_send(porter_t* self, int dest);

// Functions in codata.c
bool porter_compress(porter_t* self, buffer_t* buffer, int dest);
void porter_decompress(porter_t* self, buffer_t* buffer, int source);
bool porter_make_codata(porter_t* self);
int porter_setup_codata(porter_t* self, size_t item_size, size_t max_items);
void porter_reset_codata(porter_t* self);
void porter_free_codata(porter_t* self);

// Error codes for porter_setup_codata()
enum {
  codata_error_PLAN = -3,
  codata_error_ALIGN,
  codata_error_ALLOC,
};

// Functions in packer.c
void porter_choose_packer(porter_codata_t* codata);

#endif
