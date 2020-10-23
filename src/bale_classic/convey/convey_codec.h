// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


/** \file convey_codec.h
 * The auxiliary API for buffer compression.
 *
 * For compression of remote transfers to be active, several things
 * must happen:
 *  - The conveyor constructor must be given the \c convey_opt_COMPRESS option.
 *  - The chosen type of conveyor must support compression.  At present, only
 *    tensor conveyors support it.
 *  - If \c convey_set_codec() has been called on the conveyor, then its most
 *    recent call must have supplied a non-NULL [codec](@ref convey_codec).
 *    (If it has not been called, then the \c convey_standard_codec will be used.)
 *  - The codec's \c plan() method, called as part of \c convey_begin(), must
 *    have returned \c true and filled in a valid [layout](@ref convey_layout).
 *
 * If these conditions are satisfied, then buffer compression will be
 * attempted on non-local transfers.  But the codec's \c compress() method
 * may still refuse to compress a buffer for any reason, in which case
 * the buffer will be sent in its original form.  (The usual reason for not
 * compressing a buffer is that not enough savings can be achieved.)
 */


#ifndef CONVEY_CODEC_H
#define CONVEY_CODEC_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/// The structure that a conveyor provides to a codec to describe the items
/// and routing tags being transmitted.
typedef struct convey_cargo {
  size_t tag_size;     ///< bytes per routing tag: a multiple of 4
  size_t item_size;    ///< bytes per item
  size_t packet_size;  ///< total size of tag + item + padding
  size_t max_items;    ///< maximum number of items per buffer
  void* arg;           ///< argument given to convey_set_codec()
} convey_cargo_t;

// (In the elastic case, the plan would be to separate tags+headers
// from payload units and compress them separately.)

/// A structure that records the way a codec wants data to be laid out
/// in memory.
typedef struct convey_layout {
  /// The desired alignment of both uncompressed and compressed buffers.
  /// It must be a power of two no greater than \c CONVEY_MAX_ALIGN and
  /// no less than 8.
  size_t align;
  /// The desired stride between packets in uncompressed buffers.  It must
  /// satisfy stride >= tag_size + item_size.  Within each chunk of
  /// 'stride' bytes, the first tag_size + item_size bytes hold valid data
  /// and the remaining bytes, if any, may be garbage.
  size_t stride;
  /// Allow the compressor to read and modify, and the decompressor to
  /// write, \a slack extra packets beyond the valid packets (\a slack
  /// times \a stride additional bytes) in uncompressed buffers.
  size_t slack;
  /// The desired offset from a multiple of \c align at which each
  /// compressed buffer should begin.  This feature allows the compressor
  /// to prepend a header whose size is not a multiple of \c align.
  size_t offset;
  /// Allow the compressor to write, and the decompressor to read and
  /// modify, \a overrun bytes beyond the end of of the compressed data.
  size_t overrun;
} convey_layout_t;

/** A structure that describes compression/decompression functions.
 *
 * An application that wants to compress the data being transmitted by a
 * conveyor can supply its own codec by passing a structure of this type
 * to \c convey_set_codec.
 */
typedef struct convey_codec {
  /// Name of this codec, used only in diagnostic messages.
  const char* name;

  /// Fill in \a *layout with the appropriate values for this codec and
  /// \a cargo.  Return \c true if everything is OK or \c false if this
  /// cargo specification (in particular, the combination of item size
  /// and tag size) is not supported.
  bool (*plan)(const convey_cargo_t* cargo, convey_layout_t* layout);

  /// Create compression or decompression data for a single link.  This
  /// method may be NULL, in which case all link data will be NULL.
  void* (*link)(const convey_cargo_t* cargo, bool decompress);

  /// Compress the \a inbuf, which is an array of \a n_items slots of size
  /// \c stride, each holding a packet of \a cargo->tag_size + 
  /// \a cargo->item_size bytes.  Write the result to the (disjoint)
  /// \a outbuf, which holds \a limit bytes, and return the number of bytes
  /// written.  Modification of \a inbuf is allowed.  The special return
  /// value 0 means that compression was not done, perhaps because the
  /// result would have overrun the \a outbuf.  The layout produced by the
  /// \c plan method governs the alignment and padding of both \a inbuf and
  /// \a outbuf.
  size_t (*compress)(const convey_cargo_t* cargo, void* link_data,
                     size_t n_items, void* inbuf, size_t limit, void* outbuf);

  /// Decompress the \a inbuf, whose length is \a n_bytes bytes, and place
  /// the resulting \a n_items items in the (disjoint) \a outbuf, an array
  /// of slots of size \c stride each holding a packet of \a cargo->tag_size
  /// + \a cargo->item_size bytes.  Modification of \a inbuf is allowed.
  /// The layout produced by the \c plan method governs the alignment and
  /// padding of \a inbuf and \a outbuf.
  void (*decompress)(const convey_cargo_t* cargo, void* link_data,
                     size_t n_items, size_t n_bytes, void* inbuf, void* outbuf);

  /// Release link data; may be NULL if \c link is NULL.
  void (*unlink)(const convey_cargo_t* cargo, void* link_data);
} convey_codec_t;

/** Standard codec that squeezes out constant bits or bytes.
 *
 * In most cases this codec need not be referenced.  It suffices to pass
 * the \c convey_COMPRESS option to the conveyor constructor; if the
 * conveyor type supports compression, then the standard codec will be used
 * automatically.  But an application may wish to build its own codec on
 * top of the functions in the standard codec, or it may wish to switch
 * between the standard codec and a specialized codec.
 */
extern const convey_codec_t convey_standard_codec;

#endif
