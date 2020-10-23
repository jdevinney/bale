// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: Compression support for porters.


#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "convey.h"
#include "porter_impl.h"


/*** Compression and Decompression ***/

static size_t
codec_padding(convey_layout_t* layout)
{
  size_t align = layout->align;
  size_t offset = layout->offset;
  size_t used = sizeof(buffer_t) & (align - 1);
  return ((offset < used) ? align : 0) + offset - used;
}

bool
porter_compress(porter_t* self, buffer_t* buffer, int dest)
{
  porter_codata_t* codata = self->codata;
  size_t n_items = buffer->limit / self->packet_bytes;
  assert(n_items > 0);

  // As a first implementation, compress tags together with items

  // Step 1: Extract the items into the workspace
  if (codata->layout.stride > self->packet_bytes)
    codata->unpack(codata, n_items, buffer->data);
  else
    memcpy(codata->work, buffer->data, buffer->limit);

  // Step 2: Calculate where the compressed data should begin
  size_t skip = codec_padding(&codata->layout);

  // Step 3: Apply the compressor
  const size_t capacity = self->buffer_bytes - sizeof(buffer_t);
  const size_t reserved = skip + codata->layout.overrun;
  size_t n_bytes = 0;
  if (capacity >= reserved)
    n_bytes = self->codec->compress
      (&codata->cargo, codata->compressors[dest], n_items, codata->work,
       capacity - reserved, buffer->data + skip);

  // Step 4: Update the header if compression occurred
  if (n_bytes > 0) {
    buffer->limit = n_bytes + skip;
    buffer->n_items = n_items;
  }

  return (n_bytes > 0);
}

void
porter_decompress(porter_t* self, buffer_t* buffer, int source)
{
  porter_codata_t* codata = self->codata;
  size_t skip = codec_padding(&codata->layout);
  size_t n_items = buffer->n_items;
  size_t n_bytes = buffer->limit - skip;

  self->codec->decompress
    (&codata->cargo, codata->decompressors[source], n_items, n_bytes,
     buffer->data + skip, codata->work);

  size_t packet_bytes = self->packet_bytes;
  if (codata->layout.stride > packet_bytes)
    codata->repack(codata, n_items, buffer->data);
  else
    memcpy(buffer->data, codata->work, n_items * packet_bytes);

  buffer->limit = n_items * self->packet_bytes;
  buffer->n_items = 0;  
}


/*** Compression Setup ***/

bool
porter_make_codata(porter_t* self)
{
  porter_codata_t* codata = calloc(1, sizeof(porter_codata_t));
  bool ok = (codata != NULL);
  if (ok) {
    const int n = self->n_ranks;
    codata->work = NULL;
    codata->compressors = malloc(n * sizeof(void*));
    codata->decompressors = malloc(n * sizeof(void*));
    ok = (codata->compressors && codata->decompressors);
    if (ok)
      for (int i = 0; i < n; i++) {
        codata->compressors[i] = NULL;
        codata->decompressors[i] = NULL;
      }
  }
  self->codata = codata;
  return ok;
}

int
porter_setup_codata(porter_t* self, size_t item_size, size_t max_items)
{
  porter_codata_t* codata = self->codata;
  if (item_size == codata->cargo.item_size)
    return convey_OK;

  // If the item size changes, reset the codec
  if (codata->cargo.item_size != 0)
    porter_set_codec(self, self->codec, codata->cargo.arg);
  assert(codata->cargo.item_size == 0);

  const convey_codec_t* codec = self->codec;
  if (codec == NULL)
    return convey_FAIL;

  const size_t tag_size = self->tag_bytes;
  convey_layout_t* layout = &codata->layout;
  convey_cargo_t _cargo = {
    .tag_size = tag_size, .item_size = item_size,
    .packet_size = tag_size ? (item_size + 2*tag_size - 1) & ~(tag_size - 1) : item_size,
    .max_items = max_items, .arg = codata->cargo.arg,
  };
  bool ok = codec->plan(&_cargo, layout);
  if (!ok)
    return convey_FAIL;
  ok = (layout->stride >= tag_size + item_size) && (layout->align >= 8) &&
    (layout->align <= CONVEY_MAX_ALIGN) && !(layout->align & (layout->align - 1)) &&
    (layout->offset < layout->align);
  if (!ok)
    return codata_error_PLAN;
  if (self->buffer_align < layout->align)
    return codata_error_ALIGN;

  // The next line ensures that any allocations will be released by the
  // next call to porter_setup_codata().
  codata->cargo = _cargo;

  const int n = self->n_ranks;
  if (codec->link)
    for (int i = 0; ok && i < n; i++) {
      codata->compressors[i] = codec->link(&_cargo, false);
      codata->decompressors[i] = codec->link(&_cargo, true);
    }
  codata->work_bytes = layout->stride * (max_items + layout->slack);
  if (posix_memalign(&codata->work, layout->align, codata->work_bytes) != 0)
    return codata_error_ALLOC;
  memset(codata->work, 0, codata->work_bytes);

  porter_choose_packer(codata);
  return convey_OK;
}

void
porter_reset_codata(porter_t* self)
{
  porter_codata_t* codata = self->codata;
  if (codata->cargo.item_size > 0) {
    if (self->codec->unlink)
      for (int i = self->n_ranks - 1; i >= 0; i--) {
        self->codec->unlink(&codata->cargo, codata->compressors[i]);
        self->codec->unlink(&codata->cargo, codata->decompressors[i]);
        codata->compressors[i] = NULL;
        codata->decompressors[i] = NULL;
      }
    free(codata->work);
    codata->work = NULL;
    codata->cargo.item_size = 0;
  }
}

void
porter_free_codata(porter_t* self)
{
  porter_codata_t* codata = self->codata;
  if (codata) {
    free(codata->work);
    free(codata->decompressors);
    free(codata->compressors);
    free(codata);
  }
}
