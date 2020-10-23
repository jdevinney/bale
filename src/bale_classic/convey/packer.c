// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: Functions for reformatting buffers to satisfy a codec.


#include <string.h>
#include "porter_impl.h"

#if HAVE_SSSE3
#include <emmintrin.h>
#include <tmmintrin.h>

static void
ssse3_unpack_12to16(porter_codata_t* self, size_t n_items, const void* area)
{
  const __m128i* source = area;
  __m128i* target = self->work;
  const __m128i mask = _mm_set_epi32(0, ~0, ~0, ~0);
  const __m128i z = _mm_setzero_si128();
  
  for (size_t n = n_items / 4; n > 0; n--, source += 3, target += 4) {
    // 3 input vectors --> 4 output vectors, all aligned
    __m128i a = _mm_load_si128(source + 0);
    __m128i b = _mm_load_si128(source + 1);
    __m128i c = _mm_load_si128(source + 2);
    _mm_store_si128(target + 0, _mm_and_si128(mask, a));
    _mm_store_si128(target + 1, _mm_and_si128(mask, _mm_alignr_epi8(b, a, 12)));
    _mm_store_si128(target + 2, _mm_and_si128(mask, _mm_alignr_epi8(c, b, 8)));
    _mm_store_si128(target + 3, _mm_alignr_epi8(z, c, 4));
  }

  for (size_t i = 0; i < n_items % 4; i++) {
    char* bytes = (char*) (target + i);
    memcpy(bytes, (char*)source + 12*i, 12);
    memset(bytes + 12, 0, 4);
  }
}

static void
ssse3_repack_16to12(porter_codata_t* self, size_t n_items, void* area)
{
  const __m128i* source = self->work;
  __m128i* target = area;
  const __m128i z = _mm_setzero_si128();

  for (size_t n = n_items / 4; n > 0; n--, source += 4, target += 3) {
    // 4 input vectors --> 3 output vectors, all aligned
    __m128i a = _mm_load_si128(source + 0);
    __m128i b = _mm_load_si128(source + 1);
    __m128i c = _mm_load_si128(source + 2);
    __m128i d = _mm_load_si128(source + 3);
    _mm_store_si128(target + 0, _mm_or_si128(_mm_alignr_epi8(b, z, 4), a));
    b = _mm_alignr_epi8(z, b, 4);
    _mm_store_si128(target + 1, _mm_or_si128(_mm_alignr_epi8(c, z, 8), b));
    c = _mm_alignr_epi8(z, c, 8);
    _mm_store_si128(target + 2, _mm_or_si128(_mm_alignr_epi8(d, z, 12), c));
  }

  for (size_t i = 0; i < n_items % 4; i++)
    memcpy((char*)target + 12*i, source + i, 12);
}

static void
ssse3_unpack_20to24(porter_codata_t* self, size_t n_items, const void* area)
{
  const __m128i* source = area;
  __m128i* target = self->work;
  const __m128i m3 = _mm_set_epi32(0, ~0, ~0, ~0);
  const __m128i z = _mm_setzero_si128();

  for (size_t n = n_items / 4; n > 0; n--, source += 5, target += 6) {
    // 5 input vectors --> 6 output vectors, all aligned
    __m128i x0 = _mm_load_si128(source + 0);  // a3 a2 a1 a0
    __m128i x1 = _mm_load_si128(source + 1);  // b2 b1 b0 a4
    __m128i x2 = _mm_load_si128(source + 2);  // c1 c0 b4 b3
    __m128i x3 = _mm_load_si128(source + 3);  // d0 c4 c3 c2
    __m128i x4 = _mm_load_si128(source + 4);  // d4 d3 d2 d1
    // a3 a2 a1 a0
    // b1 b0 __ a4
    // __ b4 b3 b2
    _mm_store_si128(target + 0, x0);
    _mm_store_si128(target + 1, _mm_shuffle_epi32(_mm_and_si128(x1, m3), 0x9C));
    _mm_store_si128(target + 2, _mm_and_si128(m3, _mm_alignr_epi8(x2, x1, 12)));
    // c3 c2 c1 c0
    // d1 d0 __ c4
    // __ d4 d3 d2
    __m128i t = _mm_alignr_epi8(x4, x3, 8);
    _mm_store_si128(target + 3, _mm_alignr_epi8(x3, x2, 8));
    _mm_store_si128(target + 4, _mm_shuffle_epi32(_mm_and_si128(t, m3), 0x9C));
    _mm_store_si128(target + 5, _mm_alignr_epi8(z, x4, 4));
  }

  for (size_t i = 0; i < n_items % 4; i++) {
    char* bytes = (char*)target + 24*i;
    memcpy(bytes, (char*)source + 20*i, 20);
    memset(bytes + 20, 0, 4);
  }
}

static void
ssse3_repack_24to20(porter_codata_t* self, size_t n_items, void* area)
{
  const __m128i* source = self->work;
  __m128i* target = area;
  const __m128i m1 = _mm_set_epi32(0, 0, 0, ~0);
  const __m128i m3 = _mm_set_epi32(0, ~0, ~0, ~0);
  const __m128i z = _mm_setzero_si128();

  for (size_t n = n_items / 4; n > 0; n--, source += 6, target += 5) {
    __m128i y0 = _mm_load_si128(source + 0);
    __m128i y1 = _mm_load_si128(source + 1);
    __m128i y2 = _mm_load_si128(source + 2);
    __m128i y3 = _mm_load_si128(source + 3);
    __m128i y4 = _mm_load_si128(source + 4);
    __m128i y5 = _mm_load_si128(source + 5);
    __m128i s, t;
    _mm_store_si128(target + 0, y0);
    t = _mm_alignr_epi8(y2, y1, 4);  // b2 b1 b0 ??
    _mm_store_si128(target + 1, _mm_or_si128(_mm_and_si128(m1, y1), _mm_andnot_si128(m1, t)));
    _mm_store_si128(target + 2, _mm_alignr_epi8(y3, _mm_alignr_epi8(y2, z, 12), 8));
    t = _mm_alignr_epi8(y4, y3, 8);  // ?? c4 c3 c2
    s = _mm_alignr_epi8(y4, z, 12);  // d0 ?? ?? ??
    _mm_store_si128(target + 3, _mm_or_si128(_mm_and_si128(m3, t), _mm_andnot_si128(m3, s)));
    _mm_store_si128(target + 4, _mm_alignr_epi8(y5, y4, 12));
  }

  for (size_t i = 0; i < n_items % 4; i++)
    memcpy((char*)target + 20*i, (char*)source + 24*i, 20);
}

#endif

static void
unpack_28to32(porter_codata_t* self, size_t n_items, const void* area)
{
  const uint32_t* source = area;
  uint32_t* target = self->work;
  for (size_t n = n_items; n > 0; n--, source += 7, target += 8) {
    memcpy(target, source, 28);
    target[7] = 0;
  }
}

static void
repack_32to28(porter_codata_t* self, size_t n_items, void* area)
{
  const uint32_t* source = self->work;
  uint32_t* target = area;
  for (size_t n = n_items; n > 0; n--, source += 8, target += 7)
    memcpy(target, source, 28);
}

static void
unpack_36to40(porter_codata_t* self, size_t n_items, const void* area)
{
  const uint32_t* source = area;
  uint32_t* target = self->work;
  for (size_t n = n_items; n > 0; n--, source += 9, target += 10) {
    memcpy(target, source, 36);
    target[9] = 0;
  }
}

static void
repack_40to36(porter_codata_t* self, size_t n_items, void* area)
{
  const uint32_t* source = self->work;
  uint32_t* target = area;
  for (size_t n = n_items; n > 0; n--, source += 10, target += 9)
    memcpy(target, source, 36);
}


/*** Generic Pack/Unpack Functions ***/

static void
generic_pack(char* target, const char* source, size_t n_items,
             size_t t, size_t s, size_t p)
{
  for (; n_items > 0; n_items--, target += t, source += s)
    memcpy(target, source, p);
}

static void
generic_unpack(porter_codata_t* self, size_t n_items, const void* area)
{
  size_t t = self->layout.stride;
  memset(self->work, 0, n_items * t);
  size_t p = self->cargo.item_size + self->cargo.tag_size;
  size_t s = self->cargo.packet_size;
  generic_pack(self->work, area, n_items, t, s, p);
}

static void
generic_repack(porter_codata_t* self, size_t n_items, void* area)
{
  size_t s = self->layout.stride;
  size_t p = self->cargo.item_size + self->cargo.tag_size;
  size_t t = self->cargo.packet_size;
  generic_pack(area, self->work, n_items, t, s, p);
}


/*** Codelet Selector ***/

void
porter_choose_packer(porter_codata_t* codata)
{
  size_t stride = codata->layout.stride;
  size_t align = codata->layout.align;
  size_t packet = codata->cargo.packet_size;

#if HAVE_SSSE3
  if (stride == 24 && align >= 16 && packet == 20) {
    codata->unpack = &ssse3_unpack_20to24;
    codata->repack = &ssse3_repack_24to20;
    return;
  }
  if (stride == 16 && align >= 16 && packet == 12) {
    codata->unpack = &ssse3_unpack_12to16;
    codata->repack = &ssse3_repack_16to12;
    return;
  }
#endif
  if (stride == 32 && packet == 28) {
    codata->unpack = &unpack_28to32;
    codata->repack = &repack_32to28;
    return;
  }
  if (stride == 40 && packet == 36) {
    codata->unpack = &unpack_36to40;
    codata->repack = &repack_40to36;
    return;
  }

  codata->unpack = &generic_unpack;
  codata->repack = &generic_repack;
}
