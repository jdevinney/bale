// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: Standard compression and decompression functions.


#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "convey_codec.h"


/*** AVX2/BMI2 CODEC ***/

#if HAVE_AVX2 && HAVE_BMI2

#include <immintrin.h>
#define ALIGN32 __attribute__((aligned(32)))
#define ZERO _mm256_setzero_si256()
#define FULL _mm256_cmpeq_epi8(ZERO,ZERO)

typedef union avx2 {
  __m256i y;
  __m128i x[2];
  uint64_t w[4];
} avx2_t;

// Efficient compressor that extracts the bits under the mask[] and packs
// them.  The input consists of 4*n units, each taking up w 64-bit words;
// only 8-byte alignment is required.  The output (packed) must be 32-byte
// aligned.  The return value is the number of bytes written.
//
// The stream of units whose index is i mod 4 is packed into the ith 64-bit
// slot of the packed[] vectors.  The reason for using 32-byte vectors in
// this way, when the bit extract/deposit instructions (PEXT/PDEP) operate
// on words, is to accelerate the shifting and merging of extracted bits.

static size_t
avx2_squeeze(size_t w, const uint64_t mask[w],
             size_t n, const uint64_t data[4 * n * w],
             __m256i packed[])
{
  ALIGN32 const __m256i splat[2] = { ZERO, FULL };
  const __m128i sixty_four = _mm_cvtsi64_si128(64);
  ALIGN32 __m256i* p = packed;
  __m256i x = ZERO;
  uint64_t bits = 0;

  // Process units in groups of four
  for (; n > 0; n--, data += 4 * w)
    for (size_t i = 0; i < w; i++) {
      uint64_t m = mask[i];
      __m128i shift = _mm_cvtsi64_si128(bits);
      bits += _mm_popcnt_u64(m);
      uint64_t step = bits >> 6;
      // Somehow broadcast this bit to all 256 bits in a vector
      __m256i next = splat[step];
      // Squeeze out the unwanted bits from the data
      uint64_t y0 = _pext_u64(data[i + 0*w], m);
      uint64_t y1 = _pext_u64(data[i + 1*w], m);
      uint64_t y2 = _pext_u64(data[i + 2*w], m);
      uint64_t y3 = _pext_u64(data[i + 3*w], m);
      // Somehow build a vector from these scalars
      __m128i y01 = _mm_set_epi64x(y1, y0);
      __m128i y23 = _mm_set_epi64x(y3, y2);
      __m256i y = _mm256_inserti128_si256(_mm256_castsi128_si256(y01), y23, 1);
      // Merge in the new vector
      x = _mm256_or_si256(x, _mm256_sll_epi64(y, shift));
      *p = x;
      y = _mm256_srl_epi64(y, _mm_sub_epi16(sixty_four, shift));
      x = _mm256_or_si256(_mm256_and_si256(next, y), _mm256_andnot_si256(next, x));
      p += step;
      bits &= 63;
    }

  if (bits > 0)
    *(p++) = x;
  return (char*)p - (char*)packed;
}

// Efficient decompressor that inverts the operation performed by
// avx2_squeeze().  The bits dropped by avx2_squeeze() must be supplied, in
// their proper locations, in the array stamp[], and stamp[] must be zero
// under the mask[].  The array packed[] must be 32-byte aligned.  The
// return value is the number of bytes consumed.

static size_t
avx2_unsqueeze(size_t w, const uint64_t mask[w], const uint64_t stamp[w],
               size_t n, const __m256i packed[], uint64_t data[])
{
  ALIGN32 const __m256i splat[2] = { ZERO, FULL };
  const __m128i sixty_four = _mm_cvtsi64_si128(64);
  ALIGN32 const __m256i* p = packed;
  __m256i x = p[0];
  __m256i y = p[1];
  uint64_t bits = 0;

  // Process units in groups of four
  for (; n > 0; n--, data += 4 * w)
    for (size_t i = 0; i < w; i++) {
      uint64_t m = mask[i];
      uint64_t c = stamp[i];
      __m128i shift = _mm_cvtsi64_si128(bits);
      bits += _mm_popcnt_u64(m);
      uint64_t step = bits >> 6;
      __m256i next = splat[step];
      // Extract and pad the data bits
      __m256i z = _mm256_or_si256(_mm256_srl_epi64(x, shift),
                                  _mm256_sll_epi64(y, _mm_sub_epi16(sixty_four, shift)));
      __m128i z01 = _mm256_castsi256_si128(z);
      __m128i z23 = _mm256_extracti128_si256(z, 1);
      data[i + 0*w] = c | _pdep_u64(_mm_cvtsi128_si64(z01), m);
      data[i + 1*w] = c | _pdep_u64(_mm_extract_epi64(z01, 1), m);
      data[i + 2*w] = c | _pdep_u64(_mm_cvtsi128_si64(z23), m);
      data[i + 3*w] = c | _pdep_u64(_mm_extract_epi64(z23, 1), m);
      // Conditionally set x to y, and load y
      x = _mm256_or_si256(_mm256_and_si256(next, y), _mm256_andnot_si256(next, x));
      p += step;
      y = p[1];
      bits &= 63;
    }

  p += (bits > 0);
  return (char*)p - (char*)packed;
}

static inline size_t
standard_words(const convey_cargo_t* cargo)
{
  return (cargo->item_size + cargo->tag_size + 7) / 8;
}

static bool
avx2_plan(const convey_cargo_t* cargo, convey_layout_t* layout)
{
  size_t w = standard_words(cargo);
  // Need to store two items (mask and stamp)
  *layout = (convey_layout_t) {
    .align = 32, .stride = 8 * w, .slack = 3,
    .offset = 16 * (w & 1), .overrun = 64,
  };
  return true;
}

static void
find_pattern(size_t valid_size, size_t w, size_t n, uint64_t inbuf[],
             uint64_t mask[w], uint64_t stamp[w])
{
  const size_t n_vectors = (w + 3) / 4;
  avx2_t valid[n_vectors];
  valid[n_vectors - 1].y = ZERO;
  memset(valid, 0xFF, valid_size);

  // Bring the stride up to at least 3 words
  size_t stride = w;
  if (stride == 1) {
    if (n & 1)
      inbuf[n] = inbuf[n - 1];
    valid[0].w[1] = valid[0].w[0];
    n = (n + 1) / 2;
    stride = 2;
  }
  if (stride == 2) {
    if (n & 1) {
      inbuf[2*n + 0] = inbuf[2*n - 2];
      inbuf[2*n + 1] = inbuf[2*n - 1];
    }
    valid[0].x[1] = valid[0].x[0];
    n = (n + 1) / 2;
    stride = 4;
  }

  for (size_t j = 0; j < n_vectors; j++) {
    uint64_t* p = inbuf + 4 * j;
    __m256i x = ZERO, y = FULL;
    for (size_t i = 0; i < n; i++, p += stride) {
      __m256i z = _mm256_loadu_si256((__m256i*) p);
      x = _mm256_or_si256(x, z);
      y = _mm256_and_si256(y, z);
    }
    // Fold the vectors if necessary
    if (w <= 2) {
      x = _mm256_or_si256(x, _mm256_permute4x64_epi64(x, 0x4E));
      y = _mm256_and_si256(y, _mm256_permute4x64_epi64(y, 0x4E));
    }
    if (w == 1) {
      x = _mm256_or_si256(x, _mm256_permute4x64_epi64(x, 0xB1));
      y = _mm256_and_si256(y, _mm256_permute4x64_epi64(y, 0xB1));
    }
    // Erase invalid bits and combine results
    x = _mm256_and_si256(x, valid[j].y);
    y = _mm256_and_si256(y, valid[j].y);
    x = _mm256_xor_si256(x, y);
    // Now x is the mask and y is the stamp
    size_t words = (w - 4*j < 4) ? w - 4*j : 4;
    memcpy(mask + 4*j, &x, 8 * words);
    memcpy(stamp + 4*j, &y, 8 * words);
  }
}

static size_t
avx2_compress(const convey_cargo_t* cargo, void* link_data,
              size_t n_items, void* inbuf, size_t limit, void* outbuf)
{
  if (n_items < 3)
    return 0;

  // Find the constant and varying bits
  const size_t unit_size = cargo->item_size + cargo->tag_size;
  const size_t w = standard_words(cargo);
  uint64_t mask[2 * w];
  find_pattern(unit_size, w, n_items, inbuf, mask, mask + w);

  // Compute the compressed size
  uint64_t bits = 0;
  for (size_t i = 0; i < w; i++)
    bits += _mm_popcnt_u64(mask[i]);
  size_t n = (n_items + 3) / 4;
  size_t packed_size = ((n * bits + 63) >> 6) * 32;
  size_t n_bytes = 16 * w + packed_size;
  if (n_bytes > limit)
    return 0;
  // Don't compress unless we save more than 1/8 of the bytes
  size_t original = n_items * cargo->packet_size;
  if (n_bytes >= original - (original >> 3))
    return 0;

  // Compress into the outbuf
  memcpy(outbuf, mask, 16 * w);
  uint64_t* packed = (uint64_t*) outbuf + 2 * w;
  assert((uintptr_t)packed % 32 == 0);
  size_t b = avx2_squeeze(w, mask, n, inbuf, (__m256i*)packed);
  assert(b == packed_size);

  return n_bytes;
}

static void
avx2_decompress(const convey_cargo_t* cargo, void* link_data,
                size_t n_items, size_t n_bytes, void* inbuf, void* outbuf)
{
  const size_t w = standard_words(cargo);
  uint64_t* mask = inbuf;
  uint64_t* packed = mask + 2 * w;
  assert((uintptr_t)packed % 32 == 0);
  avx2_unsqueeze(w, mask, mask + w, (n_items + 3) / 4, (__m256i*)packed, outbuf);
}

const convey_codec_t convey_standard_codec = {
  .name = "standard-AVX2",
  .plan = &avx2_plan,
  .compress = &avx2_compress,
  .decompress = &avx2_decompress,
};

#elif HAVE_SSSE3

/*** SSE2/SSSE3 CODEC ***/

#include <immintrin.h>
#include "bolite.h"
#define ZERO _mm_setzero_si128()
#define FULL _mm_cmpeq_epi8(ZERO,ZERO)

typedef union sse2 {
  __m128i x;
  uint64_t w[2];
  uint8_t b[16];
} sse2_t;

static size_t
ssse3_squeeze(size_t stride, const __m128i extract[], const uint16_t mask[],
              size_t n, const uint8_t data[], uint8_t packed[])
{
  // Number of SSE vectors per item
  const size_t v = (stride + 15) >> 4;

  uint8_t* target = packed;
  for (size_t j = 0; j < v; j++) {
    if (mask[j] == 0)
      continue;
    const size_t count = _popcnt(mask[j]);
    const __m128i perm = extract[j];
    const uint8_t* source = data + 16 * j;

    // compiler should unroll if beneficial
    for (size_t i = n; i > 0; i--, source += stride, target += count) {
      __m128i x = _mm_loadu_si128((__m128i*)source);
      _mm_storeu_si128((__m128i*)target, _mm_shuffle_epi8(x, perm));
    }
  }

  return target - packed;
}

// stamp[] need not be properly aligned
static size_t
ssse3_unsqueeze(size_t stride, const __m128i deposit[], const __m128i stamp[],
                const uint16_t mask[], size_t n, const uint8_t packed[], uint8_t data[])
{
  size_t v = (stride + 15) >> 4;
  size_t total = 0;

  // Iterate over the vectors in the order 1, 2, ..., v-1, 0
  const uint8_t* source = packed + n * _popcnt(mask[0]);
  for (size_t j = 1; v > 0; v--, j++) {
    if (v == 1) {
      j = 0;
      total = source - packed;
      source = packed;
    }

    const size_t count = _popcnt(mask[j]);
    const __m128i perm = deposit[j];
    const __m128i fixed = _mm_loadu_si128(stamp + j);
    uint8_t* target = data + 16 * j;

    for (size_t i = n; i > 0; i--, source += count, target += stride) {
      __m128i x = _mm_loadu_si128((__m128i*)source);
      x = _mm_or_si128(fixed, _mm_shuffle_epi8(x, perm));
      _mm_storeu_si128((__m128i*)target, x);
    }
  }

  return total;
}

static inline size_t
standard_vectors(const convey_cargo_t* cargo)
{
  return (cargo->item_size + cargo->tag_size + 15) / 16;
}

static bool
ssse3_plan(const convey_cargo_t* cargo, convey_layout_t* layout)
{
  size_t v = standard_vectors(cargo);
  size_t packet_quads = (cargo->item_size + cargo->tag_size + 3) / 4;
  *layout = (convey_layout_t) {
    .align = 8, .stride = 4 * MAX(packet_quads, 2), .slack = 1,
    .offset = 2 * (v & 1), .overrun = 16,
  };
  return true;
}

// Given a packed array inbuf[] of n items, each with valid_size valid
// bytes but starting on a 4-byte boundary, determine which bytes are
// constant.  Join items in pairs if valid_size <= 8.  Write bitmasks of the
// nonconstant bytes into mask[], and write the constant bytes themselves
// into stamp[].
static void
find_pattern(size_t valid_size, size_t v, size_t n, uint32_t inbuf[],
             uint16_t mask[v], __m128i stamp[v])
{
  sse2_t valid[v];
  valid[v - 1].x = ZERO;
  memset(valid, 0xFF, valid_size);

  // Optimize small items by working with pairs of items
  size_t quads = (valid_size + 3) / 4;
  if (quads <= 2) {
    if (n & 1) {
      assert(n > 1);
      inbuf[2*n + 0] = inbuf[2*n - 4];
      inbuf[2*n + 1] = inbuf[2*n - 3];
    }
    valid[0].w[1] = valid[0].w[0];
    n = (n + 1) / 2;
    quads = 4;
  }

  for (size_t j = 0; j < v; j++) {
    uint32_t* p = inbuf + 4 * j;
    __m128i x = ZERO, y = FULL;
    for (size_t i = 0; i < n; i++, p += quads) {
      __m128i z = _mm_loadu_si128((__m128i*) p);
      x = _mm_or_si128(x, z);
      y = _mm_and_si128(y, z);
    }
    // Combine results and erase invalid bits
    __m128i z = _mm_cmpeq_epi8(x, y);
    __m128i ok = valid[j].x;
    mask[j] = _mm_movemask_epi8(_mm_andnot_si128(z, ok));
    stamp[j] = _mm_and_si128(_mm_and_si128(z, y), ok);
  }
}

static void
compute_extract(size_t v, const uint16_t mask[v], sse2_t extract[v])
{
  for (size_t j = 0; j < v; j++) {
    uint64_t m = mask[j], c = 0;
    extract[j].x = FULL;
    for (int i = 0; i < 15; i++, m >>= 1) {
      extract[j].b[c] = i;
      c += m & 1;
    }
    extract[j].b[c] = (m << 4) - 1;  // same as (m & 1) ? 15 : -1
  }
}

static size_t
ssse3_compress(const convey_cargo_t* cargo, void* link_data,
               size_t n_items, void* inbuf, size_t limit, void* outbuf)
{
  if (n_items < 3)
    return 0;

  // Find the constant and varying bits
  const size_t unit_size = cargo->item_size + cargo->tag_size;
  const size_t v = standard_vectors(cargo);
  uint16_t mask[v];
  sse2_t stamp[v];
  find_pattern(unit_size, v, n_items, inbuf, mask, &stamp[0].x);
  const size_t s = (unit_size > 8) ? (unit_size + 3) & ~(size_t)3 : 16;
  const size_t n = (unit_size > 8) ? n_items : (n_items + 1) / 2;

  // Compute the compressed size
  uint64_t bytes = 0;
  for (size_t i = 0; i < v; i++)
    bytes += _popcnt(mask[i]);
  size_t n_bytes = 18 * v + n * bytes;
  if (n_bytes > limit)
    return 0;
  // Don't compress unless we save more than 1/8 of the bytes
  size_t original = n_items * cargo->packet_size;
  if (n_bytes >= original - (original >> 3))
    return 0;

  // Compress into the outbuf
  uint8_t* q = outbuf;
  memcpy(q, mask, 2 * v);
  memcpy(q + 2 * v, stamp, 16 * v);
  compute_extract(v, mask, stamp);
  size_t c = ssse3_squeeze(s, &stamp[0].x, mask, n, inbuf, q + 18 * v);
  assert(c == n * bytes);

  return n_bytes;
}

static void
ssse3_decompress(const convey_cargo_t* cargo, void* link_data,
                 size_t n_items, size_t n_bytes, void* inbuf, void* outbuf)
{
  const size_t unit_size = cargo->item_size + cargo->tag_size;
  const size_t v = standard_vectors(cargo);
  uint8_t* p = inbuf;
  uint16_t* mask = (uint16_t*) p;
  __m128i* stamp = (__m128i*) (p + 2 * v);
  uint8_t* packed = p + 18 * v;

  // Construct the 'deposit' control vectors
  sse2_t deposit[v];
  for (size_t j = 0; j < v; j++) {
    uint64_t m = mask[j], c = 0;
    for (int i = 0; i < 16; i++, m >>= 1) {
      uint64_t bit = m & 1;
      deposit[j].b[i] = c | (bit - 1);
      c += bit;
    }
  }
  const size_t s = (unit_size > 8) ? (unit_size + 3) & ~(size_t)3 : 16;
  const size_t n = (unit_size > 8) ? n_items : (n_items + 1) / 2;
  ssse3_unsqueeze(s, &deposit[0].x, stamp, mask, n, packed, outbuf);
}

const convey_codec_t convey_standard_codec = {
  .name = "standard-SSSE3",
  .plan = &ssse3_plan,
  .compress = &ssse3_compress,
  .decompress = &ssse3_decompress,
};

#else

const convey_codec_t convey_standard_codec = {
  .name = "unavailable",
};

#endif
