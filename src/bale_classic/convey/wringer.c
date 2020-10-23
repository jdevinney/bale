// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if HAVE_CONFIG_H
# include "config.h"
#endif
#include "bolite.h"
#include "convey_codec.h"


static void*
grab(size_t n_bytes, size_t align)
{
  void* ptr = NULL;
  if (align < sizeof(void*))
    align = sizeof(void*);
  int err = posix_memalign(&ptr, align, n_bytes);
  if (err || ptr == NULL)
    return NULL;
  memset(ptr, 0, n_bytes);
  return ptr;
}

static void
print_bits(const char* label, size_t w, const uint64_t bits[w])
{
  printf("%s:", label);
  uint8_t bytes[8 * w];
  memcpy(bytes, bits, 8 * w);
  for (size_t i = 0; i < 8 * w; i++)
    printf(" %02x", bytes[i]);
  printf("\n");
}

int
main(int argc, char* argv[])
{
  // Parse command line and prepare PRNG
  if (argc < 3) {
    fprintf(stderr, "usage: %s item_size n_items [varying_bits [seed]]\n", argv[0]);
    return 2;
  }
  size_t s = atoi(argv[1]);
  size_t n = atoi(argv[2]);
  uint64_t b = (argc > 3) ? atoi(argv[3]) : 4 * (s + 4);
  if (s >> 16 || n == 0 || n * s >> 20 || b == 0 || b > 8 * (s + 4)) {
    fprintf(stderr, "%s: bad arguments\n", argv[0]);
    return 2;
  }
  uint64_t seed = (argc > 4) ? strtoul(argv[4], NULL, 10) : time(NULL);
  brand_t _prng;
  brand_init(&_prng, seed);

  // Prepare and check the layout
  const convey_codec_t* codec = &convey_standard_codec;
  convey_cargo_t cargo = { .tag_size = 4, .item_size = s,
                           .packet_size = (s + 7) & ~(size_t)3, .max_items = n };
  convey_layout_t layout;
  bool ok = codec->plan(&cargo, &layout);
  if (!ok) {
    fprintf(stderr, "%s: parameters are not supported\n", argv[0]);
    return 2;
  }
  if (layout.stride & 7) {
    fprintf(stderr, "%s: stride (%zu) is not a multiple of 8\n",
            argv[0], layout.stride);
    return 2;
  }
  printf("layout: %zu %zu %zu %zu %zu\n",
         layout.stride, layout.align, layout.offset, layout.overrun, layout.slack);
  size_t w = layout.stride >> 3;
  size_t q = 1 + ((s + 3) >> 2);
  size_t capacity = n * q * 4;
  if (layout.offset + layout.overrun > capacity) {
    printf("insufficient capacity\n");
    return 0;
  }
  size_t limit = capacity - (layout.offset + layout.overrun);

  // Allocate buffers
  uint64_t* loose = grab((n + layout.slack) * 8 * w, layout.align);
  uint64_t* restored = grab((n + layout.slack) * 8 * w, layout.align);
  uint8_t* tight = grab(capacity, layout.align);
  if (!(loose && restored && tight)) {
    fprintf(stderr, "failed to allocate %zu-byte aligned memory\n", layout.align);
    return 1;
  }

  // Prepare the pattern of fixed and varying bits
  uint64_t valid[w];
  uint64_t mask[w];
  uint64_t stamp[w];

  memset(valid, 0, 8 * w);
  memset(valid, 0xFF, s + 4);

  memset(mask, 0, 8 * w);
#if HAVE_AVX2 && HAVE_BMI2
  int64_t m = 8 * (s + 4);
  for (int64_t i = 0, k = b; i < m; i++) {
    uint64_t bit = (dbrand(&_prng) * (m - i) < k);
    mask[i >> 6] |= bit << (i & 63);
    k -= bit;
  }
#else
  int64_t m = s + 4;
  for (int64_t i = 0, k = (b + 3) / 8; i < m; i++) {
    uint64_t bit = (dbrand(&_prng) * (m - i) < k);
    mask[i >> 3] |= (-bit & UINT64_C(0xFF)) << 8 * (i & 7);
    k -= bit;
  }
#endif

  for (size_t i = 0; i < w; i++)
    stamp[i] = brand(&_prng) & valid[i] & ~mask[i];

  printf("seed = %" PRIu64 "\n", seed);
  print_bits("valid", w, valid);
  print_bits(" mask", w, mask);
  print_bits("stamp", w, stamp);

  // Set up the data, filling the extra space with garbage
  for (size_t j = 0; j < n; j++)
    for (size_t i = 0; i < w; i++)
      loose[w * j + i] = (brand(&_prng) & mask[i]) | stamp[i];
  for (size_t j = n; j < n + layout.slack; j++)
    for (size_t i = 0; i < w; i++)
      loose[w * j + i] = brand(&_prng);

  // Perform the compression and decompression
  size_t n_bytes =
    codec->compress(&cargo, NULL, n, loose, limit, tight + layout.offset);
  if (n_bytes == 0) {
    printf("compression was not done\n");
    return 0;
  }
  codec->decompress(&cargo, NULL, n, n_bytes, tight + layout.offset, restored);

  // Check the result
  for (size_t j = 0; j < n; j++)
    for (size_t i = 0; i < w; i++) {
      uint64_t result = restored[w * j + i] & valid[i];
      if (result != loose[w * j + i]) {
        printf("failure at word %zu of item %zu\n", i, j);
        print_bits("right", 1, &loose[w * j + i]);
        print_bits("wrong", 1, &result);
        return 1;
      }
    }
  printf("ok\n");
  return 0;
}
