// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// SUMMARY: This header is "bitops lite".  It provides _popcnt, _leadz,
// _trailz, and _maskr for 64-bit integers, _divbymul32 for 32-bit integers,
// and _prefetch_x.  Except for _divbymul32 and _prefetch_x, they are not
// performance-critical, though a mild effort is made to optimize them.
// We also need brand, so this is copied from bitops.
//
// *** If common.h will be included, it must be included first. ***

#ifndef CONVEY_BOLITE_H
#define CONVEY_BOLITE_H

#ifdef __GNUC__
#define BOLITE_NOWARN __attribute__ ((unused))
#else
#define BOLITE_NOWARN
#endif
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))


#if HAVE_MPP_UTIL

// In this case we need bitops.h anyway (for actimer), so we must not
// define our own versions of the standard operations.
#include "bitops.h"

#else

#include <stdint.h>


/*** Bit operations ***/

#if defined(_CRAYC) || defined(__cray__)

#include <intrinsics.h>
static inline int64_t _trailz(uint64_t x)
{
  return _popcnt((x-1) & ~x);
}

#else

#define _maskr(x) (((x) == 0) ? UINT64_C(0) : ~UINT64_C(0) >> (64-(x)))

#ifdef __GNUC__

#define _popcnt __builtin_popcountll
#define _leadz(x) (((x) == 0) ? 64 : __builtin_clzll(x))
#define _trailz(x) (((x) == 0) ? 64 : __builtin_ctzll(x))

#else

static inline int64_t _popcnt(uint64_t x)
{
  x -= (x >> 1) & UINT64_C(0x5555555555555555);
  x = ((x >> 2) & UINT64_C(0x3333333333333333)) +
    (x & UINT64_C(0x3333333333333333));
  x = ((x >> 4) + x) & UINT64_C(0x0F0F0F0F0F0F0F0F);
  x *= UINT64_C(0x0101010101010101);
  return (x >> 56);
}

static inline int64_t _leadz (uint64_t x)
{
  x |= x >>  1;
  x |= x >>  2;
  x |= x >>  4;
  x |= x >>  8;
  x |= x >> 16;
  x |= x >> 32;
  return _popcnt(~x);
}

static inline int64_t _trailz (uint64_t x)
{
  return _popcnt((x-1) & ~x);
}

#endif  // defined(__GNUC__)

#endif  // defined(_CRAYC) || defined(__cray__)


/*** Prefetch exclusive (write intent) ***/

#if defined(__GNUC__)
# define _prefetch_x(a) __builtin_prefetch(a, 1, 3)
#elif (defined(__PPC__) || defined(__POWERPC__) || defined(__ppc__) || defined(__powerpc__)) \
  && (defined(__IBMC__) || defined(__IBMCPP__))
# define _prefetch_x __dcbtst
#else
# define _prefetch_x(a)
#endif


/*** Pseudorandom number generator ***/

#define BRAND_RUNUP_      128
#define BRAND_LG_TABSZ_     7
#define BRAND_TABSZ_      (INT64_C(1)<<BRAND_LG_TABSZ_)

typedef struct {
  uint64_t hi, lo, ind;
  uint64_t tab[BRAND_TABSZ_];
} brand_t;

#define BRAND_64STEP_(H,L,A,B) {                  \
    uint64_t x;                                 \
    x = H ^ (H << A) ^ (L >> (64-A));           \
    H = L | (x >> (B-64));                      \
    L = x << (128 - B);                         \
  }

static inline uint64_t brand (brand_t *p)
{
  uint64_t hi=p->hi, lo=p->lo, i=p->ind, ret;

  ret = p->tab[i];

  // 64-step a primitive trinomial LRS:  0, 45, 118
  BRAND_64STEP_(hi,lo,45,118);

  p->tab[i] = ret + hi;

  p->hi  = hi;
  p->lo  = lo;
  p->ind = hi & _maskr(BRAND_LG_TABSZ_);

  return ret;
}

static inline double dbrand (brand_t *p)
{
#if __UPC__
  const double n = 1.0 / 9007199254740992.0;
#else
  const double n = 0x1.0p-53;  // C99 hexadecimal floating-point
#endif
  uint64_t x = brand(p) & _maskr(53);
  return (x * n);
}

static BOLITE_NOWARN void brand_init (brand_t *p, uint64_t val)
{
  int64_t i;
  uint64_t hi, lo;

  hi = UINT64_C(0x9ccae22ed2c6e578) ^ val;
  lo = UINT64_C(0xce4db5d70739bd22) & ~_maskr(10);
  // ~_maskr(10) == _maskl(118-64)

  // we 64-step 0, 33, 118 in the initialization
  for (i = 0; i < 64; i++)
    BRAND_64STEP_(hi,lo,33,118);

  for (i = 0; i < BRAND_TABSZ_; i++) {
    BRAND_64STEP_(hi,lo,33,118);
    p->tab[i] = hi;
  }
  p->ind = BRAND_TABSZ_/2;
  p->hi  = hi;
  p->lo  = lo;

  for (i = 0; i < BRAND_RUNUP_; i++)
    brand(p);
}

#endif  // HAVE_MPP_UTIL


/*** Division (unsigned 32-bit) via multiplication ***/

typedef struct {
  uint32_t shift;
  uint32_t magic;
} divbymul32_t;

static inline uint32_t _divbymul32(uint32_t numerator, divbymul32_t div)
{
  uint64_t product = (uint64_t) numerator * (uint64_t) div.magic;
  return ((product >> 32) + numerator) >> div.shift;
}

static BOLITE_NOWARN divbymul32_t _divbymul32_prep(uint32_t divisor)
{
  // division by 0 will always yield 0
  if (divisor == 0)
    return (divbymul32_t) { .shift = 32, .magic = 0 };

  // division by a power of 2 is easy
  int k = 63 - _leadz((uint64_t) divisor);
  if (divisor == (uint32_t)1 << k)
    return (divbymul32_t) { .shift = (uint32_t) k, .magic = 0 };

  // now we have 2**k < divisor < 2**(k+1)
  // the 32-bit magic value is ceil(2**(33+k) / divisor) - 2**32
  uint64_t recip = 1 + ((UINT64_C(1) << 33 << k) - 1) / divisor;
  return (divbymul32_t) { .shift = (uint32_t)(k + 1), .magic = (uint32_t) recip };
}


#endif  // !defined(CONVEY_BOLITE_H)
