// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


// The next three functions convert the PE number into a tag and the
// destination rank for pushing into porters[0].  A matrix routing tag is
// 16 bits (local).  A tensor routing tag is either 24 bits (remote) +
// 8 bits (local) or 8 bits + 8 bits.

static ROUTER_HINT route_t
vector_route(tensor_t* vector, int64_t pe)
{
  return (route_t) { .tag = 0, .next = pe };
}

static ROUTER_HINT route_t
matrix_route(tensor_t* matrix, int64_t pe)
{
  // dest is (x',y'), we are (x,y); hop to (x',y), tag is (y')
  uint32_t dest = pe;
  uint32_t upper = _divbymul32(dest, matrix->div_local);
  uint32_t lower = dest - matrix->n_local * upper;
#if MATRIX_REMOTE_HOP == 0
  return (route_t) { .tag = lower, .next = upper };
#else
  return (route_t) { .tag = upper, .next = lower };
#endif
}

static ROUTER_HINT route_t
tensor_route(tensor_t* tensor, int64_t pe)
{
  // dest is (x',y',z'), we are (x,y,z)
  // hop to (x,y,y'), tag is (x',z') [24 bits, 8 bits]
  uint32_t dest = pe;
  uint32_t upper = _divbymul32(dest, tensor->div_square);
  uint32_t middle = _divbymul32(dest, tensor->div_local);
  uint32_t lower = dest - tensor->n_local * middle;
  middle -= tensor->n_local * upper;
  uint32_t tag = (upper << 8) | lower;
  return (route_t) { .tag = tag, .next = middle };
}
