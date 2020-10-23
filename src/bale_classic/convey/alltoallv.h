// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#ifndef CONVEY_ALLTOALLV_H
#define CONVEY_ALLTOALLV_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "common.h"
#include "convey.h"
#include "bolite.h"

typedef struct checksum {
  uint64_t sent;
  uint64_t rcvd;
} checksum_t;

void conveyor_bug(convey_t* c, const char* call, int error);

checksum_t basictest(convey_t* conveyor, brand_t* prng, double load,
                     size_t size, int entropy, convey_t* tally,
                     bool elastic, double reject, bool p2p_sums);

checksum_t aligntest(convey_t* conveyor, brand_t* prng, double load,
                     size_t size, bool elastic);

uint64_t* global_table_init(size_t echo_size, size_t entries, brand_t* prng);
void global_table_free(uint64_t* table);

checksum_t indexgather(convey_t* request, convey_t* reply, size_t reply_size,
                       brand_t* prng, double load, size_t entries,
                       const uint64_t source[entries]);

checksum_t histogram(convey_t* conveyor, size_t size, brand_t* prng,
                     double load, size_t entries, uint64_t target[entries]);

#endif
