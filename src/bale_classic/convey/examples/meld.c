// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include "example.h"
#include "biconvey.h"

// Try different ways of merging data from two biconveyors.  The details
// don't matter for our purposes, as long as the lookups into the two
// structures are statistically independent.  We form a distributed list of
// random pairs (a_i, b_i) and make two distributed hash tables F and G
// such that F(a_i) = b_i and G(b_i) = a_i.  Then we iterate over i and
// check that (G(b_i), F(a_i)) == (a_i, b_i).

#if __MACH__
#include <sys/time.h>
#else
#include <time.h>
#endif

static uint64_t
now_nsec(void)
{
#if __MACH__
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return UINT64_C(1000000000) * tv.tv_sec + UINT64_C(1000) * tv.tv_usec;
#else
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC, &tp);
  return UINT64_C(1000000000) * tp.tv_sec + tp.tv_nsec;
#endif
}

// The following definitions implement a very crude hash table with
// linear probing.  It assumes that all keys are nonzero and that the
// table doesn't fill up.

// typedef __uint128_t item_t;
typedef uint64_t item_t;

typedef struct entry {
  item_t key;
  item_t value;
} entry_t;

typedef struct hashtable {
  size_t mask;
  entry_t table[];
} hashtable_t;

static item_t
ht_get(const hashtable_t* ht, item_t key)
{
  for (size_t i = key & ht->mask; ht->table[i].key; i = (i + 1) & ht->mask)
    if (ht->table[i].key == key)
      return ht->table[i].value;
  return (item_t) 0;
}

static void
ht_put(hashtable_t* ht, item_t key, item_t value)
{
  size_t i = key & ht->mask;
  while (ht->table[i].key && ht->table[i].key != key)
    i = (i + 1) & ht->mask;
  ht->table[i] = (entry_t) { key, value };
}

static hashtable_t*
ht_new(int lg_size)
{
  if (lg_size < 8)
    lg_size = 8;
  size_t mask = ~UINT64_C(0) >> (64 - lg_size);
  hashtable_t* ht = calloc(1, sizeof(hashtable_t) + (mask + 1) * sizeof(entry_t));
  ht->mask = mask;
  return ht;
}

static item_t
brand_item(brand_t* prng)
{
  item_t item = 0;
  for (int i = 0; i < sizeof(item_t) / sizeof(uint64_t); i++)
    item = ((item << 32) << 32) | brand(prng);
  return item;
}

static void
lookup(const void* query, void* reply, void* context)
{
#if 1
  item_t key, value;
  memcpy(&key, query, sizeof(item_t));
  value = ht_get(context, key);
  memcpy(reply, &value, sizeof(item_t));
#else
  // fake it! the test will still work
  memcpy(reply, query, sizeof(item_t));
#endif
}

// a version of lookup that carries along an extra item
typedef struct pair {
  item_t item;
  item_t extra;
} pair_t;

static void
lookup_x(const void* query, void* reply, void* context)
{
  pair_t pair;
  memcpy(&pair, query, sizeof(pair_t));
  pair.item = ht_get(context, pair.item);
  memcpy(reply, &pair, sizeof(pair_t));
}

static void
check(size_t k, item_t a, item_t b, item_t fa, item_t gb)
{
  if ((a ^ b) != (fa ^ gb)) {
    printf("rank %ld: failure at position %zu (%llx %llx %llx %llx)\n",
	   MY_PROC, k, (uint64_t) a, (uint64_t) b, (uint64_t) fa, (uint64_t) gb);
    exit(EXIT_FAILURE);
  }
}

int
main(int argc, char* argv[])
{
  example_start();

  // Parse command line and environment
  const long n_procs = PROCS;
  const long my_proc = MY_PROC;
  if (n_procs & (n_procs - 1)) {
    if (my_proc == 0)
      printf("number of ranks is %ld, must be a power of 2\n", n_procs);
    return 1;
  }
  int lg_size = (argc > 1) ? atoi(argv[1]) : 16;
  if (my_proc == 0)
    printf("command: %s %d [= log2(local_size)]\n", argv[0], lg_size);
  if (lg_size < 0 || lg_size > 30) {
    if (my_proc == 0)
      printf("argument is out of valid range [0,30]\n");
    return 1;
  }

  // Prepare local data
  brand_t _prng;
  brand_init(&_prng, 1 + my_proc);
  const size_t n = (size_t)1 << lg_size;
  item_t* alist = malloc(n * sizeof(item_t));
  item_t* blist = malloc(n * sizeof(item_t));
  for (size_t i = 0; i < n; i++) {
    alist[i] = brand_item(&_prng);
    blist[i] = brand_item(&_prng);
  }
  hashtable_t* ftable = ht_new(lg_size + 1);
  hashtable_t* gtable = ht_new(lg_size + 1);
  const int bits = _popcnt(ftable->mask);

#define HOME(item) (((unsigned long)(item) >> bits) & (unsigned long)(n_procs - 1))

  // Fill distributed hash tables
  convey_t* convey = convey_new(SIZE_MAX, 0, NULL, convey_opt_ALERT);
  for (int swap = 0; swap <= 1; swap++) {
    item_t* keys = swap ? blist : alist;
    item_t* values = swap ? alist : blist;
    hashtable_t* table = swap ? gtable : ftable;

    convey_begin(convey, sizeof(entry_t), alignof(entry_t));
    for (size_t i = 0; convey_advance(convey, i == n); ) {
      entry_t pair;
      for (; i < n; i++) {
	pair = (entry_t) { keys[i], values[i] };
	if (! convey_push(convey, &pair, HOME(pair.key)))
	  break;
      }
      while (convey_pull(convey, &pair, NULL))
	ht_put(table, pair.key, pair.value);
    }
    convey_reset(convey);
  }
  convey_free(convey);

  // Now do the real test
  // biconvey_t* f = biconvey_new(SIZE_MAX, 0, NULL, convey_opt_ALERT);
  // biconvey_t* g = biconvey_new(SIZE_MAX, 0, NULL, convey_opt_ALERT);
  biconvey_t* f = biconvey_new_simple(10000, NULL, NULL, convey_opt_ALERT);
  biconvey_t* g = biconvey_new_simple(10000, NULL, NULL, convey_opt_ALERT);
  uint64_t ticks = 0;
  
  for (int loop = 0; loop < 10; loop++) {
    if (loop == 1)
      ticks = now_nsec();
#if 1
    biconvey_begin(f, sizeof(item_t), sizeof(item_t), &lookup, ftable);
    biconvey_begin(g, sizeof(item_t), sizeof(item_t), &lookup, gtable);

    bool have_fx = false;
    size_t i = 0, j = 0, k = 0;
    item_t fx, gx;
    while (biconvey_advance(f, i == n) | biconvey_advance(g, j == n)) {
      for (; i < n; i++)
	if (! biconvey_push(f, &alist[i], HOME(alist[i])))
	  break;
      for (; j < n; j++)
	if (! biconvey_push(g, &blist[j], HOME(blist[j])))
	  break;
      for (; k < n; k++, have_fx = false) {
	have_fx = have_fx || biconvey_pull(f, &fx);
	if (have_fx && biconvey_pull(g, &gx))
	  check(k, alist[k], blist[k], fx, gx);
	else
	  break;
      }
    }
#else
    // worse idea, but not a lot worse on a single node?
    biconvey_begin(f, sizeof(item_t), sizeof(item_t), &lookup, ftable);
    biconvey_begin(g, sizeof(pair_t), sizeof(pair_t), &lookup_x, gtable);

    bool have_query = false;
    pair_t query, reply;
    size_t i = 0, j = 0, k = 0;
    while (biconvey_advance(f, i == n) | biconvey_advance(g, j == n)) {
      for (; i < n; i++)
	if (! biconvey_push(f, &alist[i], HOME(alist[i])))
	  break;
      for (; ; j++, have_query = false) {
	query.item = blist[j];
	have_query = have_query || biconvey_pull(f, &query.extra);
	if (!have_query || !biconvey_push(g, &query, HOME(blist[j])))
	  break;
      }
      for (; biconvey_pull(g, &reply); k++)
	check(k, alist[k], blist[k], reply.extra, reply.item);
    }

#endif
    assert(i == n && j == n && k == n);

    biconvey_reset(g);
    biconvey_reset(f);
  }
  ticks = now_nsec() - ticks;

  biconvey_free(g);
  biconvey_free(f);

  if (my_proc == 0)
    printf("%.1f nsec per pair of lookups\n", ticks * 1.0 / (9 * n));

  example_end();
  return 0;
}
