// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


//#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdalign.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "alltoallv.h"
// we need convey_imp_alloc, convey_imp_alloc_align
#include "private.h"


void
conveyor_bug(convey_t* c, const char* call, int error)
{
  mfprint(stderr, MY_PROC, 0, "%s FAILED: %s (%d)\n", call,
          convey_error_string(c, error), error);
  mpp_exit(-error);
}

static long n_allocs = 0;
static void** alloc_stack = NULL;

static void*
lifo_grab(void* alc8r, size_t size, const char* tag, uint64_t value)
{
  convey_alc8r_t* inner = alc8r;
  void* ptr = inner->grab(inner->alc8r, size, tag, value);
  if (ptr != NULL) {
    long n = n_allocs++;
    alloc_stack = realloc(alloc_stack, (n + 1) * sizeof(void*));
    alloc_stack[n] = ptr;
    // mprint(0, 0, "pushing %p onto stack (%s:%lu)\n", ptr, tag, value);
  }
  return ptr;
}

static void
lifo_free(void* alc8r, void* ptr)
{
  if (ptr != NULL) {
    // mprint(0, 0, "popping %p off stack\n", ptr);
    long n = --n_allocs;
    if (n < 0 || ptr != alloc_stack[n]) {
      mprint(MY_PROC, 0, (n < 0) ? "unmatched parallel free\n"
             : "parallel allocations are not LIFO\n");
      mpp_exit(7);  // -convey_error_ALLOC
    }
  }
  convey_alc8r_t* inner = alc8r;
  (inner->free)(inner->alc8r, ptr);
}

static convey_alc8r_t lifo_checking_alc8r = {
  .grab = &lifo_grab,
  .free = &lifo_free,
};


/*** Main Program ***/

static bool quiet = true;

static double now(void)
{
  struct timeval _t;
  gettimeofday(&_t, NULL);
  return _t.tv_sec + 1.0e-6 * _t.tv_usec;
}

static void
report_config(void)
{
#if MPP_USE_UPC
  const char* model = "UPC";
#elif MPP_RAW_MPI
  const char* model = "MPI";
#elif MPP_RAW_SHMEM
  const char* model = "SHMEM";
#elif MPP_USE_MPI
  const char* model = "MPP_UTIL/MPI";
#elif MPP_USE_SHMEM
  const char* model = "MPP_UTIL/SHMEM";
#elif HAVE_MPP_UTIL
  const char* model = "MPP_UTIL/NONE";
#else
  const char* model = "TRIVIAL";
#endif

#if ENABLE_PROFILING
  bool profile = true;
#else
  bool profile = false;
#endif

#if ENABLE_NONBLOCKING
  bool nonblock = true;
#else
  bool nonblock = false;
#endif

#if CONVEY_INTRANODE
  bool atomics = true;
#else
  bool atomics = false;
#endif

#if HAVE_AVX2 && HAVE_BMI2
  const char* codelets = "BMI2";
#elif HAVE_SSSE3
  const char* codelets = "SSSE3";
#else
  const char* codelets = "no ";
#endif

  mprint(0, 0, "libconvey version %s, %s build, profiling %s\n",
         PACKAGE_VERSION, model, (profile ? "on" : "off"));
#if MPP_USE_SHMEM || MPP_USE_UPC
  mprint(0, 0, "%sblocking puts, %splicit local transfers\n",
         (nonblock ? "non" : ""), (atomics ? "ex" : "im"));
#endif
  mprint(0, 0, "MATRIX_REMOTE_HOP = %d\n", MATRIX_REMOTE_HOP);
  mprint(0, 0, "%s compression codelets\n", codelets);
}

static void
usage(const char* prog, int status)
{
  mprint(0, 0,
         "%s [options] flavor load [item_size]\n"
         "flavor is from { simple, twohop, vector, matrix, tensor }\n"
         "options are:\n"
         "  -a #  autotune mpp_util all-to-all for # seconds { simple }\n"
         "  -b #  use # buffers per link { vector, matrix, tensor }\n"
         "  -c #  set buffer capacity to # bytes\n"
         "  -d    allocate buffers at each iteration\n"
         "  -e    send items of varying size { vector, matrix, tensor }\n"
         "  -f #  return # if the constructor fails\n"
         "  -g    use 'indexgather' benchmark\n"
         "  -h    use 'histogram' benchmark\n"
         "  -m #  access # words per PE (for -g or -h), default 100000\n"
         "  -n #  do # timed iterations, default 1\n"
         "  -p    use conveyor that guarantees steady progress\n"
         "  -q    suppress mpp_util timings\n"
         "  -r #  reject (unpull) a fraction # of pulled items\n"
         "  -s #  set pseudorandom seed to #\n"
         "  -t #  use # procs per local group { twohop, matrix, tensor }\n"
         "  -u    use mpp_util all-to-all functions { simple }\n"
         "  -v    increase verbosity\n"
         "  -w #  do # warmup iterations, default 1\n"
         "  -x    turn on prefetching for scatters\n"
         "  -y    test conveyor for steadiness\n"
         "  -z #  turn on compression if possible; vary # bits per item\n"
         "  -A    test that internal allocations are LIFO\n"
         "  -B    insist on blocking puts (if SHMEM or UPC)\n"
         "  -C    report library configuration\n"
         "  -H    report host name for each process\n"
         "  -L    force legacy (standard) buffer format\n"
	 "  -M    check alignment rather than data\n"
	 "  -P    verify point-to-point checksums\n"
         "  -R    randomize item size at each iteration\n",
         prog);
  if (status != 0)
    mpp_exit(status);
}

static void
print_args(int argc, char* argv[])
{
  mprint(0, 0, "args = ");
  for (int i = 1; i < argc; i++)
    mprint_np(0, 0, "%s%c", argv[i], ((i == argc-1) ? '\n' : ' '));
}

int
alltoallv(int argc, char* argv[])
{
  convey_alc8r_t* alloc = NULL;
  int autotune = 0;
  bool blocking = false;
  bool config = false;
  bool vendor = true;
  bool dynamic = false;
  bool elastic = false;
  bool gather = false;
  bool histo = false;
  bool hosts = false;
  bool steady = false;
  bool sorter = false;
  bool yoke = false;
  bool randomize = false;
  bool squeeze = false;
  bool standard = false;
  bool p2p_sums = false;
  bool testalign = false;
  int failcode = 1;
  size_t capacity = CONVEY_BUFFER_SIZE;
  size_t n_buffers = 1;
  size_t entries = 100000;
  int64_t n_timed = 1;
  int64_t n_warmup = 1;
  double reject = 0.0;
  uint64_t seed = 0;
  size_t row_procs = 0;
  int verbosity = 0;
  int entropy = 0;
  const char* flavor;
  double load;
  size_t item_size = 8;

  quiet = false;
  print_args(argc, argv);

  for (int opt;
       (opt = getopt(argc, argv, "a:b:c:def:ghm:n:pr:qs:t:uvw:xyz:ABCHLMPR")) != -1; ) {
    switch (opt) {
    case 'a': autotune = atoi(optarg); break;
    case 'b': n_buffers = strtoul(optarg, NULL, 10); break;
    case 'c': capacity = strtoul(optarg, NULL, 10); break;
    case 'd': dynamic = true; break;
    case 'e': elastic = true; break;
    case 'f': failcode = atoi(optarg); break;
    case 'g': gather = true; break;
    case 'h': histo = true; break;
    case 'm': entries = strtoul(optarg, NULL, 0); break;
    case 'n': n_timed = strtoul(optarg, NULL, 10); break;
    case 'p': steady = true; break;
    case 'r': reject = atof(optarg); break;
    case 'q': quiet = true; break;
    case 's': seed = strtoul(optarg, NULL, 10); break;
    case 't': row_procs = strtol(optarg, NULL, 10); break;
    case 'u': vendor = false; break;
    case 'v': verbosity++; break;
    case 'w': n_warmup = strtoul(optarg, NULL, 10); break;
    case 'x': sorter = true; break;
    case 'y': yoke = true; break;
    case 'z': squeeze = true; entropy = atoi(optarg); break;
    case 'A': alloc = &lifo_checking_alc8r; break;
    case 'B': blocking = true; break;
    case 'C': config = true; break;
    case 'H': hosts = true; break;
    case 'L': standard = true; break;
    case 'M': testalign = true; break;
    case 'P': p2p_sums = true; break;
    case 'R': randomize = true; break;
    default: usage(argv[0], 2);
    }
  }
  if (!seed)
    seed = time(NULL);
  if ((argc != optind + 2) && (argc != optind + 3))
    usage(argv[0], 2);
  if (yoke + gather + histo > 1) {
    mprint(0, 0, "%s: -g, -h, and -y are incompatible\n", argv[0]);
    return 2;
  }

  flavor = argv[optind];
  load = atof(argv[optind+1]);
  char* slash = strchr(argv[optind+1], '/');
  if (slash)
    load /= atof(slash + 1);
  if (argc > optind+2)
    item_size = strtoul(argv[optind+2], NULL, 10);
  if (!finite(load) || load <= 0.0 || item_size == 0) {
    mprint(0, 0, "%s: bogus parameters %f, %zu\n", argv[0], load, item_size);
    return 2;
  }
  if ((gather || histo) && item_size < 8) {
    mprint(0, 0, "%s: -g and -h require item_size >= 8\n", argv[0]);
    return 2;
  }

  if (config)
    report_config();
  if (row_procs == 0) {
    extern size_t convey_procs_per_node(void);
    row_procs = convey_procs_per_node();
    mprint(0, 0, "procs per group = %zu\n", row_procs);
  }
  if (hosts) {
    char name[64];
    gethostname(name, 64);
    name[63] = 0;
    mprint(MY_PROC, 0, "host = %s\n", name);
  }

  size_t n_procs = PROCS;
  convey_t* conveyor = NULL;
  convey_t* echo = NULL;
  mpp_alltoall_t* a2a = NULL;
  const char* constructor = NULL;
  const uint64_t sortopt = sorter * convey_opt_SCATTER;
  const uint64_t options = convey_opt_RECKLESS | dynamic * convey_opt_DYNAMIC |
    steady * convey_opt_PROGRESS | convey_opt_ALERT | squeeze * convey_opt_COMPRESS |
    standard * convey_opt_STANDARD | blocking * convey_opt_BLOCKING;
  const size_t echo_size = yoke ? sizeof(uint32_t) : 4 + item_size;
  const uint64_t echo_opts = options - steady * convey_opt_PROGRESS;
  if (gather)
    item_size = 8;
  if (alloc)
    alloc->alc8r = (void*)(elastic ? &convey_imp_alloc_align : &convey_imp_alloc);

  int order = 0;
  if (strcmp(flavor, "vector") == 0)
    order = 1;
  else if (strcmp(flavor, "matrix") == 0)
    order = 2;
  else if (strcmp(flavor, "tensor") == 0)
    order = 3;

  if (elastic && order) {
    constructor = ((const char* []) { "convey_new_etensor(1)",
          "convey_new_etensor(2)", "convey_new_etensor(3)" }) [order - 1];
    conveyor = convey_new_etensor(capacity, item_size, order,
                                  row_procs, n_buffers, alloc, options);
    if (gather || yoke)
      echo = convey_new_etensor(capacity, echo_size, order,
                                row_procs, n_buffers, alloc, echo_opts);
  }

  else if (order) {
    constructor = ((const char* []) { "convey_new_tensor(1)",
          "convey_new_tensor(2)", "convey_new_tensor(3)" }) [order - 1];
    conveyor = convey_new_tensor(capacity, order,
                                 row_procs, n_buffers, alloc, options);
    if (gather || yoke)
      echo = convey_new_tensor(capacity, order,
                               row_procs, n_buffers, alloc, echo_opts);
  }

  else if (strcmp(flavor, "twohop") == 0) {
    constructor = "convey_new_twohop";
    conveyor = convey_new_twohop(capacity, row_procs, alloc, options);
    if (gather || yoke)
      echo = convey_new_twohop(capacity, row_procs, alloc, echo_opts);
  }

  else if (strcmp(flavor, "simple") == 0) {

#if HAVE_MPP_UTIL
    if (!vendor || autotune) {
      a2a = mpp_alltoall_create(1);
      if (!autotune)
        mpp_alltoall_init(a2a, 0, 0, 0, 0, 0);
      else {
        size_t all_bytes = n_procs * capacity;
        void* src = mpp_alloc(all_bytes);
        void* dst = mpp_alloc(all_bytes);
        if (!src || !dst) {
          mprint(0, 0, "mpp_alloc failed (2 x %zu bytes)\n", all_bytes);
          return 1;
        }
        mpp_alltoall_tune(a2a, dst, src, capacity, 0, 0, autotune);
        mpp_free(dst);
        mpp_free(src);
      }
    }
#endif

    constructor = "convey_new_simple";
    conveyor = convey_new_simple(capacity, alloc, a2a, options | sortopt);
    // Don't use a sorter for the reply conveyor
    if (gather || yoke)
      echo = convey_new_simple(capacity, alloc, a2a, echo_opts);
  }

  else {
    mprint(0, 0, "%s: unknown conveyor flavor '%s'\n", argv[0], flavor);
    return 2;
  }

  mprint(0, 0, "%s buffer capacity: %zu bytes\n", flavor, capacity);

  // Check the conveyor
  long built = conveyor && (!(gather | yoke) || echo);
  built = mpp_and_long(built);
  if (!built) {
    mprint(0, 0, "%s FAILED\n", constructor);
    return failcode;
  }

  bool failed = false;
  if (squeeze) {
    int rc = convey_set_codec(conveyor, &convey_standard_codec, NULL);
    if (rc < 0)
      conveyor_bug(conveyor, "convey_set_codec", rc);
    bool ok = (rc == convey_OK);
    mprint(0, 0, "compression is %s\n", (ok ? "ON" : "OFF"));
    if (echo && ok) {
      rc = convey_set_codec(echo, &convey_standard_codec, NULL);
      if (rc != convey_OK)  {
        mprint(0, 0, "FAILURE: compression is OFF for the reverse conveyor!\n");
        failed = true;
      }
    }
  }

  brand_t _prng;
  brand_init(&_prng, seed ^ ((uint64_t)MY_PROC << 32));
  uint64_t* table = (gather || histo)
    ? global_table_init(echo_size, entries, &_prng) : NULL;

  // Test during warmup
  for (int64_t loop = 0; loop < n_warmup; loop++) {
    checksum_t sums;
    if (gather) {
      sums = indexgather(conveyor, echo, echo_size, &_prng, load, entries, table);
      if (sums.sent != sums.rcvd) {
        mprint(MY_PROC, 1, "LOCAL FAILURE: expected %08" PRIx64 " got %08" PRIx64 "\n",
               sums.sent, sums.rcvd);
        failed = true;
      }
    }
    else if (histo)
      sums = histogram(conveyor, item_size, &_prng, load, entries, table);
    else {
      size_t temp_size = item_size;
      if (randomize && loop > 0) {
        if (MY_PROC == 0)
          temp_size = 1 + (size_t)(dbrand(&_prng) * item_size);
        temp_size = mpp_min_long(temp_size);
        if (verbosity >= 1)
          mprint(0, 0, "testing item size = %zu\n", temp_size);
      }
      if (testalign)
	sums = aligntest(conveyor, &_prng, load, temp_size, elastic);
      else
	sums = basictest(conveyor, &_prng, load, temp_size, entropy,
			 echo, elastic, reject, p2p_sums);
    }
    sums.sent = mpp_accum_long(sums.sent);
    sums.rcvd = mpp_accum_long(sums.rcvd);
    if (sums.sent != sums.rcvd) {
      mprint(0, 0, "CHECKSUM FAILURE: sent %016" PRIx64 " rcvd %016" PRIx64 "\n",
             sums.sent, sums.rcvd);
      failed = true;
    }
  }

  if (n_warmup > 0)
    mprint(0, 0, "finished %" PRId64 " warm-up iterations\n", n_warmup);
  mpp_barrier(1);

  // Now do benchmarking
  double bw_sum = 0.0, bw_sqr = 0.0, bw_max = 0.0, seconds = 0.0;
  double start = now();
  for (int64_t loop = 0; loop < n_timed; loop++) {
    if (gather)
      indexgather(conveyor, echo, echo_size, &_prng, load, entries, table);
    else if (histo)
      histogram(conveyor, item_size, &_prng, load, entries, table);
    else
      basictest(conveyor, &_prng, load, item_size, entropy,
                echo, elastic, reject, false);
    mpp_barrier(1);
    // Compute some timing data
    double finish = now();
    double bw = (n_procs * load) / (finish - start);
    if (elastic)
      bw *= 0.5 * (item_size + 1);
    else
      bw *= (gather ? echo_size - 4 : item_size);
    bw_sum += bw;
    bw_sqr += bw * bw;
    bw_max = fmax(bw_max, bw);
    seconds += finish - start;
    start = finish;
  }

  if (n_timed > 0) {
    mprint(0, 0, "elapsed time: %g seconds\n", seconds);
    if (n_timed == 1)
      mprint(0, 0, "bandwidth per PE: %g bytes/sec\n", bw_sum / n_timed);
    else {
      double mean = bw_sum / n_timed;
      double stdv = sqrt((bw_sqr - mean * mean * n_timed) / (n_timed - 1));
      mprint(0, 0, "bandwidth per PE: %g +/- %.1f%% bytes/sec (max %g)\n",
             mean, 100.0 * stdv / mean, bw_max);
    }
  }

  long n_comms = convey_statistic(conveyor, convey_COMMS + convey_CUMULATIVE);
  long n_syncs = convey_statistic(conveyor, convey_SYNCS + convey_CUMULATIVE);
  long n_bytes = convey_statistic(conveyor, convey_BYTES + convey_CUMULATIVE);
  if (gather) {
    n_comms += convey_statistic(echo, convey_COMMS + convey_CUMULATIVE);
    n_syncs += convey_statistic(echo, convey_SYNCS + convey_CUMULATIVE);
    n_bytes += convey_statistic(echo, convey_BYTES + convey_CUMULATIVE);
  }
  n_comms = mpp_accum_long(n_comms);
  n_syncs = mpp_accum_long(n_syncs);
  n_bytes = mpp_accum_long(n_bytes);
  double terms = n_procs * (n_warmup + n_timed);
  mprint(0, 0, "average of %g sends, %g syncs, and %g bytes per PE\n",
         n_comms / terms, n_syncs / terms, n_bytes / terms);
  if (reject > 0.0) {
    long n_unpulls = convey_statistic(conveyor, convey_UNPULLS + convey_CUMULATIVE);
    n_unpulls = mpp_accum_long(n_unpulls);
    mprint(0, 0, "average of %g unpulls per PE\n", n_unpulls / terms);
  }

  if (gather || histo)
    global_table_free(table);
  if (echo)
    convey_free(echo);
  convey_free(conveyor);

  return failed;
}

int
main(int argc, char* argv[])
{
  // Initialize mpp_util without threading support
  argc = mpp_util_init(argc, argv, NULL);
  int status = 0;

  if (argc == 2 && strcmp(argv[1], "--help") == 0)
    usage(argv[0], 0);
  else if (argc == 2 && strcmp(argv[1], "--version") == 0)
    report_config();
  else if (argc == 2 && strcmp(argv[1], "--") == 0) {
    mprint(0, 0, "total PEs = %ld\n", PROCS);

    static char line[256];
    char* tokens[130];
    tokens[0] = argv[0];

    // Read, broadcast, and parse command lines from standard input,
    // until reaching an empty line or end of input.
    while (1) {
      if (MY_PROC == 0) {
        char* temp = fgets(line, 256, stdin);
        if (!temp)
          line[0] = 0;
      }
      mpp_broadcast64(line, 32, 0);

      int n = 1;
      char* saveptr = NULL;
      while ((tokens[n] = strtok_r((n == 1) ? line : NULL,
                                   " \t\n", &saveptr)) != NULL)
        n++;
      if (n == 1)
        break;
      int failure = alltoallv(n, tokens);
      if (status < failure)
        status = failure;
      fflush(stdout);
#if HAVE_OPTRESET
      optind = 1;
      optreset = 1;
#else
      optind = 0;  // correct value for glibc
#endif
    }

    mprint(0, 0, "job finished\n");
  }
  else
    status = alltoallv(argc, argv);

  if (status > 0)
    mpp_exit(status);

#if HAVE_MPP_UTIL
  if (quiet)
    mprint_set_filter(NULL);
#endif
  free(alloc_stack);
  mpp_util_end();
}
