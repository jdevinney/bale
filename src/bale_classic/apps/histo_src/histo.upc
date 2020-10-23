/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
// 
//
//  All rights reserved.
//  
//   This file is a part of Bale.  For license information see the
//   LICENSE file in the top level directory of the distribution.
//  
// 
*****************************************************************/
/*! \file histo.upc
 * \brief Demo program that computes a histogram of uint64_t's
 *  The number of histogram bins should be large enough that they need 
 *  to be spread across the whole machine
 */

#include "histo.h"
#include <std_options.h>

/*!
\page histogram_page Histogram
  
The histogram application computes the histogram of a large number
of \f$\tt int64\_t\tt \f$'s into a large number, \f$\tt M \tt \f$, of buckets.

Each processor has an array of updates which are random indicies 
in the range \f$\tt [0:M] \tt \f$.  The job of histogram is to count the number 
of occurances of each index across all PEs. 
We can accomplish this task by creating a distributed array, counts,  of length \f$\tt M \tt \f$ 
and initializing it to zero. Then each PE can run through its list and 
increments the appropriate entry in the count array.
As described, this technique could lead to a race condition 
where two or more PE's simultaneously update the same entry in count. 
Since incrementing an entry is commutative it is independent of the 
order in which the updates are done.
So, doing the updates atomically is sufficient.

Histogram is extremely order and latency tolerant.
We believe it to be representative of a number of PGAS loops are dominated by random puts.
As a buffered communication pattern, we think of it as one-sided pushes.

Run with the --help, -?, or --usage flags for run details.
*/


typedef struct args_t{
  int64_t num_ups;      /*!< number of updates for all threads */
  int64_t l_num_ups;    /*!< number of updates for each thread */
  int64_t l_tbl_size;   /*!< per thread size of the counts table */
  std_args_t std;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'N':
      args->num_ups = atol(arg); break;
    case 'n':
      args->l_num_ups = atol(arg); break;
    case 'T':
      args->l_tbl_size = atol(arg); break;
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"l_num_updates",'n', "NUM", 0, "Number of updates per PE to the histogram table"},
    {"num_updates",  'N', "NUM", 0, "Number of updates for all PE to the histogram table"},
    {"table_size",   'T', "SIZE", 0, "Number of entries per PE in the histogram table"},
    {0}
  };

static struct argp_child children_parsers[] =
  {
    {&std_options_argp, 0, "Standard Options", -2},
    {0}
  };


int main(int argc, char * argv[]) {
  int64_t i;

  /* process command line */
  args_t args;
  args.l_tbl_size = 10000;
  args.num_ups = 0;
  args.l_num_ups = 5000000;
  struct argp argp = {options, parse_opt, 0,
                      "Accumulate updates into a table.", children_parsers};

  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);

  if(args.num_ups == 0) 
    args.num_ups = args.l_num_ups * THREADS;
  else
    args.l_num_ups = (args.num_ups + THREADS - MYTHREAD - 1)/THREADS;

  if(!MYTHREAD){
    bale_app_write_int(&args.std, "num_updates_total", args.num_ups);
    bale_app_write_int(&args.std, "num_updates_per_pe", args.l_num_ups);
    bale_app_write_int(&args.std, "table_size_per_pe", args.l_tbl_size);
    write_std_options(&args.std);
  }


  /* package all the inputs up into a struct */
  histo_t data;  
  data.l_num_ups = args.l_num_ups;
  data.lnum_counts = args.l_tbl_size;
  data.num_counts = data.lnum_counts*THREADS;
  data.counts = lgp_all_alloc(data.num_counts, sizeof(int64_t));
  assert(data.counts != NULL);  
  data.lcounts = lgp_local_part(int64_t, data.counts);  
  for(i = 0; i < data.lnum_counts; i++)
    data.lcounts[i] = 0L;
  
  // index is a local array of indices into the shared counts array.
  // This is used by the _agp version. 
  // To avoid paying for computing index[i]/THREADS and index[i]%THREADS
  // when using the exstack and conveyor models
  // we also store a packed version that holds the
  // pe (= index%THREADS) and lindx (=index/THREADS)
  data.index   = calloc(data.l_num_ups, sizeof(int64_t)); assert(data.index != NULL);
  data.pckindx = calloc(data.l_num_ups, sizeof(int64_t)); assert(data.pckindx != NULL);

  int64_t indx, lindx, pe;

  lgp_rand_seed(args.std.seed);
  
  for(i = 0; i < data.l_num_ups; i++) {
    //indx = i % data.num_counts;          //might want to do this for debugging
    indx = lgp_rand_int64(data.num_counts);
    data.index[i] = indx;                 
    lindx = indx / THREADS;
    pe  = indx % THREADS;
    assert(indx < data.num_counts);
    assert(lindx < data.lnum_counts);
    assert(lindx < (1L<<44));
    assert(pe < (1L<<20));
    data.pckindx[i]  =  (lindx << 20L) | (pe & 0xfffff);
  }
  double volume_per_node = (8*data.l_num_ups*args.std.cores_per_node)*(1.0E-9);
  
  lgp_barrier();

  int use_model;
  double laptime = 0.0;
  double injection_bw = 0.0;
  int64_t num_models = 0L;               // number of models that are executed
  char model_str[32];

  for( use_model=1L; use_model < 32; use_model *=2 ) {

    switch( use_model & args.std.models_mask ) {
    case AGP_Model:
      sprintf(model_str, "AGP");
      laptime = histo_agp(&data);
      num_models++;
      lgp_barrier();
      break;
    
    case EXSTACK_Model:
      sprintf(model_str, "Exstack");
      laptime = histo_exstack(&data, args.std.buf_cnt);
      num_models++;
      lgp_barrier();
      break;

    case EXSTACK2_Model:
      sprintf(model_str, "Exstack2");
      laptime = histo_exstack2(&data, args.std.buf_cnt);
      num_models++;
      lgp_barrier();
      break;

    case CONVEYOR_Model:
      sprintf(model_str, "Conveyor");
      laptime = histo_conveyor(&data);
      num_models++;
      lgp_barrier();
      break;

    case ALTERNATE_Model:
      T0_fprintf(stderr,"There is no alternate model here!\n"); continue;      
      //T0_fprintf(stderr,"Alternate: ");
      //laptime = histo_exstack2_cyclic(pckindx, l_num_ups, lcounts, buf_cnt);
      //laptime = histo_exstack2_goto(pckindx, l_num_ups, lcounts, buf_cnt);
      //laptime = histo_exstack_function(pckindx, l_num_ups, lcounts, buf_cnt);
      //laptime = histo_exstack2_function(pckindx, l_num_ups, lcounts, buf_cnt);
      //num_models++;
      lgp_barrier();
      break;

    default:
      continue;
      
    }
    injection_bw = volume_per_node / laptime;
    if(args.std.json == 0){
      T0_fprintf(stderr,"%10s: %8.3lf seconds  %8.3lf GB/s injection bandwidth\n", model_str, laptime, injection_bw);
    }else{
      bale_app_write_time(&args.std, model_str, laptime);
    }
  }
  
  lgp_barrier();

  // Check the results
  // Assume that the atomic add version will correctly zero out the counts array
  for(i = 0; i < data.l_num_ups; i++) {
    #if __cray__ || _CRAYC
      #pragma pgas defer_sync
    #endif
    lgp_atomic_add(data.counts, data.index[i], -num_models);
  }
  lgp_barrier();

  int64_t num_errors = 0, totalerrors = 0;
  for(i = 0; i < data.lnum_counts; i++) {
    if(data.lcounts[i] != 0L) {
      num_errors++;
      if(num_errors < 5)  // print first five errors, report number of errors below
        fprintf(stderr,"ERROR: Thread %d error at %"PRId64" (= %"PRId64")\n", MYTHREAD, i, data.lcounts[i]);
    }
  }
  totalerrors = lgp_reduce_add_l(num_errors);
  if(totalerrors) {
     T0_fprintf(stderr,"FAILED!!!! total errors = %"PRId64"\n", totalerrors);   
  }
  
  lgp_all_free(data.counts);
  free(data.index);
  free(data.pckindx);
  bale_app_finish(&args.std);

  return(totalerrors);
}

