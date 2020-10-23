/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file histo.c
\brief A program that computes a histogram of uint64_t's
The intent is to histogram a large number of words into a large number of bins.

Run histo --help or --usage for insructions on running
*/

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_opts.h"


/*! \brief This routine is a generic histogram loop
\param *index array of indices into the array of counts
\param num_ups the length of the index array (number of updates)
\param *counts pointer to the count array
\param v the amount to add to the counts array for each update.
        we can use v = (-1) to check against other implementations.
\return run time
*/
double histo_generic(int64_t *index, int64_t num_ups, int64_t *counts, int64_t v) 
{
  double tm;
  int64_t i;

  tm = wall_seconds();

  for(i = 0; i < num_ups; i++)
    counts[index[i]] += v;

  tm = wall_seconds() - tm;
  return( tm );
}

#define LOG_NUM_BUFFERS 6                 /*!< parameters to play with buffering histo */
#define NUM_BUFFERS (1L<<LOG_NUM_BUFFERS) /*!< number of buffers */
#define BUFFER_SIZE 128                   /*!< size of the buffers */
/*! \brief This routine does a buffered version of histogram.
\param *index array of indices into the array of counts
\param num_ups the length of the index array (number of updates)
\param *counts pointer to the count array
\param table_size size of the histogram table
\return run time
*/
double histo_buffered(int64_t *index, int64_t num_ups, int64_t *counts, int64_t table_size) 
{
  double tm;
  int64_t i,j;
  int64_t nbits;
  int64_t sort_shift;
  int64_t s, cnts[NUM_BUFFERS]; 
  int64_t counts_idx[NUM_BUFFERS][BUFFER_SIZE];
  
  assert(table_size > 0);
  nbits = 0;
  while(table_size>>nbits){   // shift table size to find the number of bit in it
    nbits++;
  }
  // We will put indices into buffer according to their top LOG_NUM_BUFFERS 
  sort_shift = nbits - LOG_NUM_BUFFERS;
  sort_shift = (sort_shift > 0) ? sort_shift : 0;
  
  for(i = 0; i < NUM_BUFFERS; i++)
     cnts[i] = 0L; 
  
  tm = wall_seconds();
   
  for(i = 0; i < num_ups; i++){
    s = index[i] >> sort_shift;
    assert((0 <= s) && (s<NUM_BUFFERS));
    assert((0 <= cnts[s]) && (cnts[s]<BUFFER_SIZE));
    counts_idx[s][cnts[s]] = index[i];
    cnts[s]++;
    if(cnts[s] >= BUFFER_SIZE) { 
      for(j = 0; j < cnts[s]; j++)
        counts[ counts_idx[s][j] ]++;
      cnts[s] = 0;
    }
  }
  for(s = 0; s < NUM_BUFFERS; s++)            // empty any remaining elements in the buffers
    for(j = 0; j < cnts[s]; j++)
      counts[ counts_idx[s][j] ]++;
  
  tm = wall_seconds() - tm;
  return( tm );
}

/*! \brief comparison function to support histo_sorted */
static int comp(const void *a, const void *b) 
{
  return( *(uint64_t *)a - *(uint64_t *)b );
}

/*! \brief This routine does a sorting version of histogram.
\param *index array of indices into the array of counts
\param num_ups the length of the index array (number of updates)
\param *counts pointer to the count array
\return run time of just the updates NOT the sorting

Given the random index array, we first sort it so that the 
updates to the counts array occur in index order.
NB: we don't charge for the sorting.
*/
double histo_sorted(int64_t *index, int64_t num_ups,  int64_t *counts) 
{
  double tm;
  int64_t i;

  qsort( index, num_ups, sizeof(int64_t),  comp );

  tm = wall_seconds();
  for(i = 0; i < num_ups; i++) 
    counts[index[i]] += 1;
  
  tm = wall_seconds() - tm;
  return( tm );
}


/********************************  argp setup  ************************************/
typedef struct args_t {
  int64_t num_ups;
  int64_t tbl_size;
  std_args_t std;
} args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key)
  {
  case 'N':
    args->num_ups = atol(arg); break;
  case 'T':
    args->tbl_size = atol(arg); break;
  case ARGP_KEY_INIT:
    state->child_inputs[0] = &args->std;
    break;
  }
  return(0);
}

static struct argp_option options[] = {
  {"num_updates",'N', "NUM", 0, "Number of updates to the histogram table"},
  {"table_size", 'T', "SIZE", 0, "Number of entries in the histogram table"},
  {0}
};

static struct argp_child children_parsers[] = {
  {&std_options_argp, 0, "Standard Options", -2},
  {0}
};

int main(int argc, char * argv[]) 
{
  enum MODEL {GENERIC_Model=1, BUF_Model=2, SORT_Model=4, ALL_Models=8};
  args_t args;
  args.tbl_size = HISTO_TABLE_SIZE;
  args.num_ups = HISTO_NUM_UPDATES;
  args.std.models_mask = ALL_Models-1;
  struct argp argp = {options, parse_opt, 0, "Accumulate updates into a table.", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);
  write_std_options(&args.std);
  bale_app_write_int(&args.std, "num_updates", args.num_ups);
  bale_app_write_int(&args.std, "tbl_size", args.tbl_size);
    
  // index is an array of indices into the counts array.
  int64_t * index  = calloc(args.num_ups, sizeof(int64_t)); assert(index != NULL);
  int64_t * counts = calloc(args.tbl_size, sizeof(int64_t)); assert(counts != NULL);  

  rand_seed( args.std.seed );
  int64_t i;
  for(i = 0; i < args.num_ups; i++)
    index[i] = rand_int64(args.tbl_size);
  
  uint32_t use_model;
  double laptime = 0.0;
  char model_str[64];
  int32_t num_runs = 0;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & args.std.models_mask ) {
    case GENERIC_Model:
      sprintf(model_str, "Generic  histo: ");
      laptime = histo_generic(index, args.num_ups, counts, 1);
      num_runs++;
      break;
    case BUF_Model:
      sprintf(model_str, "Buffered histo: ");
      laptime = histo_buffered(index, args.num_ups, counts, args.tbl_size);
      num_runs++;
      break;
    case SORT_Model:
      sprintf(model_str, "Sorted   histo: ");
      laptime = histo_sorted(index, args.num_ups, counts);
      num_runs++;
      break;
    default:
      continue;
    }
    bale_app_write_time(&args.std, model_str, laptime);
  }
  
  // Check the buffered and sorted against the generic again
  int64_t errors = 0;
  laptime = histo_generic(index, args.num_ups,  counts, -num_runs);
  for(i = 0; i < args.tbl_size; i++) {
    if(counts[i] != 0L){
      errors++;
      if(errors < 5)  // print first five errors, report number of errors below
        fprintf(stderr,"ERROR: first five errors at %"PRId64" (= %"PRId64")\n", i, counts[i]);
    }
  }
  if(errors)
    fprintf(stderr,"FAILED!!!! total errors = %"PRId64"\n", errors);   
  
  free(index);
  free(counts);
  bale_app_finish(&args.std);
  return(errors);
}
