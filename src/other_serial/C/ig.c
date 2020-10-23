/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file ig.c
\brief A program that computes the gather of a large number of things
from a large table.

Run ig --help or --usage for insructions on running.
*/


#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_opts.h"


/*! \brief check that the indexgather worked. 
\param tgt the target array
\param index the indices that drive the gather
\param len length
\returns the number of disagreements.

THIS REQUIRES that the source array, table, is set to just be minus the index
*/
int64_t ig_check_and_zero(int64_t *tgt, int64_t *index, int64_t len)
{
  int64_t errors=0;
  int64_t i;
  for(i=0; i<len; i++){
    if( tgt[i] != -index[i] ) {
      errors ++;
      if( errors < 5 )  // print the first 5 errors and count the rest
        printf("  error tgt[%"PRId64"] = %"PRId64" != %"PRId64"\n", i, tgt[i], -index[i] );
    }
    tgt[i] = 0;
  }
  if( errors ) 
    printf(" total of %"PRId64" errors\n", errors);
  return(errors);
}


/*! \brief This is the generic serial version of indexgather
\param *tgt array of target locations for the gathered values
\param *index array of indices into the source array of counts
\param num_req the length of the index array (number of updates)
\param *table the array from which the values are gathered
\return run time
*/
double ig_generic(int64_t *tgt, int64_t *index, int64_t num_req,  int64_t *table) 
{
  int64_t i;
  double tm;

  tm = wall_seconds();

  for(i = 0; i < num_req; i++)
    tgt[i] = table[ index[i] ];

  tm = wall_seconds() - tm;
  return( tm );
}


#define LOG_NUM_BUFFERS 6                     /*!< parameters to play with buffering */
#define NUM_BUFFERS (1L<<LOG_NUM_BUFFERS)     /*!< the number of buffers */
#define BUFFER_SIZE 128                       /*!< the size of the buffers */
/*!
\brief This routine implements a buffered version of indexgather
\param tgt array of target locations for the gathered values
\param index array of indices into the source array of counts
\param num_req the length of the index array (number of updates)
\param table the array from which the values are gathered
\param table_size the size of the table 
\return runtime
*/
double ig_buffered(int64_t *tgt, int64_t *index, int64_t num_req,  int64_t *table, int64_t table_size) 
{
  double tm;
  int64_t i, j;
  int64_t nbits;
  int64_t sort_shift;

  int64_t s, cnts[NUM_BUFFERS]; 
  int64_t table_idx[NUM_BUFFERS][BUFFER_SIZE];
  int64_t tgt_idx[NUM_BUFFERS][BUFFER_SIZE];

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

  for(i = 0; i < num_req; i++){
    s = index[i] >> sort_shift;
    assert( (0 <= s) && (s<NUM_BUFFERS));
    assert( (0 <= cnts[s]) && (cnts[s]<BUFFER_SIZE));
    tgt_idx[s][cnts[s]] = i;
    table_idx[s][cnts[s]] = index[i];
    cnts[s]++;

    if( cnts[s] >= BUFFER_SIZE ) { 
      for(j = 0; j < cnts[s]; j++){
        tgt[ tgt_idx[s][j] ] = table[ table_idx[s][j] ];
      }
      cnts[s] = 0;
    }
  }
  for(s = 0; s < NUM_BUFFERS; s++){
    for(j = 0; j < cnts[s]; j++){
      tgt[ tgt_idx[s][j] ] = table[ table_idx[s][j] ];
    } 
  }

  tm = wall_seconds() - tm;
  return( tm );
}

/********************************  argp setup  ************************************/
typedef struct args_t{
  int64_t num_req;
  int64_t tbl_size;
  std_args_t std;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key)
  {
  case 'N':
    args->num_req = atol(arg); break;
  case 'T':
    args->tbl_size = atol(arg); break;
  case ARGP_KEY_INIT:
    state->child_inputs[0] = &args->std;
    break;
  }
  return(0);
}

static struct argp_option options[] = {
  {"num_requests", 'N', "NUM",  0, "Number of requests to the table"},
  {"table_size",   'T', "SIZE", 0, "Number of entries in look-up table"},
  {0}
};

static struct argp_child children_parsers[] = {
  {&std_options_argp, 0, "Standard Options", -2},
  {0}
};

int main(int argc, char * argv[])
{
  enum MODEL {GENERIC_Model=1, BUF_Model=2, ALL_Models=4};
  args_t args;
  args.tbl_size = IG_TABLE_SIZE;
  args.num_req = IG_NUM_UPDATES;
  args.std.models_mask = ALL_Models-1;
  struct argp argp = {options, parse_opt, 0, "Index gather from a table", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);
  write_std_options(&args.std);
  bale_app_write_int(&args.std, "num_updates", args.num_req);
  bale_app_write_int(&args.std, "tbl_size", args.tbl_size);


  int64_t *table   = calloc(args.tbl_size, sizeof(int64_t));
  int64_t *tgt     = calloc(args.num_req, sizeof(int64_t));
  int64_t *index   = calloc(args.num_req, sizeof(int64_t));

  //populate table array and the index array
  int64_t i;
  for(i=0; i<args.tbl_size; i++)
    table[i] = -i;   // fill the table with minus the index, so we can check it easily

  rand_seed(args.std.seed);
  for(i = 0; i < args.num_req; i++)
    index[i] = rand_int64(args.tbl_size);

  uint32_t use_model;
  double laptime = 0.0;
  char model_str[64];
  int64_t errors = 0L;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    switch( use_model & args.std.models_mask ){
    case GENERIC_Model:
      sprintf(model_str, "Generic  IG: ");
      laptime = ig_generic(tgt, index, args.num_req, table);
      errors += ig_check_and_zero(tgt, index, args.num_req);
      break;
    case BUF_Model:
      sprintf(model_str, "Buffered IG: ");
      laptime =ig_buffered(tgt, index, args.num_req,  table,  args.tbl_size); 
      errors += ig_check_and_zero(tgt, index, args.num_req);
      break;
    default:
      continue;
    }
    bale_app_write_time(&args.std, model_str, laptime);
  }

  free(table);
  free(tgt);
  free(index);

  if( errors ) 
    fprintf(stderr,"FAILED!!!!\n"); 
  bale_app_finish(&args.std);
  return(errors);
}


