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

#include <libgetput.h>
#include <spmat.h>
#include "randperm_alternates.h"
#include <std_options.h>

/*! \file randperm.upc
 * \brief Demo program that runs the variants of randperm kernel. This program
 * generates a random permutation in parallel.
 */

/*! 
\page randperm_page Random Permutation

Demo program that runs the variants of randperm kernel. This program
generates a random permutation in parallel. The algorithm used is the 
"dart throwing algorithm" found in 
 P.B.Gibbon, Y.Matias, and V.L.Ramachandran. Efficient low-contention Parallel algorithms.
 J. of Computer and System Sciences, 53:417-442, Dec 1992.

Interestingly, we discovered what looks to be a superior and simpler algorithm. This is implemented
in alternates/randperm_agp_opt.upc.

See files spmat_agp.upc, spmat_exstack.upc, spmat_exstack2.upc, and spmat_conveyor.upc
for the source for the kernels.

Run with the --help, -?, or --usage flags for run details.
 */

#include <unistd.h>

typedef struct args_t{
  int64_t num_rows;       /*!< total number of elements in the shared permutation */
  int64_t l_num_rows;     /*!< per thread number of elements in the shared permutation */
  std_args_t std;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'N':
      args->num_rows = atol(arg); break;
    case 'n':
      args->l_num_rows = atol(arg); break;
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"perm_size",'N', "NUM", 0, "Total length of permutation"},
    {"l_perm_size",'n', "NUM", 0, "Per PE length of permutation"},
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
  int ret = 0;
  args_t args;
  args.num_rows = 0;
  args.l_num_rows = 1000000;
  struct argp argp = {options, parse_opt, 0,
                      "Create a random permutation in parallel.", children_parsers};

  ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);

  if(args.num_rows == 0) 
    args.num_rows = args.l_num_rows * THREADS;
  else
    args.l_num_rows = (args.num_rows + THREADS - MYTHREAD - 1)/THREADS;
  
  if(!MYTHREAD){
    bale_app_write_int(&args.std, "num_rows_total", args.num_rows);
    bale_app_write_int(&args.std, "num_rows_per_pe", args.l_num_rows);
    write_std_options(&args.std);
  }

  double t1;
  minavgmaxD_t stat[1];
  int64_t error = 0;
  SHARED int64_t * out;
  int64_t seed = args.std.seed + MYTHREAD;
  int64_t use_model;
  char model_str[32];
  for( use_model=1L; use_model < 32; use_model *=2 ) {
    t1 = wall_seconds();
    switch( use_model & args.std.models_mask ) {

    case AGP_Model:
      sprintf(model_str, "AGP");
      out = rand_permp_agp(args.num_rows, seed);
      break;

    case EXSTACK_Model:
      sprintf(model_str, "Exstack");
      out = rand_permp_exstack(args.num_rows, seed, args.std.buf_cnt);
      break;

    case EXSTACK2_Model:
      sprintf(model_str, "Exstack2");
      out = rand_permp_exstack2(args.num_rows, seed, args.std.buf_cnt);
      break;

    case CONVEYOR_Model:
      sprintf(model_str, "Conveyor");
      out = rand_permp_conveyor(args.num_rows, seed);
      break;

    case ALTERNATE_Model:
      sprintf(model_str, "rand_permp_agp_opt");
      out = rand_permp_agp_opt(args.num_rows, seed);
      break;

    case 0:
      continue;
    }
    
    t1 = wall_seconds() - t1;    
    lgp_min_avg_max_d( stat, t1, THREADS );
    bale_app_write_time(&args.std, model_str, stat->avg);
    
    
    if(!is_perm(out, args.num_rows)){
      error++;
      T0_printf("\nERROR: rand_permp_%"PRId64" failed!\n\n", use_model & args.std.models_mask);
    }
    lgp_all_free(out);
  }
  
  if( error ) {
    T0_fprintf(stderr,"BALE FAIL!!!!\n"); 
  }
  bale_app_finish(&args.std);
  return(error);
}
