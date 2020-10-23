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
#include <getopt.h>
#include <libgetput.h>
#include <spmat.h>
#include <std_options.h>

//#include "alternates/transpose_matrix_alternates.h"

/*! \file transpose_matrix.upc
 * \brief Driver program that runs the variants of transpose_matrix kernel.
 */

/*! 
 * \page transpose_matrix_page Transpose Matrix
 *
 * Demo program that runs the variants of transpose matrix kernel. This program
 * generates a random square asymmetrix (via the Erdos-Renyi model) and then transposes
 * it in parallel.
 *
 * See files spmat_agp.upc, spmat_exstack.upc, spmat_exstack2.upc, and spmat_conveyor.upc
 * for the source for the kernels.
 * 
 * Run with the --help, -?, or --usage flags for run details.
 */

/********************************  argp setup  ************************************/
typedef struct args_t{
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_child children_parsers[] =
  {
    {&std_options_argp, 0, "Standard Options", -2},
    {&std_graph_options_argp, 0, "Standard Graph Options", -3},
    {0}
  };

int main(int argc, char * argv[])
{

  /* process command line */
  args_t args = {0}; // initialize args struct to all zero
  struct argp argp = {NULL, parse_opt, 0,
                      "Parallel sparse matrix transpose.", children_parsers};

  args.gstd.l_numrows = 500000;  
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);  
  if(ret < 0) return(ret);
  else if(ret) return(0);

  if(!MYTHREAD){
    write_std_graph_options(&args.std, &args.gstd);
    write_std_options(&args.std);
  }

  // read in a matrix or generate a random graph
  sparsemat_t * inmat = get_input_graph(&args.std, &args.gstd);
  if(!inmat){T0_fprintf(stderr, "ERROR: transpose: inmat is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(inmat, "transpose_inmat");

  int write_out = 0;
  double t1;
  minavgmaxD_t stat[1];
  int64_t error = 0;
  int64_t use_model;
  sparsemat_t * outmat;
  char model_str[32];
  for( use_model=1L; use_model < 32; use_model *=2 ) {
    t1 = wall_seconds();
    switch( use_model & args.std.models_mask ) {
    case AGP_Model:
      outmat = transpose_matrix_agp(inmat);
      sprintf(model_str, "AGP");
      break;
    case EXSTACK_Model:
      outmat = transpose_matrix_exstack(inmat, args.std.buf_cnt);
      sprintf(model_str, "Exstack");
      break;
    case EXSTACK2_Model:
      outmat = transpose_matrix_exstack2(inmat, args.std.buf_cnt);
      sprintf(model_str, "Exstack2");
      break;
    case CONVEYOR_Model:
      outmat = transpose_matrix_conveyor(inmat);
      sprintf(model_str, "Conveyor");
      break;    
    case ALTERNATE_Model:
      continue;
    case 0:
      continue;
    }
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    bale_app_write_time(&args.std, model_str, stat->avg);
    
    /* correctness check */
    sparsemat_t * outmatT = transpose_matrix(outmat);
    if(compare_matrix(outmatT, inmat)){
      T0_fprintf(stderr,"ERROR: transpose of transpose does not match!\n");
      error = 1;
    }
    
    if((write_out == 0) && (error == 0)){
      if(args.std.dump_files) write_matrix_mm(outmat, "transpose_outmat");
      write_out = 1;
    }

    lgp_barrier();
    
    clear_matrix(outmatT);    
    clear_matrix(outmat);    
  }
  
  clear_matrix(inmat);
  lgp_barrier();
  bale_app_finish(&args.std);
  return(error);
}


