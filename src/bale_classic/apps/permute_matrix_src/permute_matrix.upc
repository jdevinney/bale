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
//#include "alternates/permute_matrix_alternates.h"

/*! \file permute_matrix.upc
 * \brief Demo program that runs the variants of permute_matrix kernel. This program generates
 * a random square matrix according to the Erdos-Renyi model, two random permutations and then
 * permutes the rows and columns of the matrix according to the two random permutations.
 */

/*!
 * 
 * \page permute_matrix_page Permute Matrix
 *
 * This is a demo program that runs the variants of permute_matrix kernel. This program generates
 * a random square matrix according to the Erdos-Renyi model, two random permutations and then
 * permutes the rows and columns of the matrix according to the two random permutations.
 *
 * See files spmat_agp.upc, spmat_exstack.upc, spmat_exstack2.upc, and spmat_conveyor.upc
 * for the source for the kernels.
 * 
 * Run with the --help, -?, or --usage flags for run details.
 */

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

int main(int argc, char * argv[]) {
  
  /* process command line */
  args_t args = {0};  // initialize args struct to all zero
  struct argp argp = {NULL, parse_opt, 0,
                      "Parallel permute sparse matrix.", children_parsers};

  args.gstd.l_numrows = 1000000;
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);
  
  if(!MYTHREAD){
    write_std_graph_options(&args.std, &args.gstd);
    write_std_options(&args.std);
  }

  // read in a matrix or generate a random graph
  sparsemat_t * inmat = get_input_graph(&args.std, &args.gstd);
  if(!inmat){T0_fprintf(stderr, "ERROR: permute_matrix: inmat is NULL!\n");return(-1);}
  
  SHARED int64_t * rp = rand_permp(inmat->numrows, args.std.seed);
  SHARED int64_t * cp = rand_permp(inmat->numcols, args.std.seed + 12345);  

  if(args.std.dump_files){
    write_matrix_mm(inmat, "permute_matrix_inmat");
    if(!MYTHREAD){
      int64_t i;
      FILE * fp = fopen("permute_matrix_rperm", "w");
      for(i = 0; i < inmat->numrows; i++)
        fprintf(fp, "%ld\n", lgp_get_int64(rp, i));
      fclose(fp);
      fp = fopen("permute_matrix_cperm", "w");
      for(i = 0; i < inmat->numcols; i++)
        fprintf(fp, "%ld\n", lgp_get_int64(cp, i));
      fclose(fp);
    }
  }
  
  int64_t use_model;
  sparsemat_t * outmat;
  sparsemat_t * refmat = NULL;
  char model_str[32];
  int error = 0;
  for( use_model=1L; use_model < 32; use_model *=2 ) {
    double t1 = wall_seconds();
    switch( use_model & args.std.models_mask ) {

    case AGP_Model:
      outmat = permute_matrix_agp(inmat, rp, cp);
      sprintf(model_str, "AGP");
      break;

    case EXSTACK_Model:
      outmat = permute_matrix_exstack(inmat, rp, cp, args.std.buf_cnt);
      sprintf(model_str, "Exstack");
      break;

    case EXSTACK2_Model:
      outmat = permute_matrix_exstack2(inmat, rp, cp, args.std.buf_cnt);
      sprintf(model_str, "Exstack2");
      break;

    case CONVEYOR_Model:
      outmat = permute_matrix_conveyor(inmat, rp, cp);
      sprintf(model_str, "Conveyor");
      break;
    case ALTERNATE_Model:
      T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
      break;
    case 0:
      continue;
    }

    minavgmaxD_t stat[1];
    t1 = wall_seconds() - t1;
    lgp_min_avg_max_d( stat, t1, THREADS );
    bale_app_write_time(&args.std, model_str, stat->avg);

    /* if running more than one implmentation, save the first to check against the others*/
    if(!refmat){
      refmat = outmat;
    }else{
      if(compare_matrix(refmat, outmat)){
        T0_fprintf(stderr,"ERROR: permute_matrix does not match!\n");
        error = 1;
      }
      clear_matrix(outmat);
    }
  }
  
  if(args.std.dump_files){
    write_matrix_mm(refmat, "permute_matrix_outmat");
  }
  
  clear_matrix(inmat);
  clear_matrix(refmat);
  lgp_all_free(rp);
  lgp_all_free(cp);
  lgp_barrier();

  bale_app_finish(&args.std);
  
  return(error);
}


