/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*!  \file transpose_matrix.c
\brief Program that runs the transpose_matrix routine from the spmat library.

Run transpose_matrix --help or --usage for usage.
*/

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_opts.h"

/*!
\brief timing wrapper around the transpose_matrix routine in `spmat_utils.c`
\param A  pointer to the given matrix
\param dump_files flag to output the matrix or not
\return runtime
*/
double transpose_generic(sparsemat_t *A, int64_t dump_files)
{
  double tm;
  //write_matrix_mm(A, "trans_orig.mm");

  tm = wall_seconds();
  sparsemat_t * At = transpose_matrix(A);
  tm = wall_seconds() - tm;

  //write_matrix_mm(At, "trans_tran.mm");
  clear_matrix(At);
  return(tm);
}

/********************************  argp setup  ************************************/
typedef struct args_t {
  std_args_t std;
  std_graph_args_t gstd;
} args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
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

static struct argp_option options[] = {
  {0}
};

static struct argp_child children_parsers[] = {
  {&std_options_argp, 0, "Standard Options", -2},
  {&std_graph_options_argp, 0, "Standard Graph Options", -3},
  {0}
};

int main(int argc, char * argv[])
{
  args_t args = {{0}};
  enum FLAVOR {GENERIC=1, ALL_Models=2};
  args.std.models_mask = ALL_Models-1;
  args.gstd.numrows = TRANSPOSE_NUM_ROWS;
  struct argp argp = {options, parse_opt, 0, "Permute rows and columns of a sparse matrix", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);
  write_std_options(&args.std);
  write_std_graph_options(&args.std, &args.gstd);

  // read in a matrix or generate a random graph
  sparsemat_t * mat = get_input_graph(&args.std, &args.gstd);
  if(!mat){fprintf(stderr, "ERROR: triangle: L is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(mat, "trans_inmat");

  uint32_t use_model;
  double laptime = 0.0;
  char model_str[64];
  for (use_model=1; use_model < 2; use_model *=2) {
    switch (use_model & args.std.models_mask) {
    case GENERIC:
      sprintf(model_str, "transpose matrix : ");
      laptime = transpose_generic(mat, args.std.dump_files);
    break;
    default:
      continue;
    }
    bale_app_write_time(&args.std, model_str, laptime);
  }
  // TODO: Should we test the result here or separate unit test on spmat_utils?
  return(0);
}
