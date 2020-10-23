/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file triangle.c
\brief Program that counts the number of triangles in a graph given by its adjacency matrix

Run triangle --help or --usage for insructions on running
*/

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_opts.h"

/*!
\brief This routine counts the number of triangles in a graph
 using the matrix operation sum{A .& A * B}
\param triangles a place to write the number of triangles found
\param A the sparse matrix for A
\param B the sparse matrix for B transpose
\return run time
*/
double triangles_matrix(int64_t *triangles, sparsemat_t *A, sparsemat_t *B) 
{
  int64_t j, k, l, numtriangles;
  int64_t u,v;
  numtriangles = 0;

  double t1 = wall_seconds();
  
  // for each non-zero (i,j) in L accumulate the size 
  // of the intersection of row_i of A  and row_j of B.
  // Note: Because the matrix is a {0,1} tidy matrix,
  // (the columns in each row appear in increasing order)
  // we can find the intersection in a single pass over both rows.
  // Because  B is the transpose of the second matrix, this is the 
  // masked matrix product.

  for(u = 0; u < A->numrows; u++){ 
    for(j = A->offset[u]; j < A->offset[u+1]; j++){
      v = A->nonzero[j];
      for( l = A->offset[u], k = B->offset[v];  l < A->offset[u+1] && k < B->offset[v+1];  ) {
        if( A->nonzero[l] == B->nonzero[k] ) {
          numtriangles++;
          k++;
          l++;
        } else if( A->nonzero[l] > B->nonzero[k] ){
          k++;
        }else{ // ( A->nonzero[u] > A->nonzero[W] )
          l++;
        }
      }
    }
  }

  t1 = wall_seconds() - t1;
 
  *triangles = numtriangles;
  return(t1);
}

/********************************  argp setup  ************************************/
typedef struct args_t {
  std_args_t std;
  std_graph_args_t gstd;
} args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key) {
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
  enum FLAVOR {L_LxL=1, L_LxU=2, ALL_Models=4};
  args.std.models_mask = ALL_Models-1;
  struct argp argp = {options, parse_opt, 0, "Triangle counting", children_parsers};  
  args.gstd.numrows = TRIANGLE_NUM_ROWS;
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);
  args.gstd.loops = 0;   //overwrite these: the algorithm requires no loops, and undirected 
  args.gstd.directed = 0;
  write_std_options(&args.std);
  write_std_graph_options(&args.std, &args.gstd);
  
  // read in a matrix or generate a random graph
  sparsemat_t * L = get_input_graph(&args.std, &args.gstd);
  if(!L){fprintf(stderr, "ERROR: triangle: L is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(L, "triangle_inmat");

  /* make sure the matrix is in legal form */
  if (args.gstd.readfile) {
    ret = tril(L, -1);
    if (ret)
      fprintf(stderr,"WARNING: input graph was not lower-triangular with zero diagonal. Removing illegal nonzeros.\n");
  } else if (!is_lower_triangular(L, 0)) {
    fprintf(stderr,"ERROR: L is not lower triangular!\n");
    exit(1);  
  }  

  sparsemat_t * U = transpose_matrix(L);
  if(!U){fprintf(stderr, "ERROR: triangle: transpose failed!\n");return(-1);}

  // if KRONECKER, calculate the number of triangles 
  int64_t correct_answer = -1;
  if (args.gstd.model == KRONECKER) {
    correct_answer = calc_num_tri_kron_graph(args.gstd.kron_mode, args.gstd.kron_spec, args.gstd.kron_num);
    bale_app_write_int(&args.std, "known_num_triangles", correct_answer);
  }

  uint32_t use_model;
  double laptime = 0.0;
  int64_t tri_cnt;
  char model_str[64];
  int models_mask = (args.std.models_mask) ? args.std.models_mask : 3;
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    tri_cnt = 0;
    switch( use_model & models_mask ){
    case L_LxL:
      sprintf(model_str, "L .& (L * U)");
      laptime = triangles_matrix(&tri_cnt, L, L);
      break;
    case L_LxU:
      sprintf(model_str, "L .& (L * L)");
      laptime = triangles_matrix(&tri_cnt, L, U);
      break;
    default:
      continue;
    }
    bale_app_write_int(&args.std, "computed num tri   ", tri_cnt);
    bale_app_write_time(&args.std, model_str, laptime);
  }

  clear_matrix(L);
  clear_matrix(U);

  bale_app_finish(&args.std);
  return(0);
}

