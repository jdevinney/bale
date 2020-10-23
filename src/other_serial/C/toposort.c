/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file toposort.c
\brief Find permutations to verify that a given matrix is morally upper triangular

Run topo --help or --usage for insructions on running.
*/

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_opts.h"

/*! \brief check the result toposort 
\param mat the original matrix
\param rperminv the row permutation
\param cperminv the column permutation
\param dump_files debugging flag
\return 0 on success, 1 otherwise

Check that the permutations are in fact permutations and check that applying
them to the original matrix yields an upper triangular matrix.
*/
int check_result(sparsemat_t * mat, int64_t * rperminv, int64_t * cperminv, int64_t dump_files) 
{
  sparsemat_t * mat2;
  int ret = 0;
  
  int64_t rf = is_perm(rperminv, mat->numrows);
  int64_t cf = is_perm(cperminv, mat->numcols);
  if(!rf || !cf){
    fprintf(stderr,"ERROR: check_result is_perm(rperminv2) = %"PRId64" is_perm(cperminv2) = %"PRId64"\n",rf,cf);
    return(1);
  }
  mat2 = permute_matrix(mat, rperminv, cperminv);
  if(!is_upper_triangular(mat2, 1))
    ret = 1;
  if(dump_files) 
    dump_matrix(mat2, 20, "mat2.out");
  clear_matrix(mat2);
  free(mat2);
  return(ret);
}

/*! \brief Generate a matrix that is the a random permutation of a sparse uppper triangular matrix.
\param sargs the standard command line arguments
\param gargs the graph arguments line arguments
\return a pointer to the permuted upper triangular matrix

Either read in a matrix and force to have the right form or generate a matrix with the sparse matrix library.
Then generate and apply row and column permutations.
*/
sparsemat_t *generate_toposort_input(std_args_t *sargs, std_graph_args_t *gargs)
{
  int64_t nr = gargs->numrows;
  sparsemat_t * L = get_input_graph(sargs, gargs);
  if (!L) { fprintf(stderr, "ERROR: topo: L is NULL!\n");return(NULL); }
  sparsemat_t * U = transpose_matrix(L);
  if (!U) { fprintf(stderr, "ERROR: topo: U is NULL!\n");return(NULL); }
  if (!is_upper_triangular(U, 1)) {
    fprintf(stderr,"ERROR: generate_toposort did not start with an upper triangular\n");
    return(NULL);
  }
  clear_matrix(L); free(L);

  if(sargs->dump_files) write_matrix_mm(U, "topo_orig");
  
  // get random row and column permutations
  int64_t * rperminv = rand_perm(nr, 1234);
  int64_t * cperminv = rand_perm(nr, 5678);
  if(!rperminv || !cperminv){
    printf("ERROR: generate_toposort_input: rand_perm returned NULL!\n");
    exit(1);
  }
  if(sargs->dump_files){
    dump_array(rperminv, nr, 20, "topo_rperm.out");
    dump_array(cperminv, nr, 20, "topo_cperm.out");
  }
  
  sparsemat_t * mat = permute_matrix(U, rperminv, cperminv);
  if(!mat) {
    printf("ERROR: generate_toposort_input: permute_matrix returned NULL");
    return(NULL);
  }
  
  clear_matrix( U ); free( U );
  free(rperminv);
  free(cperminv);
  
  return(mat);
}


/*! \brief This routine implements a variant of toposort that enqueues the pivot positions as they are discovered
\param *rperm a place to return the row permutation that is found
\param *cperm a place to return the column permutation that is found
\param *mat the input sparse matrix
\param *tmat the transpose of mat
\return average run time
*/
double toposort_matrix_queue(int64_t *rperm, int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat) 
{
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
  
  int64_t * queue  = calloc(nr, sizeof(int64_t));
  int64_t * rowtrck = calloc(nr, sizeof(int64_t));
  
  int64_t start, end;
  int64_t i, k, row, col, t_row;
   
  /* initialize rowsum, rowcnt, and queue (queue holds degree one rows) */
  start = end = 0;
  for(i = 0; i < mat->numrows; i++){
    rowtrck[i] = 0L;
    for(k = mat->offset[i]; k < mat->offset[i+1]; k++)
      rowtrck[i] += (1L<<32) + mat->nonzero[k];
    if((rowtrck[i] >> 32) ==  1)
      queue[end++] = i;    
  }
  
  // we a pick a row with a single nonzero = col.
  // setting rperm[pos] = row and cprem[pos] = col
  // in a sense 'moves' the row and col to the bottom
  // right corner of matrix.
  // Next, we cross out that row and col by decrementing 
  //  the rowcnt for any row that contains that col
  // repeat
  //
  double t1 = wall_seconds();
  
  int64_t n_pivots = 0;
  while(start < end){      
    row = queue[start++];
    col = rowtrck[row] & 0xFFFF;  // see cool trick
    
    rperm[row] = nr - 1 - n_pivots;
    cperm[col] = nc - 1 - n_pivots;
    n_pivots++;
  
    // look at this column (tmat's row) to find all the rows that hit it
    for(k=tmat->offset[col]; k < tmat->offset[col+1]; k++) {
      t_row = tmat->nonzero[k];
      assert((t_row) < mat->numrows);
      rowtrck[t_row] -= (1L<<32) + col;
      if( (rowtrck[t_row] >> 32) == 1L ) {
        queue[end++] = t_row;
      }
    }
  }
  
  t1 = wall_seconds() - t1;
  
  if(n_pivots != nr){
    printf("ERROR! toposort_matrix_queue: found %"PRId64" pivots but expected %"PRId64"!\n", n_pivots, nr);
    exit(1);
  }
  free(queue);
  free(rowtrck);
  return(t1);
}

/*! \brief This routine implements toposort by continually looping over rows looking for a pivot.
\param *rperm a place to return the row permutation that is found
\param *cperm a place to return the column permutation that is found
\param *mat the input sparse matrix
\param *tmat the transpose of mat
\return average run time
 */
double toposort_matrix_loop(int64_t *rperm, int64_t *cperm, sparsemat_t *mat, sparsemat_t *tmat) 
{
  int64_t nr = mat->numrows;
  int64_t nc = mat->numcols;
  
  int64_t * rowtrck = calloc(nr, sizeof(int64_t));
  
  int64_t i, j, col, t_row;
   
  /* initialize rowtrck */
  for(i = 0; i < nr; i++){
    rowtrck[i] = 0L;
    for(j = mat->offset[i]; j < mat->offset[i+1]; j++)
      rowtrck[i] += (1L<<32) + mat->nonzero[j];
  }
  
  // we a pick a row with a single nonzero = col.
  // setting rperm[pos] = row and cprem[pos] = col
  // in a sense 'moves' the row and col to the bottom
  // right corner of matrix.
  // Next, we cross out that row and col by decrementing 
  //  the rowcnt for any row that contains that col
  // repeat
  //
  double t1 = wall_seconds();
  
  int64_t n_pivots = 0;
  while(n_pivots < nr){      
    for(i = 0; i < nr; i++){
      if( (rowtrck[i] >> 32) == 1 ){
        col = rowtrck[i] & 0xFFFF;  // see cool trick
        rperm[i] = nr - 1 - n_pivots;
        cperm[col] = nc - 1 - n_pivots;
        n_pivots++;
        rowtrck[i] = 0L;
		
        // look at this column (tmat's row) to find all the rows that hit it
        for(j=tmat->offset[col]; j < tmat->offset[col+1]; j++) {
          t_row = tmat->nonzero[j];
          assert((t_row) < mat->numrows);
          rowtrck[t_row] -= (1L<<32) + col;
        }
      }
    }
  }
  
  t1 = wall_seconds() - t1;
  
  if(n_pivots != nr){
    printf("ERROR! toposort_matrix_queue: found %"PRId64" pivots but expected %"PRId64"!\n", n_pivots, nr);
    exit(1);
  }
  free(rowtrck);
  return(t1);
}


/********************************  argp setup  ************************************/
typedef struct args_t{
  int alg;
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch (key) {
  case 'a': args->alg = atoi(arg); break;     
  case ARGP_KEY_INIT:
    state->child_inputs[0] = &args->std;
    state->child_inputs[1] = &args->gstd;
    break;
  }
  return(0);
}

static struct argp_option options[] = {
  {"toposort", 'a', "ALG", 0, "Algorithm: 0 means loops, 1 means queue"},  
  {0}
};

static struct argp_child children_parsers[] = {
  {&std_options_argp, 0, "Standard Options", -2},
  {&std_graph_options_argp, 0, "Standard Graph Options", -3},
  {0}
};


int main(int argc, char * argv[])
{
  args_t args = {0};  
  struct argp argp = {options, parse_opt, 0, "Toposort", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  args.gstd.numrows = 500;
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);
  args.gstd.loops = 1;      //override command line: user cant' pick these
  args.gstd.directed = 0;
  write_std_options(&args.std);
  write_std_graph_options(&args.std, &args.gstd);
  
  sparsemat_t * mat = generate_toposort_input(&args.std, &args.gstd);
  if(!mat){printf("ERROR: topo: generate_toposort_input failed\n"); exit(1);}
  
  sparsemat_t * tmat = transpose_matrix(mat);
  if(!tmat){printf("ERROR: topo: transpose_matrix failed\n"); exit(1);}

  if (args.std.dump_files) {
    write_matrix_mm(mat, "topo_inmat");
    write_matrix_mm(tmat, "topo_tmat.mm");
    dump_matrix(mat,20, "mat.out");
    dump_matrix(tmat,20, "trans.out");
  }

  uint32_t use_model;
  double laptime = 0.0;
  enum FLAVOR {GENERIC=1, LOOP=2, ALL=4};
  // arrays to hold the row and col permutations
  int64_t *rperminv2 = calloc(mat->numrows, sizeof(int64_t));
  int64_t *cperminv2 = calloc(mat->numrows, sizeof(int64_t));
  char model_str[64];
  int models_mask = (args.std.models_mask) ? args.std.models_mask : 3;
  for( use_model=1; use_model < ALL; use_model *=2 ) {
    switch( use_model & models_mask ) {
    case GENERIC:
      sprintf(model_str,"queue toposort");
      laptime = toposort_matrix_queue(rperminv2, cperminv2, mat, tmat);
      break;
    case LOOP:
      sprintf(model_str,"loop  toposort");
      laptime = toposort_matrix_loop(rperminv2, cperminv2, mat, tmat);
      break;
    }
    if( check_result(mat, rperminv2, cperminv2, args.std.dump_files) ) {
      fprintf(stderr,"\nERROR: After toposort_matrix_queue: mat2 is not upper-triangular!\n");
      exit(1);
    }
    bale_app_write_time(&args.std, model_str, laptime);
  }
  bale_app_finish(&args.std);
  return(0);
}

