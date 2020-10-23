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
/*! \file triangle.upc
 * \brief Demo application that counts triangles in a graph.
 */

#include "triangle.h"
#include <std_options.h>
/*!
  \page triangles_page Triangles

This uses matrix algebra approach to counting triangles in a graph.

The input is the lower-half of the adjacency matrix <b>A</b> for a
simple graph.  The matrix must be lower-triangular and have all zeros
on the diagonal. Further, it is a {0,1} matrix where the rows and cols
correspond to the vertices and \f$a_{ij} \f$ = <b>A[i][j]</b> is 1
exactly when there is a edge between vertices <i>v_i</i> and
<i>v_j</i>.

The triangle with vertices <i>{v_i, v_j, v_k}</i> has associated 
edges <i>{v_i, v_j}</i>, <i>{v_j, v_k}</i> and <i>{v_k, v_i}</i> 
which correspond to non-zero entries 
\f$a_{ij}\f$,
\f$a_{jk}\f$, and
\f$a_{ki}\f$
in the adjacency matrix.  
Hence the sum
\f$ \sum_{i,j,k} a_{ij}a_{jk}a_{ki} \f$ counts the triangles in the graph.
However, it counts each triangle 6 times according to the 6 symmetries of a triangle
and the 6 symmetric ways to choose the three nonzeros in <b>A</b>.
To count each triangle once, we compute the sum
\f[ \sum_{i=1}^{n}\sum_{j=1}^{i-1}\sum_{k=1}^{j-1} a_{ij}a_{jk}a_{ik} = 
    \sum_{i=1}^{n}\sum_{j=1}^{i-1} a_{ij} \sum_{k=1}^{j-1} a_{jk}a_{ik}. \f]

This picks out a unique labelling from the 6 possible and it means that 
all the information we need about edges is contained in the lower triangular 
part of symmetric adjacency matrix.  We call this matrix <b>L</b>.

The mathematical expression: 
for each nonzero \f$ a_{ij} \f$ compute the dot product 
of row \f$ i\f$ and row \f$ j \f$ becomes
\verbatim
  For each non-zero L[i][j] 
     compute the size of the intersection of the nonzeros in row i and row j
\endverbatim

Kronecker product graphs were implemented specifically for this app.
See "Design, Generation, and Validation of Extreme Scale Power-Law Graphs" by Kepner et. al.
for more information on Kronecker product graphs. 
 */

/*! \brief structure to hold counts of the number of pushes and pulls */
typedef struct push_pull_ctr_t{
  int64_t npushes;                //!< number of pushes
  int64_t npulls;                 //!< number of pulls
}push_pull_ctr_t;

typedef struct args_t{
  int alg;                        //!< which alg (pushing to or pulling from) remote rows
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'a': args->alg = atoi(arg); break;     
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"triangle_alg", 'a', "ALG", 0, "Algorithm: 0 means L&L*U, 1 means L&U*L"},  
    {0}
  };

static struct argp_child children_parsers[] =
  {    
    {&std_options_argp, 0, "Standard Options", -2},
    {&std_graph_options_argp, 0, "Standard Graph Options", -3},
    {0}
  };

push_pull_ctr_t calculate_num_pushes_pulls(sparsemat_t * L, args_t * args);

int main(int argc, char * argv[]) {

  double t1;
  int64_t i, j;

  /* process command line */
  args_t args = {0};  // initialize args struct to all zero
  struct argp argp = {options, parse_opt, 0,
                      "Parallel triangle counting.", children_parsers};  

  args.gstd.l_numrows = 100000;
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);
  
  //override command line (these will lead to matrices with not quite the right number of nonzeros
  // if the user also used the -z flag.
  if(args.gstd.loops == 1){
    T0_fprintf(stderr,"WARNING: triangles requires no-loops, overriding command line.\n");
    T0_fprintf(stderr,"Suggest running without this option.\n");
    args.gstd.loops = 0;
  }
  if(args.gstd.directed == 1){
    T0_fprintf(stderr,"WARNING: triangles requires undirected graph, overriding command line.\n");
    T0_fprintf(stderr,"Suggest running without this option.\n");
    args.gstd.directed = 0;
  }

  if(!MYTHREAD){
    write_std_graph_options(&args.std, &args.gstd);
    write_std_options(&args.std);
  }
  
  // read in a matrix or generate a random graph
  sparsemat_t * L = get_input_graph(&args.std, &args.gstd);
  if(!L){T0_fprintf(stderr, "ERROR: transpose: L is NULL!\n");return(-1);}

  /* make sure the matrix is in legal form */
  if(args.gstd.readfile){
    ret = tril(L, -1);
    if(ret){
      T0_fprintf(stderr,"WARNING: input graph was not lower-triangular with zero diagonal. Removing illegal nonzeros.\n");
    }
  }else if(!is_lower_triangular(L, 0)){
    T0_fprintf(stderr,"ERROR: L is not lower triangular!\n");
    lgp_global_exit(1);  
  }  

  if(args.std.dump_files) write_matrix_mm(L, "triangle_inmat");
  
  lgp_barrier();
  
  /* calculate the number of triangles */
  int64_t correct_answer = -1;
  int wrote_num_triangles = 0;
  if(args.gstd.model == KRONECKER){
    correct_answer = calc_num_tri_kron_graph(args.gstd.kron_mode, args.gstd.kron_spec, args.gstd.kron_num);
    bale_app_write_int(&args.std, "num_triangles", correct_answer);
    wrote_num_triangles = 1;
  }

  /* if we are using alg 1, we need U too */
  sparsemat_t * U;
  if(args.alg == 1) U = transpose_matrix(L);

  lgp_barrier();

  /* calculate number of pushes and pulls required for this matrix (just for kicks) */
  calculate_num_pushes_pulls(L, &args);
  
  int64_t use_model;
  double laptime = 0.0;
  char model_str[32];
  int error = 0;  
  for( use_model=1L; use_model < 32; use_model *=2 ) {

    int64_t tri_cnt = 0;           // partial count of triangles on this thread
    int64_t total_tri_cnt = 0;     // the total number of triangles on all threads
    int64_t sh_refs = 0;           // local number of shared reference or pushes
    int64_t total_sh_refs = 0;     

    switch( use_model & args.std.models_mask ) {
    case AGP_Model:
      sprintf(model_str, "AGP");
      laptime = triangle_agp(&tri_cnt, &sh_refs, L, U, args.alg); 
      break;
    
    case EXSTACK_Model:
      sprintf(model_str, "Exstack");
      laptime = triangle_exstack_push(&tri_cnt, &sh_refs, L, U, args.alg, args.std.buf_cnt);
      break;

    case EXSTACK2_Model:
      sprintf(model_str, "Exstack2");
      laptime = triangle_exstack2_push(&tri_cnt, &sh_refs, L, U, args.alg, args.std.buf_cnt);
      break;

    case CONVEYOR_Model:
      sprintf(model_str, "Conveyor");
      laptime = triangle_convey_push(&tri_cnt, &sh_refs, L, U, args.alg);
      break;

    case ALTERNATE_Model:
      sprintf(model_str, "Alternate");
      laptime = triangle_agp_iter(&tri_cnt, &sh_refs, L, U, args.alg);
      break;
    case 0:
      continue;
    }
    
    lgp_barrier();
    total_tri_cnt = lgp_reduce_add_l(tri_cnt);
    total_sh_refs = lgp_reduce_add_l(sh_refs);

    if(!wrote_num_triangles){
      bale_app_write_int(&args.std, "triangles", total_tri_cnt);
      wrote_num_triangles = 1;
    }
    
    bale_app_write_time(&args.std, model_str, laptime);    

    if(correct_answer == -1){
      correct_answer = total_tri_cnt; // set the answer we just got as the reference answer
    }else if(total_tri_cnt != correct_answer){
      T0_fprintf(stderr, "ERROR: Triangle: total_tri_cnt (%ld) != (%ld)\n", total_tri_cnt, correct_answer);
      error = 1;
    }
    
  }
  
  lgp_barrier();
  
  bale_app_finish(&args.std);
  return(error);
}

/*!
\brief calculate the number of remote push or pull needed depending on the algorithm used
\param L the lower triangular matrix
\param args  to get which algorithm is being used
\return a pointer to struct that holds both answers
*/
push_pull_ctr_t calculate_num_pushes_pulls(sparsemat_t * L, args_t * args){
  int64_t i;
  SHARED int64_t * cc = lgp_all_alloc(L->numrows, sizeof(int64_t));
  int64_t * l_cc = lgp_local_part(int64_t, cc);
  for(i = 0; i < L->lnumrows; i++)
    l_cc[i] = 0;
  lgp_barrier();
  
  /* calculate col sums */
  for(i = 0; i < L->lnnz; i++){
    lgp_fetch_and_inc(cc, L->lnonzero[i]);
  }
  
  lgp_barrier();
  
  int64_t rtimesc_calc = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t deg = L->loffset[i + 1] - L->loffset[i];        
    rtimesc_calc += deg*l_cc[i];
  }

  /* calculate sum (r_i choose 2) */
  int64_t rchoose2_calc = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t deg = L->loffset[i + 1] - L->loffset[i];
    rchoose2_calc += deg*(deg-1)/2;
  }
  
  /* calculate sum (c_i choose 2) */
  int64_t cchoose2_calc = 0;
  for(i = 0; i < L->lnumrows; i++){
    int64_t deg = l_cc[i];
    cchoose2_calc += deg*(deg-1)/2;
  }
  int64_t pulls_calc = 0;
  int64_t pushes_calc = 0;
  if(args->alg == 0){
    pulls_calc = lgp_reduce_add_l(rtimesc_calc);
    pushes_calc = lgp_reduce_add_l(rchoose2_calc);
  }else{
    pushes_calc = lgp_reduce_add_l(rtimesc_calc);
    pulls_calc = lgp_reduce_add_l(cchoose2_calc);
  }

  lgp_all_free(cc);

  if(!args->std.json)
    T0_fprintf(stderr,"Calculated: Pulls = %"PRId64"\n            Pushes = %"PRId64"\n\n",pulls_calc, pushes_calc);

  push_pull_ctr_t p;
  p.npushes = pushes_calc;
  p.npulls = pulls_calc;
  return(p);
}
