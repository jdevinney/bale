/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file sssp.c
\brief The application that calls different Single Source Shortest Path alogrithms.

Run sssp --help or --usage for insructions on running
*/

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_opts.h"

double sssp_dijsktra_linear(d_array_t * tent, sparsemat_t * mat, int64_t v0);            //!< Dijsktra's alg with linear search for smallest tent[v]
double sssp_dijsktra_heap(d_array_t * tent, sparsemat_t * mat, int64_t r0);              //!< Dijsktra's alg with heap to track the smallest tent[v]
double sssp_bellmanford_simple(d_array_t * tent, sparsemat_t *mat, int64_t r0);          //!< Bellman-Ford simple with one tentative array
double sssp_bellmanford_dynprog(d_array_t * tent, sparsemat_t *mat, int64_t r0);         //!< Bellman-Ford dynamic program with n squared storage
double sssp_bellmanford(d_array_t * tent, sparsemat_t *mat, int64_t r0);                 //!< Bellman-Ford recycling tentative arrays from dynamic program
double sssp_delta_stepping(d_array_t * tent, sparsemat_t *mat, int64_t r0, double del);  //!< Delta-Stepping algorithm

/*! \brief Compare two arrays of doubles
\param *A one array (vector)
\param *B the other
\return the l_2 norm of the given arrays
*/
static double sssp_answer_diff(d_array_t *A, d_array_t *B)
{
  int64_t i;
  double diff = 0.0;

  for(i=0; i<A->num; i++) {
    if( A->entry[i] == INFINITY && B->entry[i] == INFINITY )
      continue;
    diff += (A->entry[i] - B->entry[i]) * (A->entry[i] - B->entry[i]);
  }
  return(sqrt(diff));
}

/********************************  argp setup  ************************************/
typedef struct args_t {
  double deltaStep; 
  int64_t V0;
  int64_t compare;
  char compare_file[256];
  std_args_t std;
  std_graph_args_t gstd;
} args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key) {
  case 'S': args->deltaStep = atof(arg); break; 
  case 'V': args->V0 = atoi(arg); break;
  case 'C': args->compare = 1; strcpy(args->compare_file, arg); break;
          
  case ARGP_KEY_INIT:
    args->deltaStep = 0.0;
    args->V0 = 0;
    args->compare = 0;
    state->child_inputs[0] = &args->std;
    state->child_inputs[1] = &args->gstd;
    break;
  }
  return(0);
}

static struct argp_option options[] = {
  {"deltaStep", 'S', "STEPSIZE", 0, "user supplied delta step size"},  
  {"V0",        'V',      "NUM", 0, "initial vertex"},  
  {"compare",   'C',     "FILE", 0, "read compare weights from a file"},  
  {0}
};

static struct argp_child children_parsers[] = {    
  {&std_options_argp, 0, "Standard Options", -2},
  {&std_graph_options_argp, 0, "Standard Graph Options", -3},
  {0}
};


int main(int argc, char * argv[]) 
{
  args_t args={0};  
  enum MODEL {DIJSKTRA_HEAP=1, DELTA_STEPPING=2, BELLMAN=4, ALTERNATE=8, ALL_Models=8};
  args.std.models_mask = ALL_Models - 1;
  args.gstd.numrows = SSSP_NUM_ROWS;
  struct argp argp = {options, parse_opt, 0, "SSSP for a weighted graph.", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);
  args.gstd.loops = 0;    //override command line: the user can't pick these
  args.gstd.directed = 1;
  args.gstd.weighted = 1;
  write_std_options(&args.std);
  write_std_graph_options(&args.std, &args.gstd);
  
  // read in a matrix or generate a random graph
  sparsemat_t * mat = get_input_graph(&args.std, &args.gstd);
  if(!mat){fprintf(stderr, "ERROR: SSSP: mat is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(mat, "sssp_inmat");

  d_array_t *comp_weight=NULL;
  if (args.compare) {
    comp_weight = read_d_array(args.compare_file);
  }
  d_array_t *weight = init_d_array(mat->numrows);
  uint32_t use_model;
  double laptime = 0.0;
  char model_str[64];
  for(use_model=1; use_model < ALL_Models; use_model *=2 ){
    model_str[0] = '\0';
    switch( use_model & args.std.models_mask ){
    case DIJSKTRA_HEAP:
      strcat(model_str, "Dijsktra Heap");
      laptime = sssp_dijsktra_heap(weight, mat, args.V0);
      break;
    case DELTA_STEPPING:
      strcat(model_str, "Delta Stepping");
      laptime = sssp_delta_stepping(weight, mat, args.V0, 0.0);
      break;
    case BELLMAN:
      strcat(model_str, "Bellman Ford");
      laptime = sssp_bellmanford(weight, mat, args.V0);
      break;
    case ALTERNATE:
      strcat(model_str, "Alernate");
      //set_d_array(weight, INFINITY);
      //strcat(model_str, "Dijsktra Linear");
      //laptime = sssp_dijsktra_linear(weight, mat, args.V0);
      //strcat(model_str, "Simple Bellman");
      //laptime = sssp_bellmanford_simple(weight, mat, args.V0);
      //strcat(model_str, "Bellman Dynprog");
      //laptime = sssp_bellmanford_dynprog(weight, mat, args.V0);
      break;
    }
    if(model_str[0]) {
      if(comp_weight == NULL){
        comp_weight = copy_d_array(weight);
        strcat(model_str, " ()");
      }else{
        if( sssp_answer_diff(comp_weight, weight) < 1.0e-8)
          strcat(model_str, " (compares)");
      }
      bale_app_write_time(&args.std, model_str, laptime);
    }
  }

  if(args.std.dump_files){
    write_d_array(comp_weight, "SSSP weights to each node", "ssspout.dst");
  }
  
  clear_matrix(mat); free(mat);
  clear_d_array(weight); free(weight);
  clear_d_array(comp_weight); free(comp_weight);
  bale_app_finish(&args.std);
  return(0);
}
