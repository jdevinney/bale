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
/*! \file sssp.upc
 * \brief Implementation of Single Source Shortest Path algorithms
 */

#include "sssp.h"
#include <std_options.h>

/*!
  \page sssp_page Single Source Shortest Path
  Demo Program that runs single source shortest path algoriths. See README for more info.
*/

/*! \brief debugging rountine to dump the d_array holding the tentative weights
 * \param str a string to prefix the line of weights
 * \param tent the d_array of tentative weights
 */
void dump_tent(char *str, d_array_t *tent)
{
  int64_t i;
  if( MYTHREAD == 0 ){
    printf("%s ", str);
    for(i=0; i<tent->num; i++){
      printf(" %lg ", lgp_get_double(tent->entry, i) );
    }
    printf("\n");
  }
}

/*!
 * \brief Compare two d_arrays
 * \param *A one array (as a vector)
 * \param *B the other (as a vector)
 * \return the l_2 norm of the given arrays
 */
double sssp_answer_diff(d_array_t *A, d_array_t *B)
{
  int64_t i;
  double ldiff = 0.0;

  if(A->num != B->num)
    return(INFINITY);

  for(i=0; i<A->lnum; i++) {
    if( A->lentry[i] == INFINITY && B->lentry[i] == INFINITY )
      continue;
    ldiff += (A->lentry[i] - B->lentry[i]) * (A->lentry[i] - B->lentry[i]);
  }
  return(sqrt(lgp_reduce_add_d(ldiff)));
}

/********************************  argp setup  ************************************/
typedef struct args_t {
  double deltaStep;         //!< command line supplied delta step value
  int64_t V0;               //!< the starting vertex
  int64_t alg;              //!< which algorithm (1=Bellman-Ford | 2=Delta-stepping)
  std_args_t std;
  std_graph_args_t gstd;
} args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'a': args->alg = atoi(arg); break;
    case 'S': args->deltaStep = atof(arg); break;     
    case 'V': args->V0 = atoi(arg); break;
    case ARGP_KEY_INIT:
      args->deltaStep = 0.0;
      args->V0 = 0;
      args->alg = 3;
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_option options[] = {
  {"alg", 'a', "flag", 0, "alg: 1==bellman | 2==delta"},
  {"deltaStep", 'S', "STEPSIZE", 0, "user supplied delta step size"},  
  {"V0", 'V', "NUM", 0, "initial vertex"},  
  {0}
};

static struct argp_child children_parsers[] = {    
  {&std_options_argp, 0, "Standard Options", -2},
  {&std_graph_options_argp, 0, "Standard Graph Options", -3},
  {0}
};

// use the command line switch -a (alg) as an extension of the -M switch
#define USE_BELLMAN (32) //!< run the Bellman-Ford algorithm
#define USE_DELTA   (64) //!< run the Delta-Stepping algorithm
  

int main(int argc, char * argv[])
{
  double t1;
  int64_t i, j;

  /* process command line */
  args_t args = {0};  // initialize args struct to all zero
  struct argp argp = {options, parse_opt, 0,
                      "Parallel Single Source Shortest Path (SSSP).", children_parsers};  

  // set reasonable default for lnumrows and required defaults for graph params
  args.gstd.l_numrows = 100000;  
  args.gstd.directed = 1;
  args.gstd.weighted = 1;

  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);

  // SSSP only applies to weighted directed graphs with no loops
  if(args.gstd.loops){
    T0_fprintf(stderr,"WARNING: SSSP requires no self loops. overriding -l flag\n");
    args.gstd.loops = 0;
  }

  if(!MYTHREAD){
    write_std_graph_options(&args.std, &args.gstd);
    write_std_options(&args.std);
  }
  
  // read in a matrix or generate a random graph
  sparsemat_t * mat = get_input_graph(&args.std, &args.gstd);
  if(!mat){T0_fprintf(stderr, "ERROR: sssp: mat is NULL!\n");return(-1);}

  if(args.V0 < 0 || args.V0 >= mat->numrows){
    T0_fprintf(stderr,"Setting V0 to 0\n");
    args.V0 = 0;
  }
  if(args.alg < 1 || args.alg > 3){
    T0_fprintf(stderr,"Setting alg to 3\n");
    args.alg = 3;
  }
  
  if(args.std.dump_files) write_matrix_mm(mat, "sssp_inmat");

  lgp_barrier();
  
  d_array_t * tent = init_d_array(mat->numrows);
  d_array_t * comp_tent = NULL;

  uint64_t use_model, use_alg, alg;
  double laptime = 0.0;
  char model_str[64];

  //T0_printf("delta step = %lf\n", args.deltaStep);
  bale_app_write_double(&args.std, "delta", args.deltaStep);
  
  // To our understanding, the AGP model requires an atomic min of a double.
  // Until we get it, we haven't taken the AGP model out of the models_mask.
  lgp_barrier();
  args.std.models_mask &=0xFE;
  int64_t bz = args.std.buf_cnt;
  int64_t V0 = args.V0;
  double delta = args.deltaStep;

  for( alg = 1; alg < 3; alg *=2 ){
    use_alg = (args.alg & alg) * 32;
    for( use_model=2L; use_model < 32; use_model *=2 ) {
      model_str[0] = '\0';
      switch( (use_model & args.std.models_mask) + use_alg ) {

      case (EXSTACK_Model + USE_BELLMAN):
        strcpy(model_str, "Bellman-Ford Exstack");
        set_d_array(tent, INFINITY);
        laptime = sssp_bellman_exstack(tent, mat, bz, V0);
        break;

      case (EXSTACK_Model + USE_DELTA):
        strcpy(model_str, "Delta Exstack");
        set_d_array(tent, INFINITY);
        laptime = sssp_delta_exstack(tent, mat, bz, V0, delta);
        break;

      case (EXSTACK2_Model + USE_BELLMAN):
        strcpy(model_str, "Bellman-Ford Exstack2");
        set_d_array(tent, INFINITY);
        laptime = sssp_bellman_exstack2(tent, mat, bz, V0);
        break;

      case (EXSTACK2_Model + USE_DELTA):
        strcpy(model_str, "Delta Exstack2");
        set_d_array(tent, INFINITY);
        laptime = sssp_delta_exstack2(tent, mat, bz, V0, delta);
        break;
        
      case (CONVEYOR_Model + USE_BELLMAN):
        strcpy(model_str, "Bellman-Ford Conveyor");
        set_d_array(tent, INFINITY);
        laptime = sssp_bellman_convey(tent, mat, V0);
        break;

      case (CONVEYOR_Model + USE_DELTA):
        strcpy(model_str, "Delta Conveyor");
        set_d_array(tent, INFINITY);
        laptime = sssp_delta_convey(tent, mat, V0, delta);
        break;

      case (ALTERNATE_Model + USE_BELLMAN):
        strcpy(model_str, "Bellman-Ford AGP");
        set_d_array(tent, INFINITY);
        laptime = sssp_bellman_agp(tent, mat, V0); 
        break;
      default:
        continue;
      }
      if(model_str[0]) {
        if(comp_tent == NULL){
          comp_tent = copy_d_array(tent);
          strcat(model_str, " ()");
        }else{
          if( sssp_answer_diff(comp_tent, tent) < 1.0e-8)
            strcat(model_str, " (compares)");
        }
        lgp_barrier();
        bale_app_write_time(&args.std, model_str, laptime);
      }
    }
  }
  
  lgp_barrier();

  if(args.std.dump_files && !MYTHREAD){
    FILE * fp = fopen("sssp_dist", "w");
    for(i = 0; i < tent->num; i++)
      fprintf(fp, "%lf\n", lgp_get_double(tent->entry, i));
    fclose(fp);
  }
  lgp_barrier();
    
  clear_d_array(tent); free(tent);
  if(comp_tent) clear_d_array(comp_tent); free(comp_tent);
 
  bale_app_finish(&args.std);

  return(0);
}
