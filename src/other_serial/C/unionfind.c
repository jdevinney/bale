/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/


/*! \file unionfind.c
\brief Program that uses a unionfind data structure to find connected components in a graph

Run unionfind --help or --usage for insructions on running.
 */

#include "spmat_utils.h"
#include "std_options.h"
#include "default_app_opts.h"


/*! \brief the nodes of the disjoint union trees */
typedef struct comp_tree_t {
  int64_t parent;  //!< pointer to the nodes parent in the tree
  int64_t rank;    //!< if the node is the root, size of the component
} comp_tree_t;

/*! 
\brief Find the name (the root of the parent tree) for the component containing a given node.
\param cc the component tree
\param x a given node
\return the index of the root of the tree
 */
static int64_t comp_find( comp_tree_t *cc,  int64_t x)
{
  int64_t p;

  p = cc[x].parent;
  if( p == x )
    return(x);

  cc[x].parent = comp_find( cc, p );
  return( cc[x].parent );
};

/*! 
\brief Merge or take the union of two component trees, based on the rank of the components.
\param cc the component tree
\param r the root of one tree
\param s the root of the other
\param verbose print flag

This is the best way to join the trees.
*/
static void comp_rank_union(comp_tree_t *cc, int64_t r, int64_t s, int verbose)
{
   if( r == s )
     return;
   if( cc[r].rank < cc[s].rank ) {
     if(verbose) printf("set parent of %"PRId64" to %"PRId64"\n", r, s); 
     cc[r].parent = s;
   } else if( cc[s].rank < cc[r].rank ) {
     if(verbose) printf("set parent of %"PRId64" to %"PRId64"\n", s, r); 
     cc[s].parent = r;
   } else {
     if(verbose) printf("set parent of %"PRId64" to %"PRId64" rank++\n", r, s); 
     cc[r].parent = s;
     cc[s].rank  += 1;
   }
}

/*!
\brief Merge or take the union of two component trees, based on the rank of the components.
\param cc the component tree
\param r the root of one tree
\param s the root of the other
\param e the node that caused us to realize the components were connected,
         essentially a random node in one of the trees
\param verbose print flag

This is a bad way to join the trees.
*/
static void comp_bad_union(comp_tree_t *cc, int64_t r, int64_t s, int64_t e, int verbose)
{
   if( r == s )
     return;
   if(verbose) printf("set root %"PRId64" to limb  %"PRId64"\n", r, e); 
   cc[r].parent = e;
}

/*!
\brief debugging routine that prints the comp_tree struct
\param *prefix a string one can use to keep the output straight
\param *cc the comp_tree_t data structure 
\param numverts the number of vertices 
\param verbose print flag
*/
void dump_comp_tree( char *prefix, comp_tree_t *cc, int64_t numverts, int verbose)
{
  int64_t i;

  if(verbose) printf("tree: %s ",prefix);
  for(i=0; i<numverts; i++) {
    if(verbose) printf(" %2"PRId64"",cc[i].parent);
  }
  if(verbose) printf("\n");
}

/*!
\brief This routine finds the connected components in a graph
\param *numcomps address to return the number of components
\param *cc the comp_tree_t data structure 
\param *mat the matrix that specifies the graph
\param which union to use: the one based on rank or just them wherever you land on it
\param verbose print tracing statements
\return run time
 */
double concomp(int64_t *numcomps, comp_tree_t * cc, sparsemat_t *mat, int which_union, int verbose)
{
  int64_t i, j, k, r, s;
  
  double t1 = wall_seconds();

  *numcomps = 0;
  for( i = 0; i< mat->numrows; i++){
    cc[i].parent = i;
    cc[i].rank   = 1;
  }
   
  dump_comp_tree(" ",cc, mat->numrows, verbose);
  for(i = 0; i < mat->numrows; i++){ 
    for(k=mat->offset[i]; k < mat->offset[i+1]; k++) {
      j = mat->nonzero[k];
  
      if( verbose > 0 ) printf("looking at (%"PRId64",%"PRId64"): ", i, j);
  
      r = comp_find(cc,i);
      s = comp_find(cc,j);
  
      if( verbose > 0 ) printf("parentof %"PRId64" is %"PRId64" , parentof %"PRId64" is %"PRId64"\n", i,r,j,s);
      if( verbose > 1 ) dump_comp_tree(">", cc, mat->numrows, verbose);
  
      if( which_union == 0 ) {
        comp_bad_union(cc, r, s, j, verbose);
      } else {
        comp_rank_union(cc, r, s, verbose);
      }
      if( verbose > 1 ) dump_comp_tree("<",cc, mat->numrows, verbose);
      if( verbose > 0 ) printf("parentof %"PRId64" is %"PRId64" , parentof %"PRId64" is %"PRId64"\n", r,comp_find(cc,r),s,comp_find(cc,s));
    }
  }
  for(i = 0; i < mat->numrows; i++){ 
    cc[i].parent = comp_find(cc,i);
  }
  dump_comp_tree("-",cc, mat->numrows, verbose);
  
  int64_t *comp_size = (int64_t *) calloc(mat->numrows, sizeof(int64_t));
  for(i = 0; i < mat->numrows; i++){ 
    comp_size[cc[i].parent] += 1;
  }
  
  for(i = 0; i < mat->numrows; i++){ 
    if( comp_size[i] > 0 ){
      *numcomps += 1;
      if(verbose) printf("component %"PRId64" has size %"PRId64"\n", i, comp_size[i]);
    }
  }
  
  t1 = wall_seconds() - t1;
  return(t1);
}


/********************************  argp setup  ************************************/
typedef struct args_t{
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

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
  enum FLAVOR {BAD_UNION=1, OPT_UNION=2, ALL_Models=4};
  args.std.models_mask = ALL_Models-1;
  struct argp argp = {options, parse_opt, 0, "Union Find for connected components", children_parsers};
  args.gstd.numrows = UNIONFIND_NUM_ROWS;
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);
  write_std_options(&args.std);
  write_std_graph_options(&args.std, &args.gstd);
  
  //override command line (these will lead to matrices with not quite the right number of nonzeros
  // if the user also used the -z flag.
  if ( (args.gstd.loops == 1) || (args.gstd.directed == 1) ) {
    fprintf(stderr,"WARNING: unionfind counting requires undirected no-loops graph\n");
    args.gstd.loops = 0;
    args.gstd.directed = 0;
  }


  sparsemat_t * mat = get_input_graph(&args.std, &args.gstd);
  if(!mat){fprintf(stderr, "ERROR: unionfind: mat is NULL!\n");return(-1);}

  comp_tree_t * cc = calloc(mat->numrows, sizeof(comp_tree_t));
  
  if(args.std.dump_files) write_matrix_mm(mat, "unionfind_inmat");
  int verbose=0;

  uint32_t use_model;  
  double laptime;
  int64_t num_components = 0;
  char model_str[64];
  int models_mask = (args.std.models_mask) ? args.std.models_mask : 3;
  for( use_model=1; use_model < ALL_Models; use_model *=2 ) {
    switch( use_model & models_mask ) {
    case BAD_UNION:
      sprintf(model_str, "unionfind: bad union :");
      laptime = concomp(&num_components, cc, mat, 0, verbose);
      break;
    case OPT_UNION:
      sprintf(model_str, "unionfind: opt union :");
      laptime = concomp(&num_components, cc, mat, 1, verbose);
      break;
    default:
      continue;
    }
    bale_app_write_int(&args.std, "number of components:     ", num_components);
    bale_app_write_time(&args.std, model_str, laptime);
  }
  
  clear_matrix(mat);
  bale_app_finish(&args.std);
  return(0);
}

