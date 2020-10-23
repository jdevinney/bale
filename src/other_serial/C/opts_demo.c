/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file opts_demo.c
example programs that play with commandline options
*/
  
/*!
\brief Examples to play with modifying and adding commandline options

This is really the same toy program written four different times,
each given by it own #ifdef.  Un-comment one of the defines
at a time to play with the different versions.
*/

#include "spmat_utils.h"
#include "std_options.h"

// Only un-comment one of these at time
//#define STANDARD_OPT_ONLY
//#define STANDARD_AND_GRAPH_OPTS
//#define APP_SPECIFIC_OPTS
#define STANDARD_AND_REUSE

#ifdef STANDARD_OPT_ONLY
// Compile the demo and try these
//  ./opts_demo
//  ./opts_demo --help
//  ./opts_demo -s 131
//  ./opts_demo -s 131 -M15
// This should make a json
//  ./opts_demo -s 131 -M15 -j opts_demo.json
// Then find models_mask and uncomment the line to set your own default 
// Next try this, which should fail until we add the GRAPH options.
// ./opts_demo -F
/********************************  argp setup  ************************************/
typedef struct args_t {
  std_args_t std;
} args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key) {
  case ARGP_KEY_INIT:
    state->child_inputs[0] = &args->std;
    break;
  }
  return(0);
}

static struct argp_option options[] = {
  {0}
};

static struct argp_child children_parsers[] = {
  {&std_options_argp, 0, "Standard Options", -2},
  {0}
};

int main(int argc, char * argv[]) 
{
  args_t args = {{0}};
  //args.std.models_mask = 63;  // default value for model_mask can be overridden here
  struct argp argp = {options, parse_opt, 0, "opts_demo with only standard options", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  write_std_options(&args.std);
  
  if(args.std.dump_files) printf("dumpfile is on\n");
  
  if(args.std.json) printf("All bale init/write rountines send everything to json file %s\n", args.std.json_output);

  double laptime = 3.141592;
  bale_app_write_time(&args.std, "playing with options", laptime);
  bale_app_finish(&args.std);
  return(0);
}
#endif
#ifdef STANDARD_AND_GRAPH_OPTS
// With the graph option included, these should now work.
//  ./opts_demo --help
//  ./opts_demo -G
//  ./opts_demo -s 131 -N 17 -M15 -j opts_demo.json
//
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
  args.std.models_mask = 63;  // default value for model_mask can be overridden here
  args.gstd.numrows = 100;
  struct argp argp = {options, parse_opt, 0, "opts_demo including graph options", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  write_std_options(&args.std);
  write_std_graph_options(&args.std, &args.gstd);
  
  if(args.std.dump_files) printf("dumpfile is on\n");
  
  if(args.std.json) printf("All bale init/write rountines send everything to json file %s\n", args.std.json_output);

  double laptime = 3.141592;
  bale_app_write_time(&args.std, "playing with options", laptime);
  bale_app_finish(&args.std);
  return(0);
}
#endif
#ifdef APP_SPECIFIC_OPTS
// You can over write the default values for options like with did above.
// If you want to add a option, find an available key and do the following
// and try:
// ./opts_demo --help
// ./opts_demo -W 2
// ./opts_demo -W 2 -Y"You bet ja"
/********************************  argp setup  ************************************/

typedef struct args_t {
  int64_t Wacky;  
  char Yes_str[128];
  std_args_t std;
  std_graph_args_t gstd;
} args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key) {
  case 'W': args->Wacky = atoi(arg); break;
  case 'Y': strcpy(args->Yes_str, arg); break;
  case ARGP_KEY_INIT:
    args->Wacky = 0;
    strcpy(args->Yes_str,"Yes sir.");
    state->child_inputs[0] = &args->std;
    state->child_inputs[1] = &args->gstd;
    break;
  }
  return(0);
}

static struct argp_option options[] = {
  {"Wacky",        'W', "NUM", 0, "Wacky FLAG"},
  {"Yes_str",      'Y', "YES",  0, "Yes string"},
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
  args.std.models_mask = 63;  // default value for model_mask can be overridden here
  args.gstd.numrows = 100;
  struct argp argp = {options, parse_opt, 0, "opts_demo standard, graph, and app specific options", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  write_std_options(&args.std);
  write_std_graph_options(&args.std, &args.gstd);
  
  if(args.std.dump_files) printf("dumpfile is on\n");
  
  if(args.std.json) printf("All bale init/write rountines send everything to json file %s\n", args.std.json_output);

  printf("%s! This is Wacky times %ld\n", args.Yes_str, args.Wacky);

  double laptime = 3.141592;
  bale_app_write_time(&args.std, "playing with options", laptime);
  bale_app_finish(&args.std);
  return(0);
}
#endif
#ifdef STANDARD_AND_REUSE
// For apps without graph options, histo, ig, randperm, 
// we wanted to make -N be "the order" of the work.
// Since -N was now be used it becomes available 
// for an app specific option
// try: 
//  ./opts_demo --help
//  ./opts_demo -s 131
/********************************  argp setup  ************************************/

typedef struct args_t {
  int64_t num_things;
  std_args_t std;
} args_t;

static int parse_opt(int key, char * arg, struct argp_state * state)
{
  args_t * args = (args_t *)state->input;
  switch(key) {
  case 'N': args->num_things = atoi(arg); break;
  case ARGP_KEY_INIT:
    state->child_inputs[0] = &args->std;
    break;
  }
  return(0);
}

static struct argp_option options[] = {
  {"num_things",      'N', "NUM",  0, "numthings, not the graph numrows"},
  {0}
};

static struct argp_child children_parsers[] = {
  {&std_options_argp, 0, "Standard Options", -2},
  {0}
};


int main(int argc, char * argv[]) 
{
  args_t args = {0};
  args.std.models_mask = 63;  // default value for model_mask can be overridden here
  struct argp argp = {options, parse_opt, 0, "opts_demo reuse hack", children_parsers};
  argp_parse(&argp, argc, argv, 0, 0, &args);
  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if (ret < 0) return(ret);
  else if (ret) return(0);

  write_std_options(&args.std);
  
  if(args.std.dump_files) printf("dumpfile is on\n");
  
  if(args.std.json) printf("All bale init/write rountines send everything to json file %s\n", args.std.json_output);

    printf("option -N now means num_things = %ld\n", args.num_things);

  double laptime = 3.141592;
  bale_app_write_time(&args.std, "playing with options", laptime);
  bale_app_finish(&args.std);
  return(0);
}
#endif
