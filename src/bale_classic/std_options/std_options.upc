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

/*! \file std_options.upc
 * \brief The main source for the std_options library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include <libgetput.h>
#include "std_options.h"

// the arg_parse function for std_args_t.
static int std_parse_opt(int key, char * arg, struct argp_state * state){
  std_args_t * args = (std_args_t *)state->input;
  switch(key){
  case 'b': args->buf_cnt = atol(arg); break;
  case 'c': args->cores_per_node = atoi(arg); break;
  case 'D': args->dump_files = 1; break;
  case 'j': args->json = 1; strcpy(args->json_output,arg); break;
  case 'M': args->models_mask = atol(arg); break;
  case 's': args->seed = atol(arg); break;
  case ARGP_KEY_INIT:
    args->buf_cnt = 1024;
    args->cores_per_node = 0;
    args->seed = 122222;
    args->models_mask = ALL_Models;
    args->dump_files = 0;
    args->json = 0;
    break;
  }
  return(0);
}

static struct argp_option std_options[] =
  {
    {"buf_cnt",   'b', "BUF", 0, "Exstack or exstack2 buffer size."},
    {"cores_per_node",'c', "CPN", 0, "Specify cores per node for network injection rate statistics"},
    {"dump_files",    'D', 0,     0, "Dump files for debugging"},
    {"json_output",   'j', "FILE",0, "Output results to a json file, rather than to stderr"},
    {"models_mask",   'M', "MASK",0, "Which implementations to run. 1: AGP, 2: Exstack, 4: Exstack2, 8: Conveyor, 16: alternatives. You can add these together. For example, -M 7 means AGP, Exstack, and Exstack2."},
    {"seed",          's', "SEED",0, "Seed for RNG"},
    {0}
  };

struct argp std_options_argp = {
  std_options, std_parse_opt, 0, 0, 0
};

// This function writes some of the standard options and their values to the screen (or json file).
void write_std_options(std_args_t * sargs){
  if(sargs->json == 0){
    fprintf(stderr,"Standard options:\n");
    fprintf(stderr,"----------------------------------------------------\n");
    fprintf(stderr,"buf_cnt (buffer size)    (-b): %"PRId64"\n", sargs->buf_cnt);
    fprintf(stderr,"seed                     (-s): %"PRId64"\n", sargs->seed);
    fprintf(stderr,"cores_per_node           (-c): %d\n", sargs->cores_per_node);
    fprintf(stderr,"Models Mask              (-M): %d\n\n", sargs->models_mask);
  }else{
    FILE * jp = fopen(sargs->json_output, "a");
    fprintf(jp, "\"buf_cnt\": \"%"PRId64"\",\n", sargs->buf_cnt);
    fclose(jp);
  }
}


// the arg_parse function for std_graph_args_t.
static int graph_parse_opt(int key, char * arg, struct argp_state * state){
  std_graph_args_t * args = (std_graph_args_t *)state->input;
  char * ptr;
  int i, ptrshf;
  switch(key){
  case 'd': args->directed = 1; break;
  case 'e': args->edge_prob = atof(arg); break;
  case 'f': args->readfile=1; strcpy(args->filename, arg); break;
  case 'F': args->model = FLAT; break;
  case 'G': args->model = GEOMETRIC; break;
  case 'K':
    args->model = KRONECKER;
    strcpy(args->kron_string,arg);
    ptr = arg;
    //T0_fprintf(stderr, "%s\n", arg);
    sscanf(ptr, "%d:", &args->kron_mode);
    assert( args->kron_mode == 0 || args->kron_mode == 1 || args->kron_mode == 2);
    ptr+=2;  // move pass the ':'
    for(i=0; i<63; i++) {
      sscanf(ptr, "%d%n", &args->kron_spec[i], &ptrshf);
      ptr += ptrshf;  // move pass the digits
      if(*ptr != 'x') // better be 'x' or NULL
        break;
      ptr++;          // move pass the digits
    }
    args->kron_num = i+1;
    break;
  case 'l': args->loops = 1; break;
  case 'n': args->l_numrows = atol(arg); break;
  case 'N': args->numrows = atol(arg); break;
  case 'w': args->weighted = 1; break;
  case 'z': args->nz_per_row = atof(arg); break;
  case ARGP_KEY_INIT:
    args->edge_prob = 0.0;
    args->readfile = 0;
    args->model = FLAT;
    args->numrows = 0;
    //args->l_numrows = 10000;
    args->nz_per_row = 10.0;
    //args->directed = 0;
    //args->weighted = 0;
    //args->loops = 0;
    break;
  case ARGP_KEY_END:
    if(args->numrows == 0) args->numrows = args->l_numrows*THREADS;
    if(args->directed)
      resolve_edge_prob_and_nz_per_row(&args->edge_prob, &args->nz_per_row,
                                       args->numrows,
                                       (args->weighted ? DIRECTED_WEIGHTED : DIRECTED),
                                       (args->loops ? LOOPS : NOLOOPS));
    else{
      resolve_edge_prob_and_nz_per_row(&args->edge_prob, &args->nz_per_row,
                                       args->numrows,
                                       (args->weighted ? UNDIRECTED_WEIGHTED : UNDIRECTED),
                                       (args->loops ? LOOPS : NOLOOPS));
    }
    break;
  }
  return(0);
}

static struct argp_option graph_options[] =
  {
    {0, 0, 0, 0, "Input (as file):", 5},
    {"readfile",   'f', "FILE",  0, "Read input from a file"},
    {0, 0, 0, 0, "Input (as random graph):", 6},
    {"l_numrows",  'n', "NUM",   0, "Number of rows per PE in the matrix."},
    {"numrows",    'N', "NUM",   0, "Number of rows in the matrix. This option takes precedence if both -n and -N options are used."},
    {"directed",   'd', 0,       0, "Specify a directed graph"},
    {"edge_prob",  'e', "EDGEP", 0, "Probability that an edge appears. Use this or -z option to control the density of the graph."},
    {"flat",       'F', 0,       0, "Specify flat random graph model"},
    {"geometric",  'G', 0,       0, "Specify geometric random graph model"},
    {"kronecker",  'K', "KSTR",  0, "Specify a Kronecker product graph.\nKSTR must be a string of the form MODE:S1xS2x...Sk where MODE is 0, 1, or 2 and the Si are small integers that specify the stars whose product is the kronecker product graph. For instance -K 0:3x4x5 specifies MODE 0 and takes the product of K_{1,3}, K_{1,4}, and K_{1,5}. MODE 0 : No triangles. MODE 1: Many triangles. MODE 2: Few triangles. "},
    {"loops",      'l', 0,       0, "Specify you want to force loops into graph"},
    {"weighted",   'w', 0,       0, "Specify you want the edges to be weighted"},
    {"nz_per_row", 'z', "NZPR",  0, "Average number of nonzeros per row. Specify this or -e option to control the density of the graph. Default = 10.0"},
    {0}
  };

struct argp std_graph_options_argp = {
  graph_options, graph_parse_opt, 0, 0, 0
};


/* this function looks at the std_args_t and std_graph_args_t structs to decide
   whether to read a matrix or generate a random matrix. The function returns
   the read or generated matrix. It also writes some matrix statistics to stderr,
   or the output json file.
*/
sparsemat_t * get_input_graph(std_args_t * sargs, std_graph_args_t * gargs){

  sparsemat_t * mat;
  
  if(!gargs->readfile){
    
    if(gargs->model == KRONECKER){

      // generate a random kronecker graph from a string (mode:#x#x...#)
      mat = generate_kronecker_graph_from_spec(gargs->kron_mode,
                                               gargs->kron_spec,
                                               gargs->kron_num);
      
    }else{
      
      // Generate a random FLAT or GEOMETRIC graph
      //int64_t numrows = gargs->numrows;
      int64_t seed = MYTHREAD + sargs->seed;
      edge_type et;
      if(gargs->directed){
        et = (gargs->weighted ? DIRECTED_WEIGHTED : DIRECTED);
      }else{
        et = (gargs->weighted ? UNDIRECTED_WEIGHTED: UNDIRECTED);
      }
      self_loops loops = (gargs->loops ? LOOPS : NOLOOPS);
      
      mat = random_graph(gargs->numrows, gargs->model, et, loops, gargs->edge_prob, seed + 2);
    }
  }else{
    
    // Read a matrix from a file
    mat = read_matrix_mm_to_dist(gargs->filename);
    
  }
  if(!mat){T0_printf("ERROR: get_input_mat: mat is NULL!\n"); lgp_global_exit(1);}

  if(sargs->json == 0){
    T0_fprintf(stderr,"Input matrix:\n");
    T0_fprintf(stderr,"----------------------------------------------------\n");
    T0_fprintf(stderr,"\t%"PRId64" rows\n\t%"PRId64" columns\n\t%"PRId64" nonzeros\n\n",
               mat->numcols, mat->numrows, mat->nnz);
  }else{
    if(MYTHREAD == 0){
      FILE * jp = fopen(sargs->json_output, "a");    
      T0_fprintf(jp, "\"matrix_numrows\": \"%"PRId64"\",\n", mat->numrows);
      T0_fprintf(jp, "\"matrix_numcols\": \"%"PRId64"\",\n", mat->numcols);
      T0_fprintf(jp, "\"matrix_nnz\": \"%"PRId64"\",\n", mat->nnz);
      fclose(jp);
    }
  }

  return(mat);
}


// Writes some parameters for the input matrix before it is generated to stderr or an output json file.
void write_std_graph_options(std_args_t * sargs, std_graph_args_t * gargs){
  if(!gargs->readfile){
    char model[32];
    if(gargs->model == FLAT)
      sprintf(model, "FLAT        (-F)");
    else if(gargs->model == GEOMETRIC)
      sprintf(model, "GEOMETRIC   (-G)");
    else if(gargs->model == KRONECKER)
      sprintf(model, "KRONECKER   (-K)");
    
    if(sargs->json && !MYTHREAD){
      FILE * jp = fopen(sargs->json_output, "a");
      fprintf(jp,"\"graph_model\": \"%.*s\",\n", 9, model);
      fclose(jp);
    }else{
      T0_fprintf(stderr,"Input Graph/Matrix parameters:\n");
      T0_fprintf(stderr,"----------------------------------------------------\n");
      T0_fprintf(stderr, "Graph model: %s.\n",model);
      T0_fprintf(stderr,"%s, %s, %s\n",
                 (gargs->directed ? "Directed": "Undirected"),
                 (gargs->weighted ? "Weighted": "Unweighted"),
                 (gargs->loops ? "Loops": "No Loops"));
      
      
      T0_fprintf(stderr,"Number of rows           (-N): %"PRId64"\n", gargs->numrows);
      T0_fprintf(stderr,"Avg # nnz per row        (-z): %2.2lf\n", gargs->nz_per_row);
      T0_fprintf(stderr,"Edge probability         (-e): %lf\n\n", gargs->edge_prob);
    }
  }else{
    if(sargs->json && !MYTHREAD){
      FILE * jp = fopen(sargs->json_output, "a");
      fprintf(jp,"\"input_file\": \"%s\",\n", gargs->filename);
      fclose(jp);
    }else{
      T0_fprintf(stderr,"Reading input from %s\n\n", gargs->filename);
    }
  }
}


/* all threads look at argv for the strings "--help", "--usage", or "-?".
   If those strings are found, we know, main needs to exit after parsing the command
   line.
*/
int check_for_exit(int argc, char * argv[], int ret){
  int i;
  for(i = 0; i < argc; i++){
    if( (strcmp(argv[i], "--help") == 0) || 
        (strcmp(argv[i], "-?") == 0) ||
        (strcmp(argv[i], "--usage") == 0)){
      lgp_finalize();
      return 1;
    }
  }

  // check if PE 0 hit any reason to exit when it parsed the command line with argp.
  ret = (int)lgp_reduce_add_l((long)ret);
  if(ret){
    lgp_finalize();
    return(-1);
  }
  return(0);
}


void share_args(void * args, size_t n){
  SHARED char * temp = lgp_all_alloc(THREADS, n);
  if(!MYTHREAD)
    lgp_memput(temp, (void*)args, n, 0);
  lgp_barrier();
  lgp_memget((void*)args, temp, n, 0);
  lgp_barrier();
  lgp_all_free(temp);
}


int distribute_cmd_line(int argc, char ** argv, void * args, size_t args_len, int ret){
  
  ret = check_for_exit(argc, argv, ret);
  if(ret){
    lgp_finalize();
    return(ret);
  }

  share_args(args, args_len);
  return(0);
}



/*! \brief This function initializes the parallel environment, parses the command line, and prints
 * some basic app data to stderr or json.
 * 
 * \param argc Standard argc
 * \param argv Standard argv
 * \param args The arguments struct from the app. This struct will be different in general for each app.
 * \param arg_len The sizeof(the arguments struct)
 * \param argp The argp struct which must be initialized in main.
 * \param sargs A pointer to an std_args_t, which must be a field in args.
 * \return 0 if success, < 0 if error, > 0 if main should exit after this call, but main should return(0)
 */
int bale_app_init(int argc, char ** argv,
                  void * args, int arg_len,
                  struct argp * argp,
                  std_args_t * sargs){

  // call shmem_init / set up atomic domain variable for UPC
  lgp_init(argc, argv);

  // parse the command line with PE 0
  int ret = 0;  
  if(MYTHREAD == 0){
    ret = argp_parse(argp, argc, argv, ARGP_NO_EXIT, 0, args);
  }

  // PE 0 broadcasts command line and everyone decides if they should return to main or continue
  ret = distribute_cmd_line(argc, argv, args, arg_len, ret);
  if(ret) return(ret);
  
  // print header for every bale app
  time_t now = time(NULL);
  struct tm *date = localtime(&now);
  if(sargs->json){
  /* open the json file */ 
    if(MYTHREAD == 0){
      FILE * fp = fopen(sargs->json_output, "w");
      fprintf(fp, "{\n\"bale_version\": \"%4.2f\",\n", BALE_VERSION);
      fprintf(fp, "\"date\": \"%04d-%02d-%02d.%02d:%02d\",\n",
                date->tm_year+1990, date->tm_mon, date->tm_mday,
                date->tm_hour, date->tm_min);
      fprintf(fp, "\"app\": \"%s\",\n", argv[0]);
      fprintf(fp, "\"pes\": \"%d\",\n", THREADS);
      fclose(fp);
    }
  }else{
    T0_fprintf(stderr,"\n***************************************************************\n");
#if __UPC__
    T0_fprintf(stderr,"Bale Version %4.2f (UPC %ld): %04d-%02d-%02d.%02d:%02d\n",
               BALE_VERSION,
               __UPC_VERSION__,
               date->tm_year+1990, date->tm_mon, date->tm_mday,
               date->tm_hour, date->tm_min);
#elif USE_SHMEM
    T0_fprintf(stderr,"Bale Version %4.2f (OpenShmem version %d.%d): %04d-%02d-%02d.%02d:%02d\n",
               BALE_VERSION,
               SHMEM_MAJOR_VERSION, SHMEM_MINOR_VERSION,
               date->tm_year+1990, date->tm_mon+1, date->tm_mday,
               date->tm_hour, date->tm_min); 
#endif
    
    int i;
    
    T0_fprintf(stderr,"Running command on %d PEs:", THREADS);
    for(i=0; i<argc;i++){
      T0_fprintf(stderr," %s", argv[i]);
    }
    T0_fprintf(stderr,"\n");
    T0_fprintf(stderr,"***************************************************************\n\n");
  }
  return(0);
}


/*! \brief Finish off the json output file if appropriate. Call lgp_finalize()
 */
void bale_app_finish(std_args_t * sargs){
  if(sargs->json && !MYTHREAD){    
    FILE * jp = fopen(sargs->json_output, "r+");
    fseek(jp,-2,SEEK_END);
    fprintf(jp,"\n}");
    fclose(jp);
  }
  lgp_finalize();
}

/*! \brief Write the time of an app implmentation to stderr, or the json file.
 *
 */
void bale_app_write_time(std_args_t * sargs, char * model_str, double time){
  if(sargs->json && !MYTHREAD){    
    FILE * jp = fopen(sargs->json_output, "a");
    fprintf(jp,"\"%s\": \"%lf\",\n", model_str, time);
    fclose(jp);
  }else{
    T0_fprintf(stderr, "%14s: %8.3lf\n", model_str, time);
  }
}

void bale_app_write_int(std_args_t * sargs, char * key, int64_t val){
  if(sargs->json && !MYTHREAD){    
    FILE * jp = fopen(sargs->json_output, "a");
    fprintf(jp,"\"%s\": \"%"PRId64"\",\n", key, val);
    fclose(jp);
  }else{
    T0_fprintf(stderr, "%14s: %"PRId64"\n", key, val);
  }
}

void bale_app_write_double(std_args_t * sargs, char * key, double val){
  if(sargs->json && !MYTHREAD){    
    FILE * jp = fopen(sargs->json_output, "a");
    fprintf(jp,"\"%s\": \"%lf\",\n", key, val);
    fclose(jp);
  }else{
    T0_fprintf(stderr, "%14s: %lf\n", key, val);
  }
}
