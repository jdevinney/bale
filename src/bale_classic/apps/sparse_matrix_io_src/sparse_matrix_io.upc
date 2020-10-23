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

/*! \file sparse_matrix_io.upc
 * \brief Demo program that runs the variants of sparse_matrix_io kernel. This program
 * writes a sparsemat_t to disk in a binary format.
 */

/*! 
\page sparse_matrix_io_page Sparse Matrix I/O

Demo program that runs the variants of sparse_matrix_io kernel. It first generates 
a random matrix and then it writes this matrix to disk
in a directory called 'write_sparse_dir'. Finally, it reads the sparse matrix
dataset back in and compares the resulting matrix with the original matrix.

We define a sparse matrix dataset to be the following:
 - It lives in a directory of its own
 - It contains one ASCII file called 'metadata' which contains 5 lines of the form key=val. The keys are numrows, numcols,
   nnz, nwriters, and values. nwriters stands for the number of PEs that were involved in writing the dataset. The values
   value is 0/1 if the matrix doesn't have / or has values.
 - It contains nwriters binary 'rowcnt' files. The ith record in the jth file is the the offset in the jth nonzero file for
 where the ith row's data begins.
 - It contains nwriters binary 'nonzero' files. These files contain the nonzeros in each row and are ordered by row.
 - If the matrix has values, it contains nwriters binary 'values' files. These files contain the values for the nonzeros.

See files spmat_agp.upc, spmat_exstack.upc and spmat_io.upc
for the source for the kernels.

* Run with the --help, -?, or --usage flags for run details.
 */

/********************************  argp setup  ************************************/
typedef struct args_t{
  int64_t num_readers;
  char dir_path[128];
  std_args_t std;
  std_graph_args_t gstd;
}args_t;

static int parse_opt(int key, char * arg, struct argp_state * state){
  args_t * args = (args_t *)state->input;
  switch(key)
    {
    case 'r': args->num_readers = atoi(arg); break;
    case 'p': strcpy(args->dir_path,arg); break;
    case ARGP_KEY_INIT:
      state->child_inputs[0] = &args->std;
      state->child_inputs[1] = &args->gstd;
      break;
    }
  return(0);
}

static struct argp_option options[] =
  {
    {"num_readers", 'r', "NUM", 0, "Specify the number of PEs to read the written matrix back in."},
    {"dir_path", 'p', "PATH", 0, "Specify an alternate name for the dataset directory. Default is 'write_matrix_dir"},
    {0}
  };


static struct argp_child children_parsers[] =
  {
    {&std_options_argp, 0, "Standard Options", -2},
    {&std_graph_options_argp, 0, "Standard Graph Options", -3},
    {0}
  };

int main(int argc, char * argv[])
{

  /* process command line */
  args_t args = {0};
  struct argp argp = {options, parse_opt, 0,
                      "Parallel write and read sparse matrix.", children_parsers};

  args.gstd.l_numrows = 1000000;
  args.num_readers = -1;
  strcpy(args.dir_path, "write_matrix_dir");

  int ret = bale_app_init(argc, argv, &args, sizeof(args_t), &argp, &args.std);
  if(ret < 0) return(ret);
  else if(ret) return(0);
  
  if(args.num_readers == -1)
    args.num_readers = THREADS;
  
  if(!MYTHREAD){
    write_std_graph_options(&args.std, &args.gstd);
    write_std_options(&args.std);
  }

  // read in a matrix or generate a random graph
  sparsemat_t * inmat = get_input_graph(&args.std, &args.gstd);
  if(!inmat){T0_fprintf(stderr, "ERROR: sparse_matrix_io: inmat is NULL!\n");return(-1);}

  if(args.std.dump_files) write_matrix_mm(inmat, "inmat");
  //print_matrix(inmat);
  double t1;
  minavgmaxD_t stat[1];
  char model_str[32];
  int64_t use_model;
  int error = 0;
  sparsemat_t * readmat;
  
  for( use_model=1L; use_model < 32; use_model *=2 ) {

    switch( use_model & args.std.models_mask ) {
    case AGP_Model:
      sprintf(model_str, "AGP_write");
      t1 = wall_seconds();
      write_sparse_matrix_agp(args.dir_path, inmat);
      t1 = wall_seconds() - t1;
      lgp_min_avg_max_d( stat, t1, THREADS );
      bale_app_write_time(&args.std, model_str, stat->avg);

      sprintf(model_str, "AGP_read");
      t1 = wall_seconds();
      readmat = read_sparse_matrix_agp(args.dir_path, args.num_readers);
      t1 = wall_seconds() - t1;
      lgp_min_avg_max_d( stat, t1, THREADS );
      bale_app_write_time(&args.std, model_str, stat->avg);
      
      break;
    case EXSTACK_Model:
      
      sprintf(model_str, "Exstack_write");
      t1 = wall_seconds();
      write_sparse_matrix_exstack(args.dir_path, inmat, args.std.buf_cnt);
      t1 = wall_seconds() - t1;
      lgp_min_avg_max_d( stat, t1, THREADS );
      bale_app_write_time(&args.std, model_str, stat->avg);

      sprintf(model_str, "Exstack_read");
      t1 = wall_seconds();
      readmat = read_sparse_matrix_exstack(args.dir_path, args.num_readers, args.std.buf_cnt);
      t1 = wall_seconds() - t1;
      lgp_min_avg_max_d( stat, t1, THREADS );
      bale_app_write_time(&args.std, model_str, stat->avg);
      
      break;
    case EXSTACK2_Model:
      continue;
      //sprintf(model_str, "Exstack2");
      //break;
    case CONVEYOR_Model:
      continue;
      //sprintf(model_str, "Conveyor");
      //break;
    case ALTERNATE_Model:
      T0_fprintf(stderr,"There is no alternate model here!\n"); continue;
      continue;
    case 0:
      continue;
    }

    if(readmat == NULL){
      T0_fprintf(stderr,"ERROR: sparse_matrix_io: read failed!\n");
      error = 1;
    }else{
      if(compare_matrix(readmat, inmat)){
        T0_fprintf(stderr,"ERROR: sparse_matrix_io: read version of written matrix does not match!\n");
        error = 1; 
      }
      clear_matrix(readmat); free(readmat);
    }
    

  }
  
  clear_matrix(inmat); free(inmat);
  lgp_barrier();
  bale_app_finish(&args.std);
  return(error);
}
