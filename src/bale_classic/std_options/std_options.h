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

/*! \file std_options.h
 * \brief Header file for std_options library. 
 * The std_options library is a support library for all bale apps, mostly for command line parsing.
 */
#ifndef STD_OPTIONS_H
#define STD_OPTIONS_H

#include <stdlib.h>
#include <stdint.h>
#include <spmat.h>
#include <argp.h>


/*! \struct std_args_t
 * \brief Stuct for the arguments used by almost all apps.
 */
typedef struct std_args_t{
  int64_t buf_cnt;         /*!< the number of packages in each exstack or exstack2 buffer */
  int cores_per_node;          /*!< core per node */
  int dump_files;              /*!< flag to dump files or not */
  int json;                    /*!< flag to output json performance files or not */
  char json_output[128];       /*!< string buffer json output lines */
  int models_mask;             /*!< the OR of the implementations to run (AGP=1, EXSTACK=2, EXSTACK2=4, CONVEY=8) */
  int64_t seed;                /*!< random number generator seed */
}std_args_t;

extern struct argp std_options_argp;

/*! \struct std_graph_args_t
 * \brief Stuct for the arguments that control the graph (matrix) options
 */
typedef struct std_graph_args_t{ 
  int64_t l_numrows;          /*!< local number of rows in the matrix */
  int64_t numrows;            /*!< global number of rows in the matrix */
  int readfile;               /*!< whether or not to read from a file */
  char filename[128];         /*!< string for the filename */
  graph_model model;          /*!< graph model (currently: FLAT (Erdos-Renyi), GEOMETRIC, KRONECKER) */
  double edge_prob;           /*!< edge_prob for FLAT AND GEOMETRIC */
  double nz_per_row;          /*!< average number of non-zeros in a row (computed from edge_prob or vice-versa) */
  int directed;               /*!< whether or not the graph is directed (the matrix has non-zeros above the diagonal)*/
  int loops;                  /*!< whether or not the matrix as all ones on the diagonal */
  int weighted;               /*!< whether or not the edges have weights (doubles in internal (0,1] */
  char kron_string[128];      /*!< string to given the stars in the Kronecker product construction (format M:S1xS2x..Sn where M=mode {0,1,2} and Si's are the sizes of the stars) */
  int kron_spec[64];          /*!< array to hold the star sizes */
  int kron_num;               /*!< number of stars */
  int kron_mode;              /*!< Kronecker mode {0,1,2} resulting in no, lots and few triangles */
} std_graph_args_t;

extern struct argp std_graph_options_argp;

/* defines to support the different models of global and buffered references */
#define AGP_Model        1L /*!< the Atomic or Generic Interface (straight UPC/SHMEM) */
#define EXSTACK_Model    2L /*!< the exstack bulk synchronous buffering model */
#define EXSTACK2_Model   4L /*!< the exstack2 asynchronous buffering model */
#define CONVEYOR_Model   8L /*!< the conveyor buffering model */
#define ALTERNATE_Model 16L /*!< an alternate model (meant for user supplied code) */
#define ALL_Models      15L /*!< default for running all models */


/*! \ingroup service_functions */
sparsemat_t *  get_input_graph(std_args_t * sargs, std_graph_args_t * gargs); /*!< parses the args and calls the appropriate generator */
/*! \ingroup service_functions */
void           write_std_options(std_args_t * sargs); /*!< displays the args before the run */
/*! \ingroup service_functions */
void           write_std_graph_options(std_args_t * sargs, std_graph_args_t * gargs); /*!< displays the args before the run */

/*! \ingroup service_functions */
int  bale_app_init(int argc, char ** argv, void * args, int arg_len, struct argp * argp, std_args_t * sargs); /*!< init service structs and print args */
/*! \ingroup service_functions */
void bale_app_finish(std_args_t * sargs); /*<! clean up after run */
/*! \ingroup service_functions */
/*! write out a key, val pair */
void bale_app_write_int(std_args_t * sargs, char * key, int64_t val);
/*! write out a key, val pair */
void bale_app_write_double(std_args_t * sargs, char * key, double val);
/*! \ingroup service_functions */
void bale_app_write_time(std_args_t * sargs, char * model_str, double time); /*!< write out a simple wall clock timer */

#endif


