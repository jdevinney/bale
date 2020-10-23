/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*!  \file std_options.h
\brief Header file for std_options library. 
The std_options library is a support library for parsing command line options.
*/

#ifndef STD_OPTIONS_H
#define STD_OPTIONS_H

#include <stdlib.h>
#include <stdint.h>
#include "spmat_utils.h"
#include <argp.h>
#include "default_app_opts.h"

/*! \brief the version (equals the version of bale_classic) */
#define C_BALE_VERSION 3.0

/*!
\struct std_args_t
\brief Stuct for the arguments used by almost all apps
*/
typedef struct std_args_t{
  int dump_files;              /*!< flag to dump files or not */
  int json;                    /*!< flag to output json performance files or not */
  char json_output[128];       /*!< string buffer json output lines */
  int models_mask;             /*!< the OR of different implementations in each app. Gives command line control, but is different for each app. */
  int64_t seed;                /*!< random number generator seed */
}std_args_t;

extern struct argp std_options_argp;

/*!
\struct std_graph_args_t
\brief Stuct for the arguments that control the graph (matrix) options
*/
typedef struct std_graph_args_t{
  int64_t numrows;            /*!< number of rows in the matrix */
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

// TODO should this be in spmat_utils
sparsemat_t *  get_input_graph(std_args_t * sargs, std_graph_args_t * gargs); //!< parses the args and calls the appropriate generator

void           write_std_options(std_args_t * sargs); //!< displays the args before the run 
void           write_std_graph_options(std_args_t * sargs, std_graph_args_t * gargs); //!< displays the args before the run

int  bale_app_init(int argc, char ** argv, void * args, int arg_len, struct argp * argp, std_args_t * sargs); //!< init service structs and print args
void bale_app_finish(std_args_t * sargs); //<! clean up after run 
void bale_app_write_int(std_args_t * sargs, char * key, int64_t val); //!< write out a key, val pair
void bale_app_write_time(std_args_t * sargs, char * model_str, double time); //!< write out a simple wall clock timer 

#endif
