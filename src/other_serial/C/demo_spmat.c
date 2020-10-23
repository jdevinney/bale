/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file demo_spmat.c
 * \brief Program that demonstrates (checks) some of the rountines in
 * spmat_util.c
 */

#include "spmat_utils.h"

int main(int argc, char * argv[]) 
{
  int i;

  int64_t *p = rand_perm(10, 0);
  printf("Identity perm = ");
  for(i=0; i<10; i++)
    printf(" %"PRId64,p[i]);
  printf("\n");
  free(p);

  int64_t *q = rand_perm(12, 1);
  printf("Random perm = ");
  for(i=0; i<12; i++)
    printf(" %"PRId64,q[i]);
  printf("\n");

  printf("Is q a perm?  %s\n", (is_perm(q,12))?"yes":"no");
  q[0] = 12;
  printf("Is q a perm?  %s\n", (is_perm(q,12))?"yes":"no");
  free(q);



  sparsemat_t *graph;
  graph = read_matrix_mm("../../../example_matrices/undirected_flat_100.mm");
  spmat_stats(graph);

  dump_matrix(graph, 4, "dump_4.out");
  dump_matrix(graph, 0, "dump_0.out");

  write_matrix_mm(graph, "tidy_demo.mm.out");
  clear_matrix(graph);

#if 0
  printf("Generate a Kronecker Product of Stars\n");
  kron_args_t * kron_args = kron_args_init("2: 2 2");
  printf("-- input %s\n", kron_args->str);
  printf("-- mode %"PRId32"\n", kron_args->mode);
  printf("-- num_stars %"PRId32"\n", kron_args->num_stars);
  for(i=0; i<kron_args->num_stars; i++)
    printf("%"PRId32" ", kron_args->star_size[i]); 
  printf("\n-- numrows %"PRId64"\n", kron_args->numrows);
  printf("Known number of triangles = %"PRId64"\n", calc_num_tri_kron_graph(kron_args));

  sparsemat_t *Kron = kronecker_product_graph(kron_args);
  spmat_stats(Kron);
  dump_matrix(Kron, 0, "kron.out");

  clear_matrix(Kron);
  free(kron_args);
#endif
  return(0);

}
