/*****************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//
//  All rights reserved.
//  
//  This file is a part of Bale.  For license information see the
//  LICENSE file in the top level directory of the distribution.
//  
 *****************************************************************/
/*! \file spmat_utils.upc
 * \brief Utilities to support spmat and wrappers to make it easier to switch 
 * between routines written in the various models.
 */
#include <math.h>
#include <spmat.h>
#include <exstack.h>

/*!
\brief Produce a global array the holds a uniform random permutation.
\param N the global length of the permutaion
\param seed the seed for the random number generator
\return a pointer to the matrix that has been produced or NULL if the model can't be used

This is wrapper for implementations written in the different models.
It is interest enough to be its own apps, one should experiment with it
within the apps framework. 

This is a collective call.

this implements the random dart algorithm to generate the permutation.
Each thread throws its elements of the perm array randomly at large target array.
Each element claims a unique entry in the large array using compare_and_swap.
This gives a random permutation with spaces in it, then you squeeze out the spaces.
\ingroup spmatgrp
*/
SHARED int64_t * rand_permp(int64_t N, int seed) {
  SHARED int64_t * p;
  //p = rand_permp_agp(N, seed);
  p = rand_permp_exstack(N, seed, 128);
  //p = rand_permp_exstack2(N, seed, 1024);
  
  if(!is_perm(p, N)){
    T0_printf("ERROR: rand_permp: not a permutation!\n");fflush(0);
    return(NULL);
  }
  return(p);
}


/*!
\brief apply row and column permutations to a sparse matrix using straight UPC
\param omat pointer to the original matrix
\param rperminv pointer to the global array holding the inverse of the row permutation
       rperminv[i] = j means that row i of A goes to row j in matrix Ap
\param cperminv pointer to the global array holding the inverse of the column permutation
       cperminv[i] = j means that col i of A goes to col j in matrix Ap

\return a pointer to the matrix that has be computed or NULL on failure

This is wrapper for implementations written in the different models.
It is interest enough to be its own apps, one should experiment with it
within the apps framework. 

\ingroup spmatgrp
 */
sparsemat_t * permute_matrix(sparsemat_t *omat, SHARED int64_t *rperminv, SHARED int64_t *cperminv) {
  //return( permute_matrix_agp(omat, rperminv, cperminv) );
    return( permute_matrix_exstack(omat, rperminv, cperminv, 512) );
  //return( permute_matrix_exstack2(omat, rperminv, cperminv, 512) );
  //return( permute_matrix_conveyor(omat, rperminv, cperminv) );
}

/*! \brief produce the transpose of a sparse matrix 
 * \param omat  pointer to the original matrix
 * \return a pointer to the matrix that has be computed or NULL on failure
 *
 * This is wrapper for implementations written in the different models.
 * It is interest enough to be its own apps, one should experiment with it within the apps framework. 
 * This is provided as library function. 
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix(sparsemat_t *omat) {
  sparsemat_t * A;
  //A = transpose_matrix_agp(omat);
  A = transpose_matrix_exstack(omat, 512);
  //A = transpose_matrix_exstack2(omat, 512);
  //A = transpose_matrix_conveyor(omat);
  if(!A){return(NULL);}
  
  sort_nonzeros(A);
  return(A);
}


/*! \brief checks that the sparse matrix is lower triangluar
 * \param A pointer to the sparse matrix
 * \param unit_diagonal set to 1 to make sure all pivots are nonzero
 * \return 0 on success, non-0 on error.
 * kind of a speciality routine to check that toposort might of worked
 * \ingroup spmatgrp
 */
int is_lower_triangular(sparsemat_t *A, int64_t unit_diagonal) {
  int64_t i,j, row, * ltri_rperm, * ltri_cperm;
  int64_t err = 0, err2 = 0;

  lgp_barrier();

  /* we are hoping for an lower triangular matrix here */
  for(i=0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j = A->loffset[i]; j < A->loffset[i+1]; j++){
      if( A->lnonzero[j] > global_row ) {
        err++;
        //if(global_row < 10)
        //printf("row %ld : col %ld (j = %ld)\n", global_row, A->lnonzero[j], j);
      }else if( A->lnonzero[j] == global_row ){
        pivot = 1;
      }
    }
    if(!pivot)
      err2++;
  }

  err = lgp_reduce_add_l(err);
  err2 = (unit_diagonal ? lgp_reduce_add_l(err2) : 0);
  if( err || err2 ){
    if(!MYTHREAD)printf("\nThere are %"PRId64" nz above diag. and %"PRId64" missing pivots in lower.\n", err, err2);
    fflush(0);
  }

  lgp_barrier();

  return(!(err || err2));
}

/*! \brief checks that the sparse matrix is upper triangluar
 * \param A pointer to the sparse matrix
 * \param unit_diagonal set to 1 to make sure all pivots are nonzero
 * \return 0 on success, non-0 on error.
 *
 * \ingroup spmatgrp
 */
int is_upper_triangular(sparsemat_t *A, int64_t unit_diagonal) {
  int64_t i,j, row, * ltri_rperm, * ltri_cperm;
  int64_t err = 0, err2 = 0;
  lgp_barrier();

  /* we are hoping for an upper triangular matrix here */
  for(i=0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j=A->loffset[i]; j < A->loffset[i+1];j++){
      if( A->lnonzero[j] < global_row ) {
        err++;
      }else if( A->lnonzero[j] == global_row){
        pivot = 1;
      }
    }
    if(!pivot)
      err2++;
  }

  err = lgp_reduce_add_l(err);
  err2 = (unit_diagonal ? lgp_reduce_add_l(err2) : 0);
  if( err || err2 ){
    //if(!MYTHREAD)printf("\nThere are %"PRId64" nz below diag. and %"PRId64" missing pivots in upper.\n", err, err2);
    fflush(0);
  }

  lgp_barrier();

  return(!(err || err2));
}

/*! \brief Sets all entries above the kth diagonal to zero.
 * \param A A pointer to a sparse matrix
 * \param k Anything above the kth diagonal will be zero (k > 0 is above the main diagonal k < 0 is below)
 * \return The number of nonzeros that were removed.
 * \ingroup spmatgrp
 */
int64_t tril(sparsemat_t * A, int64_t k) {
  // remove entries below the diagonal
  int64_t i, j, col, pos = 0, start = 0;
  int weighted = (A->value != NULL);
  for(i = 0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j = start; j < A->loffset[i+1]; j++){
      col = A->lnonzero[j];
      if(col - global_row <= k){
        A->lnonzero[pos] = col;
        if(weighted) A->lvalue[pos] = A->lvalue[j];
        pos++;
      }
    }
    start = A->loffset[i+1];
    A->loffset[i+1] = pos;
  }
  int64_t new_nnz = lgp_reduce_add_l(pos);
  int64_t ret = A->nnz - new_nnz;
  A->nnz = new_nnz;
  A->lnnz = pos;
  return(ret);
}

/*! \brief Sets all entries below the kth diagonal to zero.
 * \param A A pointer to a sparse matrix
 * \param k Anything below the kth diagonal will be zero (k > 0 is above the main diagonal k < 0 is below)
 * \return The number of nonzeros that were removed.
 * \ingroup spmatgrp
 */
int64_t triu(sparsemat_t * A, int64_t k) {
  // remove entries below the diagonal
  int64_t i, j, col, pos = 0, start = 0;
  int weighted = (A->value != NULL);
  for(i = 0; i < A->lnumrows; i++){
    int64_t global_row = i*THREADS + MYTHREAD;
    int pivot = 0;
    for(j = start; j < A->loffset[i+1]; j++){
      col = A->lnonzero[j];
      if(col - global_row >= k){
        A->lnonzero[pos] = col;
        if(weighted) A->lvalue[pos] = A->lvalue[j];
        pos++;
      }
    }
    start = A->loffset[i+1];
    A->loffset[i+1] = pos;
  }
  int64_t new_nnz = lgp_reduce_add_l(pos);
  int64_t ret = A->nnz - new_nnz;
  A->nnz = new_nnz;
  A->lnnz = pos;
  return(ret);
}


/*! \brief checks that a global array is in fact a permutation
 * \param perm SHARED pointer to the global array
 * \param N the length of the global array
 * \return 1 if it is permutation
 * \ingroup spmatgrp
 */
int is_perm(SHARED int64_t * perm, int64_t N) {
  int64_t i;
  int64_t * lperm = lgp_local_part(int64_t, perm);
  SHARED int64_t * flag = lgp_all_alloc(N, sizeof(int64_t));
  if( flag == NULL ) return(0);
  int64_t * lflag = lgp_local_part(int64_t, flag);
  int64_t l_N = (N + THREADS - MYTHREAD - 1)/THREADS;

  for(i = 0; i < l_N; i++) 
    lflag[i] = 0;
  lgp_barrier();
  for(i = 0; i < l_N; i++) 
   lgp_put_int64(flag, lperm[i], 1);
  lgp_barrier();
  int64_t err = 0L;
  for(i = 0; i < l_N; i++) 
    if(lflag[i] == 0) 
      err++;
  lgp_all_free(flag);
  err = lgp_reduce_add_l(err);
  return(err == 0L);
}


 /*************************************************************************************/
 /*                               RANDOM MATRICES                                     */
 /*************************************************************************************/

/*! 
\brief A routine to generate the adjacency matrix of a random graph.
\param n The number of vertices in the graph.
\param model FLAT: Erdos-Renyi random, GEOMETRIC: geometric random graph
\param edgetype See enum edge_type enum. Directed, or not, Weighted or not.
\param loops see self_loops enum. Does every node have a loop or not.
\param edge_density: d in [0, 1), target fraction of edges present.
\param seed: RNG seed.

If the graph is undirected, this routine only returns a lower-triangular
adjancency matrix (since the adjancency matrix would be symmetric and we don't need
the redundant entries).
*/

sparsemat_t * random_graph(int64_t n, graph_model model, edge_type edgetype,
                           self_loops loops,
                           double edge_density, int64_t seed){

  if(model == FLAT){

    return(erdos_renyi_random_graph(n, edge_density, edgetype, loops, seed));
    
  }else if(model == GEOMETRIC){
    edge_type et = edgetype;
    if(et == DIRECTED)
      et = UNDIRECTED;
    else if(et == DIRECTED_WEIGHTED)
      et = UNDIRECTED_WEIGHTED;

    SHARED point_t * op;
    
    double r;

    // determine the r that will lead to the desired edge density
    // Expected degree = n*pi*r^2
    // The expected number of edges E = n^2*pi*r^2/2
    // for undirected density = E/(n choose 2)
    // for directed   density = E/(n^2 - n)
    if(edgetype == UNDIRECTED || edgetype == UNDIRECTED_WEIGHTED)
      r = sqrt((n-1)*edge_density/(M_PI*n));
    else
      r = sqrt(2*(n-1)*edge_density/(M_PI*n));
    sparsemat_t * L = geometric_random_graph(n, r, et, loops, seed, NULL);// &op);
    if(L == NULL){
      T0_fprintf(stderr,"Error: random_graph: geometric_r_g returned NULL!\n");
      return(NULL);
    }
    lgp_barrier();
    // hack to print points
    if(!MYTHREAD && 0){
      point_t pt;
      for(int i = 0; i < n; i++){
        lgp_memget(&pt, op, sizeof(point_t), i);
        //printf("%ld: %lf %lf\n", i, pt.x, pt.y);
      }
    }
    if(edgetype == DIRECTED || edgetype == DIRECTED_WEIGHTED){
      sparsemat_t * A = direct_undirected_graph(L);
      clear_matrix(L);
      //print_matrix(A);
      return(A);
    }
    
    return(L);
    
  }else{
    T0_fprintf(stderr,"ERROR: random_graph: illegal model (%d) for random_graph!\n", model);
    return(NULL);
  }
}

/*! \brief Generates the lower half of the adjacency matrix (non-local) for an Erdos-Renyi random
 * graph. This subroutine uses ALG1 from the paper "Efficient Generation of Large Random Networks" 
 * by Batageli and Brandes appearing in Physical Review 2005. 
 * Instead of flipping a coin for each potential edge this algorithm generates a sequence of 
 * "gaps" between 1s in the upper or lower triangular portion of the adjacency matrix using 
 * a geometric random variable.
 *
 * We parallelized this algorithm by noting that the geometric random variable is memory-less. This means, we 
 * can start the first row on each PE independently. In fact, we could start each row from scratch if we
 * we wanted to! This makes this routine embarrassingly parallel.
 *
 * \param n The total number of vertices in the graph.
 * \param p The probability that each non-loop edge is present.
 * \param edgetype See enum edge_type enum. DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED
 * \param loops See self_loops enum. Does all or no vertices have self loops.
 * \param seed A random seed. This should be a single across all PEs (it will be modified by each PE individually).
 * \return A distributed sparsemat_t that holds the graph.  It is either a lower triangular  matrix
           or square matrix, weighted or not with or without a diagonal.
 */
  sparsemat_t * erdos_renyi_random_graph(int64_t n, double p, edge_type edgetype, self_loops loops, uint64_t seed){
  
  int64_t row, col, i;
  int64_t ln = (n + THREADS - MYTHREAD - 1)/THREADS;
  int64_t lnnz = 0, lnnz_orig=0;
  double D = log(1 - p);
  int64_t r;
  int64_t end = n;
  int64_t ndiag = ln;
  
  lgp_rand_seed(seed);

  /* count lnnz so we can allocate A correctly */   
  row = MYTHREAD;
  col = 1 + floor(log(1 - lgp_rand_double()) / D);
  while(row < n){
    if(edgetype == UNDIRECTED || edgetype == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      // if we just hit a diagonal entry (we don't have to generate this one later)
      if(col == row) ndiag--; 
      lnnz_orig++;
      col += 1 + floor(log(1 - lgp_rand_double()) / D);
    }

    row += THREADS;
    col -= end;
  }
  if(loops == LOOPS) lnnz_orig += ndiag;

  lgp_barrier();

  int weighted = (edgetype == UNDIRECTED_WEIGHTED || edgetype == DIRECTED_WEIGHTED);
  sparsemat_t * A = init_matrix(n, n, lnnz_orig, weighted);
  if(!A){T0_printf("ERROR: erdos_renyi_random_graph: init_matrix failed!\n"); return(NULL);}

  /* reset the seed so we get the same sequence of coin flips */
  lgp_rand_seed(seed);

  /* now go through the same sequence of random events and fill in nonzeros */
  A->loffset[0] = 0;
  row = MYTHREAD;
  col = 1 + floor(log(1 - lgp_rand_double()) / D);
  while(row < n){
    int need_diag = (loops == LOOPS);
    if(edgetype == UNDIRECTED || edgetype == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      if(col == row) need_diag = 0;
      A->lnonzero[lnnz++] = col;
      col += 1 + floor(log(1 - lgp_rand_double()) / D);
    }
    if(need_diag) {
      A->lnonzero[lnnz++] = row;
    }
    row+=THREADS;
    col -= end;
    A->loffset[row/THREADS] = lnnz;
  }
  
  if(lnnz != lnnz_orig){
    printf("ERROR: lnnz (%"PRId64") != lnnz_orig (%"PRId64")\n", lnnz, lnnz_orig);
    return(NULL);
  }

  // fill in the weights
  if(weighted){
    for(i = 0; i < lnnz; i++){
      A->lvalue[i] = lgp_rand_double();
    }
  }

  if(loops == LOOPS && (edgetype == DIRECTED || edgetype == DIRECTED_WEIGHTED))
    sort_nonzeros(A); // to get the diagonal entry sorted correctly
  
  return(A);
}

/*! \brief This is included because it is the definition of an Erdos-Renyi random matrix.
 * This is the naive O(n^2) algorithm. It flips an unfair coin for each possible edge.
 * \param n The total number of vertices in the graph.
 * \param p The probability that each non-loop edge is present.
 * \param edgetype See enum edge_type enum. DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED
 * \param loops See self_loops enum. Do all or no vertices have self loops.
 * \param seed A random seed. This should be a single across all PEs (it will be modified by each PE 
 individually).
 * \return A distributed sparsemat_t that holds the graph.  It is either a lower triangular  matrix
           or square matrix, weighted or not with or without a diagonal.
 */
sparsemat_t * erdos_renyi_random_graph_naive(int64_t n, double p, edge_type edgetype, self_loops loops, uint64_t seed){
  
  int64_t row, col, i, j;
  int64_t ln = (n + THREADS - MYTHREAD - 1)/THREADS;
  int64_t lnnz, lnnz_orig;
  //int64_t P = p*RAND_MAX;
  int64_t end = n;
  
  /* count lnnz so we can allocate A correctly */
  lgp_rand_seed(seed);
  lnnz_orig = 0;
  for(row = MYTHREAD; row < n; row += THREADS){
    if(edgetype == UNDIRECTED || edgetype == UNDIRECTED_WEIGHTED)
      end = row;
    for(col = 0; col < end; col++){
      if(col == row)//skip diagonal (we handle loops separately)
        continue;
      if(lgp_rand_double() < p)
        lnnz_orig++;
    }    
  }
  if(loops == LOOPS) lnnz_orig += n;
  
  lgp_barrier();

  int weighted = (edgetype == UNDIRECTED_WEIGHTED || edgetype == DIRECTED_WEIGHTED);
  sparsemat_t * A = init_matrix(n, n, lnnz_orig, weighted);
  if(!A){T0_printf("ERROR: erdos_renyi_random_graph_naive: init_matrix failed!\n"); return(NULL);}

  /* reset the seed so we get the same sequence of coin flips */
  lgp_rand_seed(seed);

  /* now go through the same sequence of random events and fill in nonzeros */
  A->loffset[0] = 0;  
  lnnz = 0;
  for(row = MYTHREAD; row < n; row += THREADS){
    if(edgetype == UNDIRECTED || edgetype == UNDIRECTED_WEIGHTED)
      end = row + (loops == LOOPS);
    for(col = 0; col < row; col++){
      if(col == row && loops == LOOPS){
        A->nonzero[lnnz++] = row;
        continue;
      }
      if(lgp_rand_double() < p){
        A->lnonzero[lnnz++] = col;
      }      
    }
    A->loffset[row/THREADS + 1] = lnnz;
  }

  if(lnnz != lnnz_orig){
    printf("ERROR: erdos_renyi_random_graph_naive: lnnz (%"PRId64") != lnnz_orig (%"PRId64")\n", lnnz, lnnz_orig);
    return(NULL);
  }

  //fill in the weights
  if(weighted){
    int64_t i;
    for(i = 0; i < lnnz; i++){
      A->lvalue[i] = lgp_rand_double();
    }
  }

  return(A);
}

/*! \brief Generate a distributed graph that is the product of a
 collection of star graphs. This is done in two stages. In the first
 stage, the list of stars (parameterized by an integer m, K_{1,m}) is
 split in half and each half-list forms a local adjacency matrix for
 the Kronecker product of the stars (matrices B and C). Then the two
 local matrices are combined to form a distributed adjacency matrix
 for the Kronecker product of B and C.
 *
 * \param mode Mode 0 graphs have no triangles, mode 1 graphs have lots of triangles 
 *             and mode 2 graphs have few triangles.
 * \param spec An array of star sizes (if spec[i] = m, the ith element in the Kronecker product
 *             will be K_{1,m}.
 * \param num The number of stars in the spec list.
 *
 * See "Design, Generation, and Validation of Extreme Scale Power-Law Graphs" by Kepner et. al.
 * for more information on Kronecker product graphs. 
 *
 * \return A distributed matrix which represents the adjacency matrix for the 
 * Kronecker product of all the stars (B and C lists).
 */

sparsemat_t * generate_kronecker_graph_from_spec(int mode, int * spec, int num){
  int i;
  //T0_fprintf(stderr,"Generating Mode %d Kronecker Product graph (A = B X C) with parameters:  ", mode);
  if(num <=2){
    T0_fprintf(stderr,"ERROR: generate_kronecker: spec must contain more than 2 products\n");
    return(NULL);
  }
  int64_t *B_spec = calloc(num, sizeof(int64_t));
  int64_t *C_spec = calloc(num, sizeof(int64_t));
  int64_t B_num = num/2;
  for(i = 0; i < B_num; i++) B_spec[i] = spec[i];
  for(; i <  num; i++) C_spec[i - B_num] = spec[i];
  
  sparsemat_t * B = kronecker_product_of_stars(B_num, B_spec, mode);
  sparsemat_t * C = kronecker_product_of_stars(num - B_num, C_spec, mode);
  if(!B || !C){
    T0_fprintf(stderr,"ERROR: triangles: error generating input!\n"); lgp_global_exit(1);
  }
  
  //T0_fprintf(stderr,"B has %"PRId64" rows/cols and %"PRId64" nnz\n", B->numrows, B->lnnz);
  //T0_fprintf(stderr,"C has %"PRId64" rows/cols and %"PRId64" nnz\n", C->numrows, C->lnnz);
  
  sparsemat_t * A = kronecker_product_graph_dist(B, C);
  
  return(A);
}


/*! \brief Generate kron(B,C) in a distributed matrix (B and C are local matrices and all PEs have the same B and C).
 * \param B A local sparsemat_t (perhaps generated with kron_prod) that is the adjacency matrix of a graph.
 * \param C Another local sparsemat_t that is the adjacency matrix of a graph.
 * \return The adjacency matrix the kronecker product of B and C.
 * \ingroup spmatgrp
 */
sparsemat_t * kronecker_product_graph_dist(sparsemat_t * B, sparsemat_t * C) {

  int64_t lrow, row, rowB, rowC, col, pos, i, j, k;
  int64_t n = B->lnumrows*C->lnumrows;
  int64_t lnnz = 0;

  /* we first calculate the number of nonzeros on each thread */
  row = MYTHREAD;
  rowB = row / C->lnumrows;
  rowC = row % C->lnumrows;
  while(row < n){
    lnnz += (B->loffset[rowB + 1] - B->loffset[rowB])*(C->loffset[rowC + 1] - C->loffset[rowC]);
    row += THREADS;
    rowB = row / C->lnumrows;
    rowC = row % C->lnumrows;
  }
  
  sparsemat_t * A = init_matrix(n, n, lnnz, 0);
  if(A == NULL){
    printf("ERROR: kron_prod_dist: init_matrix returned NULL with inputs %"PRId64" %"PRId64" %"PRId64"\n", n, n, lnnz);
    return(NULL);
  }
  row = MYTHREAD;
  lrow = 0;
  rowB = row / C->lnumrows;
  rowC = row % C->lnumrows;
  pos = 0;
  A->loffset[0] = 0;
  
  while(row < n){
    
    for(j = B->loffset[rowB]; j < B->loffset[rowB + 1]; j++){
      int64_t Bcol = B->lnonzero[j];
      for(k = C->loffset[rowC]; k < C->loffset[rowC+1]; k++){
        col = Bcol*C->numcols + C->lnonzero[k];
        if(col < row)
          A->lnonzero[pos++] = col;
        else
          break;
      }      
    }
    A->loffset[lrow+1] = pos;
    
    row += THREADS;
    lrow++;
    rowB = row / C->lnumrows;
    rowC = row % C->lnumrows;
  }
  //if(pos != lnnz){
  //printf("ERROR: kron_prod_dist: pos (%"PRId64") != lnnz (%"PRId64")\n", pos, lnnz);
  //return(NULL);
  //}
  lgp_barrier();
  A->lnnz = pos;
  A->nnz = lgp_reduce_add_l(pos);

#if 0
  SHARED int64_t * sh_data = lgp_all_alloc(THREADS, sizeof(int64_t));
  int64_t * lsh_data = lgp_local_part(int64_t, sh_data);
  lsh_data[0] = 0;
  lgp_barrier();
  
  /* remove the self loop */
  if(!MYTHREAD){
    fprintf(stderr,"%"PRId64" %"PRId64" %"PRId64"\n", A->lnonzero[A->loffset[0]], A->loffset[0], A->loffset[1]);
    if(A->lnonzero[A->loffset[0]] == 0 && A->loffset[1] != A->loffset[0]){
      fprintf(stderr,"Removing A[0,0]!\n");
      /* set A(0,0) = 0*/
      for(i = 1; i < A->lnumrows+1; i++)
        A->loffset[i]--;
      for(i = 0; i < A->lnnz - 1; i++){
        A->lnonzero[i] = A->lnonzero[i+1];
      }
      A->lnnz--;      
      for(i = 0; i < THREADS; i++)
        lgp_atomic_add(sh_data, i, 1);
    }
  }
  if(MYTHREAD == ((A->numrows-1) % THREADS)){
    if(A->lnonzero[A->lnnz - 1] == (A->numrows-1)){
      fprintf(stderr,"Removing A[m,m]!\n");
      /* set A(m1*m2 - 1, m1*m2 - 1) = 0 */
      A->loffset[A->lnumrows]--;
      A->lnnz--;
      for(i = 0; i < THREADS; i++)
        lgp_atomic_add(sh_data, i, 1);
    }
  }
  lgp_barrier();
  if(lsh_data[0] == 1){
    A->nnz--;
  }
  if(lsh_data[0] > 1){
    T0_fprintf(stderr,"ERROR: kron_prod_dist: two diagonal elements!\n");
    return(NULL);
  }
#endif
  
  return(A);
}


/*! \brief Get the Kronecker product of two local matrices and store the result in a local matrix.
 * \param B A local sparse matrix (which is the adjacency matrix of a graph).
 * \param C A local sparse matrix (which is the adjacency matrix of a graph).
 * \return The adjacency matrix of the Kronecker product of B and C.
 * \ingroup spmatgrp
*/
sparsemat_t * kronecker_product_graph_local(sparsemat_t * B, sparsemat_t * C) {
  int64_t row, col, i, j, k;
  int64_t n = B->lnumrows*C->lnumrows;
  int64_t nnz = B->lnnz*C->lnnz;

  sparsemat_t * A = init_local_matrix(n, n, nnz);

  /* get the number of nonzeros in each row */
  int64_t * tmp_offset = calloc(n + 1, sizeof(int64_t));
  for(row = 0; row < B->lnumrows; row++){
    int64_t d1 = B->loffset[row + 1] - B->loffset[row];
    for(i = 0; i < C->lnumrows; i++){
      int64_t d2 = C->loffset[i + 1] - C->loffset[i]; 
      tmp_offset[1 + row*C->lnumrows + i] = d1*d2;
    }
  }

  /* get the cumulative sum */
  for(i = 0; i < n; i++){
    tmp_offset[i + 1] = tmp_offset[i + 1] + tmp_offset[i];
    A->loffset[i+1] = tmp_offset[i+1];
  }
  
  for(row = 0; row < B->lnumrows; row++){
    for(j = B->loffset[row]; j < B->loffset[row+1]; j++){
      col = B->lnonzero[j];
      for(i = 0; i < C->lnumrows; i++){
        for(k = C->loffset[i]; k < C->loffset[i+1]; k++){
          A->lnonzero[tmp_offset[row*C->lnumrows + i]++] = col*C->numcols + C->lnonzero[k];
        }
      }
    }
  }
  
  free(tmp_offset);

  return(A);
}

/*! \brief Generate the adjacency matrix for K_{1,m} (with 0 or 1 loop edges) 
 * \param m The parameter m in K_{1,m} (this is the number of non-hub vertices).
 * \param mode
 *  - mode == 0: default
 *  - mode == 1: add a self loop to center vertex of star
 *  - mode == 2: add a self loop to an outer vertex
 * \return The adjacency matrix of the graph.
 * \ingroup spmatgrp
*/
sparsemat_t * gen_star_graph(int64_t m, int mode) {
  sparsemat_t * G = init_local_matrix(m + 1, m + 1, 2*m + (mode > 0 ? 1 : 0));
  
  int64_t i, pos = 0;

  if(mode == 1) /* add a self loop to center vertex */    
    G->lnonzero[pos++] = 0;

  /* Add the top row */
  for(i = 0; i < m; i++){
    G->lnonzero[pos++] = i+1;
  }  
  G->loffset[1] = pos;

  /* rest of the rows */
  for(i = 1; i < m+1; i++){
    G->lnonzero[pos++] = 0;
    G->loffset[i+1] = pos;
  }
  
  if(mode == 2){/* add a self loop to last outer vertex */
    G->lnonzero[pos++] = m;
    G->loffset[m+1] = pos;
  }
  
  return(G);
}

/*! \brief Generate the kroncker product of a collection of star graphs. The result is a
 * local submatrix.
 * 
 * \param M The number of elements in the list. M must be at least 2 
 * \param m The list of parameters: each element m[i] represents a K_{1,m[i]} graph.
 * \param mode Add 0 or 1 loop edges.
 * - mode == 0: default
 * - mode == 1: add a self loop to center vertex of star
 * - mode == 2: add a self loop to an outer vertex
 * \return The local adjacency matrix for the kronecker product of the list of star graphs.
 * \ingroup spmatgrp
*/
sparsemat_t * kronecker_product_of_stars(int64_t M, int64_t * m, int mode) {
  int64_t i;

  if(M < 1){
    T0_printf("ERROR: gen_local_mat_from_stars: M < 1!\n");
    return(NULL);
  }

  if(M == 1){
    return(gen_star_graph(m[0], mode));
  }
  
  sparsemat_t ** Aarr = calloc(2*M, sizeof(sparsemat_t *));
  for(i = 0; i < M; i++)
    Aarr[i] = gen_star_graph(m[i], mode);

  Aarr[M] = kronecker_product_graph_local(Aarr[0], Aarr[1]);
  for(i = 0; i < M-2; i++)
    Aarr[M + 1 + i] = kronecker_product_graph_local(Aarr[M+i], Aarr[2+i]);
 
  sparsemat_t * A = Aarr[2*M-2];
  Aarr[2*M - 2] = NULL;

  for(i = 0; i < 2*M-2; i++){
    clear_matrix(Aarr[i]); free(Aarr[i]);
  }
  return(A);
}

/*! \brief sort the non-zeros in each row of a sparse matrix (make it tidy)
 * \param mat pointer to the sparse matrix
 * \ingroup spmatgrp
 */
int sort_nonzeros( sparsemat_t *mat) {
  int i,j;
  if(mat->value){
    // we have to sort the column indicies, but we also have to permute the value array accordingly
    // this is annoying in C
    // we have to create an array of stucts that holds col,val pairs for a row
    // sort that array according to the col keys
    // and then overwrite the row data
    int64_t max = 0;
    for(i = 0; i < mat->lnumrows; i++)
      if(mat->loffset[i+1] - mat->loffset[i] > max)
        max = mat->loffset[i+1] - mat->loffset[i];
    
    // allocate a temporary array to hold a row's worth of col, value pairs
    col_val_t * tmparr = calloc(max, sizeof(col_val_t));
    
    for(i = 0; i < mat->lnumrows; i++){
      int64_t pos = 0;
      for(j = mat->loffset[i]; j < mat->loffset[i+1]; j++){
        tmparr[pos].col = mat->lnonzero[j];
        tmparr[pos++].value = mat->lvalue[j];
        //fprintf(stderr, "value[%ld] = %lf\n", j , mat->lvalue[j]);
      }
      qsort(tmparr, mat->loffset[i+1] - mat->loffset[i], sizeof(col_val_t), col_val_comp );
      pos = 0;
      for(j = mat->loffset[i]; j < mat->loffset[i+1]; j++){
        mat->lnonzero[j] = tmparr[pos].col;
        mat->lvalue[j] = tmparr[pos++].value;
      }
    }
    free(tmparr);
    
  }else{

    for(i = 0; i < mat->lnumrows; i++){
      qsort( &(mat->lnonzero[mat->loffset[i]]), mat->loffset[i+1] - mat->loffset[i], sizeof(int64_t), nz_comp );
    }
  }
  lgp_barrier();
  return(0);
}

/*!
\brief compute the number of triangles in a Kronecker Product of stars graph
\param kron_mode the mode 0,1,2
\param kron_spec the list of sizes of the stars (note the +1)
\param kron_num the number of stars

See:
On Large-Scale Graph Generation with Validation of Diverse Triangle Statistics at Edges and Vertices
by Geoffrey Sanders, Roger Pearce, Timothy La Fond, Jeremy Kepner
*/

int64_t calc_num_tri_kron_graph(int kron_mode, int * kron_spec, int kron_num)
{
  int i;
  double correct_answer = 0;
  if(kron_mode == 0){
    correct_answer = 0.0;
  }else if(kron_mode == 1){
    correct_answer = 1;
    for(i = 0; i < kron_num; i++)
      correct_answer *= (3*kron_spec[i] + 1);
    correct_answer *= 1.0/6.0;
    double x = 1;
    for(i = 0; i < kron_num; i++)
      x *= (kron_spec[i] + 1);
    correct_answer = correct_answer - 0.5*x + 1.0/3.0;
  }else if(kron_mode == 2){
    correct_answer = (1.0/6.0)*pow(4,kron_num) - pow(2.0,(kron_num - 1)) + 1.0/3.0;
  }
  
  return((int64_t)round(correct_answer));
}


/*! \brief compare the structs that hold two sparse matrices
 * \param lmat pointer to the left sparse matrix
 * \param rmat pointer to the right sparse matrix
 * \return 0 on success
 * \ingroup spmatgrp
 */
int compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat) {
  int i,j;

  if( lmat->numrows != rmat->numrows ){
    fprintf(stderr,"PE %d: lmat->numrows = %"PRId64" != rmat->numrows = %"PRId64"\n",
           MYTHREAD, lmat->numrows, rmat->numrows );
    return(1);
  }
  if( lmat->lnumrows != rmat->lnumrows ){
    fprintf(stderr,"PE %d: lmat->lnumrows = %"PRId64"  !=rmat->lnumrows = %"PRId64"\n", 
            MYTHREAD, lmat->lnumrows, rmat->lnumrows );
    return(1);
  }
  if( lmat->numcols != rmat->numcols ){
    fprintf(stderr,"lmat->numcols = %"PRId64" != rmat->numcols = %"PRId64"\n", lmat->numcols, rmat->numcols );
    return(1);
  }
  if( lmat->nnz != rmat->nnz ){
    fprintf(stderr,"lmat->nnz = %"PRId64" != rmat->nnz = %"PRId64"\n",
            lmat->nnz, rmat->nnz );
    return(1);
  }
  if( lmat->lnnz != rmat->lnnz ){
    fprintf(stderr,"PE %d: lmat->lnnz = %"PRId64" != (rmat->lnnz = %"PRId64"\n", 
            MYTHREAD, lmat->lnnz, rmat->lnnz );
    return(1);
  }

  if( lmat->loffset[0] != 0 || rmat->loffset[0] != 0 
    || (lmat->loffset[0] != rmat->loffset[0] ) ){
    fprintf(stderr,"PE %d: lmat->loffset[0] = %"PRId64" != rmat->loffset[0] = %"PRId64"\n", 
            MYTHREAD, lmat->loffset[0], rmat->loffset[0] );
    return(1);
  }

  
  for(i = 0; i < lmat->lnumrows; i++){
    if( lmat->loffset[i+1] != rmat->loffset[i+1] ){
      fprintf(stderr,"PE %d: lmat->loffset[%d] = %"PRId64" != rmat->loffset[%d] = %"PRId64"\n", 
          MYTHREAD, i+1, lmat->loffset[i+1], i+1, rmat->loffset[i+1] );
       return(1);
    }
  }
  
  for(j=0; j< lmat->lnnz; j++) {
    if( lmat->lnonzero[j] != rmat->lnonzero[j] ){
      fprintf(stderr,"PE %d:lmat->lnonzero[%d] = %"PRId64" != rmat->lnonzero[%d] = %"PRId64"\n", 
                MYTHREAD, j, lmat->lnonzero[j], j, rmat->lnonzero[j] );
      return(1);
    }
  }
  if((lmat->value != 0) != (rmat->value != 0)){
    fprintf(stderr,"ERROR: compare_matrix: both matrices need to be weighted or unweighted\n");
    return(1);
  }
  if(lmat->value){
    for(j=0; j< lmat->lnnz; j++){
      if( spmat_are_eq_doubles(lmat->lvalue[j],rmat->lvalue[j])) {
        fprintf(stderr,"PE %d: lmat->lvalue[%d] = %lf != rmat->lvalue[%d] = %lf\n", 
                MYTHREAD, j, lmat->lvalue[j], j, rmat->lvalue[j] );
        return(1);
      }
    }
  }
  
  return(0);
}


/*! \brief produce a directed graph from a undirected graph
 * \param L undirected graph, a lower triangular matrix
 * \return the directed graph
 * \ingroup spmatgrp
 * We randomly orient each edge and build edge lists on the appropriate threads.
 * Then convert the edge lists to a sparse matrix.
 */
sparsemat_t * direct_undirected_graph(sparsemat_t * L){
  int64_t i, j;
  int weighted = (L->value != NULL);
  edge_list_t * el = init_edge_list(L->lnnz, weighted);

  w_edge_t edge;
  int64_t pe, col, row = 0;

  i = 0;
  exstack_t * ex = exstack_init(128, sizeof(w_edge_t));
  if(ex == NULL){return(NULL);}

  while(exstack_proceed(ex, (i == L->lnnz))){
    for(; i < L->lnnz; i++){
      while((row < L->lnumrows) && (i >= L->loffset[row+1]))
        row++;
      if(lgp_rand_double() < 0.5){
        if(weighted)
          append_weighted_edge(el, row*THREADS + MYTHREAD, L->lnonzero[i], L->lvalue[i]);
        else
          append_edge(el, row*THREADS + MYTHREAD, L->lnonzero[i]);
      }else{
        edge.row = L->lnonzero[i];
        edge.col = row*THREADS + MYTHREAD;
        if(weighted) edge.val = L->lvalue[i];
        pe = edge.row % THREADS;
        if(exstack_push(ex, &edge, pe) == 0L)
          break;
      }
    }
    
    exstack_exchange(ex);
      
    while(exstack_pop(ex, &edge, NULL)){
      //printf("B Appending edge %ld %ld\n", edge.row, edge.col);
      if(weighted)
        append_weighted_edge(el, edge.row, edge.col, edge.val);
      else
        append_edge(el, edge.row, edge.col);
    }
  }

  lgp_barrier();
  exstack_clear(ex);
  
  sparsemat_t * A = init_matrix(L->numrows, L->numcols, el->num, weighted);
  if(!A){
    T0_fprintf(stderr,"ERROR: direct_undirected_graph: Could not initialize A\n");
    return(NULL);
  }

  int64_t * lrowcounts = calloc(A->lnumrows, sizeof(int64_t));  
  for(i = 0; i < el->num; i++){
    if(weighted)
      lrowcounts[el->wedges[i].row/THREADS]++;
    else
      lrowcounts[el->edges[i].row/THREADS]++;
  }
  
  // initialize the offsets for the sparse matrix
  A->loffset[0] = 0;
  for(i = 1; i <= A->lnumrows; i++){
    A->loffset[i] = A->loffset[i - 1] + lrowcounts[i-1];
    lrowcounts[i-1] = 0;
  }
  assert(A->loffset[A->lnumrows] == el->num);

  // populate A
  for(i = 0; i < el->num; i++){
    if(weighted){
      row = el->wedges[i].row/THREADS;
      int64_t pos = A->loffset[row] + lrowcounts[row]++;
      A->lnonzero[pos] = el->wedges[i].col;
      A->lvalue[pos] = el->wedges[i].val;
    }else{
      row = el->edges[i].row/THREADS;
      int64_t pos = A->loffset[row] + lrowcounts[row]++;
      A->lnonzero[pos] = el->edges[i].col;
    }
  }

  sort_nonzeros(A);
  free(lrowcounts);
  clear_edge_list(el);
  
  return(A);
  
}

/*! 
 * \brief produce an exact copy of a given sparse matrices
 * \param srcmat pointer to the original sparse matrix
 * \return A pointer to the cloned sparse matrix
 */
sparsemat_t * copy_matrix(sparsemat_t *srcmat) {
  int i,j;
  int64_t numrows, numcols, lnumrows, lnumcols;
  
  sparsemat_t * destmat = init_matrix(srcmat->numrows, srcmat->numcols, srcmat->lnnz, (srcmat->value != NULL));
  if(!destmat) return(NULL);

  for(i = 0; i < (srcmat->lnumrows)+1; i++){
     destmat->loffset[i] = srcmat->loffset[i];
  }
  for(j=0; j < srcmat->lnnz; j++) {    
    destmat->lnonzero[j] = srcmat->lnonzero[j];
    if(srcmat->value) destmat->lvalue[j] = srcmat->lvalue[j];
  }

  lgp_barrier();
  return(destmat);
}

// IF we are given nz_per_row (z), we calculate edge_prob (e)
// or if we are given edge_prob, we calculate nz_per_row
// using with the following formulas
// z*n = e*(n*(n-1)/2)
// or
// (z-1)*n = e*(n*(n-1)/2) (if we are forcing loops into the graph)
//

/*! 
 * \brief computes the number of non-zeros per row from the edge probability or vice-versa.
 * \param edge_prob given edge probability (or place to put the computed edge_prob)
 * \param nz_per_row  given number of nonzeros per row (or place to put the computed nz_per_row)
 * \param numrows = numcols global order of the matrix 
 * \param edgetype enum edge_type weighted or unweighted
 * \param loops whether the diagonal is all zeros or all ones
 */
void resolve_edge_prob_and_nz_per_row(double * edge_prob, double * nz_per_row, int64_t numrows, edge_type edgetype, self_loops loops){
  if(*edge_prob == 0.0){ // use nz_per_row to get erdos_renyi_prob
    if(loops == LOOPS)
      *edge_prob = (*nz_per_row - 1)/(numrows-1);
    else
      *edge_prob = (*nz_per_row)/(numrows-1);
    
    if (edgetype == UNDIRECTED || edgetype == UNDIRECTED_WEIGHTED)
      *edge_prob = *edge_prob*2;
        
    if(*edge_prob > 1.0)
      *edge_prob = 1.0;
  } else {    // use erdos_renyi_prob to get nz_per_row
    *nz_per_row = *edge_prob * ((numrows - 1)/2.0);
  }
  assert(*edge_prob <= 1.0);
}


// return 1 if the two values are different "enough" relative to their size.
/*! 
 * \brief compare to doubles
 * \param a first one
 * \param b the other one
 * \return int 0, if they are close enough 
 */
int spmat_are_eq_doubles(double a, double b) {
  if(a == 0.0){
    if(fabs(b) < 1e-8) return(0);
    return(1);
  }
  if((fabs(a - b)  / fabs(a)) < FLT_EPSILON) return(0);
  return(1);
}

/*! \brief initializes the struct that holds a sparse matrix
 *    given the total number of rows and columns and the local number of non-zeros
 * \param numrows total (global) number of rows
 * \param numcols total (global) number of columns  (always equals numrows in bale)
 * \param nnz_this_thread number of nonzero on this thread
 * \param weighted whether the nonzeros have weights or not
 * \return An initialized sparsemat_t or NULL on error.
 * \ingroup spmatgrp
 */
sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz_this_thread, int weighted) 
{
  sparsemat_t * mat = calloc(1, sizeof(sparsemat_t));
  mat->local = 0;
  mat->numrows  = numrows;
  mat->lnumrows = (numrows + THREADS - MYTHREAD - 1)/THREADS;
  mat->numcols  = numcols;  
  mat->offset   = lgp_all_alloc(mat->numrows + THREADS, sizeof(int64_t));
  if(mat->offset == NULL){
    T0_printf("ERROR: init_matrix: could not allocate %"PRId64" bytes for offset array\n", mat->numrows*8);
    return(NULL);
  }
  mat->loffset  =  lgp_local_part(int64_t, mat->offset);
  int64_t max = lgp_reduce_max_l(nnz_this_thread);
  int64_t total = lgp_reduce_add_l(nnz_this_thread);
  mat->nonzero = lgp_all_alloc(max*THREADS, sizeof(int64_t));
  if(mat->nonzero == NULL){
    T0_printf("ERROR: init_matrix: could not allocate %"PRId64" bytes for nonzero array (max = %"PRId64")\n",
              max*THREADS*8, max);
    return(NULL);
  }
  mat->lnonzero = lgp_local_part(int64_t, mat->nonzero);

  if(weighted){
    mat->value = lgp_all_alloc(max*THREADS, sizeof(double));
    if(mat->value == NULL){
      T0_printf("ERROR: init_matrix: could not allocate %"PRId64" bytes for value array (max = %"PRId64")\n",
                max*THREADS*8, max);
      return(NULL);
    }
    mat->lvalue = lgp_local_part(double, mat->value);
  }else{
    mat->value = NULL;
    mat->lvalue = NULL;
  }
  
  mat->nnz = total;
  mat->lnnz = nnz_this_thread;

  lgp_barrier();
  
  return(mat);
}

/*! \brief initialize the local parts of a sparsemat_t
 * \param numrows the number of rows
 * \param numcols the number of cols
 * \param nnz the number of nonzeros
 * \ingroup spmatgrp
 */
sparsemat_t * init_local_matrix(int64_t numrows, int64_t numcols, int64_t nnz) {
  sparsemat_t * mat = calloc(1, sizeof(sparsemat_t));
  mat->local = 1;
  mat->numrows  = numrows;
  mat->lnumrows = numrows;
  mat->numcols  = numcols;  
  mat->offset   = NULL;
  mat->loffset  = calloc(mat->lnumrows + 1, sizeof(int64_t));
  if(mat->loffset == NULL) 
    return(NULL);
  
  mat->nonzero = NULL;
  mat->lnonzero = calloc(nnz, sizeof(int64_t));
  if(mat->lnonzero == NULL) 
    return(NULL);

  mat->nnz = nnz;
  mat->lnnz = nnz;

  return(mat);
}

/*! 
 * \brief initialize an edge list
 * \param nalloc the number tuples or triples to allocate
 * \param weighted whether or not we are using weights
 * \return a pointer to the list
 * \ingroup spmatgrp
 */
edge_list_t * init_edge_list(int64_t nalloc, int weighted){
  edge_list_t * el = calloc(1, sizeof(edge_list_t));
  el->num = 0;
  el->nalloc = nalloc;
  if(!weighted){
    el->edges = calloc(nalloc, sizeof(edge_t));
    el->wedges = NULL;
  }else{
    el->edges = NULL;
    el->wedges = calloc(nalloc, sizeof(w_edge_t));
  }
  return(el);
}

/*! 
 * \brief free the memory using by the edge list
 * \param el the edge list
 * \ingroup spmatgrp
 */
void clear_edge_list(edge_list_t * el){
  free(el->edges);
  free(el);
}


/*! 
 * \brief initialize a triples_t struct to hold a sparse matrix
 * \param numrows the number of rows in the matrix
 * \param numcols the number of column in the matrix (same as numrows for bale)
 * \param nalloc number of triples to allocate
 * \param weighted whether or not we are using weights
 * \return a pointer to the triples_t 
 * \ingroup spmatgrp
 */
triples_t * init_triples(int64_t numrows, int64_t numcols, int64_t nalloc, int weighted){
  triples_t * T = calloc(1, sizeof(triples_t));
  T->numrows = numrows;
  T->numcols = numcols;
  T->lnnz = 0;
  T->nalloc = nalloc;
  T->row = calloc(T->nalloc, sizeof(int64_t));
  T->col = calloc(T->nalloc, sizeof(int64_t));
  if(weighted) T->val = calloc(T->nalloc, sizeof(double));
  else T->val = NULL;
  return(T);
}

/*! \brief frees the space allocated for a sparse matrix
 * \param mat pointer to the sparse matrix
 * \ingroup spmatgrp
 */
void clear_matrix(sparsemat_t * mat) {
  if(mat->local){
    free(mat->lnonzero);
    free(mat->loffset);
  }else{
    lgp_all_free(mat->nonzero);
    lgp_all_free(mat->offset);
  }
}

/*! \brief frees the space allocated for a triples_t
 * \param T pointer to the triples_t
 * \ingroup spmatgrp
 */
void clear_triples(triples_t * T) {
  free(T->val);
  free(T->row);
  free(T->col);
}

/*! \brief frees the space allocated for an spmat_dataset_t
 * \param spd pointer to the spmat_dataset_t
 * \ingroup spmatgrp
 */
void clear_spmat_dataset(spmat_dataset_t * spd){
  free(spd->dirname);
  free(spd->nrows_in_file);
  free(spd->nrows_read_in_file);
  free(spd->rowcnt);
  free(spd->global_first_row_to_me);
}

/*! \brief Append a triple to a triples_t struct and expand the storage if necessary
 *
 * \param T The triples_t struct
 * \param row The row index of the triple
 * \param col The column index of the triples
 * \param val If T->val is not NULL, the value of the triple.
 * \return The number of triples after this triple was appended.
 */
int64_t append_triple(triples_t * T, int64_t row, int64_t col, double val){
  if(T->nalloc == T->lnnz){
    // we need to expand our allocations!
    if(T->nalloc < 10000)
      T->nalloc = 2*T->nalloc;
    else
      T->nalloc = T->nalloc*1.25;
    T->row = realloc(T->row, T->nalloc*sizeof(int64_t));
    T->col = realloc(T->col, T->nalloc*sizeof(int64_t));
    if(T->val)T->val = realloc(T->val, T->nalloc*sizeof(double));
  }
  T->row[T->lnnz] = row;
  T->col[T->lnnz] = col;
  if(T->val)T->val[T->lnnz] = val;
  T->lnnz++;
  return(T->lnnz);
}

/*! 
 * \brief Append an edge to a edge_list and expand the storage if necessary
 * \param el the edge_list
 * \param row the row index of the edge
 * \param col the column index of the edge
 * \return the number of edges after this edge was appended.
 */
int64_t append_edge(edge_list_t * el, int64_t row, int64_t col){
  if(el->nalloc == el->num){
    //printf("PE %d: out of space! nalloc = %ld\n", MYTHREAD, el->nalloc);
    // we need to expand our allocations!
    if(el->nalloc < 10000)
      el->nalloc = 2*el->nalloc;
    else
      el->nalloc = el->nalloc*1.25;
    //printf("PE %d: new space! nalloc = %ld\n", MYTHREAD, el->nalloc);
    el->edges = realloc(el->edges, el->nalloc*sizeof(edge_t));
  }
  el->edges[el->num].row = row;
  el->edges[el->num].col = col;
  el->num++;
  return(el->num);
}

/*! 
 * \brief Append a weighted edge to a edge_list and expand the storage if necessary
 * \param el the edge_list
 * \param row the row index of the edge
 * \param col the column index of the edge
 * \param val the value 
 * \return the number of edges after this edge was appended.
 */
int64_t append_weighted_edge(edge_list_t * el, int64_t row, int64_t col, double val){
    if(el->nalloc == el->num){
    //printf("PE %d: out of space! nalloc = %ld\n", MYTHREAD, el->nalloc);
    // we need to expand our allocations!
    if(el->nalloc < 10000)
      el->nalloc = 2*el->nalloc;
    else
      el->nalloc = el->nalloc*1.25;
    //printf("PE %d: new space! nalloc = %ld\n", MYTHREAD, el->nalloc);
    el->wedges = realloc(el->wedges, el->nalloc*sizeof(w_edge_t));
  }
  el->wedges[el->num].row = row;
  el->wedges[el->num].col = col;
  el->wedges[el->num].val = val;
  el->num++;
  return(el->num);
}

/*!
 * \brief Convert a triples_t to a sparsemat_t.
 * \param T A triples_t struct.
 * \return The sparsemat_t version of T-> Or NULL on error.
 */
sparsemat_t * triples_to_sparsemat(triples_t * T){

  int64_t i, j;
  // get rowcounts (to create tmp offset array)
  int64_t * tmp = calloc(T->lnumrows + 1, sizeof(int64_t));
  for(i = 0; i < T->lnnz; i++){
    tmp[T->row[i] + 1]++;
  }
  for(i = 0; i < T->lnumrows; i++)
    tmp[i+1] += tmp[i];

  sparsemat_t * A = init_matrix(T->numrows, T->numcols, T->lnnz, (T->val != NULL));
  if(!A){printf("ERROR: triples_to_sparsemat: init failed.\n"); return(NULL);}

  // copy tmp array to A
  for(i = 0; i <= T->lnumrows; i++)
    A->loffset[i] = tmp[i];
  
  for(i = 0; i < T->lnnz; i++){
    j = tmp[T->row[i]];
    A->lnonzero[j] = T->col[i];
    if(T->val) A->value[j] = T->val[i];
  }
  free(tmp);
  
  return(A);
}

/*!
 * \brief debugging routine to print a matrix to stdout
 * \param A the sparsemat
 */
void print_matrix(sparsemat_t * A){
  int64_t i, j;
  
  for(i = 0; i < A->numrows; i++){
    T0_printf("row %ld: ",i);
    int64_t start = lgp_get_int64(A->offset, i);
    int64_t end = lgp_get_int64(A->offset, i+THREADS);
    for(j = start; j < end; j++){
      if( A->value != NULL ){
        T0_printf("(%ld, %lg) ", lgp_get_int64(A->nonzero, j*THREADS + i%THREADS),
                  lgp_get_double(A->value, j*THREADS + i%THREADS));
      } else {
        T0_printf("%ld ", lgp_get_int64(A->nonzero, j*THREADS + i%THREADS));
      }
    }
    T0_printf("\n");
  }
}

/*!
 * \brief returns the number of nonzeros in a row of the localize part of a sparse matrix
 * \param *mat pointer to the sparse matrix 
 * \param l_row row index into local version (loffset) of the offset array
 * \return the number of nonzero in the request row
 */
int64_t rowcount_l( sparsemat_t *mat, int64_t l_row ) {
   return( mat->loffset[l_row+1] - mat->loffset[l_row] );
}

/*!
 * \brief returns the number of nonzeros in a row of a sparse matrix
 * \param *mat pointer to the sparse matrix 
 * \param S_row global row index into the shared offset array
 * \return the number of nonzero in the request row
 */
int64_t rowcount_S( sparsemat_t *mat, int64_t S_row ) {
   return( mat->offset[S_row+THREADS] - mat->offset[S_row] );
}

/*!
 * \brief Initialize an nxnz_t struct 
 * \param *mat pointer to the sparse matrix 
 * \return a pointer to the nxnz_t struct associated with the matrix, rest of the state is set to meaningless values.
 */
nxnz_t * init_nxnz( sparsemat_t *mat ) {
  nxnz_t *ret = calloc(1,sizeof(nxnz_t));
  ret->mat = mat;
  ret->row = -1;
  ret->first = 0x7FFFFFFFFFFFFFFF;
  ret->idx = -1;
  ret->stop = -1;
  ret->col = -1;
  return(ret);
};

/*!
 * \brief Sets the state of an nxnz_t struct for local first touch of the given row
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param l_row the local version of said row
 * \return void
 */
void first_l_nxnz( nxnz_t *nxz, int64_t l_row ) {
  assert( l_row < nxz->mat->lnumrows );
  nxz->row   = l_row;
  nxz->first = nxz->mat->loffset[l_row];
  nxz->idx   = nxz->mat->loffset[l_row];
  nxz->stop  = nxz->mat->loffset[l_row+1];
  nxz->col   = nxz->mat->lnonzero[nxz->idx];
}

/*!
 * \brief The condition part of a for loop across the local view of the given row
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param l_row the local version of said row
 * \return true if we are at a valid index, false otherwise
 */
bool has_l_nxnz( nxnz_t *nxz, int64_t l_row ) {
  assert( l_row == nxz->row );
  if( nxz->idx < nxz->stop )
    return(true);
  return(false);
}

/*!
 * \brief Attempts to move to the next nonzero in the local view of the given row,
 *   by incrementing the index and setting the col field on the nxnz_t.
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param l_row the local version of said row
 * \return void
 */
void incr_l_nxnz( nxnz_t *nxz, int64_t l_row ) {
  assert( l_row == nxz->row );
  nxz->idx += 1;
  nxz->col   = nxz->mat->lnonzero[nxz->idx];
}

/*!
 * \brief Sets the state of an nxnz_t struct for first touch of the given row,
 *   from the shared view of the matrix
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param S_row the shared version of said row
 * \return void
 */
void first_S_nxnz( nxnz_t *nxz, int64_t S_row ) {
  assert( S_row < nxz->mat->numrows );
  nxz->row   = S_row;
  nxz->first = lgp_get_int64(nxz->mat->offset, S_row);
  nxz->idx   = nxz->first;
  nxz->stop  = lgp_get_int64(nxz->mat->offset, S_row+THREADS);
  nxz->col   = lgp_get_int64(nxz->mat->nonzero, (nxz->idx) * THREADS +  (nxz->row)%THREADS );
}

/*!
 * \brief The condition part of a for loop across the shared view of the given row
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param S_row the shared version of said row
 * \return true if we are at a valid index, false otherwise
 */
bool has_S_nxnz( nxnz_t *nxz, int64_t S_row ) {
  assert( S_row == nxz->row );
  if( nxz->idx < nxz->stop )
    return(true);
  return(false);
}

/*!
 * \brief Attempts to move to the next nonzero in the shared view of the given row,
 *   by incrementing the index and setting the col field of the nxnz_t.
 * \param *nxz pointer the nxnz_t to hold the state from the specified row
 * \param S_row the shared version of said row
 * \return void
 */
void incr_S_nxnz( nxnz_t *nxz, int64_t S_row ) {
  assert( S_row == nxz->row );
  nxz->idx += 1; 
  nxz->col   = lgp_get_int64(nxz->mat->nonzero, (nxz->idx) * THREADS +  (nxz->row)%THREADS );
}

/*! \brief comparison function to support qsort */
int nz_comp(const void *a, const void *b) {
  return( *(uint64_t *)a - *(uint64_t *)b );
}

/*! \brief comparison function to support qsort */
int point_comp(const void *a, const void *b) {
  point_t *p1 = (point_t *)a;
  point_t *p2 = (point_t *)b;
  if(p1->x > p2->x) return 1;
  else if(p1->x == p2->x)
    if (p1->y > p2->y)
      return 1;
  return(-1);
}

/*! \brief comparison function to support qsort in sort_nonzeros */
int col_val_comp(const void *a, const void *b) {
  col_val_t * A = (col_val_t*)a;
  col_val_t * B = (col_val_t*)b;
  return((int)(A->col - B->col));
}

/*! \brief the compare function for qsort called while reading a MatrixMarket format.
 * NB. We sort on the rows so that we can fill the offset array sequentially in one pass. 
 * We sort on the columns so that the matrix will be "tidy"
 */
int edge_comp(const void *a, const void *b) 
{
  edge_t * eltA = (edge_t *)a;
  edge_t * eltB = (edge_t *)b;
  if( (eltA->row - eltB->row) == 0 )
    return( eltA->col - eltB->col );
  return( eltA->row - eltB->row );
}

/*! \brief the compare function for qsort of w_edge_t structs.
 * NB. We sort on the rows so that we can fill the offset array sequentially in one pass. 
 * We sort on the columns so that the matrix will be "tidy"
 */
int w_edge_comp(const void *a, const void *b) 
{
  w_edge_t * A = (w_edge_t *)a;
  w_edge_t * B = (w_edge_t *)b;
  if( (A->row - B->row) == 0 )
    return( A->col - B->col );
  return( A->row - B->row );
}



/*! 
 * \brief initializes the struct that holds a distributed array of doubles
 * \param num total number of entries
 * \return an allocated d_array_t or NULL on error
 * \ingroup spmatgrp
 */
d_array_t * init_d_array(int64_t num) 
{
  d_array_t * array = calloc(1, sizeof(d_array_t));
  array->num  = num;
  array->lnum = (num + THREADS - MYTHREAD - 1)/THREADS;
  array->entry   = lgp_all_alloc(num + THREADS, sizeof(double));
  if(array->entry == NULL){
    T0_printf("ERROR: init_d_array: could not allocate %"PRId64" bytes for array\n", num*sizeof(double));
    return(NULL);
  }
  array->lentry  =  lgp_local_part(double, array->entry);

  return(array);
}

/*! 
 * \brief sets all the entries of d_array_t to a value
 * \param A the array
 * \param v the value
 * \ingroup spmatgrp
 */
void set_d_array(d_array_t * A, double v) 
{
  int64_t i;
  for(i=0; i<A->lnum; i++) {
    A->lentry[i] = v;
  }
}

/*! 
 * \brief produces a copy of a source array
 * \param S the source array
 * \ingroup spmatgrp
 * this is a collective operation with a barrier before and aft
 */
d_array_t * copy_d_array(d_array_t * S) 
{
  int64_t i;
  lgp_barrier(); 

  d_array_t * ret = init_d_array(S->num);

  ret->num  = S->num;
  ret->lnum = S->lnum;
  for(i=0; i<S->lnum; i++) {
    ret->lentry[i] = S->lentry[i];
  }
  lgp_barrier(); 
  return(ret);
}


/*! 
 * \brief replaces the destination array with a copy of the source array
 * \param D the destination array
 * \param S the source array
 * \ingroup spmatgrp
 * this is a collective operation with a barrier before and aft
 * Note: The destination must be allocated and of the right size.
 */
int64_t replace_d_array(d_array_t * D, d_array_t * S) 
{

  int64_t i;

  lgp_barrier(); 
  if(D->num != S->num || D->lnum != S->lnum){
    T0_printf("ERROR: replace_d_array: arrays lengths don't match\n");
    return(0);
  }
  for(i=0; i<S->lnum; i++) {
    D->lentry[i] = S->lentry[i];
  }
  lgp_barrier(); 
  return(1);
}

/*! 
 * \brief frees the space used by the d_array_t
 * \param A the d_array
 */
void clear_d_array(d_array_t * A)
{
  lgp_all_free(A->entry);
}


