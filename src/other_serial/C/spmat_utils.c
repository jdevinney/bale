/*******************************************************************/
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
/*******************************************************************/

/*! \file spmat_utils.c
 * \brief Routines to implement basic spmat functions and some
 * utilities for debugging and timing.
 */

#include "spmat_utils.h"


/*! \brief  
A routine to give access to the wall clock timer on most UNIX-like systems.
\return the number of seconds since the beginning of time on UNIX machines.
Uses gettimeofday.
*/
double wall_seconds() 
{
  struct timeval tp;
  int retVal = gettimeofday(&tp,NULL);
  if (retVal == -1){ perror("gettimeofday"); fflush(stderr); }
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );

}



/*! 
\brief set the seed for the random number generator
\param seed the seed
*/
void rand_seed(int64_t seed){
#ifdef USE_KNUTH
  ranf_start(seed);
#else
  srand48(seed);
#endif
}

/*!
\brief return the next random integer mod N
\param N the modulus
\return a random number from 0 to N-1
*/
int64_t rand_int64(int64_t N){
  assert(N < CBALE_RAND_MAX);
#ifdef USE_KNUTH
  return((int64_t)(ranf_arr_next()*N));
#else
  return((int64_t)(drand48()*N));
#endif
}


/*!
\brief return the next random double
\return a random number from 0 to N-1
*/
double rand_double(){
#ifdef USE_KNUTH
  return(ranf_arr_next());
#else
  return(drand48());
#endif
}


/*!
\brief create an int64_t array which holds a uniform random permutation
\param N the length of the global array
\param seed seed for the random number generator
\return the permutation (return identity perm if seed is 0)

This implements the standard serial algorithm, known at least as
Fisher-Yates or Knuth shuffle, to generate the uniform permutation.

Start with an array holding the identity permutation,
  - swap the last entry with a random entry in the array,
  - this determines the last entry,
  - shorten the array to be all but the last entry 
  -  and repeat
*/
int64_t * rand_perm(int64_t N, int64_t seed) 
{
  int64_t r, i, L, s;

  if( seed != 0 ) 
    rand_seed( seed );

  int64_t * perm = malloc(N * sizeof(int64_t));
  if( perm == NULL ) return(NULL);
  for(i=0; i<N; i++)
    perm[i] = i;
  if( seed == 0 )
    return(perm);
  for(L=N-1; L>0; L--){
    r = rand_int64(L+1);
    s = perm[L];
    perm[L] = perm[r];
    perm[r] = s;
  }
  return(perm);
}


/*! \brief checks that an array is in fact a permutation
 * \param perm pointer to the array
 * \param N the length of the array
 * \return 1 if it is permutation
 *
 * Every element in the flag array will equal 1 iff perm is a permutation.
 */
int64_t is_perm(int64_t * perm, int64_t N) 
{
  int64_t i, err=0L;
  int64_t * flag = calloc(N, sizeof(int64_t));
  if( flag == NULL ) return(0);

  for(i = 0; i < N; i++){
    if( 0 > perm[i] || perm[i] >= N )
      err++;
    else 
      flag[perm[i]]++;
  }
  for(i = 0; i < N; i++) 
    if(flag[i] != 1L) 
      err++;
  free(flag);
  return(err == 0L);
}

/*!
\brief writes the first and last part of an array of int64_t's to the specified file
\param A pointer to the array
\param len the length of the array
\param maxdisp the number of entries that are written, 0 means everything, 
  otherwise write the first and last maxdisplay/2 entries
\param name the filename to written to
\return 0 on success, non-0 on error.
*/
int64_t dump_array(int64_t *A, int64_t len, int64_t maxdisp, char * name) 
{
  int64_t i, stoprow, startrow;
  stoprow = 0;
  startrow = 0;
  if((maxdisp > 0) && (maxdisp < len)){
    stoprow = maxdisp / 2;
    startrow = len - maxdisp/2;
  }

  FILE * fp = fopen(name, "w");
  for(i=0; i<len; i = (i<stoprow || i>=startrow)? i+1: startrow){
    fprintf(fp, "%3"PRId64"\n", A[i]);
  }
  fclose(fp);
  return(0);
}


/*!
\brief dumps a sparse matrix to a file in a ASCII format
\param A pointer to the sparse matrix
\param maxrows the number of rows that are written, 0 means everything, 
  otherwise write the first and last maxrows/2 rows
\param name the filename to written to
\return 0 on success, non-0 on error.
*/
int64_t dump_matrix(sparsemat_t *A, int64_t maxrows, char * name) 
{
  int64_t i,j, off, nxtoff;
  int64_t stoprow, startrow;

  if(!A){
    fprintf(stderr,"ERROR: trying to dump an empty matrix");
    return(0);
  }

  stoprow = A->numrows;
  startrow = A->numrows;
  if( (maxrows > 0) && (maxrows < A->numrows) ){
    stoprow = maxrows / 2;
    startrow = A->numrows - maxrows/2;
  }

  FILE * fp = fopen(name, "w");
  fprintf(fp,"\n----- offsets:\n");
  for(i=0; i<stoprow; i++) 
    fprintf(fp, "%3"PRId64" ", A->offset[i]);
  if( stoprow < startrow)
    fprintf(fp, " ... ");
  for(i=startrow; i<A->numrows; i++) 
    fprintf(fp, "%3"PRId64" ", A->offset[i]);
  fprintf(fp, "%3"PRId64" ", A->offset[A->numrows]);

  fprintf(fp,"\n--------- nonzeros:\n");
  for(i=0; i<stoprow; i++){
    off    = A->offset[i];
    nxtoff = A->offset[i+1];
    for(j=off; j < nxtoff;j++){
      fprintf(fp, "%9"PRId64" %9"PRId64"",i, A->nonzero[j] );
      if(A->value)
        fprintf(fp," %9.5f\n",A->value[j]);
      else
        fprintf(fp,"\n");
    }
  }
  if( stoprow < startrow)
    fprintf(fp, ".\n.\n.\n");
  for(i=startrow; i<A->numrows; i++){
    off    = A->offset[i];
    nxtoff = A->offset[i+1];
    for(j=off; j < nxtoff;j++){
      fprintf(fp, "%9"PRId64" %9"PRId64"",i, A->nonzero[j] );
      if(A->value)
        fprintf(fp," %9.5f\n",A->value[j]);
      else
        fprintf(fp,"\n");
    }
  }
  fclose(fp);
  return(1);
}

/*! \brief writes a sparse matrix to a file in a MatrixMarket ASCII format
 * \param A pointer to the sparse matrix
 * \param name the filename to be written to
 * \return 0 on success, non-0 on error.
 */
int64_t write_matrix_mm(sparsemat_t *A, char * name){
  int64_t i,j;

  FILE * fp = fopen(name, "w");
  if( fp == NULL ){
    fprintf(stderr,"write_matrix_mm: can't open file %s \n", name);
    exit(1);
  }
  if(A->value)
    fprintf(fp,"%%%%MatrixMarket matrix coordinate real\n");
  else
    fprintf(fp,"%%%%MatrixMarket matrix coordinate pattern\n");
  fprintf(fp, "%"PRId64" %"PRId64" %"PRId64"\n", A->numrows, A->numcols, A->nnz);

  for(i=0; i<A->numrows; i++){
    for(j=A->offset[i]; j < A->offset[i+1];j++)
      if(A->value)
        fprintf(fp, "%"PRId64" %"PRId64" %lf\n",i+1, A->nonzero[j]+1, A->value[j]);
      else
        fprintf(fp, "%"PRId64" %"PRId64"\n",i+1, A->nonzero[j]+1);
  }
  fclose(fp);
  return(0);
}


/*! \brief the compare function for edges for qsort called while reading 
 * a MatrixMarket format (and other places).
 * NB. We sort on the rows so that we can fill the offset array
 * sequentially in one pass. We sort on the columns so that
 * the matrix will be "tidy"
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
 * NB. We sort on the rows so that we can fill the offset array
 * sequentially in one pass. We sort on the columns so that
 * the matrix will be "tidy"
 */
int w_edge_comp(const void *a, const void *b) 
{
  w_edge_t * A = (w_edge_t *)a;
  w_edge_t * B = (w_edge_t *)b;
  if( (A->row - B->row) == 0 )
    return( A->col - B->col );
  return( A->row - B->row );
}

/*! \brief read a sparse matrix from a file in a MatrixMarket ASCII format
 * \param name the filename to be read
 * \return a pointer to the sparse matrix or NULL on failure
 */
sparsemat_t  *read_matrix_mm(char * name) 
{
  int64_t i;
  int64_t nr, nc, nnz;
  int fscanfret;

  // Read the header line of the MatrixMarket format 
  FILE * fp = fopen(name, "r");
  if ( fp == NULL ) {
    fprintf(stderr,"read_matrix_mm: can't open file %s \n", name);
    exit(1);
  }
    
  char object[24], format[24], field[24];
  fscanfret = fscanf(fp,"%%%%MatrixMarket %s %s %s\n", object, format, field);
  if ( (fscanfret != 3 ) || strncmp(object,"matrix",24) || strncmp(format,"coordinate",24) ) {
    fprintf(stderr,"read_matrix_mm: Incompatible matrix market format.\n");
    fprintf(stderr,"                First line should be either:\n");
    fprintf(stderr,"                matrix coordinate pattern\n");
    fprintf(stderr,"                OR\n");
    fprintf(stderr,"                matrix coordinate real\n");
    fprintf(stderr,"                OR\n");
    fprintf(stderr,"                matrix coordinate integer\n");
    exit(1);
  }
  
  if (strncmp(field,"pattern",24) && strncmp(field,"real",24) && strncmp(field,"integer",24) ) {
    fprintf(stderr,"read_matrix_mm: Incompatible matrix market field.\n");
    fprintf(stderr,"                Last entry on first line should be pattern, real, or integer\n");
    exit(1);
  }
  int values;
  if (strncmp(field,"pattern",24) == 0) {
    values = 0; // no values
  } else if (strncmp(field,"real",24) == 0) {
    values = 1; // real values
  } else {
    values = 2; // integer values
  }

  // Read the header (nr, nc, nnz)
  fscanfret = fscanf(fp,"%"SCNd64" %"SCNd64" %"SCNd64"\n", &nr, &nc, &nnz);
  if ( (fscanfret != 3 ) || (nr<=0) || (nc<=0) || (nnz<=0) ) {
    fprintf(stderr,"read_matrix_mm: reading nr, nc, nnz\n");
    exit(1);
  }
  
  sparsemat_t * ret_mat;
  edge_t *edges;
  w_edge_t *w_edges;
  if (values == 0L) { // no values
    
    // read all the nonzeros into the edges array of (row,col,[values])
    // and sort them before building the sparsemat format
    edges = calloc(nnz, sizeof(edge_t));
    assert((edges != NULL) && "Memory Error");

    for (i=0; i<nnz; i++) {
      fscanfret = fscanf(fp,"%"SCNd64" %"SCNd64"\n", &(edges[i].row), &(edges[i].col));
      assert (fscanfret == 2);
      //fprintf(stderr,"--- %"PRId64" %"PRId64"\n",  edges[i].row, edges[i].col);
      assert ( 0<edges[i].row && edges[i].row <=nr);
      assert ( 0<edges[i].col && edges[i].col <=nc);
      edges[i].row -=1;    // MatrixMarket format is 1-up, not 0-up
      edges[i].col -=1;
    }
    qsort( edges, nnz, sizeof(edge_t), edge_comp);
  } else { // real 
    assert(values == 1);
    w_edges = calloc(nnz, sizeof(w_edge_t));
    assert((w_edges != NULL) && "Memory Error");

    for (i=0; i<nnz; i++) {
      fscanfret = fscanf(fp,"%"SCNd64" %"SCNd64" %lf\n", &(w_edges[i].row), &(w_edges[i].col), &(w_edges[i].val));
      assert (fscanfret == 3);
      assert ( 0<w_edges[i].row && w_edges[i].row <=nr);
      assert ( 0<w_edges[i].col && w_edges[i].col <=nc);
      w_edges[i].row -= 1;    // MatrixMarket format is 1-up, not 0-up
      w_edges[i].col -= 1;
    }
    qsort( w_edges, nnz, sizeof(w_edge_t), w_edge_comp);
  }
  fclose(fp);

  ret_mat = init_matrix( nr, nc, nnz, (values > 0));
  assert((ret_mat != NULL) && "Memory Error");

  int64_t pos = 0;
  int64_t row = 0;
  ret_mat->offset[row] = 0;
  if (values == 0){
    while (pos<nnz) {
      if (edges[pos].row == row) {
        ret_mat->nonzero[pos] = edges[pos].col;
        pos++;
        continue;
      }
      while (row < edges[pos].row) {
        row++;
        ret_mat->offset[row] = pos;
      }
    }
    free(edges);
  } else {
    assert(values == 1);
    while (pos<nnz) {
      if (w_edges[pos].row == row) {
        ret_mat->nonzero[pos] = w_edges[pos].col;
        ret_mat->value[pos]   = w_edges[pos].val;
        pos++;
        continue;
      }
      while (row < w_edges[pos].row) {
        row++;
        ret_mat->offset[row] = pos;
      }
    }
    free(w_edges);
  }
  ret_mat->offset[row+1] = pos;

  return(ret_mat);
}



/*! \brief apply row and column permutations to a sparse matrix 
 * \param A pointer to the original matrix
 * \param rperminv pointer to the global array holding the inverse of the row permutation
 * \param cperminv pointer to the global array holding the inverse of the column permutation
 * rperminv[i] = j means that row i of A goes to row j in matrix Ap
 * cperminv[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced 
 */
sparsemat_t * permute_matrix(sparsemat_t *A, int64_t *rperminv, int64_t *cperminv) 
{
  int64_t i, j, row, pos;

  sparsemat_t * Ap = init_matrix(A->numrows, A->numcols, A->nnz, (A->value!=NULL));
  if(!Ap){printf("ERROR: permute_matrix: init_matrix failed!\n");return(NULL);}

  int64_t * rperm = calloc(A->numrows, sizeof(int64_t));
  if( rperm == NULL ) return(NULL);
  //compute rperm from rperminv 
  for(i=0; i < A->numrows; i++)
    rperm[rperminv[i]] = i;
  
  // fill in permuted rows with permuted columns
  Ap->offset[0] = pos = 0;
  for(i = 0; i < Ap->numrows; i++){
    row = rperm[i];
    for(j = A->offset[row]; j < A->offset[row+1]; j++){
      Ap->nonzero[pos] = cperminv[A->nonzero[j]];
      if(A->value) Ap->value[pos] = A->value[j];
      pos++;      
    }
    Ap->offset[i+1] = pos;
  }
  
  free(rperm);
  return(Ap);
}

/*! 
\brief produce the transpose of a sparse matrix
\param A  pointer to the original matrix
\return a pointer to the matrix that has been produced or NULL on error.
 
Note: Given a tidy matrix, this will produce a tidy matrix because we
process the rows of the given matrix in order. Thus, we fill the 
columns of the transpose in order.
*/
sparsemat_t * transpose_matrix(sparsemat_t *A) 
{
  int64_t i, j, row, col;
  sparsemat_t * At;
  
  // get column counts
  int64_t * colcnt = calloc( A->numcols, sizeof(int64_t));
  if( colcnt == NULL ) return(NULL);
  int64_t * tmpoffset = colcnt;
  
  // histogram the column counts of A into colcnt
  for( i=0; i< A->nnz; i++){  
    assert( A->nonzero[i] < A->numcols );
    assert( A->nonzero[i] >= 0 ); 
    colcnt[A->nonzero[i]]++;
  }

  At = init_matrix(A->numcols, A->numrows, A->nnz, (A->value!=NULL));
  if(!At){printf("ERROR: transpose_matrix: init_matrix failed!\n");return(NULL);}

  //use the colcnt array to build the offset array and
  //reuse it as we fill in the nonzeros 
  At->offset[0] = 0;
  for(i = 0; i < At->numrows; i++){
    At->offset[i+1] = At->offset[i] + colcnt[i];
    tmpoffset[i] = At->offset[i];
  }

  //redistribute the nonzeros 
  for(row=0; row<A->numrows; row++){
    for(j=A->offset[row]; j<A->offset[row+1]; j++){
      col = A->nonzero[j];
      At->nonzero[ tmpoffset[col] ] = row;
      if(A->value)  // expect (or hope) that the compiler would pull this out of the loop
				At->value[ tmpoffset[col] ] = A->value[j];
      tmpoffset[col]++;
    }
  }

  free(colcnt);
  return(At);
}


/*!
 * \brief builds a full symmetric matrix from a lower triangular matrix
 * \param *L sparsemat_t that holds a lower triangular matrix (no diagonal elements) 
 * \return the full matrix
 */
sparsemat_t * make_symmetric_from_lower(sparsemat_t * L)
{
	int64_t i, l, k, pos;
	sparsemat_t *U  = transpose_matrix(L);
  int hasvalues = (L->value == NULL) ? 0 : 1;
	sparsemat_t *retmat = init_matrix(L->numrows, L->numrows, 2*(L->nnz), hasvalues);

  pos = 0;
  retmat->offset[0] = 0;
  for(i=0; i<L->numrows; i++){
    for(l=L->offset[i]; l<L->offset[i+1]; l++){
      retmat->nonzero[pos] = L->nonzero[l];
      pos++;
    }
    for(k=U->offset[i]; k<U->offset[i+1]; k++){
      retmat->nonzero[pos] = U->nonzero[k];
      pos++;
    }
    retmat->offset[i+1] = pos;
  }

  if(hasvalues) {
    pos = 0;
    for(i=0; i<L->numrows; i++){
      for(l=L->offset[i]; l<L->offset[i+1]; l++){
        retmat->value[pos] = L->value[l];
        pos++;
      }
      for(k=U->offset[i]; k<U->offset[i+1]; k++){
        retmat->value[pos] = U->value[k];
        pos++;
      }
    }
  } 

	clear_matrix(U);
	free(U);
  return(retmat);	
}


/*! \brief checks that the sparse matrix is triangular with or without a diagonal
 * \param A pointer to the sparse matrix
 * \param upper is it upper=1 or lower upper=0
 * \param ondiag whether or not to include diagonal: 0=no-diagonal, 1=all-diagonal, else don't care
 * \return 0 on success, non-0 on error.
 * Only called by is_(upper|lower)_triangular
 */
static int64_t is_triangular(sparsemat_t *A, int64_t upper, int64_t ondiag) 
{
  int64_t row, k;
  int64_t err = 0;
  int64_t diag_cnt = 0;

  for (row=0; row < A->numrows; row++){
    for (k=A->offset[row]; k < A->offset[row+1];k++) {
      if (A->nonzero[k] == row) {
        diag_cnt++;
      }
      err += ((upper==1) && (A->nonzero[k] < row)) ? 1 : 0;
      err += ((upper==0) && (A->nonzero[k] > row)) ? 1 : 0;
    }
  }
  if(err > 0)
    return(0);
  if ((ondiag == 0) && (diag_cnt !=0))
    return(0);
  if ((ondiag == 1) && (diag_cnt != A->numrows))
    return(0);
  return(1);
}

/*! \brief checks that the sparse matrix is (strictly, i.e. zero on diagonal) lower triangluar
 * \param A pointer to the sparse matrix
 * \param ondiag whether or not to include diagonal: 0=no-diagonal, 1=all-diagonal, else don't care
 * \return 0 on success, non-0 on error.
 */
int64_t is_lower_triangular(sparsemat_t *A, int64_t ondiag) 
{
  return(is_triangular(A, 0, ondiag));
}

/*! \brief checks that the sparse matrix is (strictly) upper triangluar
 * \param A pointer to the sparse matrix
 * \param ondiag whether or not to include diagonal: 0=no-diagonal, 1=all-diagonal, else don't care
 * \return 0 on success, non-0 on error.
 */
int64_t is_upper_triangular(sparsemat_t *A, int64_t ondiag) 
{
  return(is_triangular(A, 1, ondiag));
}

/*! \brief Sets all entries above the kth diagonal to zero.
 * \param A A pointer to a sparse matrix
 * \param kth the kth diagonal (kth > 0 is above the main diagonal kth < 0 is below)
 * \return The number of nonzeros that were removed.
 * NB: this packs the nonzero and values into the front of the arrays;
 * we don't reallocate.
 * \ingroup spmatgrp
 */
int64_t tril(sparsemat_t * A, int64_t kth)
{
  int64_t i, k, col, pos, start;
  int weighted = (A->value != NULL);
  pos = 0; 
  start = 0;
  for (i = 0; i < A->numrows; i++) {
    for (k = start; k < A->offset[i+1]; k++) {
      col = A->nonzero[k];
      if ((col - i) <= kth) {
        A->nonzero[pos] = col;
        if (weighted) A->value[pos] = A->value[k];
        pos++;
      }
    }
    start = A->offset[i+1];
    A->offset[i+1] = pos;
  }
  int64_t ret = A->nnz - pos;
  A->nnz = pos;
  return(ret);
}

/*! \brief Sets all entries below the kth diagonal to zero.
 * \param A A pointer to a sparse matrix
 * \param kth the kth diagonal will be zero (kth > 0 is above the main diagonal kth < 0 is below)
 * \return The number of nonzeros that were removed.
 * NB: this packs the nonzero and values into the front of the arrays;
 * we don't reallocate.
 * \ingroup spmatgrp
 */
int64_t triu(sparsemat_t * A, int64_t kth) 
{
  int64_t i, k, col, pos, start;
  int weighted = (A->value != NULL);
  pos = 0; 
  start = 0;
  for (i = 0; i < A->numrows; i++) {
    for (k = start; k < A->offset[i+1]; k++) {
      col = A->nonzero[k];
      if ((col - i) >= kth) {
        A->nonzero[pos] = col;
        if(weighted) A->value[pos] = A->value[k];
        pos++;
      }
    }
    start = A->offset[i+1];
    A->offset[i+1] = pos;
  }
  int64_t ret = A->nnz - pos;
  A->nnz = pos;
  return(ret);
}

/*! \brief comparison function for qsort in sort_nonzeros
 */
int nz_comp(const void *a, const void *b) 
{
  return( *(uint64_t *)a - *(uint64_t *)b );
}

/*! \brief comparison function for qsort in sort_nonzeros in a matrix with values
 */
int cv_comp(const void *a, const void *b)
 {
   return(((int64_t)((col_val_t*)a)->col) - ((int64_t)((col_val_t*)b)->col));
 }
 
/*! \brief sort the  non-zeros in each row of a sparse matrix so that their column indices
 * are in ascending order.
 * \param mat pointer to the sparse matrix
 */
int64_t sort_nonzeros( sparsemat_t *mat)
{
  int64_t i, j;
  if(mat->value){
    // we have to sort the column indicies, but we also have to permute the value array accordingly
    // this is annoying in C
    // we have to create an array of stucts that holds col,val pairs for a row
    // sort that array according to the col keys
    // and then overwrite the row data
    int64_t max = 0;
    for(i = 0; i < mat->numrows; i++)
      if(mat->offset[i+1] - mat->offset[i] > max) max = mat->offset[i+1] - mat->offset[i];

    // allocate a temporary array to hold a row's worth of col, value pairs
    col_val_t * tmparr = calloc(max, sizeof(col_val_t));

    for(i = 0; i < mat->numrows; i++){
      int64_t pos = 0;
      for(j = mat->offset[i]; j < mat->offset[i+1]; j++){
        tmparr[pos].col = mat->nonzero[j];
        tmparr[pos++].value = mat->value[j];
      }
      qsort(tmparr, mat->offset[i+1] - mat->offset[i], sizeof(col_val_t), cv_comp );
      pos = 0;
      for(j = mat->offset[i]; j < mat->offset[i+1]; j++){
        mat->nonzero[j] = tmparr[pos].col;
        mat->value[j] = tmparr[pos++].value;
      }
    }
    free(tmparr);
  }else{
    // no values, just sort the column indicies
    for(i = 0; i < mat->numrows; i++){
      qsort( &(mat->nonzero[mat->offset[i]]), mat->offset[i+1] - mat->offset[i], sizeof(int64_t), nz_comp );
    }
  }
  return(0);
}

/*! \brief compare the structs that hold two sparse matrices
 * \param lmat pointer to the left sparse matrix
 * \param rmat pointer to the right sparse matrix
 * \return 0 on success
 */
int64_t compare_matrix(sparsemat_t *lmat, sparsemat_t *rmat) 
{
  int64_t i,j;
  int values;

  if( lmat->value == NULL && rmat->value == NULL){
    values = 0;
  }else if((lmat->value && !rmat->value) || (rmat->value && !lmat->value)){
    printf("Only one matrix has values!\n");
    return(1);    
  }else{
    values = 1;
  }
  
  if( lmat->numrows != rmat->numrows ){
    printf("(lmat->numrows = %"PRId64")  != (rmat->numrows = %"PRId64")", lmat->numrows, rmat->numrows );
    return(1);
  }
  if( lmat->numrows != rmat->numrows ){
    printf("(lmat->numrows = %"PRId64")  != (rmat->numrows = %"PRId64")", lmat->numrows, rmat->numrows );
    return(1);
  }
  if( lmat->numcols != rmat->numcols ){
    printf("(lmat->numcols = %"PRId64")  != (rmat->numcols = %"PRId64")", lmat->numcols, rmat->numcols );
    return(1);
  }
  if( lmat->nnz != rmat->nnz ){
    printf("(lmat->nnz = %"PRId64")  != (rmat->nnz = %"PRId64")", lmat->nnz, rmat->nnz );
    return(1);
  }
  if( lmat->nnz != rmat->nnz ){
    printf("(lmat->lnnz = %"PRId64") != (rmat->lnnz = %"PRId64")", lmat->nnz, rmat->nnz );
    return(1);
  }

  if( lmat->offset[0] != 0 || rmat->offset[0] != 0 || (lmat->offset[0] != rmat->offset[0] ) ){
    printf("(lmat->offset[0] = %"PRId64")  != (rmat->offset[0] = %"PRId64")", lmat->offset[0], rmat->offset[0] );
    return(1);
  }
  
  for(i = 0; i < lmat->numrows; i++){
    if( lmat->offset[i+1] != rmat->offset[i+1] ){
      printf("(lmat->offset[%"PRId64"] = %"PRId64")  != (rmat->offset[%"PRId64"] = %"PRId64")", i+1, lmat->offset[i+1], i+1, rmat->offset[i+1] );
      return(1);
    }
  }
  
  for(j=0; j< lmat->nnz; j++){
    if( lmat->nonzero[j] != rmat->nonzero[j] ){
      printf("(lmat->nonzero[%"PRId64"] = %"PRId64")  != (rmat->nonzero[%"PRId64"] = %"PRId64")", j,
             lmat->nonzero[j], j, rmat->nonzero[j] );      
      return(1);
    }
    if(values){
      if( lmat->value[j] != rmat->value[j] ){
        printf("(lmat->value[%"PRId64"] = %lf)  != (rmat->value[%"PRId64"] = %lf)", j,
               lmat->value[j], j, rmat->value[j] );
        return(1);
      }
    }
  }

  return(0);
}

/*! \brief makes an exact copy of a given sparse matrices
 * \param srcmat pointer to the original sparse matrix
 * \return A pointer to the cloned sparse matrix
 */
sparsemat_t * copy_matrix(sparsemat_t *srcmat) 
{
  int64_t i,j;
  
  sparsemat_t * destmat = init_matrix(srcmat->numrows, srcmat->numcols, srcmat->nnz, (srcmat->value != NULL));
  if(!destmat) return(NULL);

  for(i = 0; i < (srcmat->numrows)+1; i++){
    destmat->offset[i] = srcmat->offset[i];
  }
  for(j=0; j < srcmat->nnz; j++){
    destmat->nonzero[j] = srcmat->nonzero[j];
    if(srcmat->value) destmat->value[j] = srcmat->value[j];
  }
  return(destmat);
}

/*! \brief initializes the struct that holds a sparse matrix
 *    given the total number of rows and columns and the local number of non-zeros
 * \param numrows total number of rows
 * \param numcols total number of columns
 * \param nnz number of nonzero
 * \param values flag: 0 = no values, 1 = doubles
 * \return An initialized sparsemat_t (numrows, numcols, nnz are set) or NULL on error.
 */
sparsemat_t * init_matrix(int64_t numrows, int64_t numcols, int64_t nnz, int values) 
{
  sparsemat_t * mat = calloc(1, sizeof(sparsemat_t));
  mat->numrows  = numrows;
  mat->numcols  = numcols;
  mat->nnz  = nnz;  
  mat->offset   = calloc(mat->numrows+1, sizeof(int64_t));
  if(mat->offset == NULL) 
    return(NULL);
  mat->nonzero = calloc(mat->nnz, sizeof(int64_t));
  if(mat->nonzero == NULL) 
    return(NULL);
  if(values){
    mat->value = calloc(mat->nnz, sizeof(double));
    if(mat->value == NULL)
      return(NULL);
  }else
    mat->value = NULL;
  return(mat);
}

 /*************************************************************************************/
 /*                               RANDOM MATRICES                                     */
 /*************************************************************************************/

 /* 
  * What kinds of methods should we have to generate "random" sparse matrices?
  * 1) erdos-renyi - square upper or lower-triangular matrices 
  *    can be used to create symmetric or nonsymmetric square matrices
  * 2) uniform sparse: use ER technology
  * 3) kron_prod_graph (square symmetric or lower-triangular)
  * 4) random geometric? random planar graphs? Social network (preferential attachment, etc)
  * geometric graph - k-nearest neighbor graph, or epsilon ball graph
  *   - has different properties than ER
  *
  */
 // apps that need input matrices:
 // topo: square, upper-triangular, unit-diagonal, random or read in
 // permute_matrix and transpose: any random matrix, or read in matrix
 // triangle: kron_graph (special lower triangular), or any random lower triangular, or read in
 // SSSP: random non-symmetric square with values

/*! \brief DPRT whether or not to do debug printing */
#define DPRT 1

/*!
\brief routine to generate the adjacency matrix of a random graph.
\param n the number of vertices in the graph.
\param model FLAT: Erdos-Renyi random, GEOMETRIC: geometric random graph
\param edge_type See edge_type enum. directed, or not, weighted or not.
\param loops see self_loops enum. Does every node have a loop or not.
\param edge_density: d in [0, 1), target fraction of edges present.
\param seed: RNG seed.

If the graph is undirected, this routine only returns a lower-triangular
adjancency matrix (since the adjancency matrix would be symmetric and we don't need
the redundant entries).
*/
sparsemat_t * random_graph(int64_t n, graph_model model, edge_type edge_type, self_loops loops,
                            double edge_density, int64_t seed)
{
   if(DPRT){
     printf("making a random graph:\n");
     printf("n = %ld, prob = %lf, seed %ld\n", n, edge_density, seed);
   }
   if(model == FLAT){
     return(erdos_renyi_random_graph(n, edge_density, edge_type, loops, seed));
   }else if(model == GEOMETRIC){
     if(DPRT){printf("making a geometric graph\n");}
     double r;
     // determine the r that will lead to the desired edge density
     // Expected degree = n*pi*r^2
     // The expected number of edges E = n^2*pi*r^2/2
     // We solve for r by setting edge_prob = E / total possible edges
     // for undirected edge_prob = E/(n choose 2)
     // for directed   edge_prob = E/(n^2 - n)
     if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
       r = sqrt((n-1)*edge_density/(M_PI*n));
     else
       r = sqrt(2*(n-1)*edge_density/(M_PI*n));
     return(geometric_random_graph(n, r, edge_type, loops, seed));
     
   }else{
     printf("ERROR: random_graph: Unknown type!\n");
     return(NULL);
   }
   
}


/****************** Kronecker Products of Star  Graphs *******************************/
/*! \brief Generate the adjacency matrix for the star K_{1,m} (with or without loop edges)
 * \param m  the number of non-hub vertices
 * \param mode
 *  - mode == 0: default, no self loops
 *  - mode == 1: add a self loop to center vertex of star (row zero)
 *  - mode == 2: add a self loop to an outer vertex (the last one)
 * \return the adjacency matrix of the graph
*/
sparsemat_t * gen_star_graph(int64_t m, int mode) 
{
  sparsemat_t * G = init_matrix(m + 1, m + 1, 2*m + (mode > 0 ? 1 : 0), 0);
  
  int64_t i,j;
  int64_t pos = 0;  // counter thru the nonzeros

  G->offset[0] = 0;
  if(mode == 1)     // add a self loop to center vertex 
    G->nonzero[pos++] = 0;

  // Add the top row 
  for(j = 1; j < m+1; j++){
    G->nonzero[pos++] = j;
  }  
  G->offset[1] = pos;

  // rest of the rows 
  for(i = 1; i < m+1; i++){
    G->nonzero[pos++] = 0;
    G->offset[i+1] = pos;
  }
  
  if(mode == 2){  // add a self loop to last vertex and fix the last offset
    G->nonzero[pos++] = m;
    G->offset[m+1] = pos;
  }
  
  return(G);
}

/*! \brief Produce the Kronecker matrix product of two square {0,1}-matrices
 * \param A sparse matrix, A is A_n x A_n
 * \param B sparse matrix, B is B_n x B_n
 * \return \f$C =  A \otimes B\f$
 * 
 * Key fact is that element wise, using the notation C[row,col], A[r,s] and B[u,v],
 * 
 *  C[row,col] = C[r*B_n + u, s*B_n + v] = A[r,s] * B[u,v]
 *
*/
sparsemat_t * kronecker_mat_product(sparsemat_t * A, sparsemat_t * B) 
{
  int64_t r, s, u, v, j, k;
  int64_t A_n = A->numrows;
  int64_t B_n = B->numrows;
  int64_t C_n = A_n * B_n;
  int64_t C_nnz = A->nnz * B->nnz;

  sparsemat_t * C = init_matrix(C_n, C_n, C_nnz, 0);  // values = 0, for a {0,1}-matrix

  // get the number of nonzeros in each row 
  int64_t * C_rowtmp = calloc(C_n + 1, sizeof(int64_t));
  for(r = 0; r < A_n; r++){
    int64_t da = A->offset[r + 1] - A->offset[r];
    for(s = 0; s < B_n; s++){
      int64_t db = B->offset[s + 1] - B->offset[s]; 
      C_rowtmp[r*B_n + s] = da * db;
    }
  }

  // set the offsets in C 
  // we will reuse C_rowtmp[] to count from C->offset[r] to C->offset[r+1] 
  // as we fill in the nonzeros in row r of C
  C->offset[0] = 0;
  for(r = 0; r < C_n; r++){
    C->offset[r+1] = C->offset[r] + C_rowtmp[r];
    C_rowtmp[r] = C->offset[r];
  }
  

  // fill in the nonzeros of C, we fill them in a "block" order.
  // ie, we write all of B for each nonzero in A 
  for(r = 0; r < A_n; r++){
    for(j = A->offset[r]; j < A->offset[r+1]; j++){
      s = A->nonzero[j];
      for(u = 0; u < B_n; u++){
        for(k = B->offset[u]; k < B->offset[u+1]; k++){
          v = B->nonzero[k];
          C->nonzero[C_rowtmp[r*B_n + u]++] = s * B_n + v;
        }
      }
    }
  }
  
  free(C_rowtmp);

  return(C);
}


/*! \brief Generate the kroncker product of a collection of star graphs.
 * \param mode 
 *  - mode == 0: default, no self loops
 *  - mode == 1: add a self loop to center vertex of each star 
 *  - mode == 2: add a self loop to an outer vertex (the last vertex) of each star
 * \param spec an array the holds the sizes of the stars
 * \param num number of stars
 * \param weighted whether or not to apply weights to the edges
 * \return the adjacency matrix for graph
*/
sparsemat_t * generate_kronecker_graph_from_spec(int mode, int * spec, int num, int weighted)
{
  int64_t i, k;
  int64_t G_n, G_nnz;
  sparsemat_t * star;

  if (num < 2) {
    fprintf(stderr,"ERROR: generate_kronecker_graph_from_spec requires at least two stars\n");
    return(NULL);
  }

  sparsemat_t ** mats = calloc(2*num, sizeof(sparsemat_t *));

  mats[0] = gen_star_graph(spec[0], mode);

  for (i = 1; i < num; i++) {
    star = gen_star_graph(spec[i], mode);
    mats[i] = kronecker_mat_product(mats[i-1] , star);
    clear_matrix(star);
  }
  sparsemat_t * R = mats[num-1];  // what we want but is fully symmetric and may have loops
  
  // count the number of nonzeros in the lower triangle part of the matrix
  G_n = R->numrows;
  G_nnz = 0;
  for (i=0; i<R->numrows; i++) {
    for (k=R->offset[i]; (k<R->offset[i+1]) && (R->nonzero[k] < i ); k++) {
      G_nnz++;
    }
  }

  // copy the lower triangle of the matrix (excluding the diagonal) to the returned matrix
  sparsemat_t * G = init_matrix(G_n, G_n, G_nnz, weighted);
  assert( (G != NULL) && "Memory Error");
  G->offset[0] = 0;
  int64_t pos = 0;
  for (i=0; i<R->numrows; i++) {
    for (k=R->offset[i]; (k<R->offset[i+1]) && (R->nonzero[k] < i ); k++) {
      G->nonzero[pos++] = R->nonzero[k];
    }
    G->offset[i+1] = pos;
  }
   // fill in the random values
  if(G->value){
    for(i = 0; i < G->nnz; i++){
      G->value[i] = rand_double();
    }
  }
  for (i = 1; i < num; i++) 
    free(mats[i]);
  free(mats);

  return(G);
}

/*! \brief compute the number of triangles in the the Kroncker Product graph
 * \param mode 
 *  - mode == 0: default, no self loops
 *  - mode == 1: add a self loop to center vertex of each star 
 *  - mode == 2: add a self loop to an outer vertex (the last vertex) of each star
 * \param spec an array the holds the sizes of the stars
 * \param num number of stars
 * \return the number of triangles
*/
int64_t calc_num_tri_kron_graph(int mode, int * spec, int num)
{
   double approx, ns, nr;
   int i;

   if ( mode == 0 ) {
     return 0;
   } else if ( mode == 1 ) {
     approx = 1.0;
     nr = 1.0;
     for (i = 0; i < num; i++) {
       approx *= (3*spec[i] + 1);
       nr *= spec[i] + 1.0;
     }
     return ( round((approx / 6.0) - 0.5 * nr + 1.0/3.0) );
   } //else if ( mode == 2 ) 
   assert(mode==2 && "Parameter error");
   ns = (double)num;
   return( round( (1.0/6.0)*pow(4,ns) - pow(2.0,(ns - 1.0)) + 1.0/3.0) );
}

// IF we are given nz_per_row (z), we calculate edge_prob (e)
// or if we are given edge_prob, we calculate nz_per_row
// using with the following formulas:
// WITH LOOPS:
// z*n = e*(n*(n-1)/2) (for UNDIRECTED*) or z*n = e*(n(n-1)) for DIRECTED
// WITHOUT LOOPS:
// (z-1)*n = e*(n*(n-1)/2) (for UNDIRECTED*) or (z-1)*n = e*(n(n-1)) for DIRECTED
//

/*! 
 * \brief computes the number of non-zeros per row from the edge probability or vice-versa.
 * \param edge_prob given edge probability (or place to put the computed edge_prob)
 * \param nz_per_row  given number of nonzeros per row (or place to put the computed nz_per_row)
 * \param numrows = numcols global order of the matrix 
 * \param edge_type weighted or unweighted
 * \param loops whether the diagonal is all zeros or all ones
 */
void resolve_edge_prob_and_nz_per_row(double * edge_prob, double * nz_per_row, int64_t numrows,
                                      edge_type edge_type, self_loops loops){
  if(*edge_prob == 0.0){ // use nz_per_row to get erdos_renyi_prob
    if(loops == LOOPS)
      *edge_prob = (*nz_per_row - 1)/(numrows - 1);
    else
      *edge_prob = (*nz_per_row)/(numrows-1);    

    if (edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      *edge_prob = *edge_prob*2;
      
    if(*edge_prob > 1.0)
      *edge_prob = 1.0;
  } else {    // use erdos_renyi_prob to get nz_per_row
    *nz_per_row = *edge_prob * ((numrows - 1)/2.0);
    if (edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      *nz_per_row = *nz_per_row * 2;
  }
  assert(*edge_prob <= 1.0);
}

//****************************** Geometric Graphs *********************************************/
/*!
\brief struct to hold a point in the unit square and index (name) for the point
*/
typedef struct points_t{
   double x;           //!< x coordinate
   double y;           //!< y coordinate
   int64_t index;      //!< index for the point
 }points_t;
 
/*!
\brief array of points in a sector
*/
typedef struct sector_t{
  points_t * points;     //!< array of points
  int64_t numpoints;     //!< number of points in the sector
}sector_t;

/*! 
\brief compute the square of the distance between two points
\param a the first point
\param b the other point
\return the square 
*/
double dist(points_t a, points_t b){
  return((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

/*!
\brief initialize an edge list
\param nalloc the number of edges the list will initially hold
*/
edge_list_t * init_edge_list(int64_t nalloc){
  edge_list_t * el = calloc(1, sizeof(edge_list_t));
  el->num = 0;
  el->nalloc = nalloc;
  el->edges = calloc(nalloc, sizeof(edge_t));
  return(el);
}


/*!
\brief add an edge to an edge list
\param el the list
\param row the row (tail) of the edge
\param col the column (head) of the edge

NB: note additional space is allocated when needed
*/
int64_t append_edge(edge_list_t * el, int64_t row, int64_t col){
  if(el->nalloc == el->num){
    //printf("out of space! nalloc = %ld\n", el->nalloc);
    // we need to expand our allocations!
    if(el->nalloc < 10000)
      el->nalloc = 2*el->nalloc;
    else
      el->nalloc = el->nalloc*1.25;
    //printf("new space! nalloc = %ld\n", el->nalloc);
    el->edges = realloc(el->edges, el->nalloc*sizeof(edge_t));
  }
  //printf("appending %ld %ld\n", row, col);
  el->edges[el->num].row = row;
  el->edges[el->num].col = col;
  el->num++;
  return(el->num);
}

/*!
\brief clear the memory used by the edge list
\param el the list
*/
void clear_edge_list(edge_list_t * el){
  free(el->edges);
  free(el);
}

  
//****************************** Geometric Graphs *********************************************/

/*!
\brief Generates the adjacency matrix for a random geometric graph. 
\param n the number of vertices
\param r the size of the neighborhoods that determine edges.
\param edge_type if undirected, this routine returns a lower triangular matrix
\param loops whether or not the diagonal is all 1's or all 0's
\param seed A seed for the RNG.
\return An adjacency matrix (or lower portion of in the undirected case).

See https://en.wikipedia.org/wiki/Random_geometric_graph
Each vertex corresponds to a point randomly placed in the unit square. Two vertices
are adjancent if their corresponding points are within distance r of each other.
*/
sparsemat_t * geometric_random_graph(int64_t n, double r, edge_type edge_type, self_loops loops, uint64_t seed){
  int64_t i, j, k, l;
  double r2 = r*r;
  // We break up the unit square into chunks that are rxr grid.
  // We generate the points uniformly at random over the unit square.
  // We calculate edges by comparing distances between each point and every other point in
  // its own sector and in neighboring sectors.
  int weighted = (edge_type == UNDIRECTED_WEIGHTED || edge_type == DIRECTED_WEIGHTED);
  int directed = (edge_type == DIRECTED || edge_type == DIRECTED_WEIGHTED);
  int64_t nsectors_across = floor(1.0/r);
  double sector_width = 1.0/nsectors_across;
  int64_t nsectors = nsectors_across*nsectors_across;

  printf("GEOMETRIC GRAPH: r = %lf number of sectors = %"PRId64" sector_width = %lf\n", r, nsectors,sector_width);
  printf("                 edge_type = %d ", edge_type);
  if(loops == NOLOOPS)
    printf("no loops\n");
  else
    printf("with loops\n");
  
  sector_t ** sectors = calloc(nsectors_across, sizeof(sector_t*));
  int64_t ** first_index_this_sector = calloc(nsectors_across, sizeof(int64_t*));
  for(i = 0; i < nsectors_across; i++){
    sectors[i] = calloc(nsectors_across, sizeof(sector_t));
    first_index_this_sector[i] = calloc(nsectors_across, sizeof(int64_t));
    for(j = 0; j < nsectors_across; j++)
      sectors[i][j].numpoints = 0;
  }


  // First pass, generate the points and count how many will fall in each sector.
  rand_seed(seed); 
  for(i = 0; i < n; i++){
    //double x = (double)rand()/RAND_MAX;
    //double y = (double)rand()/RAND_MAX;
    double x = rand_double();
    double y = rand_double();
    int64_t row = floor(y/sector_width);
    int64_t col = floor(x/sector_width);
    assert(row < nsectors_across);
    assert(col < nsectors_across);
    sectors[row][col].numpoints++;
  }

  // initialize the struct to hold the points
  for(i = 0; i < nsectors_across; i++){
    for(j = 0; j < nsectors_across; j++){
      if(j > 0)
        first_index_this_sector[i][j] += first_index_this_sector[i][j-1] + sectors[i][j-1].numpoints;
      else if(i > 0)
        first_index_this_sector[i][j] += first_index_this_sector[i-1][nsectors_across-1] + sectors[i-1][nsectors_across-1].numpoints;
      sectors[i][j].points = calloc(sectors[i][j].numpoints, sizeof(points_t));
    }
  }

  // reset numpoints
  for(i = 0; i < nsectors_across; i++){
    for(j = 0; j < nsectors_across; j++){
      sectors[i][j].numpoints = 0;
    }
  }
  // Second pass: generate the points and insert them into the struct
  rand_seed(seed);
  for(i = 0; i < n; i++){
    double x = rand_double();
    double y = rand_double();
    int64_t row = floor(y/sector_width);
    int64_t col = floor(x/sector_width);
    int64_t li = sectors[row][col].numpoints;
    sectors[row][col].points[li].x = x;
    sectors[row][col].points[li].y = y;
    sectors[row][col].points[li].index = first_index_this_sector[row][col] + li;
    sectors[row][col].numpoints++;
  }
  
  // next we calculate the edges that appear and put them in a list
  int64_t this_node, other_node;
  int64_t space = ceil(1.1*n*n*M_PI*r*r/2.0);
  //printf("Initial allocation: %ld\n", space);
  edge_list_t * el = init_edge_list(space);
  
  for(i = 0; i < nsectors_across; i++){
    for(j = 0; j < nsectors_across; j++){
      
      sector_t * sec = &sectors[i][j];
      int64_t m = sec->numpoints;
      for(k = 0; k < m; k++){
        this_node = sec->points[k].index;

        if(loops == LOOPS){
          append_edge(el, this_node, this_node);
        }
        // count the edges to lower-indexed nodes within this sector
        for(l = 0; l < k; l++){
          other_node = sec->points[l].index;
          if(dist(sec->points[k], sec->points[l]) < r2){
            if(directed && rand_int64(2))
              append_edge(el, other_node, this_node);
            else
              append_edge(el, this_node, other_node);
          }
        }

        // count the edges to lower-indexed nodes outside the sector
        // to do this, we need to look at sectors to the W, NW, N, and NE.
        // W
        if(j > 0){
          sector_t * sec2 = &sectors[i][j-1];
          for(l = 0; l < sec2->numpoints; l++){
            other_node = sec2->points[l].index;
            if(dist(sec->points[k], sec2->points[l]) < r2){
              if(directed && rand_int64(2))
                append_edge(el, other_node, this_node);
              else
                append_edge(el, this_node, other_node);
            }
          }
        } 
        // NW
        if(i > 0 && j > 0){
          sector_t * sec2 = &sectors[i-1][j-1];
          for(l = 0; l < sec2->numpoints; l++){
            other_node = sec2->points[l].index;
            if(dist(sec->points[k], sec2->points[l]) < r2){
              if(directed && rand_int64(2))
                append_edge(el, other_node, this_node);
              else
                append_edge(el, this_node, other_node);
            }
          }
        } 
        // N
        if(i > 0){
          sector_t * sec2 = &sectors[i-1][j];
          for(l = 0; l < sec2->numpoints; l++){
            other_node = sec2->points[l].index;
            if(dist(sec->points[k], sec2->points[l]) < r2){
              if(directed && rand_int64(2))
                append_edge(el, other_node, this_node);
              else
                append_edge(el, this_node, other_node);
            }
          }
        }
        // NE
        if(i > 0 && j < (nsectors_across-1)){
          sector_t * sec2 = &sectors[i-1][j+1];
          for(l = 0; l < sec2->numpoints; l++){
            other_node = sec2->points[l].index;
            if(dist(sec->points[k], sec2->points[l]) < r2){
              if(directed && rand_int64(2))
                append_edge(el, other_node, this_node);
              else
                append_edge(el, this_node, other_node);
            }
          }
        } 
      }
    }
  }
  
  //printf("nedges = %"PRId64"\n", el->num);
  
  // Free up points and sectors
  for(i = 0; i < nsectors_across; i++){
    for(j = 0; j < nsectors_across; j++)
      free(sectors[i][j].points);
    free(first_index_this_sector[i]);
    free(sectors[i]);
  }
  free(sectors);
  free(first_index_this_sector);
  
  sparsemat_t * A = init_matrix(n, n, el->num, weighted);

  // sort the edges in row / col order
  qsort(el->edges, el->num, sizeof(edge_t), edge_comp);

  int64_t row = 0;
  for(i = 0; i < el->num; i++){
    while(row != el->edges[i].row){
      row++;
      assert(row < n);
      A->offset[row] = i;
    }
    A->nonzero[i] = el->edges[i].col;
    if(weighted) A->value[i] = rand_double();
  }
  
  while(row < n){
    row++;
    A->offset[row] = el->num;
  }

  clear_edge_list(el);
  
  return(A);
  
 }

 
//****************************** Erdos Renyi Graphs *******************************************/

/*!
\brief generate the adjacency matrix for an Erdos-Renyi random graph. 
\param n The total number of vertices in the graph (rows in the matrix).
\param p The probability that each non-loop edge is present.
\param edge_type See edge_type. DIRECTED, UNDIRECTED, DIRECTED_WEIGHTED, UNDIRECTED_WEIGHTED
\param loops See self_loops.
\param seed A random seed.
\return A sparsemat_t

 * This subroutine uses ALG1 from the paper "Efficient Generation of Large Random Networks" 
 * by Batageli and Brandes appearing in Physical Review 2005. Instead of flipping a coin for each potential edge
 * this algorithm generates a sequence of "gaps" between 1s in the upper or lower triangular portion of the 
 * adjancecy matrix using a geometric random variable.
*/
sparsemat_t * erdos_renyi_random_graph(int64_t n, double p, edge_type edge_type, self_loops loops, int64_t seed)
{
  int64_t row, col;
  double D  = log(1 - p);
  int64_t nnz;
  int64_t end = n;
  int64_t ndiag = n;

  if(0){printf("making a erdos_renyi random graph:\n");}
  if(0){printf("n = %ld, p=%f, edges = %d, loops = %d, seed %ld\n", n, p, edge_type, loops, seed);}
  assert(n > 0);

  // first loop to count the number of nonzeros
  rand_seed(seed);
  row = 0;
  nnz = 0;
  col = 1 + floor(log(1 - rand_double()) / D);
  while (row < n) {
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      if(col == row) // we hit a diagonal for free!
        ndiag--;
      nnz++;
      col += 1 + floor(log(1 - rand_double()) / D);
    }
    row++;
    col -= end;
  }
  if(loops == LOOPS) nnz += ndiag;

  int weighted = (edge_type == UNDIRECTED_WEIGHTED || edge_type == DIRECTED_WEIGHTED);
  sparsemat_t * mat = init_matrix(n, n, nnz, weighted);

  if(!mat){ printf("ERROR: erdos_renyi_random_graph: init_matrix failed!\n"); return(NULL); }

  // fill in the nonzeros
  rand_seed(seed);
  nnz = 0;
  row = 0;
  mat->offset[0] = 0;
  col = 1 + floor(log(1 - rand_double()) / D);
  while(row < n){
    int need_diag = (loops==LOOPS);
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row;
    while(col < end){
      if(col == row) need_diag = 0;
      mat->nonzero[nnz++] = col;
      col += 1 + floor(log(1 - rand_double()) / D);
    }
    
    if(need_diag){
      mat->nonzero[nnz++] = row;
    }
    row++;
    mat->offset[row] = nnz;
    col -= end;
  }

  // fill in weights
  if(weighted){
    int64_t i;
    for(i = 0; i < nnz; i++){
      mat->value[i] = rand_double();
    }
  }
  assert(mat->nnz == nnz);
  if(loops == LOOPS && (edge_type == DIRECTED || edge_type == DIRECTED_WEIGHTED))
    sort_nonzeros(mat); // to get the diagonal entry sorted correctly
  return(mat);
}

/*!
\brief generates the adjacency matrix for an Erdos-Renyi random graph by flipping a coin for each possible edge.
\param n The total number of vertices in the graph (rows in the matrix).
\param p The probability that each non-loop edge is present.
\param edge_type See edge_type.
\param loops See self_loops.
\param seed A random seed.
\return A sparsemat_t

This is here mostly as an illustration. 
It is the definition of the graph, but computing it this way
is so slow that it can only be used for small examples.
*/
sparsemat_t * erdos_renyi_random_graph_naive(int64_t n, double p, edge_type edge_type, self_loops loops, int64_t seed)
{
  int64_t row, col;
  int64_t pos;
  
  assert(n > 0);
  assert(0.0<p && p<1.0);
  
  rand_seed(seed);
  int64_t nnz = 0;
  int64_t end = n;
  for(row = 0; row < n; row++){
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row;    
    for(col = 0; col < end; col++){
      if(col == row) 
        continue;
      if(rand_double() < p)
        nnz++;
    }
  }
  if(loops == LOOPS) nnz += n;
  
  int weighted = (edge_type == UNDIRECTED_WEIGHTED || edge_type == DIRECTED_WEIGHTED);
  sparsemat_t * mat = init_matrix(n, n, nnz, weighted);
  if(!mat){ printf("ERROR: erdos_renyi_random_graph_naive: init_matrix failed!\n"); return(NULL); }

  // fill in the nonzeros
  rand_seed(seed);
  pos = 0;
  mat->offset[0] = 0;
  for(row = 0; row < n; row++){
    if(edge_type == UNDIRECTED || edge_type == UNDIRECTED_WEIGHTED)
      end = row + (loops == LOOPS);
    for(col = 0; col < end; col++){
      if(col == row && loops == LOOPS){
        mat->nonzero[pos++] = row;
        continue;
      }
      if(rand_double() < p){
        mat->nonzero[pos++] = col;
      }
    }
    mat->offset[row+1] = pos;
  }
  // fill in the weights
  if(weighted){
    int64_t i;
    for(i = 0; i < pos; i++){
      mat->value[i] = rand_double();
    }
  }
  return( mat);
}


/*! \brief prints some stats of a sparse matrix
 * \param mat the sparse matrix
 */
void spmat_stats(sparsemat_t *mat) 
{
  printf("   mat->numrows  = %12"PRId64"\n", mat->numrows);
  printf("   mat->numcols  = %12"PRId64"\n", mat->numcols);
  printf("   mat->nnz      = %12"PRId64"\n", mat->nnz);
  
  int64_t i, d, mindeg, maxdeg, cntdeg;
  double avgdeg;
  mindeg = mat->numcols;
  maxdeg = 0;
  cntdeg = 0;
  for(i=0; i<mat->numrows; i++){
    d = mat->offset[i+1]-mat->offset[i];
    cntdeg += d;
    mindeg = (d < mindeg) ? d : mindeg;
    maxdeg = (d > maxdeg) ? d : maxdeg;
  }
  avgdeg = (double)cntdeg/(double)mat->numrows;
  
  printf("  min, avg, max degree = %"PRId64", %g, %"PRId64"\n", mindeg, avgdeg, maxdeg);
}


/*!
\brief frees the space allocated for a sparse matrix
\param mat pointer to the sparse matrix
*/
void clear_matrix(sparsemat_t * mat)
{
  free(mat->nonzero);
  free(mat->offset);
  free(mat->value);
}


/*!
\brief initializes the struct that holds an array of doubles
\param num total number of entries
\return an allocated d_array_t or NULL on error
*/
d_array_t * init_d_array(int64_t num) 
{
  d_array_t * array = calloc(1, sizeof(d_array_t));
  array->num  = num;
  array->entry   = malloc(num * sizeof(double));
  return(array);
}

/*!
\brief sets all the entries of d_array_t to a value
\param A the array
\param v the value
*/
void set_d_array(d_array_t * A, double v) 
{
  int64_t i;
  for(i=0; i<A->num; i++) {
    A->entry[i] = v;
  }
}

/*!
\brief produces a copy of a source array
\param src the sourcearray
*/
d_array_t *copy_d_array(d_array_t * src) 
{
  int64_t i;
  d_array_t * ret = init_d_array(src->num);
 
  for(i=0; i<src->num; i++)
    ret->entry[i] = src->entry[i];
  return(ret);
}

/*!
\brief replaces the destination array by overwriting the source array
\param dest the destination array
\param src the source array
\return 0 or 1  on an error
NB: The destination must be allocated and of the right size.
*/
int64_t replace_d_array(d_array_t *dest, d_array_t *src) 
{
  int64_t i;
  if(dest->num != src->num){
    fprintf(stderr,"ERROR: replace_d_array: arrays lengths don't match\n");
    return(0);
  }
  for(i=0; i<src->num; i++) {
    dest->entry[i] = src->entry[i];
  }
  return(1);
}


/*!
\brief read a double array from a file
\param name the filename to be read
\return a pointer to the double array or NULL on failure

Note: file format is: 
 - lines that begin with '#' are comment lines
 - first non-comment line is the num of entries in the array,
 - that many lines of one double per line in the format "%lf"
*/
d_array_t * read_d_array(char *name)
{
  int64_t i, num;
  double t;
  char inputstr[1024];
  d_array_t *retarr;

  FILE * fp = fopen(name, "r");
  if ( fp == NULL ){
    fprintf(stderr,"read_d_array: can't open file %s \n", name);
    return(NULL);
  }

  i = 0;
  num = 0; //also flag for not started yet
  while (fgets(inputstr, sizeof inputstr, fp)) {
    if (inputstr[0] == '#') continue;  // skip comment lines
    if (num) {
      if(sscanf(inputstr, "%lf\n", &t) != 1) {
        fprintf(stderr,"read_d_array: file format error, line %ld\n", i);
        return(NULL);
      }
      retarr->entry[i++] = t;
    } else {   // read number of rows
      if(sscanf(inputstr, "%"SCNd64"\n", &num) != 1) {
        fprintf(stderr,"read_d_array: file format error reading num\n");
        return(NULL);
      }
      num = t;
      retarr = init_d_array(num);
    }
  }
  if ( i != num) {
    fprintf(stderr,"read_d_array: read %ld of expected %ld lines\n", i, num);
    return(NULL);
  }
  return(retarr);
}

/*!
\brief writes a double array to a file
\param A pointer to the double array
\param name the filename to be written to
\return 0 on success, non-0 on error.
*/
int64_t write_d_array(d_array_t *A, char *comment,  char * name)
{
  int64_t i;

  FILE * fp = fopen(name, "w");
  if( fp == NULL ){
    fprintf(stderr,"write_d_array: can't open file %s \n", name);
    exit(1);
  }

  fprintf(fp, "#%s\n", comment);
  fprintf(fp, "%"PRId64"\n", A->num);

  for(i=0; i<A->num; i++){
    fprintf(fp, "%lf\n", A->entry[i]);
  }
  fclose(fp);
  return(0);
}


/*!
\brief clears the space for d_array_t
\param A the d_array
*/
void clear_d_array(d_array_t * A)
{
  free(A->entry);
}


