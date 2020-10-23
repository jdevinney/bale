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
/*! \file spmat_agp.upc
 * \brief Sparse matrix support functions implemented with global addresses and atomics
 */
#include <spmat.h>
#include <sys/stat.h>   // for mkdir()
#include <fcntl.h>

/******************************************************************************/
/*! \brief create a global int64_t array with a uniform random permutation
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \return the permutation
 * 
 * This is a collective call.
 * this implements the random dart algorithm to generate the permutation.
 * Each thread throws its elements of the perm array randomly at large target array.
 * Each element claims a unique entry in the large array using compare_and_swap.
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 * \ingroup spmatgrp
 */
SHARED int64_t * rand_permp_agp(int64_t N, int seed) {  
  int64_t * ltarget, *lperm;
  int64_t r, i, j;
  int64_t pos, numdarts, numtargets, lnumtargets;

  lgp_rand_seed(seed);

  //T0_printf("Entering rand_permp_atomic...");fflush(0);

  SHARED int64_t * perm = lgp_all_alloc(N, sizeof(int64_t));
  if( perm == NULL ) return(NULL);
  lperm = lgp_local_part(int64_t, perm);

  int64_t l_N = (N + THREADS - MYTHREAD - 1)/THREADS;
  int64_t M = 2*N;
  int64_t l_M = (M + THREADS - MYTHREAD - 1)/THREADS;

  SHARED int64_t * target = lgp_all_alloc(M, sizeof(int64_t));
  if( target == NULL ) return(NULL);
  ltarget = lgp_local_part(int64_t, target);
  
  for(i=0; i<l_M; i++)
    ltarget[i] = -1L;
  lgp_barrier();

  i=0;
  while(i < l_N){                // throw the darts until you get l_N hits
    r = lgp_rand_int64(M);
    if( lgp_cmp_and_swap(target, r, -1L, (i*THREADS + MYTHREAD)) == (-1L) ){
      i++;
    }
  }
  lgp_barrier();

  numdarts = 0;
  for(i = 0; i < l_M; i++)    // count how many darts I got
    numdarts += (ltarget[i] != -1L );

  pos = lgp_prior_add_l(numdarts);    // my first index in the perm array is the number 
                                      // of elements produce by the smaller threads
  for(i = 0; i < l_M; i++){
    if(ltarget[i] != -1L ) {
       lgp_put_int64(perm, pos, ltarget[i]);
       pos++;
    }
  }

  lgp_all_free(target);
  lgp_barrier();
  //T0_printf("done!\n");
  return(perm);
}


/*! \brief apply row and column permutations to a sparse matrix using straight UPC
 * \param A pointer to the original matrix
 * \param rperm pointer to the global array holding the row permutation
 * \param cperm pointer to the global array holding the column permutation
 * rperm[i] = j means that row i of A goes to row j in matrix Ap
 * cperm[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix_agp(sparsemat_t *A, SHARED int64_t *rperm, SHARED int64_t *cperm) {
  //T0_printf("Permuting matrix with single puts\n");
  int weighted = (A->value != NULL);
  int64_t i, j, col, row, pos;
  int64_t * lrperm = lgp_local_part(int64_t, rperm);
  SHARED int64_t * rperminv = lgp_all_alloc(A->numrows, sizeof(int64_t));
  if( rperminv == NULL ) return(NULL);
  int64_t *lrperminv = lgp_local_part(int64_t, rperminv);

  //compute rperminv from rperm 
  for(i=0; i < A->lnumrows; i++){
    lgp_put_int64(rperminv, lrperm[i], i*THREADS + MYTHREAD);
  }

  lgp_barrier();
  
  int64_t cnt = 0, off, nxtoff;
  for(i = 0; i < A->lnumrows; i++){
    row = lrperminv[i];
    off    = lgp_get_int64(A->offset, row);
    nxtoff = lgp_get_int64(A->offset, row + THREADS);
    cnt += nxtoff - off;
  }
  lgp_barrier();

  sparsemat_t * Ap = init_matrix(A->numrows, A->numcols, cnt, weighted);
  
  // fill in permuted rows
  Ap->loffset[0] = pos = 0;
  for(i = 0; i < Ap->lnumrows; i++){
    row = lrperminv[i];
    off    = lgp_get_int64(A->offset, row);
    nxtoff = lgp_get_int64(A->offset, row + THREADS);
    lgp_memget(&Ap->lnonzero[pos], A->nonzero,
               (nxtoff-off)*sizeof(int64_t), off*THREADS + row%THREADS);    
    if(weighted)
      lgp_memget(&Ap->lvalue[pos], A->value,
                 (nxtoff-off)*sizeof(double), off*THREADS + row%THREADS);
    pos += nxtoff - off;
    //for(j = off; j < nxtoff; j++){
    //Ap->lnonzero[pos++] = lgp_get_int64(A->nonzero, j*THREADS + row%THREADS);
    //}
    Ap->loffset[i+1] = pos;
  }
  
  assert(pos == cnt);
  
  lgp_barrier();
  
  // finally permute column indices
  for(i = 0; i < Ap->lnumrows; i++){
    for(j = Ap->loffset[i]; j < Ap->loffset[i+1]; j++){
      Ap->lnonzero[j] = lgp_get_int64(cperm, Ap->lnonzero[j]);      
    }
  }
  lgp_barrier();

  lgp_all_free(rperminv);
  sort_nonzeros(Ap);

  return(Ap);
}

/*! \brief produce the transpose of a sparse matrix using UPC
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix_agp(sparsemat_t *A) {
  int64_t counted_nnz_At;
  int64_t lnnz, i, j, col, row, fromth, idx;
  int64_t pos;
  sparsemat_t * At;
  
  //T0_printf("UPC version of matrix transpose...");
  
  // find the number of nnz.s per thread

  SHARED int64_t * shtmp = lgp_all_alloc( A->numcols + THREADS, sizeof(int64_t));
  if( shtmp == NULL ) return(NULL);
  int64_t * l_shtmp = lgp_local_part(int64_t, shtmp);
  int64_t lnc = (A->numcols + THREADS - MYTHREAD - 1)/THREADS;
  for(i=0; i < lnc; i++)
    l_shtmp[i] = 0;
  lgp_barrier();

  for( i=0; i< A->lnnz; i++) {                   // histogram the column counts of A
    assert( A->lnonzero[i] < A->numcols );
    assert( A->lnonzero[i] >= 0 ); 
    pos = lgp_fetch_and_inc(shtmp, A->lnonzero[i]);
  }
  lgp_barrier();


  lnnz = 0;
  for( i = 0; i < lnc; i++) {
    lnnz += l_shtmp[i];
  }
  int weighted = (A->value != NULL);
  At = init_matrix(A->numcols, A->numrows, lnnz, weighted);
  if(!At){printf("ERROR: transpose_matrix_upc: init_matrix failed!\n");return(NULL);}

  int64_t sum = lgp_reduce_add_l(lnnz);      // check the histogram counted everything
  assert( A->nnz == sum ); 

  // compute the local offsets
  At->loffset[0] = 0;
  for(i = 1; i < At->lnumrows+1; i++)
    At->loffset[i] = At->loffset[i-1] + l_shtmp[i-1];

  // get the global indices of the start of each row of At
  for(i = 0; i < At->lnumrows; i++)
    l_shtmp[i] = MYTHREAD + THREADS * (At->loffset[i]);
    
  lgp_barrier();

  //redistribute the nonzeros 
  for(row=0; row<A->lnumrows; row++) {
    for(j=A->loffset[row]; j<A->loffset[row+1]; j++){
      pos = lgp_fetch_and_add(shtmp, A->lnonzero[j], (int64_t) THREADS);
      lgp_put_int64(At->nonzero, pos, row*THREADS + MYTHREAD);
      if(weighted) lgp_put_double(At->value, pos, A->lvalue[j]);
    }
  }

  lgp_barrier();
  //if(!MYTHREAD)printf("done\n");
  lgp_all_free(shtmp);

  return(At);
}





/*! \brief Reads a distributed sparse matrix in parallel using an arbitrary number of 
 * PEs (between 1 and THREADS) to do the actual I/O.
 *
 * \param datadir The directory that holds the dataset.
 * \param nreaders Number of PEs that will actually do the I/O (all other PEs will mostly
 * be idle.
 * \return The sparsemat_t that was read from the matrix dataset.
 */
sparsemat_t * read_sparse_matrix_agp(char * datadir, int64_t nreaders){
  int64_t i;
  
  spmat_dataset_t * spd = open_sparse_matrix_read(datadir, nreaders);

  // read row_info to get the counts of the rows we will read
  read_row_info(spd);
  
  /* distribute the row counts in CYCLIC fashion */
  SHARED int64_t * rowcnt = lgp_all_alloc(spd->numrows, sizeof(int64_t));
  assert(rowcnt != NULL);
  int64_t * lrowcnt = lgp_local_part(int64_t, rowcnt);

  for(i = 0; i < spd->lnumrows; i++){
    lgp_put_int64(rowcnt, spd->global_first_row + i, spd->rowcnt[i]);
    //fprintf(stderr,"copying rc %ld to %ld\n", spd->rowcnt[i], spd->global_first_row + i);
  }
  
  lgp_barrier();
  
  /* calculate nnz that will land on each PE */
  int64_t lnnz = 0;
  int64_t lnr = (spd->numrows + THREADS - MYTHREAD - 1)/THREADS;
  for(i = 0; i < lnr; i++){
    //fprintf(stderr, "PE %d: rc[%ld] = %ld\n", MYTHREAD, i, lrowcnt[i]);
    lnnz += lrowcnt[i];
  }
  //fprintf(stderr,"PE %d gets %ld nonzeros\n", MYTHREAD, lnnz);fflush(stderr);
  
  /* initialize sparse matrix */
  sparsemat_t * A = init_matrix(spd->numrows, spd->numcols, lnnz, spd->values);
  if(!A){
    T0_fprintf(stderr,"ERROR: read_sparse_matrix_agp: init_matrix failed!\n");return(NULL);
  }

  /* create A->offset array */
  A->loffset[0] = 0;
  for(i = 0; i < lnr; i++){
    A->loffset[i+1] = A->loffset[i] + lrowcnt[i];
  }

  lgp_barrier();
  lgp_all_free(rowcnt);

  /* read the nonzeros into the matrix data structure */  
  int64_t buf_size = 512*512;
  int64_t * buf = calloc(buf_size, sizeof(int64_t));
  double * vbuf = calloc(buf_size, sizeof(double));
  int64_t tot_rows_read = 0;  
  while(tot_rows_read < spd->lnumrows){
    /* read a buffer of nonzeros */
    int64_t pos = 0;
    int64_t num_rows_read = read_nonzeros_buffer(spd, buf, vbuf, buf_size);
    //fprintf(stderr, "PE %d: numrowsread = %ld\n", MYTHREAD, num_rows_read);
    assert(num_rows_read > 0);

    /* place the read nonzeros into the sparse matrix data structure */
    for(i = 0; i < num_rows_read; i++){
      int64_t rc = spd->rowcnt[tot_rows_read + i];
      int64_t global_row = spd->global_first_row + tot_rows_read + i;
      int64_t rs = lgp_get_int64(A->offset,global_row)*THREADS + global_row%THREADS;
      //for(int j = 0; j < rc; j++)
      //fprintf(stderr,"PE %d: putting %ld for row %ld (at pos %ld)\n", MYTHREAD, buf[pos + j], global_row, rs + j*THREADS);
      lgp_memput(A->nonzero, &buf[pos], rc*sizeof(int64_t), rs);
      if(spd->values)
        lgp_memput(A->value, &vbuf[pos], rc*sizeof(double), rs);
      pos += rc;
    }
    tot_rows_read += num_rows_read;
  }

  free(buf);
  if(spd->values) free(vbuf);
  clear_spmat_dataset(spd);
  lgp_barrier();
  //print_matrix(A);
  return(A);
}




/* \brief Write a sparsemat_t to disk in parallel. 
 * This is different than write_matrix_mm since a) it is done in parallel
 * and b) we are not using matrix market format.
 *
 * Every THREAD just writes a slice of the matrix (in BLOCK fashion).
 * This requires data to be transferred since the matrix is stored
 * in CYCLIC layout in memory. 
 *
 * We write the dataset into PEs*2 files plus an ASCII metadata
 * file. Each PE writes a nonzero file which contains all the nonzeros
 * in their slice and a row_offset file (which says the offsets for
 * each row in the nonzero file). If there are values in the matrix,
 * there will be 3*PEs files.
 *
 * The metadata file just contains some high level info about the matrix
 * dataset.
 * \param dirname The directory to write the matrix to.
 * \param A The matrix to write.
 * \return 0 on SUCCESS 1 on error.
 */
int write_sparse_matrix_agp(char * dirname, sparsemat_t * A){
  int64_t i;
  
  spmat_dataset_t * spd = open_sparse_matrix_write(dirname, A);

  write_row_info(spd, A);
  
  /* allocate write buffer */
  int64_t buf_size = 512*512;
  int64_t * buf = calloc(buf_size, sizeof(int64_t));
  double * vbuf;
  if(spd->values)
    vbuf = calloc(buf_size, sizeof(double));

  char fname[128];
  sprintf(fname, "%s/nonzero_%d", spd->dirname, MYTHREAD);
  spd->nnzfp = fopen(fname, "wb");
  if(spd->values){
    sprintf(fname, "%s/value_%d", spd->dirname, MYTHREAD);
    spd->valfp = fopen(fname, "wb");
  }
  
  /* write out the nonzeros in your block */
  int64_t pos = 0;
  for(i = spd->global_first_row; i < spd->global_first_row + spd->lnumrows; i++){
    int64_t row_start = lgp_get_int64(A->offset, i);
    int64_t row_cnt = lgp_get_int64(A->offset, i + THREADS) - row_start;
    int64_t plc = row_start*THREADS + i%THREADS;
    if(pos + row_cnt >= buf_size){
      fwrite(buf, sizeof(int64_t), pos, spd->nnzfp); 
      if(spd->values) fwrite(vbuf, sizeof(double), pos, spd->valfp);
      pos = 0;
    }

    lgp_memget(&buf[pos], A->nonzero, row_cnt*sizeof(int64_t), plc);
    if(spd->values) lgp_memget(&vbuf[pos], A->value, row_cnt*sizeof(int64_t), plc);
    pos += row_cnt;
  }
  fwrite(buf, sizeof(int64_t), pos, spd->nnzfp);
  if(spd->values) fwrite(vbuf, sizeof(int64_t), pos, spd->valfp);

  lgp_barrier();

  fclose(spd->nnzfp);  
  free(buf);

  if(spd->values){
    fclose(spd->valfp);
    free(vbuf);
  }

  clear_spmat_dataset(spd);
  
  return(0);
  
}

