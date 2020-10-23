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
#include <spmat.h>
#include <exstack.h>
#include <sys/stat.h> // for mkdir()

/*! \file spmat_exstack.upc
 * \brief spmat routines written using exstack
 */ 

/******************************************************************************/
/*! \brief create a global int64_t array with a uniform random permutation
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \param buf_cnt number of items in an exstack buffer
 * \return the permutation
 * 
 * This is a collective call.
 * This implements the random dart algorithm to generate the permutation using exstack.
 * Each thread throws its elements of the perm array randomly at large target array by
 * pushing a throw on to the appropriate exstack.
 * The popping pe then claims the requested slot or returns the dart to the sender.
 * 
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 *  
 * \ingroup spmatgrp
 */
SHARED int64_t * rand_permp_exstack(int64_t N, int seed, int64_t buf_cnt) {
  int ret;
  int64_t i, j, cnt, pe, pos, fromth, iend;
  int64_t val;
  
  int64_t lN = (N + THREADS - MYTHREAD - 1)/THREADS;
  int64_t M = N * 2;
  int64_t lM = (M + THREADS - MYTHREAD - 1)/THREADS;

  typedef struct pkg_t{
    int64_t idx; 
    int64_t val;
  }pkg_t;

  SHARED int64_t * perm = lgp_all_alloc(N, sizeof(int64_t));
  if( perm == NULL ) return(NULL);
  int64_t * lperm = lgp_local_part(int64_t, perm);

  SHARED int64_t * target = lgp_all_alloc(M, sizeof(int64_t));
  if( target == NULL ) return(NULL);
  int64_t * ltarget = lgp_local_part(int64_t, target);
  
  /* initialize perm[i] = i,  the darts*/
  for(i = 0; i < lN; i++)
    lperm[i] = i*THREADS + MYTHREAD;

  /* initialize target[i] = -1 */
  for(i = 0; i < lM; i++)
    ltarget[i] = -1L;

  lgp_rand_seed(seed);

  lgp_barrier();
  
  double t1 = wall_seconds();
  int64_t rejects = 0;
  pkg_t pkg;
  exstack_t * ex = exstack_init(buf_cnt, sizeof(pkg_t));
  if(ex == NULL){return(NULL);}
  
  iend = 0L;
  while(exstack_proceed(ex, (iend == lN))){
    i = iend;
    while(i < lN){
      int64_t r = lgp_rand_int64(M);

      pe = r % THREADS;
      pkg.idx = r/THREADS;
      pkg.val = lperm[i];
      if(exstack_push(ex, &pkg, pe) == 0L)
        break;
      i++;
    }
    iend = i;
    
    exstack_exchange(ex);
    
    while(exstack_pop(ex, &pkg, &fromth)){
      if(ltarget[pkg.idx] == -1L){
        ltarget[pkg.idx] = pkg.val;
      }else{
        /* push this pkg back to the sender */
        exstack_push(ex, &pkg, fromth);
        //rejects++;
      }
    }

    lgp_barrier();
    exstack_exchange(ex);
    
    while(exstack_pop(ex, &pkg, &fromth)){
      lperm[--iend] = pkg.val;
    }
    lgp_barrier();
  }
  
  lgp_barrier();
  t1 = wall_seconds() - t1;
  //T0_printf("phase 1 t1 = %8.3lf\n", t1);
  //rejects = lgp_reduce_add_l(rejects);
  //T0_printf("rejects = %"PRId64"\n", rejects);

  exstack_clear(ex);
  free(ex);

  /* now locally pack the values you have in target */
  cnt = 0;
  for(i = 0; i < lM; i++){
    if( ltarget[i] != -1L ) {
      ltarget[cnt++] = ltarget[i];
    }
  }
  lgp_barrier();
  
  /* sanity check */
  int64_t total = lgp_reduce_add_l(cnt);
  if(total != N){
    T0_printf("ERROR: rand_permp_exstack: total = %"PRId64" should be %"PRId64"\n", total, N);
    return(NULL);
  }

  exstack_t * ex1 = exstack_init(buf_cnt, sizeof(int64_t));
  if(ex1 == NULL){return(NULL);}

  int64_t offset = lgp_prior_add_l(cnt);
  pe = offset % THREADS;
  i = pos = 0;
  while(exstack_proceed(ex1, (i==cnt))) {
    while(i < cnt){
      if(exstack_push(ex1, &ltarget[i], pe) == 0)
        break;
      i++;
      pe++;
      if(pe == THREADS) pe = 0;
    }
    
    exstack_exchange(ex1);

    while(exstack_pop(ex1, &val, &fromth)){
      lperm[pos++] = val;
    }
  }
  
  pos = lgp_reduce_add_l(pos);
  if(pos != N){
    printf("ERROR! in rand_permp! sum of pos = %"PRId64" lN = %"PRId64"\n", pos, N);
    return(NULL);
  }
  lgp_barrier();

  exstack_clear(ex1);
  free(ex1);

  return(perm);
}

/*! \brief apply row and column permutations to a sparse matrix using exstack2 
 * \param A pointer to the original matrix
 * \param rperm pointer to the global array holding the row permutation
 * \param cperm pointer to the global array holding the column permutation
 * rperm[i] = j means that row i of A goes to row j in matrix Ap
 * cperm[i] = j means that col i of A goes to col j in matrix Ap
 * \param buf_cnt number of items in an exstack buffer
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * permute_matrix_exstack(sparsemat_t * A, SHARED int64_t * rperm, SHARED int64_t * cperm, int64_t buf_cnt) {
  typedef struct pkg_rowcnt_t{
    int64_t row;
    int64_t cnt;
  }pkg_rowcnt_t;
  typedef struct pkg_inonz_t{
    int64_t i;    
    int64_t nonz;
  }pkg_inonz_t;

  sparsemat_t * Ap;
  
  int64_t i, fromth, fromth2, pe, row, lnnz;
  pkg_rowcnt_t pkg_rc;
  int64_t * lrperm = lgp_local_part(int64_t, rperm);
  int64_t * lcperm = lgp_local_part(int64_t, cperm);
  
  //if(!MYTHREAD)printf("Permuting matrix with exstack...");
  
  /****************************************************************/
  // distribute row counts to the permuted matrix and count the number of nonzeros per thread
  // in the permuted matrix. tmprowcnts holds the post-rperm rowcounts 
  /****************************************************************/
  int64_t * tmprowcnts = calloc(A->lnumrows + 1, sizeof(int64_t));
  
  exstack_t * ex = exstack_init( buf_cnt, sizeof(pkg_rowcnt_t));
  if( ex == NULL ) return(NULL);
  lnnz = row = 0;
  while(exstack_proceed(ex, (row == A->lnumrows))) {
    while(row < A->lnumrows){
      pe = lrperm[row] % THREADS;
      pkg_rc.row = lrperm[row] / THREADS;
      pkg_rc.cnt = A->loffset[row+1] - A->loffset[row];
      if(exstack_push(ex, &pkg_rc, pe) == 0L)
        break;
      row++;
    }

    exstack_exchange(ex);

    while(exstack_pop(ex, &pkg_rc, &fromth)){
      tmprowcnts[pkg_rc.row] = pkg_rc.cnt;
      lnnz += pkg_rc.cnt;
    }
  }
  lgp_barrier();
  exstack_clear(ex);
  free(ex);

  assert(A->nnz == lgp_reduce_add_l(lnnz));

  // allocate pmat to the max of the new number of nonzeros per thread
  int weighted = (A->value != NULL);

  Ap = init_matrix(A->numrows, A->numcols, lnnz, weighted);
  if(Ap == NULL) return(NULL);

  lgp_barrier();

  // convert row counts to offsets 
  Ap->loffset[0] = 0;
  for(i = 1; i < Ap->lnumrows+1; i++)
    Ap->loffset[i] = Ap->loffset[i-1] + tmprowcnts[i-1];

  /****************************************************************/
  // re-distribute nonzeros
  // working offset: wrkoff[row] gives the first empty spot on row row  
  /****************************************************************/
  int64_t * wrkoff = calloc(A->lnumrows, sizeof(int64_t)); 

  exstack_t * ex1;
  if(weighted)
    ex1 = exstack_init( buf_cnt, sizeof(w_edge_t));
  else
    ex1 = exstack_init( buf_cnt, sizeof(edge_t));  
  if( ex1 == NULL )return(NULL);

  edge_t edge;
  w_edge_t wedge;
  i = row = 0;
  while(exstack_proceed(ex1, (i == A->lnnz))){
    while(i < A->lnnz){
      while( i == A->loffset[row+1] ) 
        row++;
      pe = lrperm[row] % THREADS;
      if(weighted){
        wedge.row = lrperm[row] / THREADS;
        wedge.col = A->lnonzero[i];
        wedge.val = A->lvalue[i];
        if( !exstack_push(ex1, &wedge, pe ))
          break;
        i++;

      }else{
        edge.row = lrperm[row] / THREADS;
        edge.col = A->lnonzero[i];
        if( !exstack_push(ex1, &edge, pe ))        
          break;
        i++;
      }
    }

    exstack_exchange(ex1);

    if(weighted){
      while(exstack_pop(ex1, &wedge, &fromth)) {
        Ap->lnonzero[ Ap->loffset[wedge.row] + wrkoff[wedge.row] ] = wedge.col;
        Ap->lvalue[ Ap->loffset[wedge.row] + wrkoff[wedge.row] ] = wedge.val;
        wrkoff[wedge.row]++;
      }
    }else{
      while(exstack_pop(ex1, &edge, &fromth)) {
        Ap->lnonzero[ Ap->loffset[edge.row] + wrkoff[edge.row] ] = edge.col;
        wrkoff[edge.row]++;
      }
    }
  }
  lgp_barrier();
  
  /* sanity check */
  int64_t error = 0L;
  for(i = 0; i < Ap->lnumrows; i++){
    if(wrkoff[i] != tmprowcnts[i]){printf("T%d: w[%"PRId64"] = %"PRId64" trc[%"PRId64"] = %"PRId64"\n", MYTHREAD, i, wrkoff[i], i, tmprowcnts[i]);error++;}
  }
  if(error){printf("ERROR! permute_matrix_exstack: error = %"PRId64"\n", error);}

  free(wrkoff);
  free(tmprowcnts);
  exstack_clear(ex1);
  free(ex1);


  /****************************************************************/
  /* do column permutation ... this is essentially an indexgather */
  /****************************************************************/
  pkg_inonz_t pkg3;
  exstack_t * ex2 = exstack_init( buf_cnt, sizeof(pkg_inonz_t));
  if( ex2 == NULL ) return(NULL);
  i=0;
  while(exstack_proceed(ex2,(i == Ap->lnnz))){
    while(i < Ap->lnnz){     // request the new name for this non-zero
      pkg3.i = i;
      pkg3.nonz = Ap->lnonzero[i] / THREADS;
      pe = Ap->lnonzero[i] % THREADS;
      if( !exstack_push(ex2, &pkg3, pe) )
        break;
      i++;
    }

    exstack_exchange(ex2);

    while(exstack_pop(ex2, &pkg3, &fromth)){ // echo the requests back to requester 
      pkg3.nonz = lcperm[pkg3.nonz];
      exstack_push(ex2, &pkg3, fromth);
    }

    lgp_barrier();

    exstack_exchange(ex2);

    while(exstack_pop(ex2, &pkg3, &fromth)){ //process the echo to your requests
      Ap->lnonzero[pkg3.i] = pkg3.nonz;
    }
  }

  exstack_clear(ex2);
  free(ex2);

  lgp_barrier();
  
  sort_nonzeros(Ap);

  return(Ap);
}

/*! \brief produce the transpose of a sparse matrix using exstack2
 * \param A  pointer to the original matrix
 * \param buf_cnt number of items in an exstack buffer
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 * \ingroup spmatgrp
 */
sparsemat_t * transpose_matrix_exstack(sparsemat_t * A, int64_t buf_cnt) {

  int64_t ret, pe;
  int64_t lnnz, i, col, row, fromth; 
  int64_t idx, *idxp;

  
  /* get the colcnts */
  int64_t lnumcols = (A->numrows + THREADS - MYTHREAD - 1)/THREADS;  
  int64_t * lcounts = calloc(lnumcols, sizeof(int64_t));
  lgp_barrier();
  int weighted = (A->value != NULL);
  
  exstack_t * ex = exstack_init( buf_cnt, sizeof(int64_t));
  if( ex == NULL ) return(NULL);
  
  lnnz = i = 0;
  while(exstack_proceed(ex, (i == A->lnnz))){
    while(i < A->lnnz){
      col = A->lnonzero[i] / THREADS;
      pe = A->lnonzero[i] % THREADS;
      if( !exstack_push(ex, &col, pe) )
        break;
      i++;        
    }
    
    exstack_exchange(ex);
    
    while(exstack_pop(ex, &idx, &fromth)){
      lcounts[idx]++; 
      lnnz++;
    }
  }
  exstack_clear(ex);
  free(ex);

  int64_t sum = lgp_reduce_add_l(lnnz);
  assert( A->nnz == sum ); 
  
  sparsemat_t * At = init_matrix(A->numcols, A->numrows, lnnz, weighted);
  if(!At){printf("ERROR: transpose_matrix: init_matrix failed!\n");return(NULL);}

  /* convert colcounts to offsets */
  At->loffset[0] = 0;  
  for(i = 1; i < At->lnumrows+1; i++)
    At->loffset[i] = At->loffset[i-1] + lcounts[i-1];
    
  lgp_barrier();
  
  /* redistribute the nonzeros */
  int64_t *wrkoff = calloc(At->lnumrows, sizeof(int64_t));
  if(!wrkoff) {printf("ERROR: transpose_matrix: wrkoff alloc fail!\n"); return(NULL);}

  exstack_t * ex1;
  if(weighted){
    ex1 = exstack_init( buf_cnt, sizeof(w_edge_t));
  }else{
    ex1 = exstack_init( buf_cnt, sizeof(edge_t));
  }
  
  if( ex1 == NULL ) return(NULL);
  w_edge_t wedge;
  edge_t edge;
  uint64_t numtimespop=0;
  i = row = 0;
  while(exstack_proceed(ex1, (i == A->lnnz))){
    while(i < A->lnnz){
      while( i == A->loffset[row+1] ){
        row++;
        assert(row < A->lnumrows);
      }
      pe = A->lnonzero[i]%THREADS;
      if(weighted){
        wedge.row = row * THREADS + MYTHREAD;
        wedge.col = A->lnonzero[i] / THREADS;
        wedge.val = A->lvalue[i];        
        if( exstack_push(ex1, &wedge, pe) == 0 )
          break;
        i++;  
      }else{
        edge.row = row * THREADS + MYTHREAD;
        edge.col = A->lnonzero[i] / THREADS;
        if( exstack_push(ex1, &edge, pe) == 0 )
          break;
        i++;
      }
    }

    exstack_exchange(ex1);

    if(weighted){
      while(exstack_pop(ex1, &wedge, &fromth)){
        numtimespop++;
        At->lnonzero[ At->loffset[wedge.col] + wrkoff[wedge.col] ] = wedge.row;
        At->lvalue[ At->loffset[wedge.col] + wrkoff[wedge.col] ] = wedge.val;
        wrkoff[wedge.col]++;
      }
    }else{
      while(exstack_pop(ex1, &edge, &fromth)){
        numtimespop++;
        At->lnonzero[ At->loffset[edge.col] + wrkoff[edge.col] ] = edge.row;
        wrkoff[edge.col]++;
      }
    }
  }

  lgp_barrier();
  exstack_clear(ex1);
  free(ex1);

  //if(!MYTHREAD)printf("done\n");

  numtimespop = lgp_reduce_add_l(numtimespop);
  if(numtimespop != A->nnz ){
    printf("ERROR: numtimespop %"PRId64" \n", numtimespop);
    printf("%d wrkoff %"PRId64"\n", MYTHREAD, wrkoff[0]);
    return(NULL);
  }

  for(i = 0; i < At->lnumrows; i++){
    if(wrkoff[i] != lcounts[i] ) {
      printf("ERROR: %d wrkoff[%"PRId64"] = %"PRId64" !=  %"PRId64" = lcounts[%"PRId64"]\n", MYTHREAD, i, wrkoff[i],lcounts[i],i);
      return(NULL);
    }
  }
  
  free(wrkoff);
  free(lcounts);
  return(At);
}


/*! \brief Reads a distributed sparse matrix in parallel using an arbitrary number of 
 * PEs (between 1 and THREADS) to do the actual I/O.
 *
 * \param datadir The directory that holds the dataset.
 * \param nreaders Number of PEs that will actually do the I/O (all other PEs will mostly be idle).
 * \param buf_cnt number of packages in an exstack buffer
 * \return The sparsemat_t that was read from the matrix dataset.
 */
sparsemat_t * read_sparse_matrix_exstack(char * datadir, int64_t nreaders, int64_t buf_cnt){
  int64_t i;
  exstack_t * ex, * vex;
  spmat_dataset_t * spd = open_sparse_matrix_read(datadir, nreaders);

  // Read row_info to get the counts of the rows we will read.
  // These land in the spd->rowcnt array.
  read_row_info(spd);

  /*****************************************************************/
  /* Distribute the row counts from spd->rowcnt (in BLOCK layout)  */  
  /* to rowcount array (in CYCLIC layout)                          */
  /*****************************************************************/
  SHARED int64_t * rowcnt = lgp_all_alloc(spd->numrows, sizeof(int64_t));
  assert(rowcnt != NULL);
  int64_t * lrowcnt = lgp_local_part(int64_t, rowcnt);

  /* offset has the first row (local index) I will see from each reader */
  int64_t * offset = calloc(THREADS, sizeof(int64_t));
  for(i = 0; i < THREADS; i++)
    offset[i] = spd->global_first_row_to_me[i]/THREADS; 
  
  ex = exstack_init(buf_cnt, sizeof(int64_t));
  assert(ex != NULL);
  
  i = 0;
  int64_t lnnz = 0, count, from_pe;  
  while(exstack_proceed(ex, (i==spd->lnumrows))){
    for(; i < spd->lnumrows; i++){
      if(exstack_push(ex, &spd->rowcnt[i], (spd->global_first_row + i) % THREADS) == 0)
        break;
    }

    exstack_exchange(ex);
    
    while(exstack_pop(ex, &count, &from_pe)){
      assert(offset[from_pe] < (spd->numrows + THREADS - MYTHREAD - 1)/THREADS);
      lrowcnt[offset[from_pe]] = count;
      lnnz += count;
      offset[from_pe]++;      
    }
  }
  
  exstack_reset(ex);
  lgp_barrier();
  
  //fprintf(stderr,"PE %d gets %ld nonzeros\n", MYTHREAD, lnnz);fflush(stderr);

  /****************************/
  /* initialize sparse matrix */
  /****************************/
  sparsemat_t * A = init_matrix(spd->numrows, spd->numcols, lnnz, spd->values);
  if(!A){
    T0_fprintf(stderr,"ERROR: read_sparse_matrix_exstack: init_matrix failed!\n");return(NULL);
  }

  /* create A->offset array */
  A->loffset[0] = 0;
  for(i = 0; i < A->lnumrows; i++){
    A->loffset[i+1] = A->loffset[i] + lrowcnt[i];
  }

  /* set up a pointer for each reader (this is where you will be putting nonzeros you receive from them) */
  for(i = 0; i < THREADS; i++){
    offset[i] = A->loffset[spd->global_first_row_to_me[i]/THREADS];
  }
  
  if(spd->values)
    vex = exstack_init(buf_cnt, sizeof(double));
  
  lgp_barrier();
  lgp_all_free(rowcnt);
  
  /****************************************************/
  /* read the nonzeros, send via exstack, and stick   */
  /* them into the sparse matrix data structure.      */
  /****************************************************/
  int64_t buf_size = 512*512;
  int64_t * buf = calloc(buf_size, sizeof(int64_t));
  double * vbuf;
  if(spd->values)
    vbuf = calloc(buf_size, sizeof(double));  
  int64_t tot_rows_read = 0;
  int64_t global_row = spd->global_first_row;
  int64_t pos = 0;
  int64_t j = 0;
  int64_t row = 0;
  int64_t num_rows_read = 0;
  int64_t nnz_read = 0;
  int64_t rc = 0;
  int64_t pe;
  
  while(exstack_proceed(ex, (tot_rows_read == spd->lnumrows))){
    int64_t col;
    double val;
    int loop_break = (tot_rows_read == spd->lnumrows);
    
    while(!loop_break){

      if(row == num_rows_read){
        /* read a buffer of nonzeros if needed */
        pos = 0;
        num_rows_read = read_nonzeros_buffer(spd, buf, vbuf, buf_size);
        assert(num_rows_read > 0);
        row = 0;
        nnz_read = 0;
        for(i = 0; i < num_rows_read; i++) nnz_read += spd->rowcnt[tot_rows_read + i];
        pe = global_row % THREADS;
        //fprintf(stderr,"PE %d: read %ld rows %ld %ld %ld\n", MYTHREAD, num_rows_read, tot_rows_read, pos, nnz_read);
      }

      /* push the nonzero data for the row you are working on */
      for(; j < spd->rowcnt[tot_rows_read]; j++){
        if(spd->values) exstack_push(vex, &vbuf[pos], pe);
        if(exstack_push(ex, &buf[pos], pe) == 0){
          loop_break = 1;
          break;
        }
        pos++;
      }

      /* if you are done pushing the nonzeros of this row, tee up the next one. */
      if(j == spd->rowcnt[tot_rows_read]){
        tot_rows_read++;
        if(tot_rows_read < spd->lnumrows){
          global_row++;
          row++;
          j = 0;
          pe = (pe == (THREADS - 1)) ? 0 : pe+1;
        }else
          loop_break = 1;
      }
    }
    
    exstack_exchange(ex);
    if(spd->values) exstack_exchange(vex);

    // pop the nonzeros that were sent to you
    while(exstack_pop(ex, &col, &from_pe)){
      A->lnonzero[offset[from_pe]] = col;
      if(spd->values){
        exstack_pop(vex, &val, &from_pe);
        A->lvalue[offset[from_pe]] = val;
      }
      offset[from_pe]++;
    }
  }
  
  free(buf);
  if(spd->values) free(vbuf);

  exstack_clear(ex);
  if(spd->values) exstack_clear(vex);
  clear_spmat_dataset(spd);
  
  lgp_barrier();
  return(A);
}




/****************************************************************************/
/*! \brief Write a sparsemat_t to disk as a sparse matrix dataset.
 *
 * This routine is collective and it's return is single-valued.
 *
 * PE i will collect the ith chunk of rows (and their nonzeros) to write. 
 * 
 * \param dirname The directory where the sparsemat will be written (must be single-valued).
 * \param mat The sparsemat_t to be written.
 * \param buf_cnt number of packets in an exstack buffer
 * \return 0 on SUCCESS nonzero on ERROR.
 */
/****************************************************************************/
int64_t write_sparse_matrix_exstack( char * dirname, sparsemat_t * mat, int64_t buf_cnt){

  int64_t *rowcnt_buf;
  int64_t *nz_buf;
  uint64_t row, col, cnt, max, lrow, pe;
  uint64_t * current_row_to_th, * last_row;

  int64_t i, j, k, rows_per_pass, nr = mat->numrows;
  int64_t pass, num_passes, first_pe;
  int64_t ret, error;
  exstack_t * ex;

  spmat_dataset_t * spd = open_sparse_matrix_write(dirname, mat);

  write_row_info(spd, mat);
  
  /* get max row density */
  max = 0;
  for(row = 0; row < mat->lnumrows; row++){
    cnt = mat->loffset[row + 1] - mat->loffset[row];
    max = (cnt > max ? cnt : max);
  }

  max = lgp_reduce_max_l(max);
  if(max==0) max = 1;

  rows_per_pass = buf_cnt/(max+1);
  while(rows_per_pass == 0){
    buf_cnt *= 2;
    rows_per_pass = buf_cnt/(max+1);
  }
  //T0_fprintf("rows_per_pass = %"PRId64", max = %lu\n", rows_per_pass, max);
  
  /* allocate space */
  rowcnt_buf         = calloc(THREADS, sizeof(int64_t));  
  current_row_to_th  = calloc(THREADS, sizeof(uint64_t));
  last_row           = calloc(THREADS, sizeof(uint64_t));
  nz_buf             = calloc(max*THREADS, sizeof(int64_t));
  
  ex = exstack_init((max + 1)*rows_per_pass, sizeof(uint64_t));
  if(!ex){
    fprintf(stderr,"write_sparse_matrix_exstack: exstack_init failed\n");lgp_barrier();
    return(-1);
  }
  
  /* first figure out the first row each PE will write. */
  current_row_to_th[0] = 0L;
  for(i = 1; i < THREADS; i++)
    current_row_to_th[i] =  current_row_to_th[i-1] + (nr + THREADS - (i-1) - 1)/THREADS;

  /* now figure out which row on this PE should get sent first to each other PE */
  /* and determine which PE owns the first row you will be writing */
  for(i = 0; i < THREADS; i++){
    if(i == MYTHREAD)
      first_pe = current_row_to_th[i] % THREADS; /* this is the pe who will send your first row to write */
    while((current_row_to_th[i] % THREADS) != MYTHREAD)
      current_row_to_th[i]++;
    if(i > 0)
      last_row[i-1] = current_row_to_th[i];
  }
  last_row[THREADS - 1] = mat->numrows;
  
  /* open nonzero file */
  char fname[128];
  sprintf(fname, "%s/nonzero_%d", spd->dirname, MYTHREAD);
  FILE * nzfp = fopen(fname, "wb");
  
  /*********************************/
  /*        WRITE OUT DATASET      */
  /*********************************/
  /* For each PE, we keep track of what row we are currently on for that PE (current_row_to_th).
     Each exstack pass, we add rows_per_pass rows worth of nonzero data to the buffer for that PE.
     This is a very different kind of exstack loop, where we don't rely on exstack to tell us when we
     are done, but rather carefully push an amount of data in each pass that will not overflow the 
     buffers.
  */
  
  int64_t room;
  num_passes = (nr + THREADS - 1)/THREADS;
  num_passes /= (rows_per_pass * THREADS);
  //T0_fprintf(stderr,"num_passes = %"PRId64"\n", num_passes);
  
  error = 0;
  int imdone = 0;
  int64_t recs_written, nnz;
  int64_t rc_pos = 0;
  int64_t nz_pos = 0;
  while(exstack_proceed(ex, imdone)){

    imdone = 1;

    /* add rows_per_pass rows to everyone's buffers */
    for(i = 0; i < THREADS; i++){
      row = current_row_to_th[i];
      lrow = row/THREADS;
      for(k = 0; (row < last_row[i]) && (k < rows_per_pass); row+=THREADS, lrow++,k++){
        cnt = mat->loffset[lrow+1] - mat->loffset[lrow];
        room = exstack_push(ex, &cnt, i);
        if(room < 0L){
          fprintf(stderr, "ERROR: write_sparse_matrix: Trying to push cnt onto full exstack!\n");
          error = 1;
        }
        for(j = mat->loffset[lrow]; j < mat->loffset[lrow+1]; j++){
          col = mat->lnonzero[j];
          room = exstack_push(ex, &col, i);
          if(room < 0L){
            fprintf(stderr,"ERROR: write_sparse_matrix: Trying to push cols onto full exstack!\n");
            error = 1;
          }
        }
      }
      current_row_to_th[i] = row;
      /* you are not done if you did not finish pushing to some PE */
      if(current_row_to_th[i] < last_row[i])
        imdone = 0;
    }
      
    exstack_exchange(ex);

    for(k = 0; k < rows_per_pass; k++){
      pe = first_pe;
      nz_pos = rc_pos = 0L;
      for(i = first_pe; i < first_pe + THREADS; i++, pe++){
        pe = (pe == THREADS ? 0 : pe);
        if(exstack_pop_thread(ex, &cnt, pe)){
          for(j = 0; j < cnt; j++){
            exstack_pop_thread(ex, &col, pe);
            nz_buf[nz_pos++] = col;
          }
          rowcnt_buf[rc_pos++] = cnt;
        }
      }

      /* count up the number of nonzeros to be written */
      nnz = 0;
      for(i = 0; i < rc_pos; i++)
        nnz += rowcnt_buf[i];
        
      /* write the nonzeros */
      recs_written = fwrite(nz_buf, sizeof(int64_t), nz_pos, nzfp);
      if(recs_written != nz_pos){
        error = 1;
        fprintf(stderr, "write_sparse_matrix_exstack: recs_written != nnz %"PRId64" %"PRId64" on PE %d\n", recs_written, nz_pos, MYTHREAD);
      }
    }
    lgp_barrier();
  }
  
  fclose(nzfp);
  
  /* make sure all threads wrote cleanly */
  ret = lgp_reduce_add_l(error);
  if(ret){
    T0_fprintf(stderr,"ERROR: write_sparse_matrix_exstack: error in main loop\n");
    lgp_barrier();
    return(-1);
  }

  /* finished writing */
  lgp_barrier();

  free(rowcnt_buf);
  free(nz_buf);
  free(current_row_to_th);
  free(last_row);
  exstack_clear(ex);

  lgp_barrier();

  return(0);
}
