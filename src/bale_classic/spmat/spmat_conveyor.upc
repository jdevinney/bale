/*****************************************************************
//
//
//  Copyright(C) 2019-2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
//
//  All rights reserved.
//  
//  This file is a part of Bale.  For license information see the
//  LICENSE file in the top level directory of the distribution.
 *****************************************************************/
#include <spmat.h>

/*! \file spmat_conveyor.upc
 * \brief spmat routines written with conveyors
 */ 

/*! \brief create a global int64_t array with a uniform random permutation using conveyors
 * \param N the length of the global array
 * \param seed seed for the random number generator
 * \return the permutation
 * 
 * This is a collective call.
 * this implements the random dart algorithm to generate the permutation.
 * Each thread throws its elements of the perm array randomly at large target array.
 *  
 *
 * This gives a random permutation with spaces in it, then you squeeze out the spaces.
 *
 * This implementation tried to mimic the exstack implementation but it is
 * a little more complicated because we are assuming the conveyors are asynchronous.
 * That makes it harder to figure out when you are done throwing darts. I ended up
 * returning the results of ALL throws (rather than just the failures like in the exstack version).
 * 
 */
SHARED int64_t * rand_permp_conveyor(int64_t N, int seed) {
  int ret;
  int64_t i, j, cnt, pe, pos, fromth, istart, iend;
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
  convey_t * conv_throw = convey_new(SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
  convey_t * conv_reply = convey_new(SIZE_MAX, 0, NULL, convey_opt_PROGRESS);
  if(conv_throw == NULL){return(NULL);}
  if(conv_reply == NULL){return(NULL);}

  convey_begin(conv_throw, sizeof(pkg_t), 0);
  convey_begin(conv_reply, sizeof(int64_t), 0);
  
  bool more;
  int64_t hits = 0;
  iend = 0;
  while(more = convey_advance(conv_throw, (hits == lN)),
        more | convey_advance(conv_reply, !more)){
    i = iend;
    while(i < lN){
      int64_t r = lgp_rand_int64(M);
      pe = r % THREADS;
      pkg.idx = r/THREADS;
      pkg.val = lperm[i];
      if(!convey_push(conv_throw, &pkg, pe))
        break;
      i++;
    }
    iend = i;
    
    while(convey_pull(conv_throw, &pkg, &fromth) == convey_OK){
      if(ltarget[pkg.idx] == -1L){
        val = pkg.val;
        if(!convey_push(conv_reply, &val, fromth)){
          convey_unpull(conv_throw);
          break;
        } 
        ltarget[pkg.idx] = pkg.val;
      }else{
        val = -(pkg.val + 1);
        if(!convey_push(conv_reply, &val, fromth)){
          convey_unpull(conv_throw);
          rejects++;
          break;
        }
      }
    }

    while(convey_pull(conv_reply, &val, NULL) == convey_OK){
      if(val < 0L)
        lperm[--iend] = -val - 1;
      else
        hits++;
    }
  }
  
  lgp_barrier();
  t1 = wall_seconds() - t1;
  //T0_printf("phase 1 t1 = %8.3lf\n", t1);
  //rejects = lgp_reduce_add_l(rejects);
  //T0_printf("rejects = %"PRId64"\n", rejects);

  convey_free(conv_throw);
  convey_reset(conv_reply);
  convey_begin(conv_reply, sizeof(int64_t), 0);
  
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
    T0_printf("ERROR: rand_permp_convey: total = %"PRId64" should be %"PRId64"\n", total, N);
    return(NULL);
  }

  int64_t offset = lgp_prior_add_l(cnt);
  pe = offset % THREADS;
  i = pos = 0;
  while(convey_advance(conv_reply, (i==cnt))) {
    while(i < cnt){
      if(!convey_push(conv_reply, &ltarget[i], pe))
        break;
      i++;
      pe++;
      if(pe == THREADS) pe = 0;
    }
    
    while(convey_pull(conv_reply, &val, &fromth)){
      lperm[pos++] = val;
    }
  }

  pos = lgp_reduce_add_l(pos);
  if(pos != N){
    printf("ERROR! in rand_permp_convey! sum of pos = %"PRId64" lN = %"PRId64"\n", pos, N);
    return(NULL);
  }
  lgp_barrier();
  
  convey_free(conv_reply);

  return(perm);
}


/*! \brief apply row and column permutations to a sparse matrix using conveyors 
 * \param A pointer to the original matrix
 * \param rperm pointer to the global array holding the row permutation
 * \param cperm pointer to the global array holding the column permutation
 * rperm[i] = j means that row i of A goes to row j in matrix Ap
 * cperm[i] = j means that col i of A goes to col j in matrix Ap
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 */
sparsemat_t * permute_matrix_conveyor(sparsemat_t * A, SHARED int64_t * rperm, SHARED int64_t * cperm) {
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
  int64_t * lrperm = lgp_local_part(int64_t, rperm);
  int64_t * lcperm = lgp_local_part(int64_t, cperm);

  //T0_printf("Permuting matrix with conveyors\n");
  
  /****************************************************************/
  // distribute row counts to the permuted matrix and count the number of nonzeros per thread
  // in the permuted matrix. tmprowcnts holds the post-rperm rowcounts 
  /****************************************************************/
  int64_t * tmprowcnts = calloc(A->lnumrows + 1, sizeof(int64_t));

  pkg_rowcnt_t pkg_rc;
  pkg_rowcnt_t pkgrc_p;
  convey_t* cnv_rc = convey_new(SIZE_MAX, 0, NULL, convey_opt_SCATTER);

  convey_begin(cnv_rc, sizeof(pkg_rowcnt_t), 0);
  lnnz = row = 0;
  while(convey_advance(cnv_rc, (row == A->lnumrows))) {
    for( ;row < A->lnumrows; row++){
      pe = lrperm[row] % THREADS;
      pkg_rc.row = lrperm[row] / THREADS;
      pkg_rc.cnt = A->loffset[row+1] - A->loffset[row];
      if( !convey_push(cnv_rc, &pkg_rc, pe) )
        break;
    }
    while( convey_pull(cnv_rc, &pkgrc_p, NULL) == convey_OK ){
      tmprowcnts[pkgrc_p.row] = pkgrc_p.cnt;
      lnnz += pkgrc_p.cnt;
    }
  }
  lgp_barrier();
  convey_free(cnv_rc);  

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
  edge_t edge;
  w_edge_t wedge;
  convey_t* cnv_nz = convey_new(SIZE_MAX, 0, NULL, convey_opt_SCATTER);
  if(weighted)
    convey_begin(cnv_nz, sizeof(w_edge_t), 0);
  else
    convey_begin(cnv_nz, sizeof(edge_t), 0);

  i = row = 0;
  while(convey_advance(cnv_nz, (i == A->lnnz))){
    for( ;i < A->lnnz; i++){
      while( i == A->loffset[row+1] ) // skip empty rows 
        row++;
      pe = lrperm[row] % THREADS;
      if(weighted){
        wedge.row = lrperm[row] / THREADS;
        wedge.col = A->lnonzero[i];
        wedge.val = A->lvalue[i];
        if( !convey_push(cnv_nz, &wedge, pe) )
          break;

      }else{
        edge.row = lrperm[row] / THREADS;
        edge.col = A->lnonzero[i];
        if( !convey_push(cnv_nz, &edge, pe) )
          break;
      }
    }

    if(weighted){
      while( convey_pull(cnv_nz, &wedge, NULL) == convey_OK) {
        Ap->lnonzero[ Ap->loffset[wedge.row] + wrkoff[wedge.row] ] = wedge.col;
        Ap->lvalue[ Ap->loffset[wedge.row] + wrkoff[wedge.row] ] = wedge.val;
        wrkoff[wedge.row]++;
      }
    }else{
      while( convey_pull(cnv_nz, &edge, NULL) == convey_OK) {
        Ap->lnonzero[ Ap->loffset[edge.row] + wrkoff[edge.row] ] = edge.col;
        wrkoff[edge.row]++;
      }
    }
  }
  lgp_barrier();
  convey_free(cnv_nz);

  /* sanity check */
  int64_t error = 0L;
  for(i = 0; i < Ap->lnumrows; i++)
    if(wrkoff[i] != tmprowcnts[i]){printf("w[%"PRId64"] = %"PRId64" trc[%"PRId64"] = %"PRId64"\n", i, wrkoff[i], i, tmprowcnts[i]);error++;}
  if(error){printf("ERROR! permute_matrix_conveyor: error = %"PRId64"\n", error);}

  free(wrkoff);
  free(tmprowcnts);

  /****************************************************************/
  /* do column permutation ... this is essentially an indexgather */
  /****************************************************************/
  pkg_inonz_t pkg_r, pkg_e, pkg_p;
  convey_t* cnv_r = convey_new(SIZE_MAX, 0, NULL, convey_opt_SCATTER);
  assert( cnv_r != NULL );
  convey_t* cnv_e = convey_new(SIZE_MAX, 0, NULL, 0);
  assert( cnv_e != NULL );
  convey_begin(cnv_r, sizeof(pkg_inonz_t), 0);
  convey_begin(cnv_e, sizeof(pkg_inonz_t), 0);
  bool more;
  i=0;
  while( more = convey_advance(cnv_r,(i == Ap->lnnz)), more | convey_advance(cnv_e, !more) ){
    for( ; i < Ap->lnnz; i++){
      pkg_r.i = i;
      pkg_r.nonz = Ap->lnonzero[i] / THREADS;
      pe = Ap->lnonzero[i] % THREADS;
      if( !convey_push(cnv_r, &pkg_r, pe) )
        break;
    }
    while( convey_pull(cnv_r, &pkg_p, &fromth2) == convey_OK ){ 
      pkg_r.i = pkg_p.i;
      pkg_r.nonz = lcperm[pkg_p.nonz];
      if( !convey_push(cnv_e, &pkg_r, fromth2) ){
        convey_unpull(cnv_r);
        break;
      }
    }
    while( convey_pull(cnv_e, &pkg_p, NULL) == convey_OK ){ 
      Ap->lnonzero[pkg_p.i] = pkg_p.nonz;
    }
  }

  lgp_barrier();
  convey_free(cnv_e);
  convey_free(cnv_r);
  sort_nonzeros(Ap);
  return(Ap);
}

/*! \brief produce the transpose of a sparse matrix using conveyors
 * \param A  pointer to the original matrix
 * \return a pointer to the matrix that has been produced or NULL if the model can't be used
 */
sparsemat_t * transpose_matrix_conveyor(sparsemat_t * A) {
  int64_t ret, pe;
  int64_t lnnz, i, col, row, fromth; 
  int64_t idx, idxp;

  /* get the colcnts */
  int64_t lnumcols = (A->numrows + THREADS - MYTHREAD - 1)/THREADS;  
  int64_t * lcounts = calloc(lnumcols, sizeof(int64_t));
  lgp_barrier();

  int weighted = (A->value != NULL);
  convey_t* cnv_cnt = convey_new(SIZE_MAX, 0, NULL, convey_opt_SCATTER);
  convey_begin(cnv_cnt, sizeof(int64_t), 0);

  lnnz = i = 0;
  while(convey_advance(cnv_cnt, (i == A->lnnz))){
    for( ;i < A->lnnz; i++){
      col = A->lnonzero[i] / THREADS;
      pe = A->lnonzero[i] % THREADS;
      if( !convey_push(cnv_cnt, &col, pe) )
        break;
    }
    while( convey_pull(cnv_cnt, &idxp, NULL) == convey_OK ) {
      lcounts[idxp]++; 
      lnnz++;
    }
  }
  convey_free(cnv_cnt);

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

  
  convey_t* cnv_rd = convey_new(SIZE_MAX, 0, NULL, convey_opt_SCATTER);
  if(weighted)
    convey_begin(cnv_rd, sizeof(w_edge_t), 0);
  else
    convey_begin(cnv_rd, sizeof(edge_t), 0);

  uint64_t numtimespop=0;
  edge_t edge;
  w_edge_t wedge;  
  i = row = 0;
  while(convey_advance(cnv_rd, (i == A->lnnz))){
    for( ; i < A->lnnz; i++){
      while( i == A->loffset[row+1] ) 
        row++;
      pe = A->lnonzero[i]%THREADS;
      if(weighted){
        wedge.row = row * THREADS + MYTHREAD;
        wedge.col = A->lnonzero[i] / THREADS;
        wedge.val = A->lvalue[i];
        if( !convey_push(cnv_rd, &wedge, pe) )
          break;
      }else{
        edge.row = row * THREADS + MYTHREAD;
        edge.col = A->lnonzero[i] / THREADS;
        if( !convey_push(cnv_rd, &edge, pe) )
          break;
      }
    }

    if(weighted){
      while(convey_pull(cnv_rd, &wedge, NULL ) == convey_OK){
        numtimespop++;
        At->lnonzero[ At->loffset[wedge.col] + wrkoff[wedge.col] ] = wedge.row;
        At->lvalue[ At->loffset[wedge.col] + wrkoff[wedge.col] ] = wedge.val;
        wrkoff[wedge.col]++;
      }
    }else{
      while(convey_pull(cnv_rd, &edge, NULL ) == convey_OK){
        numtimespop++;
        At->lnonzero[ At->loffset[edge.col] + wrkoff[edge.col] ] = edge.row;
        wrkoff[edge.col]++;
      }
    }
  }
  
  lgp_barrier();
  convey_free(cnv_rd);
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

