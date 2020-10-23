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

/*! \file spmat_io.upc
 * \brief Utilities to support reading and writing distributed sparse matrices in parallel.
 */
#include <spmat.h>
#include <sys/stat.h>   // for mkdir()
#include <fcntl.h>
#include <exstack.h>


/*! \brief This function reads the metadata file for a sparse matrix dataset.
 * The metadata file is called "metadata" and contains 5 lines
 * - numrows
 * - numcols
 * - number of nonzeros
 * - number of PEs that are writing this dataset.
 * - 0/1 if the matrix has values
 *
 * \param dirname The name of the directory where the dataset will be written.
 * \return A spmat_dataset_t struct.
 * \ingroup spmatgrp
 */

spmat_dataset_t * read_sparse_matrix_metadata(char * dirname){

  spmat_dataset_t * spd = calloc(1, sizeof(spmat_dataset_t));

  spd->dirname = calloc(strlen(dirname)+1, sizeof(char));
  strcpy(spd->dirname, dirname);
  
  if(MYTHREAD == 0){
    int ret = 0;
    char fname[64];
    sprintf(fname, "%s/metadata", dirname);
    FILE * fp = fopen(fname, "r");
    ret += fscanf(fp, "numrows: %"PRId64"\n", &spd->numrows);
    ret += fscanf(fp, "numcols: %"PRId64"\n", &spd->numcols);
    ret += fscanf(fp, "nnz: %"PRId64"\n", &spd->nnz);
    ret += fscanf(fp, "nwriters: %"PRId64"\n", &spd->nwriters);
    ret += fscanf(fp, "values: %d\n", &spd->values);
    if(ret != 5){
      fprintf(stderr,"ERROR: read_sparse_matrix_metadata\n");
      spd->numrows = -1;
    }
    fclose(fp);
  }  

  lgp_barrier();

  spd->numrows = lgp_reduce_add_l(spd->numrows);
  spd->numcols = lgp_reduce_add_l(spd->numcols);
  spd->nnz     = lgp_reduce_add_l(spd->nnz);
  spd->nwriters= lgp_reduce_add_l(spd->nwriters);
  spd->values  = lgp_reduce_add_l(spd->values);  

  if(spd->numrows < 0) return(NULL);
  return(spd);

}

/*! \brief This function writes out the metadata file for a sparse matrix dataset.
 * The metadata file is called "metadata" and contains 5 lines
 * - numrows
 * - numcols
 * - number of nonzeros
 * - number of PEs that are writing this dataset.
 * - 0/1 if the matrix has values
 *
 * \param spd The spmat_dataset_t obtained from open_sparse_matrix_write.
 * \return 0 on success
 * \ingroup spmatgrp
 */

int write_sparse_matrix_metadata(spmat_dataset_t * spd){
  if(MYTHREAD == 0){
    /* open params file for writing */
    char fname[64];
    sprintf(fname, "%s/metadata", spd->dirname);
    FILE * fp = fopen(fname, "w");
    fprintf(fp, "numrows: %"PRId64"\n", spd->numrows);
    fprintf(fp, "numcols: %"PRId64"\n", spd->numcols);
    fprintf(fp, "nnz: %"PRId64"\n", spd->nnz);
    fprintf(fp, "nwriters: %"PRId64"\n", spd->nwriters);
    fprintf(fp, "values: %d\n", spd->values);
    fclose(fp);
  }
  return(0);
  
}

/* \brief Open a sparse matrix dataset for writing.
 * 
 * \param dirname The datadir where the dataset will be written.
 * \param A The matrix in memory that will be written.
 * \return A spmat_dataset_t struct.
 * \ingroup spmatgrp
 */
spmat_dataset_t * open_sparse_matrix_write(char * dirname, sparsemat_t * A){

  spmat_dataset_t * spd = calloc(1, sizeof(spmat_dataset_t));

  spd->nwriters = THREADS; // maybe this should be a parameter like in reads?
  spd->lnumrows = A->lnumrows;
  spd->numrows = A->numrows;
  spd->numcols = A->numcols;
  spd->nnz = A->nnz;
  spd->values = (A->value != NULL);
  spd->io_rank = MYTHREAD;
  spd->global_first_row =  MYTHREAD * spd->lnumrows;
  if(MYTHREAD >= A->numrows % spd->nwriters)
    spd->global_first_row +=  A->numrows % spd->nwriters;
  spd->dirname = calloc(strlen(dirname)+1, sizeof(char));
  strcpy(spd->dirname, dirname);

  /* create the directory  */
  mkdir(dirname, 0770);
  
  /* write metadata file */
  int ret = write_sparse_matrix_metadata(spd);
  assert(ret == 0);
  
  return(spd);
  
}

/* \brief Open a sparse matrix dataset for reading.
 * 
 * \param The datadir where the dataset lives.
 * \param nreaders The number of PEs that will actually read from disk.
 * \return A spmat_dataset_t struct.
 * \ingroup spmatgrp
 */
spmat_dataset_t * open_sparse_matrix_read(char * dirname, int64_t nreaders){
  int64_t i;

  /* read the metadata file */
  spmat_dataset_t * spd = read_sparse_matrix_metadata(dirname);
  assert(spd != NULL);

  /* finish initializing spd */
  /* figure out which PEs are going to participate in reading data
     how many rows they will read and which row that start with. 
  */
  int64_t d;
  if(nreaders > 1)
    d = (THREADS-1)/(nreaders - 1);
  else
    d = 0;  
  spd->io_rank = -1;
  spd->lnumrows = 0;
  spd->global_first_row = -1L;
  spd->global_first_row_to_me = calloc(THREADS, sizeof(int64_t));
  for(i = 0; i < nreaders; i++){

    /* figure out this reader's first row they will read. */
    spd->global_first_row_to_me[d*i] = i*((spd->numrows + nreaders - i - 1)/nreaders);
    if(i >= spd->numrows % nreaders)
      spd->global_first_row_to_me[d*i] += spd->numrows % nreaders;
    /* incrememnt until I get to the first row this reader will send me. */
    while((spd->global_first_row_to_me[d*i] % THREADS) != MYTHREAD)
      spd->global_first_row_to_me[d*i]++;

    /* if I am a reader populate io_rank, lnumrows, and global_first_row */
    if(MYTHREAD == d*i){
      spd->io_rank = i;
      spd->lnumrows = (spd->numrows + nreaders - i - 1)/nreaders;
      spd->global_first_row = i*spd->lnumrows;
      if(i >= spd->numrows % nreaders)
        spd->global_first_row += spd->numrows % nreaders;
    }
    
  }
  //fprintf(stderr,"PE %d: nreaders %ld io_rank = %ld globalfr %ld lnumrows = %ld\n", MYTHREAD, nreaders, spd->io_rank, spd->global_first_row, spd->lnumrows);

  return(spd);
}


/*! \brief This function reads the row_info files for the dataset.
 * 
 * After this function we will have the row counts of every row we are
 * responsible for reading. We will also have the number of files we
 * are going to read from, the first file we are going to read from
 * and then place to seek to in that first file.
 */
void read_row_info(spmat_dataset_t * spd){
  int64_t i, j;
  
  if(spd->io_rank >= 0){
    
    /* figure out which file to start reading your rowcnts from 
        how many files you will need to read from 
        and where to seek to in the first file
    */
    int64_t current_nrows = 0;
    int64_t num_recs_this_file;
    int64_t current_file;
    int64_t stop_row = spd->global_first_row + spd->lnumrows;
    spd->first_file = -1;
    spd->nfiles = 0;
    for(current_file = 0; current_file < spd->nwriters; current_file++){

      num_recs_this_file = (spd->numrows + spd->nwriters - current_file - 1)/spd->nwriters;

      if(current_nrows + num_recs_this_file > spd->global_first_row){

        if(spd->first_file < 0){
          // we just found our first file
          spd->first_file = current_file;
          spd->start_row_first_file = spd->global_first_row - current_nrows;
        }
        spd->nfiles++;
        if( current_nrows + num_recs_this_file >= stop_row){
          // this is our last file
          break;
        }
      }
      
      current_nrows += num_recs_this_file;
      
    }    
    
    spd->current_file = spd->first_file; // this sets us up to start reading nonzeros

    //fprintf(stderr,"PE%d: first_file %ld sr %ld gfr %ld\n", MYTHREAD, spd->first_file,
    //spd->start_row_first_file, spd->global_first_row);

    /* figure out how many rows in each file this reader will need to read */
    spd->nrows_in_file      = calloc(spd->nfiles, sizeof(int64_t));
    spd->nrows_read_in_file = calloc(spd->nfiles, sizeof(int64_t));
    spd->rowcnt             = calloc(spd->lnumrows, sizeof(int64_t));
    
    if(spd->nfiles == 1)
      spd->nrows_in_file[0] = spd->lnumrows;
    else{
      num_recs_this_file = (spd->numrows + spd->nwriters - spd->first_file - 1)/spd->nwriters;
      spd->nrows_in_file[0] = num_recs_this_file - spd->start_row_first_file;
      int64_t total = spd->nrows_in_file[0];
      for(i = spd->first_file + 1; i < spd->first_file + spd->nfiles; i++){
        num_recs_this_file = (spd->numrows + spd->nwriters - i - 1)/spd->nwriters;
        spd->nrows_in_file[i - spd->first_file] = num_recs_this_file;
        //fprintf(stderr,"PE %d: nrif %ld\n", MYTHREAD, spd->nrows_in_file[i - spd->first_file]);
        total += spd->nrows_in_file[i - spd->first_file];
      }
      // we might not need to read all of last file
      spd->nrows_in_file[spd->nfiles - 1] -= total - spd->lnumrows;

    }
    
  }

  /* For each file this PE is going to read from, read all of the data
   * we are supposed to read (nrows_in_file) and convert that data
   * into rowcounts. 
   *
   * Also record the offset of our first nonzero we are to read.
  */
  int64_t num_to_read, rows_read = 0;
  for(i = 0; i < spd->nfiles; i++){

    int64_t current_file = spd->first_file + i;
    //fprintf(stderr,"PE %d will read %ld rows from file %ld\n", MYTHREAD, spd->nrows_in_file[i], current_file);

    /* open the first rowcnt file you will read from */
    char fname[128];
    sprintf(fname, "%s/rowinfo_%" PRId64, spd->dirname, current_file);

    FILE * fp = fopen(fname, "rb");
    assert(fp != NULL);

    if(i == 0 && spd->start_row_first_file){
      /* seek to the right place in the first file */
      fseek(fp, spd->start_row_first_file*sizeof(int64_t), SEEK_SET);
    }

    /* read the rowinfo into the spd->rowcnt array */
    int64_t num_read = fread(&spd->rowcnt[rows_read], sizeof(uint64_t),
                             spd->nrows_in_file[i], fp);
    assert(num_read == spd->nrows_in_file[i]);

    
    if(i == 0){  // remember the offset of the first nonzero in first file
      spd->first_file_offset = spd->rowcnt[rows_read];
    }
    
    // turn nonzero offsets into rowcounts.
    for(j = 0; j < num_read - 1; j++){
      //fprintf(stderr,"PE%d before: rc[%ld] = %ld\n", MYTHREAD, rows_read+j, spd->rowcnt[rows_read + j]);
      spd->rowcnt[rows_read+j] = spd->rowcnt[rows_read + j + 1] - spd->rowcnt[rows_read + j];
      //T0_fprintf(stderr,"PE %d after: rc[%ld] = %ld\n", MYTHREAD, rows_read+j, spd->rowcnt[rows_read + j]);
    }

    // handle the last row that you read. If it is the not the last
    // row in the file, you can just read the next offset, otherwise
    // you will have to figure out its rowcnt from the file size.
    int64_t tmp, last_num_read;
    last_num_read = fread(&tmp, sizeof(uint64_t), 1, fp);
    if(last_num_read){
      // there are still more rows in this file, read one more so that we can get the rowcnt of our last row
      spd->rowcnt[rows_read + num_read - 1] = tmp - spd->rowcnt[rows_read + num_read - 1];
      fclose(fp);
    }else{
      // else... we need to know the size of the nonzero file to get last rowcount
      sprintf(fname, "%s/nonzero_%" PRId64, spd->dirname, current_file);
      int fd = open(fname, O_RDONLY);
      struct stat buf;
      fstat(fd, &buf);
      spd->rowcnt[rows_read + num_read - 1] = (buf.st_size/8) - spd->rowcnt[rows_read + num_read - 1];
      close(fd);
    }
    //T0_fprintf(stderr,"PE %d after: rc[%ld] = %ld\n",MYTHREAD,rows_read+num_read-1,spd->rowcnt[rows_read + num_read - 1]);
    rows_read += num_read;

  }
  
  assert(rows_read == spd->lnumrows);
  
  lgp_barrier();
  
}


/*! \brief Write row_info files. These files collectively
 * have one line for each row in the matrix. Each writer writes a
 * unique file containing the offsets into the nonzero file for 
 * the rows this PE is responsible for.
 * \param spd An initialized spmat_dataset_t.
 * \param A  A pointer to the matrix
 * \ingroup spmatgrp
 */

void write_row_info(spmat_dataset_t * spd, sparsemat_t * A){

  int64_t i;

  /* get the rowcounts of the matrix into a distributed array */
  SHARED int64_t * rowcnt = lgp_all_alloc(A->numrows, sizeof(int64_t));
  int64_t * lrowcnt = lgp_local_part(int64_t, rowcnt);

  for(i = 0; i < A->lnumrows; i++)
    lrowcnt[i] = A->loffset[i + 1] - A->loffset[i];
  lgp_barrier();

  int64_t error = 0;  
  if(spd->lnumrows){

    char fname[128];
    sprintf(fname, "%s/rowinfo_%d", spd->dirname, MYTHREAD);
    FILE * fp = fopen(fname, "wb");
    assert(fp != NULL);
    
    int64_t offset = 0;
    for(i = 0; i < spd->lnumrows; i++){
      if( fwrite(&offset, sizeof(int64_t), 1,  fp) != 1 ){
        error = 1;
        break;
      }
      int64_t rc = lgp_get_int64(rowcnt, spd->global_first_row + i);
      offset += rc;
      
    }
    
    fclose(fp);
  }
  error = lgp_reduce_add_l(error);
  assert(error == 0);

  lgp_all_free(rowcnt);

}




/*! \brief This routine reads the next block of nonzeros of complete rows from the current file.
 * 
 * It will read up to buf_size nonzeros. It returns the total number of rows for which the
 * nonzeros were read and puts the nonzeros in buf (and values in vbuf in the case of
 * a matrix with values). Both buf and vbuf must have been previously allocated.
 *
 * \param spd The spmat_dataset_t obtained from calling open_sparse_matrix_read().
 * \param buf A buffer to hold the read nonzeros.
 * \param vbuf A buffer to hold the read values (if the matrix has values).
 * \param The length of these buffers.
 * \return The number of rows that were read.
 */
int64_t read_nonzeros_buffer(spmat_dataset_t * spd, int64_t * buf,
                             double * vbuf, int64_t buf_size){

  //fprintf(stderr,"PE %d: cf %ld first %ld nfiles %ld\n",
  // MYTHREAD, spd->current_file, spd->first_file, spd->nfiles);
  if(spd->current_file >= spd->first_file + spd->nfiles)
    return(-1);
  
  if(spd->file_open == 0){
    // we need to open the current file
    spd->current_row_in_file = 0;
    char fname[128];
    sprintf(fname, "%s/nonzero_%" PRId64, spd->dirname, spd->current_file);
    spd->nnzfp = fopen(fname,"rb");
    if(spd->values){
      sprintf(fname, "%s/value_%" PRId64, spd->dirname, spd->current_file);
      spd->valfp = fopen(fname,"rb");
    }
    spd->file_open = 1;
    if(spd->current_file == spd->first_file){
      spd->current_global_row = 0;
      fseek(spd->nnzfp, spd->first_file_offset*sizeof(int64_t), SEEK_SET);
      //fprintf(stderr,"PE %d: seeking to %ld in nonzero file %ld\n", MYTHREAD, spd->first_file_offset, spd->current_file);
      if(spd->values) fseek(spd->valfp, spd->first_file_offset*sizeof(double), SEEK_SET);
    }
  }
  
  /* figure out how many nonzeros you can read from this file into this buffer */
  int64_t current_buf_cnt = 0;
  int64_t num_rows = 0;
  int64_t row;
  //fprintf(stderr,"PE %d: current_row = %ld rows in file %ld\n", MYTHREAD, spd->current_row_in_file, spd->nrows_in_file[spd->current_file - spd->first_file]);  
  for(; spd->current_row_in_file < spd->nrows_in_file[spd->current_file - spd->first_file]; spd->current_row_in_file++){
    if(spd->rowcnt[spd->current_global_row] + current_buf_cnt > buf_size)
      break;
    current_buf_cnt += spd->rowcnt[spd->current_global_row];
    spd->current_global_row++;
    num_rows++;
  }
  //fprintf(stderr,"PE %d: numrows to read %ld current_buf_cnt = %ld current row = %ld crif %ld current_file %ld\n", MYTHREAD, num_rows, current_buf_cnt, spd->current_global_row, spd->current_row_in_file, spd->current_file);
  
  /* read nonzero and value data into buffers */
  int64_t num = fread(buf, sizeof(int64_t), current_buf_cnt, spd->nnzfp);
  assert(num == current_buf_cnt);
  if(spd->values){
    num = fread(vbuf, sizeof(double), current_buf_cnt, spd->valfp);
    assert(num == current_buf_cnt);
  }

  if(spd->current_row_in_file == spd->nrows_in_file[spd->current_file - spd->first_file]){
    // we just finished this file
    fclose(spd->nnzfp);
    if(spd->values) fclose(spd->valfp);
    spd->file_open = 0;
    spd->current_file++;
  }

  return(num_rows);
}



/*! \brief writes a sparse matrix to a file in Matrix Market format
 * \param A pointer to the sparse matrix
 * \param name the filename to written to
 * \return 0 on success, non-0 on error.
 * \ingroup spmatgrp
 */
int write_matrix_mm(sparsemat_t *A, char * name) {
  if(!MYTHREAD){
    FILE * fp = fopen(name, "w");
    /* write the banner */
    if(A->value){
      fprintf(fp,"%%%%MatrixMarket matrix coordinate real\n");
    }else{
      fprintf(fp,"%%%%MatrixMarket matrix coordinate pattern\n");
    }
    fprintf(fp,"%"PRId64" %"PRId64" %"PRId64"\n", A->numrows, A->numcols, A->nnz);
    fclose(fp);
  }
  
  lgp_barrier();
  
  int64_t i, j,k, row;
  for(k = 0; k < THREADS; k++){
    if(k == MYTHREAD){
      FILE * fp = fopen(name, "a");
      for(i = 0; i < A->lnumrows; i++){
        row = i*THREADS + MYTHREAD;
        for(j = A->loffset[i]; j < A->loffset[i+1]; j++){
          if(A->value)
            fprintf(fp, "%"PRId64" %"PRId64" %lf\n",
                    row + 1,
                    A->lnonzero[j] + 1,
                    A->value[j]);
          else
            fprintf(fp, "%"PRId64" %"PRId64"\n",
                    row + 1,
                    A->lnonzero[j] + 1);
        }
      }
      fclose(fp);
    }
    lgp_barrier();
  }
  
  return(0);
}


/*! \brief Read a sparse matrix in matrix market format on one PE and create a distributed matrix
  from that.
  * Only PE 0 reads the matrix file.
  * 
  * \param name The name of the file.
  * \return The sparsemat_t struct.
  * \ingroup spmatgrp
  */
sparsemat_t * read_matrix_mm_to_dist(char * name) {
  typedef struct pkg_rowcol_t{
    int64_t row;    
    int64_t col;
  }pkg_rowcol_t;

  int64_t nr, nc, nnz = 0, i, pe;
  SHARED int64_t * sh_data;
  sh_data = lgp_all_alloc (THREADS*4, sizeof(int64_t));

  int64_t * rowcount;
  edge_t * edges;
  w_edge_t * tri;
  if(!MYTHREAD){
    int fscanfret;
    int64_t * nnz_per_th = calloc(THREADS, sizeof(int64_t));
    
    FILE * fp = fopen(name, "r");
    if( fp == NULL ) {
      fprintf(stderr,"read_matrix_mm: can't open file %s \n", name);
      lgp_global_exit(1);
    }
    
    // Read the header line of the MatrixMarket format 
    char * object = calloc(64, sizeof(char));
    char * format = calloc(64, sizeof(char));
    char * field = calloc(64, sizeof(char));;
    fscanfret = fscanf(fp,"%%%%MatrixMarket %s %s %s\n", object, format, field);
    if( (fscanfret != 3 ) || strncmp(object,"matrix",24) || strncmp(format,"coordinate",24) ){
      fprintf(stderr,"read_matrix_mm: Incompatible matrix market format.\n");
      fprintf(stderr,"                First line should be either:\n");
      fprintf(stderr,"                matrix coordinate pattern\n");
      fprintf(stderr,"                OR\n");
      fprintf(stderr,"                matrix coordinate real\n");
      fprintf(stderr,"                OR\n");
      fprintf(stderr,"                matrix coordinate integer\n");
      lgp_global_exit(1);
    }

    // Make sure that this is a format we support
    if(strncmp(field,"pattern",24) && strncmp(field,"real",24) && strncmp(field,"integer",24) ){
      fprintf(stderr,"read_matrix_mm: Incompatible matrix market field.\n");
      fprintf(stderr,"                Last entry on first line should be pattern, real, or integer\n");
      lgp_global_exit(1);
    }
    int64_t values;
    if(strncmp(field,"pattern",7) == 0){
      values = 0L; // no values
    }else if(strncmp(field,"real",4) == 0){
      values = 1L; // real values
    }else{
      values = 2L; // integer values
    }
    
    // Read the header (nr, nc, nnz)
    fscanfret = fscanf(fp,"%"PRId64" %"PRId64" %"PRId64"\n", &nr, &nc, &nnz);
    if( (fscanfret != 3 ) || (nr<=0) || (nc<=0) || (nnz<=0) ) {
      fprintf(stderr,"read_matrix_mm: reading nr, nc, nnz\n");
      lgp_global_exit(1);
    }

    // allocate space to store the matrix data    
    rowcount = calloc(nr, sizeof(int64_t));
    if(!rowcount){
      T0_printf("ERROR: read_matrix_mm_to_dist: could not allocate arrays\n");
      for(i = 0; i < THREADS; i++) lgp_put_int64(sh_data, i, -1);
    }
    
    // read the data
    int64_t row, col, val, pos = 0;
    if(values == 0){
      edges = calloc(nnz, sizeof(edge_t));
      while(fscanf(fp,"%"PRId64" %"PRId64"\n", &row, &col) != EOF){
        row--;//MM format is 1-up
        col--;
        edges[pos].row   = row;
        edges[pos++].col = col;
        nnz_per_th[row % THREADS]++;
        rowcount[row]++;
      }
      qsort( edges, nnz, sizeof(edge_t), edge_comp);
    }else{
      tri = calloc(nnz, sizeof(w_edge_t));    
      while(fscanf(fp,"%"PRId64" %"PRId64" %"PRId64"\n", &row, &col, &val) != EOF){
        tri[pos].row = row - 1;
        tri[pos].col = col - 1;
        tri[pos++].val = val;
        nnz_per_th[row % THREADS]++;
        rowcount[row]++;
      }
      qsort( tri, nnz, sizeof(w_edge_t), w_edge_comp);
    }
    
    fclose(fp);
    if(nnz != pos){
      T0_printf("ERROR: read_matrix_mm_to_dist: nnz (%"PRId64") != pos (%"PRId64")\n", nnz, pos);
      for(i = 0; i < THREADS; i++) lgp_put_int64(sh_data, i, -1);
    }
    for(i = 0; i < THREADS; i++){
      lgp_put_int64(sh_data, i, nnz_per_th[i]);
      lgp_put_int64(sh_data, i+THREADS, nr);
      lgp_put_int64(sh_data, i+2*THREADS, nc);
      lgp_put_int64(sh_data, i+3*THREADS, values);
    }
    free(nnz_per_th);

  }
  
  lgp_barrier();

  int64_t * lsh_data = lgp_local_part(int64_t, sh_data);
  if(lsh_data[0] == -1)
    return(NULL);
  
  int64_t lnnz = lsh_data[0];
  nr = lsh_data[1];
  nc = lsh_data[2];
  int value = (lsh_data[3] != 0L);
  
  sparsemat_t * A = init_matrix(nr, nc, lnnz, value);
  SHARED int64_t * tmp_offset = lgp_all_alloc(nr + THREADS, sizeof(int64_t));
  if(!A || !tmp_offset){
    T0_printf("ERROR: read_matrix_mm_to_dist: failed to init matrix or tmp_offset!\n");
    return(NULL);
  }

  /* set up offset array and tmp_offset */
  lgp_barrier();
  lgp_all_free(sh_data);
  
  if(!MYTHREAD){
    for(i = 0; i < nr; i++)
      lgp_put_int64(tmp_offset, i, rowcount[i]);
    free(rowcount);
  }

  lgp_barrier();

  int64_t * ltmp_offset = lgp_local_part(int64_t, tmp_offset);
  A->loffset[0] = 0;
  for(i = 1; i <= A->lnumrows; i++){
    A->loffset[i] = A->loffset[i-1] + ltmp_offset[i-1];
    ltmp_offset[i-1] = 0;
  }

  int64_t fromth;
  w_edge_t pkg;
  exstack_t * ex = exstack_init(256, sizeof(w_edge_t));
  if( ex == NULL ) return(NULL);
  
  /* distribute the matrix to all other PEs */
  /* pass around the nonzeros */
  /* this is a strange exstack loop since only PE0 has data to push */
  i = 0;
  while(exstack_proceed(ex, (i == nnz))){
    while(i < nnz){
      if(value == 0){
        pkg.row = edges[i].row;
        pkg.col = edges[i].col;
      }else{
        pkg.row = tri[i].row;
        pkg.col = tri[i].col;
        pkg.val = tri[i].val;
      }
      pe = pkg.row % THREADS;
      if(!exstack_push(ex, &pkg, pe))
        break;
      i++;
    }
    exstack_exchange(ex);

    while(exstack_pop(ex, &pkg, &fromth)){
      int64_t row = pkg.row/THREADS;
      int64_t pos = A->loffset[row] + ltmp_offset[row];
      //printf("pos = %ld row = %ld col = %ld\n", pos, row, pkg.col);fflush(0);
      A->lnonzero[pos] = pkg.col;
      if(value) A->lvalue[pos] = pkg.val;
      ltmp_offset[row]++;
    }
  }

  lgp_barrier();
  if(!MYTHREAD){
    if(value == 0)
      free(edges);
    else
      free(tri);
  }

  lgp_all_free(tmp_offset);
  exstack_clear(ex);
  sort_nonzeros(A);
  return(A);
}
