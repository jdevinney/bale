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

/*! \file geometric.upc
 * \brief Functions to support generating random geometric graphs.
 */
#include <spmat.h>
#include <exstack.h>

/*! \brief maximum number of sectors across all threads */
int64_t sector_max;

/*! 
\brief compute the square of the distance between two points
\param a the first point
\param b the other point
\return the square 
*/
double dist(point_t a, point_t b){
  return((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y));
}

/*! 
\brief This function computes all edges between two sectors and added them to the edge list
\param this_sec_idx index of my sector
\param other_sec_idx index of the other sector
\param my_first_sector the first sector on my thread
\param nsectors number of sectors
\param points array of points
\param first_point_in_sector the first point in the sectors
\param counts the number of
\param r2 square of the radius
\param loops whether or not we have loops
\param el the edge list
*/
void append_edges_between_sectors(uint64_t this_sec_idx,
                                  uint64_t other_sec_idx,
                                  int64_t my_first_sector,
                                  int64_t nsectors,
                                  SHARED point_t * points,
                                  SHARED int64_t * first_point_in_sector,
                                  SHARED int64_t * counts,
                                  double r2,
                                  self_loops loops,
                                  edge_list_t * el){
  
  assert(this_sec_idx < nsectors);
  assert(other_sec_idx < nsectors);
  edge_t e;
  int64_t j, l;
  point_t tmp;
  int64_t lnsectors = (nsectors + THREADS - MYTHREAD - 1)/THREADS;
  int64_t * lcounts = lgp_local_part(int64_t, counts);
  int64_t * lfirst_point_in_sector = lgp_local_part(int64_t, first_point_in_sector);
  point_t * lpoints = lgp_local_part(point_t, points);
  
  //int64_t other_sec_pe = other_sec_idx % THREADS;
  int64_t this_sec_local_idx = this_sec_idx - my_first_sector;
  point_t * lpoints_this = &lpoints[this_sec_local_idx*sector_max];
  int64_t first_point_this = lfirst_point_in_sector[this_sec_local_idx];
  int64_t npts_this = lcounts[this_sec_local_idx];

  if((other_sec_idx >= my_first_sector) &&
     (other_sec_idx < my_first_sector + lnsectors)){

    // We are finding edges between points within two local sectors.
    int64_t other_sec_local_idx = other_sec_idx - my_first_sector;
    point_t * lpoints_other = &lpoints[other_sec_local_idx*sector_max];
    int64_t first_point_other = lfirst_point_in_sector[other_sec_local_idx];
    int64_t npts_other = lcounts[other_sec_local_idx];
    int64_t stop_pt = npts_other;
    
    for(j = 0; j < npts_this; j++){
      int64_t point_indexA = first_point_this + j;
      // If we are looking within the same sector, only compare
      // points where l < j.
      for(l = 0; l < npts_other; l++){

        // skip self loops if not specified
        if((other_sec_idx == this_sec_idx)
           && (loops != LOOPS)
           && (l == j))
          continue;
        
        int64_t point_indexB = first_point_other + l;
        // assert(point_indexB < 10);
        if(dist(lpoints_this[j], lpoints_other[l]) < r2){
          append_edge(el, point_indexA, point_indexB);          
        }
      }
    }
  }else{
    /* the other sector is not on our PE */
    int64_t other_sec_pe;
    int64_t other_sec_local_idx;
    //int64_t other_sec_local_idx = global_index_to_pe_and_offset(other_sec_idx, nsectors, &other_sec_pe, BLOCK);
    global_index_to_pe_and_offset(&other_sec_pe, &other_sec_local_idx, other_sec_idx, nsectors, BLOCK);
    int64_t first_point_other = lgp_get_int64(first_point_in_sector, other_sec_local_idx*THREADS + other_sec_pe);
    int64_t num_pts_other = lgp_get_int64(counts, pe_and_offset_to_global_index(other_sec_pe, other_sec_local_idx, nsectors, CYCLIC));
    if(num_pts_other == 0)
      return;
    /* get all the points from the other sector in one transfer */
    point_t * other_sec_pts = calloc(num_pts_other, sizeof(point_t));
    int64_t sector_start = sector_max*other_sec_local_idx*THREADS + other_sec_pe;
    lgp_memget(other_sec_pts, points, num_pts_other*sizeof(point_t), sector_start);

    for(j = 0; j < npts_this; j++){
      int64_t point_indexA = first_point_this + j;
      for(l = 0; l < num_pts_other; l++){
        int64_t point_indexB = first_point_other + l;
        if(dist(lpoints_this[j], other_sec_pts[l]) < r2){
          append_edge(el, point_indexA, point_indexB);        
        }
      }
    }
    free(other_sec_pts);
  }
}


/*!
\brief This routine appends to the edge_list (el) all edges for points in a given sector.
\param sector the sector
\param my_first_sector the first sector on my thread
\param nsectors number of sectors
\param points array of points
\param first_point_in_sector the first point in the sectors
\param counts the number of
\param r the radius
\param loops whether or not we have loops
\param el the edge list

It does this by looking for point pairs that are closer than r, where one point lies
in "sector" and the other point lies in any of the 8 neighboring sectors or "sector" itself.
*/ 
void append_edges_for_sector(int64_t sector,
                             int64_t my_first_sector,
                             int64_t nsectors,
                             SHARED point_t * points,
                             SHARED int64_t * first_point_in_sector,
                             SHARED int64_t * counts,
                             double r,
                             self_loops loops,
                             edge_list_t * el){

  int64_t nsectors_across = floor(1.0/r);
  double r2 = r*r;
  int64_t sector_row = sector % nsectors_across;
  int64_t sector_col = sector / nsectors_across;
  int64_t i, j;

  // Append edges for between this sector and other sectors
  // We only need to look at all neighboring sectors
  // (including diagonals) AND the sector itself.
  for(i = -1; i < 2; i++){
    for(j = -1; j < 2; j++){
      int64_t nbrsector_col = sector_col + i;
      int64_t nbrsector_row = sector_row + j;
      if(nbrsector_col < 0 || nbrsector_col == nsectors_across)
        continue;
      if(nbrsector_row < 0 || nbrsector_row == nsectors_across)
        continue;
      int64_t nbr_sector_index = (nbrsector_col)*nsectors_across + (nbrsector_row);
      append_edges_between_sectors(sector, nbr_sector_index,
                                   my_first_sector, nsectors,
                                   points, first_point_in_sector,
                                   counts, r2, loops, el);      
    }
  }
  
}
                             


/*!
\brief produces a shared array containing the prefix scan of the given array
\param counts the given array
\param n the length of the array
\return a pointer to the new array

This is the Hillis and Steele algorithm as presented on wikipedia.
*/
SHARED int64_t * prefix_scan(SHARED int64_t * counts, int64_t n){
  int64_t i, j;
  SHARED int64_t * cumsum = lgp_all_alloc(n, sizeof(int64_t));
  int64_t * lcumsum = lgp_local_part(int64_t, cumsum);
  for(i = 0; i < (n + THREADS - MYTHREAD - 1)/THREADS; i++)
    lcumsum[i] = 0;

  lgp_barrier();
  
  int64_t tn = n;
  int64_t lgn = 0;
  while(tn > 0){
    tn = tn >> 1;
    lgn++;
  }
  for(i = 0; i < lgn; i++){
    for(j = 0; j < n; j++){
      int64_t twoj = 1L<<j;
      if((j >= twoj) && ((j % n) == MYTHREAD)){
        lcumsum[j] += lgp_get_int64(cumsum, j-twoj);
      }
    }
    lgp_barrier();
  }
  return(cumsum);
}

/*!
\brief Generates the adjacency matrix for a random geometric graph. 
\param n The number of vertices
\param r The size of the neighborhoods that determine edges.
\param edgetype See enum edge_type. Must be UNDIRECTED or UNDIRECTED_WEIGHTED.
\param loops See self_loops. Are there self loops?
\param seed A seed for the RNG. This should be a single across all PEs (it will be modified by each PE individually).
\param out_points (Optional) If you supply this pointer, the routine will populate it with the points associated with the vertices in the graph. 
\return An adjacency matrix (or lower portion of it in the undirected case).

See https://en.wikipedia.org/wiki/Random_geometric_graph
Each vertex corresponds to a point randomly placed in the unit square. Two vertices
are adjancent if their corresponding points are within distance r of each other.
*/
sparsemat_t * geometric_random_graph(int64_t n, double r, edge_type edgetype, self_loops loops, uint64_t seed, SHARED point_t ** out_points){
  
  // We generate n points (p_i) uniformly at random over the unit square.
  // Each point corresponds to a vertex (v_i).
  // There is an edge between v_A and v_B if the distance between the p_A and p_B is less than r.
  // To do this in parallel, we break up the unit square into chunks that are rxr grid.  
  // Each processor is responsible for some number of chunks.
  // We calculate edges by comparing distances between each point and every other point in
  // its own sector and in neighboring sectors.

  
  int64_t i, j, k, l;
  double r2 = r*r;
  int64_t nsectors_across = floor(1.0/r);
  double sector_width = 1.0/nsectors_across;
  int64_t nsectors = nsectors_across*nsectors_across;
  int64_t lnsectors = (nsectors + THREADS - MYTHREAD - 1)/THREADS;
  int64_t ln = (n + THREADS - MYTHREAD - 1)/THREADS;
  int weighted = (edgetype == DIRECTED_WEIGHTED || edgetype == UNDIRECTED_WEIGHTED);

  if((edgetype == DIRECTED) || (edgetype == DIRECTED_WEIGHTED)){
    T0_fprintf(stderr,"Error!: geometric random graphs can only be undirected!");
    return(NULL);
  }
  
  //T0_printf("GEOMETRIC GRAPH: r = %lf number of sectors = %ld sector_width = %lf\n",
  //          r, nsectors, sector_width);
  //T0_printf("                 edgetype = %d loops = %d\n", edgetype, loops);

  
  lgp_rand_seed(seed);

  // FUTURE: permute matrix at end to get Zmorton order (which would improve locality)
  //         or round-robin point order (which would destroy locality)

  // Step 1. Figure out how many points land in every sector. We do
  // this by having each PE generate n/THREADS random sector indices,
  // incrementing the count for each.
  SHARED int64_t * counts = lgp_all_alloc(nsectors, sizeof(int64_t));
  int64_t * lcounts = lgp_local_part(int64_t, counts);

  for(i = 0; i < lnsectors; i++) lcounts[i] = 0;

  lgp_barrier();    
  
  for(i = 0; i < ln; i++) lgp_atomic_add(counts, lgp_rand_int64(nsectors), 1L);
    
  lgp_barrier();

  // Step 2. Figure out how many points landed in your sectors.
  // Note: We are not attaching sectors to anywhere in particular in
  // the unit square yet.
  int64_t my_total_points = 0;
  int64_t my_sector_max = 0;
  for(i = 0; i < lnsectors; i++){
    my_total_points += lcounts[i];
    my_sector_max = (my_sector_max < lcounts[i] ? lcounts[i] : my_sector_max);
  }
  int64_t sum = lgp_reduce_add_l(my_total_points);
  assert(sum == n);
  
  // Step 3. Allocate space for global array of points.
  // sector_max is a global variable ...yuck. 
  sector_max = lgp_reduce_max_l(my_sector_max);
  int64_t max_lnsectors = (nsectors + THREADS - 1)/THREADS; 
  SHARED point_t * points = lgp_all_alloc(sector_max*max_lnsectors*THREADS, sizeof(point_t));
  point_t * lpoints = lgp_local_part(point_t, points);
  
  lgp_barrier();

  // Step 4. We compute the global index of the first point in each sector.
  //
  // We assign sectors to PEs and put an ordering on sectors here. Sectors are assigned to PEs
  // in BLOCK (i.e not CYCLIC) fashion. The sectors will be ordered lexiographically by column and then row
  // in the unit square starting with the sector whose lower left index is (0,0).
  //
  // Points within a sector are given a set of contiguous indices and all
  // points in a sector have lower/higher indices than sectors with higher/lower indices.

  int64_t my_first_sector = 0;
  for(i = 0; i < MYTHREAD; i++)
    my_first_sector += (nsectors + THREADS - i - 1)/THREADS;

  int64_t points_before_my_block = lgp_prior_add_l(my_total_points);
  SHARED int64_t * first_point_in_sector = lgp_all_alloc(max_lnsectors*THREADS, sizeof(int64_t));
  int64_t * lfirst_point_in_sector = lgp_local_part(int64_t, first_point_in_sector);

  lfirst_point_in_sector[0] = points_before_my_block;
  for(i = 1; i < lnsectors; i++){
    lfirst_point_in_sector[i] = lfirst_point_in_sector[i-1] + lcounts[i-1];
  }
  lgp_barrier();


  // Step 5. Each PE generates the points for their own sectors.
  int64_t pt = 0;
  for(i = 0; i < lnsectors; i++){
    int64_t sector = my_first_sector + i;
    double x_off = (sector / nsectors_across) * sector_width;
    double y_off = (sector % nsectors_across) * sector_width;
    pt = i*sector_max;
    for(j = 0; j < lcounts[i]; j++){    
      lpoints[pt].x = lgp_rand_double()*sector_width + x_off;
      lpoints[pt].y = lgp_rand_double()*sector_width + y_off;
      pt++;
    }
    //Sort the points in each sector lexiographically
    //qsort(&lpoints[i*sector_max], lcounts[i], sizeof(point_t), point_comp);
  }

  
  // Step 6. Determine which edges are present
  int64_t space = ceil(1.1*my_total_points*(n*M_PI*r*r));
  // we don't worry about values at this point (thus we are creating an edge list
  // without values.
  edge_list_t * el = init_edge_list(space, 0);
  if(el == NULL){
    printf("ERROR: geometric graph: el is NULL\n");
    return(NULL);
  }

  /* For each of our local sectors we compute all possible edges */ 
  for(i = 0; i < lnsectors; i++){
    append_edges_for_sector(my_first_sector + i, my_first_sector,
                            nsectors, points, first_point_in_sector,
                            counts, r, loops, el);
    
  }

  lgp_barrier();
  
  //printf("PE %d: After comparing inter-point distances: lnnz = %ld\n",MYTHREAD, el->num);
  int64_t nnz = lgp_reduce_add_l(el->num);
  //T0_printf("Symmetric matrix has %ld nonzeros total\n", nnz);
  

  // Step 7. We need to redistribute points to PEs to flatten
  // the number of points per PE. We do this by assigning points
  // to PEs in a BLOCK fashion. Note this is not a huge change from before
  // where we assigned sectors to PEs in block fashion. But sectors
  // do not all have the same number of points, so we need to smooth things out.
  //
  // At the same time, we are preparing to create the adjacency matrix for our graph.
  // Bale matrices are only stored in CYCLIC row layout, so we need to reinterpret the
  // labels on points. So the points that PE i owns will become points i, i + THREADS, i + 2*THREADS, etc.
  //
  // Once we relabel points, we can get rid of the unnecessary edges (those that represent nonzeros
  // above the diagonal).


  int64_t pe;
  int64_t toss = 0, keep = 0;
  for(i = 0; i < el->num; i++){
    // Calling global_index_to_pe_and_offset gives us the PE and offset of a point if distribute the points
    // to PEs evenly and in BLOCK fashion.
    // We then convert that pe and offset into a new global index based on CYCLIC row layout.
    int64_t roffset;
    global_index_to_pe_and_offset(&pe, &roffset, el->edges[i].row, n, BLOCK);
    int64_t row_index = pe_and_offset_to_global_index(pe, roffset, n, CYCLIC);
    int64_t coffset;
    global_index_to_pe_and_offset(&pe, &coffset, el->edges[i].col, n, BLOCK);
    int64_t col_index = pe_and_offset_to_global_index(pe, coffset, n, CYCLIC);
    if ( col_index > row_index) {
      el->edges[i].row = -1; // mark as NULL since this edge represents a 1 above the diagonal in the adjacency matrix.
      toss++;
    }else{
      keep++;
      el->edges[i].row = row_index;
      el->edges[i].col = col_index;
    }    
  }

  // If the calling function wanted the geometric points too, we fill in that array now.
  if(out_points){
    printf("Hi!\n");
    int64_t point_index = points_before_my_block;
    SHARED point_t * op = lgp_all_alloc(n, sizeof(point_t));
    lgp_barrier();
    
    for(i = 0; i < lnsectors; i++){
      for(j = 0; j < lcounts[i]; j++){
        int64_t off;
        global_index_to_pe_and_offset(&pe, &off, point_index, n, BLOCK);
        int64_t new_index = pe_and_offset_to_global_index(pe, off, n, CYCLIC);
        lgp_memput(op, &lpoints[i*sector_max + j], sizeof(point_t), new_index);
        point_index++;
      }
    }
    *out_points = op;
  }
  lgp_barrier();
  lgp_all_free(points);
  
  // Step 8. Do a histogram to get vertex degrees (aka row_counts)
  SHARED int64_t * row_counts = lgp_all_alloc(n, sizeof(int64_t));
  int64_t * lrow_counts = lgp_local_part(int64_t, row_counts);
  for(i = 0; i < ln; i++) lrow_counts[i] = 0;

  lgp_barrier();
  
  for(i = 0; i < el->num; i++){
    if(el->edges[i].row < 0) continue;
    assert(el->edges[i].row < n);
    lgp_atomic_add(row_counts, el->edges[i].row, 1L);    
  }

  lgp_barrier();
  
  lgp_barrier();
  // this second barrier seems to be necessary in OpenMPI's OSHMEM 1.4 on SMP
  // We were seeing occasional errors without it.
  
  // Step 9. Sum the vertex degrees for the points you will receive
  int64_t lnnz = 0;
  for(i = 0; i < ln; i++) lnnz += lrow_counts[i];

  lgp_barrier();
  lgp_all_free(first_point_in_sector);
  
  // Step 10. Initialilze the sparse matrix
  //T0_printf("Trying to init matrix with %ld nonzeros on PE 0\n", lnnz);fflush(0);
  sparsemat_t * A = init_matrix(n, n, lnnz, weighted);
  if(!A){
    T0_printf("ERROR: geometric_random_graph: init_matrix failed.\n");
    return(NULL);
  }
  
  // initialize the offsets for the sparse matrix
  A->loffset[0] = 0;
  for(i = 1; i <= ln; i++){
    A->loffset[i] = A->loffset[i - 1] + lrow_counts[i-1];
  }
  //if(A->loffset[ln] != lnnz){
  //fprintf(stderr,"ERROR: PE %d: %ld %ld\n", MYTHREAD, A->loffset[ln], lnnz);
  //}
  assert(A->loffset[ln] == lnnz);
  
  for(i = 0; i < ln; i++) lrow_counts[i] = 0;

  // Step 11. Distribute the edges to the proper PEs
  // and populate the matrix struct.
  exstack_t * ex = exstack_init(128, sizeof(edge_t));
  edge_t edge;
  if(ex == NULL){return(NULL);}
  i = 0;
  while(exstack_proceed(ex, (i == el->num))){
    for(; i < el->num; i++){
      edge_t e = el->edges[i];
      if(e.row < 0) continue;
      pe = e.row % THREADS;
      if(pe == MYTHREAD){
        int64_t row = e.row/THREADS;
        int64_t pos = A->loffset[row] + lrow_counts[row]++;
        A->lnonzero[pos] = e.col;
        if(weighted){
          A->lvalue[pos] = lgp_rand_double();
        }
      }else{
        //printf("Pushing %ld %ld to pe %ld\n", e.row, e.col, pe);
        if(exstack_push(ex, &e, pe) == 0L)
          break;
      }
    }

    exstack_exchange(ex);

    while(exstack_pop(ex, &edge, NULL)){
      int64_t row = edge.row/THREADS;
      assert(row < ln);
      int64_t pos = A->loffset[row] + lrow_counts[row]++;
      assert(pos < A->lnnz);
      assert(edge.col < n);
      A->lnonzero[pos] = edge.col;
      if(weighted) {
        A->lvalue[pos] = lgp_rand_double();
      }
    }
  }

  lgp_barrier();
  
  for(i = 0; i < A->lnumrows; i++){
    assert(lrow_counts[i] == (A->loffset[i+1] - A->loffset[i]));
  }
  
  exstack_clear(ex);
  free(el);
  lgp_all_free(row_counts);

  sort_nonzeros(A);
  
#if 0
  print_matrix(A);
#endif
  
  lgp_barrier();

  return(A);
}


