# transpose_matrix

## Definition

In this application we transpose a distributed sparse matrix. See the [spmat](../../spmat/README.md) library for this data structure. The implementations of transpose_matrix are contained in the spmat library also.

## Discussion

#### Parallel Implementation Considerations

The algorithm we use to transpose a matrix has two phases. Both phases of the transpose function are well-suited for aggregation.

##### Phase 1

PEs use a histogram pattern to calculate the number of nonzeros in each column of the matrix. This allows us to allocate the exact space needed to store the transpose matrix. This phase is not strictly necessary if memory constraints are not a concern. We include it (and other similar measures in other bale code) in an effort to demonstrate efficiency when it comes to resource utilization.

##### Phase 2

We again use the histogram pattern, but this time we send the nonzeros of the transpose matrix to the correct PE. For example, if *A[i,j]* lives on PE k in the original matrix, PE *k* sends a *(i,j)* to PE *m* to create *AT[j,i]* (where *AT* is the transpose of *A*). 

#### Why is in bale?

An interesting difference between the AGP and aggregated versions comes in the second phase. In the AGP version, PEs are placing nonzeros in AT directly into the distributed sparse matrix data structure via remote writes. To do that in parallel, PEs must be able to atomically reserve a spot in the "nonzero" array for their writes. Since we don't know in what order these nonzeros will arrive. the PEs are contending for these writes. In the aggregated versions, PEs are being sent nonzeros (via an aggregation library) and process them in serial. So there is no need for atomic operations. This phenomenon occurs in other bale apps (histogram for example) and is fairly common when going from AGP style paralell code to aggregated code.

#### From the Book?

Similar to permute_matrix, the AGP version of permute_matrix (see spmat/spmat_agp.upc) is fairly concise, but not beautiful. The aggregated versions are OK, but could be improved vastly (we think) with a more modern language.

