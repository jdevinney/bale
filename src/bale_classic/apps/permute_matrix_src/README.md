# Permute Matrix

## Definition

This application takes a distributed sparse matrix (see [spmat](../../spmat/README.md) library), a permutation for the rows and a permutation for the columns of the matrix. The goal is to apply the permutations and return the permuted matrix. (A note on the permutations: If row_perm is a permutation of the rows, we interpret rowperm[i] = j to mean that row i in the original matrix should go to row j in the permuted matrix.). The implementations of permute_matrix are contained in the spmat library.

## Discussion

#### Parallel Implementation Considerations

Given our interpretation of the permutation arrays, processors know where all of their local rows are destined in the permuted matrix. Also, as a consequence of the locality properties of our distributed sparse matrix data structure, the nonzeros in a row all go to the same PE (permuting the rows and columns of a matrix does not change the number of nonzeros in each row).

This application is done in two phases. 

##### Phase 1 

In phase 1, we shuffle the row data so that the nonzeros for each row land in the correct place in the permuted matrix data structure. This phase looks a little like the histogram pattern, for each local row, PEs look up the destination of the row and send/write its data to the correct PE. 

There are a couple of ways of going about sending the nonzeros for a row in phase 1. In the AGP code, we can send it in bulk if we know where it is destined on the remote side (which we can get by remotely reading the offset array). In an aggregated code, we could send (row, col) pairs, which would make the code a little easier, but sends extra data, or we could send a small header (row, row_count) and then the row_count nonzeros. 

##### Phase 2

In the second phase, we relabel the column indices on every nonzero in the permuted sparse matrix. This phase looks exactly like indexgather. PEs look at each column index they have locally and request/read its permuted index in the distributed colperm array.

#### Why is it in bale?

permute_matrix an interesting application because it presents an example where the payload for the sends  (see phase 1) are not necessarily a fixed size.

#### From the Book?

The AGP version of permute_matrix (see spmat/spmat_agp.upc) is fairly concise, but still suffers from some pretty low level address arithmetic (that could be improved upon). The aggregated codes are pretty horrible.