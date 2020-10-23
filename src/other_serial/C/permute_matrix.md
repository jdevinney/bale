## permute_matrix

### Definition
We apply given row and column permutations to a sparse matrix.

### Algorithm
We produce a new sparse matrix data structure by copying the nonzeros and value entries
of the original matrix to their new positions in the permuted matrix.
To permute the rows, we compute the new offsets, based on the new order for the 
original rows.   
As we are copying the `nonzero[ ]` and `value[ ]` entries from the original matrix
to their new position in the permuted matrix,
we replace the nonzeros (the column indices) with the new column indices 
given by the column permutation.

### Discussion
This app is really just a wrapper that calls the routine in the sparse matrix library. 

This is C_bale to shadow the app in bale_classic.
It is here because we need it in the library.
It is an app in bale_classic because the communication pattern is interesting.

### References

