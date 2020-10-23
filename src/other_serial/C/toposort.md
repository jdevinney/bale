## toposort
### Definition
We are given a matrix that is a random row and column permutation 
of an upper triangular matrix (with ones on the diagonal).
Such a matrix has been called a morally upper triangular matrix.
This algorithm finds a row and column permutation that, when applied,
returns it to an upper triangular form.

### Algorithms
We generate the morally upper triangular matrix by 
randomly permuting the rows and columns of a triangular matrix.
In the figure below, we have marked the nonzeros in the matrix
with letter to help follow the permutations.

To find row and column permutations that would return the matrix 
to upper triangular form, we recursively find pivot positions.
A pivot is a nonzero (really a (row,col) pair) that is the single
nonzero in a row. If we permute this row and column to the last
row and column of our new matrix and delete the row and column
from the original matrix, we can recursively construct the new
matrix from the bottom right corner to the top left corner.

<img src="../../../../images/toposort.png" alt="" align=center style="height: 400px;"/>
<img src="../../../images/toposort.png" alt="" align=center style="height: 400px;"/>

A more detailed description of the algorithm is given in 
bale_classic toposort documentation.

#### enqueuing pivots
In this version when the pivots are found they are placed in a queue.
The algorithm runs until the queue is empty.

#### loop to find pivots
In the loop version we simply continue to loop over the rows 
until we have found all the pivots.  This simplifies the 
flow of the algorithm but does redundant checking of 
rows which have already been processed. 

### Discussion
In the serial case, the use of a queue seem like an obvious win.
In the parallel case, the queue has to be managed with remote operations. 
Whether or not it is a win in this case is an interesting discussion.

### References
