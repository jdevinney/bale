## triangle
### Definition
Find the number of triangles in a given simple unweighted graph. 

A triangle is a set of three vertices {u,w,v} where edges {u,w}, {w,v} and {u,v} are in the graph.

### Algorithm
This uses matrix algebra approach to counting triangles in a graph.
Given, L,  the (strictly) lower triangular matrix that holds the undirected graph,
we compute \sum_ij{L .& (L * L)} and \sum_ij{L .& (L * U)}.
Where U is the upper triangular matrix from the full adjacency matrix.
Recall that  U = L transpose. 

### Discussion
This is here to shadow the algorithms in bale_classic. 

The amount of work done in these formulations depends on the matrix.
This is interesting and well studied even in the serial case, but not here.

In the parallel case, it is even more interesting because it depends 
on ones ability to push or pull information remotely as well as the 
row densities in the matrix.

### References
See the book, "Graph Algorithms in the Language of Linear Algebra",
edited by Gilbert, and Kepner for more details on our approach to this problem.

