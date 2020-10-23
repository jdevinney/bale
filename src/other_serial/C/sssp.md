## sssp Single Source Shortest Path

### Definition
We are given the adjacency matrix for a graph with non-negative edge weights *c(v,w)*
and a given source vertex, *v_0*.
We wish to find the lightest weighted path from *v_0* to all other vertices,
where the weight of a path is the sum of the weights of the edges in the path.
Note, if the graph is undirected we work with the full (symmetric) adjacency matrix.

### Algorithms
We consider three algorithms: Dijsktra's, Delta-Stepping, and Bellman-Ford.
Dijsktra's algorithm is not in bale_classic because it is a serial algorithm.

Delta-Stepping and Bellman-Ford are here as shadows of the parallel versions.
The algorithms here are surprisingly similar to those in bale_classic.
These may be slightly easier to read because we don't have the communication layer.

### Discussion
The priority queue version of Dijsktra's algorithm is a favorite example
of the use of data structures in irregular algorithms.
We discuss this issue in the bale_classic app and in the serial unionfind app.

### References
"Delta-stepping: a parallelizable shortest path algorithm" by U. Meyer and P. Sanders.
