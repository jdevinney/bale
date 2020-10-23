## unionfind

### Definition
We use the unionfind data structure to implement the union of disjoint sets
approach to finding the connect components in a simple graph.

### Algorithms
This algorithm relies on the notion of disjoint subsets.  Given a collection
of disjoint subsets that covers a space, the union of any of the members of
the collection will result in another collection of disjoint subsets that
covers the space.  The algorithm starts with each vertex in its own subset.
Then it looks at each edge in the graph and forms the union of the subsets
that contain the incident vertices.  The resulting subsets are the connected
components.

The key to algorithm is a data structure that forms a tree for each subset.
The union of the subsets is formed by connecting the root of one tree to the
other tree.

We have implemented two versions (there are a number of versions): 
the first joins the trees by connecting
the root of the tree corresponding to one vertex of the edge to the node
corresponding to the other vertex of the connecting edge. 
This is referred to as the "bad algorithm" because it is not 
nearly as efficient as the second version. The most efficient
version connects the roots of the two trees according to the rank of the
trees, where the rank is the length of the longest branch in the tree.

### Discussion
This algorithm is in the C cousin of bale because we
don't have a parallel version.

Like the use of the priority queue in Dijsktra's algorithm, this algorithm is
a favorite example of how important data structures are.  In these algorithms,
the data structure is more than a place to storage the state.  Manipulating
the data structure *is* the computation.

In the most efficient version of the algorithm are manipulating structures
that hold pointers that encode the tree and the rank of the tree.  To do this
in an AGP model, the whole operation must be done atomically.  Doing this in a
lock-free way is currently beyond our capability.

Unlike Dijsktra's or the Fisher-Yates algorithm, the use of this forest
of trees is not necessarily serial.  There is plenty of opportunity for
asynchronous parallelism, but keeping the data structure consistent while
making parallel changes to it seems overwhelming.  

There are other algorithms
to find the connected components in a graph.  We present this algorithm
because we are interested in a discussion about parallel data structures.

### References
