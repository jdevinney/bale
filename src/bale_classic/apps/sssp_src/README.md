# Single Source Shortest Path  (sssp)

## Definition
We are given a directed graph with non-negative edge weights *c(v,w)*
and an given source vertex, *v_0*.
We wish to find the lightest weighted path from *v_0* to all other vertices,
where the weight of a path is the sum of the weights of the edges in the path.

## Discussion
This is well studied problem.  
It is used in classes on algorithms, data structures, and graph theory.
It is also used as a benchmark in Graph500 benchmark suite.
We are using this problem to evaluate what we have learned
about our programming models.
Our approach is that of someone given the assignment
to get an SSSP application running reasonable well,
at scale, using our programming models.

We consider three algorithms: Dijsktra's, Delta-Stepping, and Bellman-Ford.
Good references for these uses are easily found online.
The explanations given here are meant to explain the code and
our understanding of the parallelism, not to explain the algorithms.
To study the Delta-Stepping algorithm we are using the original Meyer/Sanders paper,
[delta-stepping algorithm](https://www.sciencedirect.com/science/article/pii/S0196677403000762).

These algorithm work by assign a tentative weight, *tent(v)*, to each vertex. 
Initially the weight of the given vertex is zero and all other weights are set to infinity.
The algorithms proceed by changing the state of the vertices
 from *unreached* (with infinite weight)
 to *unsettled* (with a tentative weight)
 to *settled* (meaning that the lightest weight path to that vertex is known).
The vertices change state via the process of *relaxing edges*.  Relaxing an edge *(v,w)*
 replaces *tent(w)* with the *min(tent(w), tent(v)+c(v,w))*.

There is subtle relationship between the length of a path
(the number of edge in the path) and the weight of the path.
Clearly, it is possible to have a longer path that is lighter than a shorter path.
If a path from *v0* to *w* contains the vertex *v*, then the weight of the path
from *v0* to *v* is less than the weight of the path from *v0* thru *v* to *w*.
If *v* is an intermediate vertex on a lightest path from *v0* to *w*,
then that path from $v0$ to $v$ is also a lightest path.
The different algorithms play off these relationships in different ways.

Dijsktra's algorithm works by considering the weight *tent(v)* of the "unsettled" vertices.
One proves that lightest such vertex does, in fact, have the correct weight.
So, that vertex can be "settled" and we need only relax edges from that vertex.
Furthermore, we never have to relax those edges again.
This gives the most efficient algorithm in the sense that edges are relaxed exactly once.
The problem of finding the lightest unsettled vertex is addressed below.

Bellman-Ford relaxes edges based on the length of paths 
to the "unsettled" and "unreached" vertices.
Basically the algorithm simply relaxes all of the edges in the graph over and over again 
until none of the tentative weights *tent(v)* change. 
The algorithm finds the lightest tentative weights to all vertices on paths of length *i*,
for *i* equal zero to the length of a longest path.

The Meyer/Sanders delta-stepping algorithm uses ideas from both of the previous algorithms.
Unsettled vertices are kept in "buckets" based on their by *tent(v)* weights; bucket *i*
contains vertices with *tent(v)* at least _i\*delta_ and less than _(i+1)\*delta_, 
where *delta* is a parameter.  The active bucket is the *i*th bucket (with the smallest *i*)
that has unsettled vertices.  Edges in the graph are considered light if their weight 
is less than or equal to *delta* and heavy otherwise.
The algorithm "settles" the vertices in the active bucket with a Bellman-Ford approach
using only the light edges. Then the algorithm relaxes the vertices that were in the active
bucket using the heavy edges. This uses an extension of Dijsktra's approach because the 
heavy edges cannot put vertices into the active bucket.  The efficiency of the algorithm
comes at the price of flow control and data-structure manipulations
to maintain the buckets.  Here is psuedo code for the algorithm:
```
// We write an edge with tail v, head w and weight (cost) as {v,w,c}.
// We will say the weight of the lightest tentative path 
// to a vertex is the price of the vertex.  The relax routine is 
// more involved than it is in other sssp algorithms.  In addition 
// to (possibly) reducing the price of the head of an edge,
// it can move the head from one bucket to another, possibly to the 
// bucket of the tail.

def relax(w, p, B):
  if tent[w] > p:
    remove.bucket(w, B(tent[w]/delta))
    add.bucket(w, B(p/delta))
    tent[w] = x


program sssp:

  set tent[v] = inf for all v

  relax(s, 0, B)  // sets tent[s] = 0  and puts s in B[0]

  while "there is a non-empty bucket" :
    let B[i] "be the first (smallest i) non-empty bucket"
    set R = NULL   // set of vertices that we will retire

    while B[i] is not empty:
      let v = a vertex in B[i]
      for all light edges {v,w,c}:
        p = tent[v] + c  // possible new price of w
        relax (w, p, B)  // possibly adding w to B[i]
        add v to R
        remove v from B[i]
    
    for all v in R:
      for all heavy edges {v,w,c}
        p = tent[v] + c  // possible new price of w
        relax (w, p, B)  // won't add w to B[i]
```        

We have implementations of Delta-Stepping and Bellman-Ford in this directory.
We also have implementations of all three algorithms in the other_serial/C "cousin" 
directory and the delta-stepping algorithm in the other_serial/Rust and other_parallel/Rust directories.

#### Parallel Considerations
Dijsktra's algorithm is a serial algorithm 
with no opportunity for our static thread based parallelism. 
One must find **the** single lightest *tent(v)* and then relax the edges from that *v*.
In addition, serial versions of the algorithm can use a priority queue to efficient find
this next unsettled vertex. In our C implementation we use a simple binary heap.
Depending on the average degree of a vertex, one could imagine relaxing the edges
in parallel with a fork-join threading model, but that is limited parallelism and
introduces difficulties in maintaining the consistency of the heap.

Bellman-Ford does a lot of redundant relaxations. We attempt to limit the amount of 
redundant work with a dynamic programming approach, but it still seems to be over-whelming.
In its favor, it has trivial flow control, no data structure manipulations 
and massive amounts (equal to the number of edges) of parallelism.

The delta-stepping algorithm has a tunable amount of parallelism 
that depends on the parameter *delta*.
The more vertices that are in the active bucket the more relaxations,
possibly redundant relaxations, that can be done in a parallel.
Note that the heavy edges are relaxed only once and they are done as a batch
of all the heavy edges from any vertex that was in the active bucket.

The parallelism available in the Bellman-Ford and Delta-stepping algorithms
is a natural fit to our buffered communication.
A thread with affinity to the tail of an edge knows the *tent(v)* 
of a vertex and the weight of the outgoing edges from *v*. It is responsible
for picking which edges to relax.  A thread with affinity to the head
of an edge knows whether or not the relaxation improves *tent(w)* of the head
and is responsible for deciding that and updating the data structures accordingly.

If our AGP model had an atomic min of a `double`, we could reasonable write
an AGP version of Bellman-Ford.  Our buggy version of the AGP is currently 
in the alternate sub-directory for sssp. An AGP version of the Delta-Stepping
algorithm does not seem feasible at this time.

#### Why it is in bale?
This problem and the different algorithms give a new use case to study programming models.
The Bellman-Ford is simple algorithm, but it requires an atomic min of a `double` to 
perform the relaxations in the AGP model.
Most of the apps in bale do lots of communication and very little work at the remote
side of the communication.
The Delta-stepping algorithm provides a new data point 
because the amount of work done "remotely" is at least as much as the work done "locally".
The analysis of the SSSP algorithms is on-going.

#### From the Book?

We don't think there's an FTB version of any of these algorithms using
the AGP model because of the amount of work that must be done
remotely and atomically for every edge relaxation. For this reason the aggregated versions
of delta-stepping are actually easier to code than an AGP version and are not far off of the serial implementation of SSSP.

See apps/sssp_src/ for the all implementations.

