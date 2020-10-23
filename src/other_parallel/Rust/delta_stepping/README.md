# Delta-stepping for single-source shortest path (SSSP)


This application finds shortest path lengths from a single source in
a directed graph, using the Meyer/Sanders [delta-stepping algorithm](
https://www.sciencedirect.com/science/article/pii/S0196677403000762).

* [How to build](#Building-the-application)
* [How to run](#Running-the-application)
* [Command line options](#Command-line-options)
* [Correctness-test data](#Sample-data)
* [Outline of the algorithm](#Outline-of-delta-stepping)
* [Discussion](#Discussion)
* [Test Matrices](#Test-Matrices)

### Building the application

Build the application from the parallel Rust directory (the parent of this directory)
according to the [instructions](../README.md) there.

### Running the application

To run on 4 threads from the parallel Rust directory (the parent of this directory) say:

```
oshrun -n 4 target/release/delta_stepping
```

### Command line options

By default, the app generates its input as a flat random 
(Erdos-Renyi) directed graph with edge lengths uniformly random in (0,1).
It can also take as input a square sparse matrix, whose values are edge lengths,
in MatrixMarket format.

* -d : dump files "sssp.mm" (the input matrix) and "sssp.dst" (the output vertex distances)
* -e real-number : edge probability for Erdos-Renyi random graph, if no input file
* -f real-number : value for Delta, overriding the algorithm's choice
* -i filename : input file, MatrixMarket format
* -n integer : number of vertices in test graph, if no input file
* -q : quietly run without printing progress to stdout
* -s integer : source vertex, default 0
* -t : trace execution, appending to a file trace.#.out from each thread

### Sample data

The directory erdosrenyi/ contains several small graphs as .mm files, 
each accompanied by a .dst file that gives the correct distances.
The source vertex varies; see the [README](erdosrenyi/README.md) for details.
The files sparse100.mm and sparse100.dst are used by the unit tests.

### Outline of delta-stepping

The Meyer/Sanders [delta-stepping algorithm](
https://www.sciencedirect.com/science/article/pii/S0196677403000762)
finds shortest paths from a single source in a directed graph with
non-negative edge lengths c(v,w).

Each vertex has a "tentative distance" during the algorithm. The source has tentative
distance 0, and all other vertices initially have tentative distance inf. We
proceed by "relaxing" edges (in a clever order): relaxing edge (v,w) changes the
tentative distance of w to min(tent(w), tent(v) + c(v,w)). Eventually each vertex's
tent() distance becomes final, or "settled"; initially only the source is settled.

Unsettled vertices with tent() < inf are kept in "buckets" by tent() value; bucket i
contains vertices with tent() at least i\*delta and less than (i+1)\*delta, where
delta is a parameter.

The algorithm has three nested loops.

The outer (serial) loop is over buckets; an iteration processes vertices in the lowest
nonempty bucket until it is empty.

The middle (serial) loop is over "phases"; a phase consists of removing all the vertices
in the active bucket and relaxing all the "light" edges out of them (an edge is "light"
if it has cost at most delta, "heavy" otherwise). The edge relaxations in a phase may
cause vertices to enter the active bucket; phases continue until the bucket is empty.
At that point all the vertices that were in that bucket are settled.  Following the
light-edge phases, one more phase relaxes all the heavy edges from vertices deleted
from the active bucket; this cannot cause any vertices to go into that bucket.

The inner (parallel) loop implements the edge relaxations in a single phase.
Those relaxations can be done in any order, provided they are done atomically.

### Discussion

#### Parallel Considerations

The parallelism in delta-stepping is within a single phase, 
when edges are relaxed in parallel.
Every rank looks makes a list of edges to relax whose tails it owns;
then it conveys each relaxation request to the rank that owns the head of the edge,
which relaxes it.
An edge relaxation happens atomically because, while a vertex may be the head or
more than one relaxation request in a phase, the rank that owns that vertex is the
only one that can relax it.

There is no shared data in the algorithm. 
The graph/matrix is partitioned by rows, that is by edge tail vertices.
Each processor has a copy of all the buckets, and links vertices it owns into its buckets.
There is a parallel reduction at each iteration of the outer loop (over buckets) to
find the next bucket that is nonempty on some rank,
and there is a parallel reduction at each iteration of the middle loop (over phases)
to determine whether the active bucket is empty on every rank.

#### Why is it in bale?	

This is a more complicated graph algorithm than toposort, 
but its parallelism is simpler; 
indeed, the single one-way conveyors
have essentially the same structure as in histogram.
The novel feature is that the "atomic" operation to be done at the far end 
of the conveyor is not just an increment or an add, 
but a complex user-defined operation (relax) that does comparisons
and modifies linked lists.
This comes out looking very nice (imho) in Rust Conveyors, where the
closure in the conveyor just calls a "relax" routine that is identical
to the one in the serial Rust delta-stepping code.


#### From The Book?

We don't think there's an FTB version based on gets, puts, and hardware atomics.
The Book version probably requires user-defined atomics,
and then probably looks pretty close to the half-page pseudocode in the Meyer/Sanders paper.
It is interesting that the parallel Rust Conveyors code looks almost identical to the serial Rust code.

### Test Matrices

The following are provided in the top Bale repository example_matrices directory:

Matrices:

- er_nn_prob.mm  reasonably dense ER graphs
- sparsennn.mm   sparser ER graphs

Distances:

- er_nn_prob.dst distances with source_vtx 0
- sparsennn.dst  distances with source_vtx 2
