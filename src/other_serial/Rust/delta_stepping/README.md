# Delta-stepping for single-source shortest path (SSSP)


This application finds shortest path lengths from a single source in
a directed graph, using the Meyer/Sanders [delta-stepping algorithm](
https://www.sciencedirect.com/science/article/pii/S0196677403000762).

It is a serial rust program.

### Test Matrices

The following are provided in the top Bale repository example_matrices directory:

Matrices:

- er_nn_prob.mm  reasonably dense ER graphs
- sparsennn.mm   sparser ER graphs

Distances:

- er_nn_prob.dst distances with source_vtx 0
- sparsennn.dst  distances with source_vtx 2
