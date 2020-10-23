# bale

The bale effort is, first and foremost, *a vehicle for discussion* for parallel programming productivity.

The bale effort attempts to:

- demonstrate some challenges of implementing interesting (i.e. irregular) scalable distributed parallel applications.

- demonstrate an approach (aggregation libraries) to achieving high performance for the internode communication in such applications 

- explore concepts that make it easier to write, maintain, and get top performance from such applications

bale does not claim to have the answers to better parallel programming.

bale is not a collection of benchmarks.

We use bale to evolve our thinking on parallel programming in the effort to make parallel programming easier, more productive, and more fun. Yes, we think making it fun is a worthy goal!

The original bale (now called "[bale_classic](src/bale_classic/README.md)") was implemented in UPC or C and SHMEM. If you are new to bale, you should start by exploring there. New in bale 3.0, the bale repository now contains other variants of bale (sequential, parallel, and custom parallel). Each variant is a separate project and comes with its own build and run instructions. All source code is found in the src directory.

### Contents

**[bale_classic](src/bale_classic/README.md)**: Built with UPC or C/SHMEM. Contains three aggregation APIs (exstack, exstack2, and conveyors) and the widest collection of [apps](src/bale_classic/apps/README.md) and implementations of those apps.  Location: src/bale_classic



**[C](src/other_serial/C/README.md)**: Sequential versions of the bale_classic apps, meant to be a gentle introduction to some of the more complicated algorithms in bale_classic. Location: src/other_serial/C



**[Serial Rust](src/other_serial/Rust/README.md)**: Similar to the C version, only in Rust. Location: src/other_serial/Rust



**[Parallel Rust](src/other_parallel/Rust/README.md)**: Parallel versions of the bale apps in Rust on top of Rust conveyors (which are currently implemented on SHMEM). Location: src/other_parallel/Rust



**[Chapel](src/other_parallel/Chapel/README.md)**: A collection of some of the bale apps written in Chapel. Location: src/other_parallel/Chapel



### Version History

* May 2018: Initial Release version 1.0.0 

* Dec. 2018: version 2.0.0 
  * New apps: transpose, randperm, permute_matrix, write_sparse_matrix

* Sep. 2019: version 2.1.0
  * update conveyors to version 0.5.0

* Nov. 2020: version 3.0.0
  * Added cousins: Rust, Parallel Rust, Serial C, and Chapel to bale
  * new bale_classic features:
    * New graph model: Geometric graphs
    * New app: SSSP
    * replaced write_sparse_matrix with sparse_matrix_io
    * convey-0.6.0 updates
    * arg_parse replaced getopt
    * unit tests with pytest
    * new make_bale script
    * new run_apps script
    * docker files
    * AGP (Atomics, Gets, and Puts) replaces AGI: simple PGAS style programming
    * FTB (From the Book) replaces As God Intended.

### Contact

bale@super.org