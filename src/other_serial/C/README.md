<!---
/* Copyright (c) 2020, Institute for Defense Analyses              */
/* 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500 */
/*                                                                 */
/* All rights reserved.                                            */
/*                                                                 */
/* This file is part of Bale.   For license information see the    */
/* LICENSE file in the top level dirctory of the distribution.     */
--->

# C_bale: the bale apps written in serial C serial
## In one sentence
A textbook, C, implementation of the apps in bale.

## The elevator pitch

The bale effort is, first and foremost, 
a vehicle for discussion for parallel programming productivity.  
We have included a simple C version of the apps in 
bale_classic as a concrete description of the apps written in a familiar language.
Some of the apps (like histo and ig) are trivial as serial apps.
It might be useful to view the serial version of the more complicated apps
(like toposort and sssp) before dealing with the parallelism and buffered communication 
of the bale_classic apps.  The C version of our implementation of a 
compressed row format data structure for a sparse matrix is also simpler than
the parallel version.  Finally, we also have examples of applications that
are efficient as serial codes, but have no known parallel implementations.

This is a self contained directory that does not depend 
on build process for the rest of bale.

### Where does it run?
This is written in generic C, hopefully it will run in any C environment.
It does depend on the ``argp`` library for command line argument parsing,
doxygen for documentation and python3 for ``pytest`` unit testing.

### What is included here:

- Makefile - simple explicit makefile for all the apps.

- APPS:
  - [histo](histo.md) -- creates a large histogram (random stores) (see doxygen: histo.c)
  - [ig](ig.md) -- a large gather (random loads) (see doxygen: ig.c)
  - [randperm](randperm.md) -- creates a random permutation (see doxygen: randperm.c)
  - [transpose_matrix](transpose_matrix.md) -- computes the transpose of a sparse matrix (see doxygen: transpose_matrix.c)
  - [permute_matrix](permute_matrix.md) -- applies row and column permutations to a sparse matrix (see doxygen: permute_matrix.c)
  - [triangle](triangle.md) -- counts the number of triangles in a graph (see doxygen: triangle.c)
  - [toposort](toposort.md) -- performs a topological sort of a morally upper triangular matrix (see doxygen: toposort.c)
  - [sssp](sssp.md) -- solves the single source shortest path problem on a graph (see doxygen: sssp.c)
  - [unionfind](unionfind.md) -- uses the union-find data structure to find connected components in a graph (see doxygen: unionfind.c)

- Other:
- [spmat_utils](spmat_utils.md) -- the sparse matrix library and some support functions (see doxygen: spmat_utils.h, spmat_utils.c)
- [std_options](std_options.md)  -- the command line parsing routines (see doxygen: std_options.h, std_options.c)
- [opts_demo](opts_demo.md) -- programs to play with modifying the command line options (see doxygen: opts_demo.c) 

## Build Instructions
This is meant to be basic C, with a simple Makefile, so hopefully just typing 'make' will work.
```
    make  # make the apps
    make test # run the unit tests (if python3 is installed)
    make doc # make the doxygen documentation
```
However, we do use the argp library from the GNU standard library. 
This is usually present by default on most linux systems. 
On Mac, you will probably have to install it by hand 
(one way to do this is 'brew install argp-standalone') and then mess with your
LD_LIBRARY_PATH and LDFLAGS variables.

## Testing
For bale 3.0 we have started using pytest for the unit testing. This requires python3.
In addition to using `make test` as part of the build process,
one can run the test specified in file ``tests/test_all.py`` by hand with the command:

```
    pytest -s
```
It is also easy to edit the file `tests/test_all.py` to add or modify the tests that are run.


## Documentation
serial_C is documented using Doxygen. 
If doxygen is installed on your system, you should be able to type ``make doc``
and then load ``./html/index.html`` into a browser.

