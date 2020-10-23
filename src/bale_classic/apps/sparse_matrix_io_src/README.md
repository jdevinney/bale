# Sparse Matrix IO

## Description

The goal of this app to write to disk and then read from disk a distributed sparse matrix dataset. It sounds like a very simple task; however, we throw in a some twists.

## Discussion

bale sparse matrices are stored so that the rows are distributed to PEs in CYCLIC order. That means rows i, i + THREADS, i + 2*THREADS... are local to (or "have affinity to") PE i. This is a contrast to BLOCK layout where PEs get continguous sets of rows. There are pros and cons to both CYCLIC and BLOCK layouts. One of the nice features of CYCLIC layout is that nonzeros tend to be better load balanced to PEs. 

Goals of our sparse matrix dataset:

* We would like to be able to read and write sparse matrices with no constraints on the number of PEs.
* We would also like the data to be laid out in row-major order on disk. That would make it more convenient to explore the data on disk if necessary.
* We would like to be able to read/write a sparse matrix without allocating a lot of extra memory (for example, a whole matrix worth of memory is unacceptable).

All of these goals are automatically and easily achieved using BLOCK layout. However, as we mentioned, CYCLIC layout has its perks too. Also, that would not be that interesting. So for bale, we are going to concentrate on reading matrices from a sparse matrix dataset (with the properties above) into memory in CYCLIC layout, and writing a matrix in BLOCK layout onto disk as a sparse matrix dataset.

Another interesting idea for the write is to always just dump the data in CYCLIC layout. Then when reading back in, the PEs need to figure out which rows they read and where those rows go. This can be tricky since the rows are laid out in a way that depends on the number of writing PEs and the communication pattern during the read depends on the number of PEs used to write and how many are reading. We do not pursue this idea in bale.

#### Sparse Matrix Dataset

The bale sparse matrix dataset contains:

* An ASCII metadata file with some basic matrix details
* N "nonzero" files. These files contain the column indicies for the nonzeros in the matrix in row-major order.
* N "values" files (optional). These files are only present if the matrix has "values", in which case these are the values of the nonzeros in row-major order. The values are floating point numbers.
* N "row_info" files. The ith record in file j tells where the nonzero (and value) data for the ith row in file j begins in the jth nonzero (and value) files. These files function much like the "offset" array in the sparse matrix data structure.

### Writing

One way to do this task would be to shuffle rows from CYCLIC layout to BLOCK layout and then just dump the data. As described this sounds just like the [permute_matrix](../permute_matrix_src/README.md) app. Since we do not want to allocate the space to store an entire second copy of the matrix, we need to send the rows in batches to the PE that would own them if the matrix were in BLOCK layout and write small buffers of data. This is where an additional challenge of this application comes in. We need to send just the right data so that each PE is able to fill a write buffer with the data it is supposed to write. This makes this code challenging for aggregation, especially asynchronous aggregation. In fact, we only have an AGP and exstack version of write_sparse_matrix currently. We don't think an exstack2 version would be reasonable. We hope to soon write an efficient conveyor version.

### Reading

Our read_sparse_matrix routines (again, we currently only have AGP and exstack versions) have an extra parameter: the number of readers. We have seen evidence that as the number of cores on multi-core processors grow, I/O performance can suffer when the number of cores concurrently reading or writing gets too large. For this reason, it seems logical that we should limit the number of PEs taking part in the reading or writing. We have elected to implement this twist in the read_sparse_matrix apps. The functions that read a sparse matrix dataset have a lot of complicated glue code that is not that interesting in the bale sense. This code is all included in spmat_io.upc. The heart of the application has two interesting communication sections. The first one is where we distribute the row counts for the rows the PEs read to their place in CYCLIC layout. The second one is where we distribute the nonzeros we read to the PEs where they will land in CYCLIC layout. Both of these are quite simple with the AGP model and fairly simple implemented in exstack as well.

### Why is this in bale?

Large dataset IO is important in HPC so we felt that this was a good reason to include the read and write sparse matrix functions as an app. Writing a sparse dataset is also more challenging than a dense dataset and requires a lot of uninteresting code that we would rather not write over and over. Writing a "dataset library" for general HPC sparse and dense datasets may be a future goal of bale.
