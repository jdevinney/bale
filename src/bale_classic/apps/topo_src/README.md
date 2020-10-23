# toposort

#### Glossary

*unit-upper-triangular matrix*: A matrix that is upper-triangular and where every diagonal entry is equal to 1.

*morally unit-upper-triangular*: A matrix **M** is morally unit-upper-triagular 
if there are square permutation matrices **P** and **Q** and a unit-upper-triangular matrix **T**, such that **P M Q** = **T**. 

## Definition

The input to toposort is a morally unit-upper-triangular sparse matrix **M**.
To obtain **M** we create a unit-upper-triangular matrix
and apply random permutations to the rows and columns of that matrix.
Note that we don't care about the values of the non-zeros
in any matrix in toposort, only their position.
The goal of toposort is to create row and column permutations such
that when these permutations are applied to **M**,
the result is an unit-upper triangular matrix.
The answer need not be unique.
A solution always exists since **M** is a row and column permutation of **T**.


#### Base Algorithm

If you observe the set of the column labels of nonzeros
in a particular row; a column permutation might change the labels of
the elements in the set, but doesn't change the cardinality of the
set. Likewise for columns. Hence, there must be a row in **M** with a
single non-zero.  If we remove that row and the column it intersects,
we are left with smaller, morally unit upper triangular matrix. This is the
motivation behind a simple algorithm (and the outline of an induction
proof of its correctness).  

In the following graphic for the algorithm we labeled the nonzeros with 
letters to help follow the permutations.

<img src="../../../../images/toposort.png" alt="picture of toposort" align=center style="height: 400px;"/>

An outline of an algorithm is as follows:

Let **M** be an N by N morally unit-upper-triangular matrix. 
```
   For all rows r with a single non-zero in a column c, put the pair (r,c) onto a queue.
   pos = N-1
   While the queue is not empty: 
     pop an (r,c) pair from the queue
     rperm[r] = pos; cperm[c] = pos; pos--;
     remove all the non-zeros in column c
     if any row now has a single non-zero, enqueue that (r,c) pair
```
Rather than changing the matrix by deleting rows and column and then searching the 
new matrix for the next row.  We do the obvious thing of keeping an array of row counts.
`rowcnt[i]` is the number of undeleted non-zeros in *row i*. 
We use a cool trick to find the last surviving column of a row. 
We initialize an array, `rowsum[i]`, to be the sum of the column indices in *row i*.
When we "delete" a column we decrement `rowcnt[i]` by 1 and `rowsum[i]` by that column index.
Hence, when the `rowcnt[i]` = 1, `rowsum[i]` contains the column that is left. 

## Discussion
This is an interesting algorithm because it is one of the simplest algorithms that has both
irregular data layout and irregular program flow. 
It enjoys some, but not unlimited, parallelism because all the pivots in the queue could
be worked on simultaneously. Because the queue only grows after pivots have been "removed"
it has some, but not unlimited, latency tolerance.

#### Parallel Considerations

In parallel there are three race conditions or synchronization issues to address.

1. The first is reading and writing the queue of rows to be processed.
   One way to handle it is to introduce the notion of a levels.
   Within a level all threads process the all the rows on their queues 
   and by doing so create new degree one rows. These rows are placed on the 
   appropriate queues for the next level. There is a barrier between levels.

2. Threads race to pick their position in *rperm* and *cperm*. 
   One could handle this race for the pivots with a fetch_and_inc for each new pivot. 
   An improvement on this idea is to use one fetch_and_add to claim enough room 
   for all the local pivots in the current level on each thread then assign them in order per thread.

3. Threads race to update the *rowcnt* and *rowsum* arrays. 

#### Why is it in bale?	

The toposort application is a major step up in complexity 
from [histogram](../histo_src/README.md) and [indexgather](../ig_src/README.md). 
Successfully implementing toposort (including all of the pre-computation of its input) 
requires 3 other bale apps (transpose_matrix, randperm, permute_matrix). 
In this way, toposort represents a crucible for any new parallel programming model.

Depending on its input, toposort typically enjoys a significant 
amount of parallelism, but it is not completely order and latency tolerant. 

It was not obvious that toposort would be so well suited for communication aggregation. 
This discovery has made toposort a poster-child for the wide applicability of aggregation. 

Over several rounds of brainstorming, we have found new ways to think 
about implementing toposort that make it even more amenable to aggregation. 
For example, see the toposort_cooler.upc code in the alternates directory. 

#### From the Book?

We are not really satisfied with any of the versions of toposort in bale_classic. 
A more modern language (like Chapel or Rust) help make toposort more readable. 
But we still haven't seen the one FTB.

