## histo
### Definition
We form the histogram of a large number of `int64_t`'s into a large table.
The loop is as simple as:
```
foreach idx in index[ ]
  counts[idx]++
```
On a serial thread this app shows the difference between random stores
and streaming stores.  On a large parallel machine this is the simplest
case of managing the latency and bandwidth of the interconnection network
and the race condition of multiple threads updating the same entry.

### Algorithms
We start by filling the index array with random numbers between 0 and the `table size`.

We have the generic algorithm (as simple as the loop above).

We also have a buffered version where we sort the indices
into buffers, based on their high bits.  When a buffer gets full
we perform all the updates from that buffer's contents all at once.
Studying this version with different number of buffers and buffer sizes 
might reveal properties of the memory hierarchy, like page sizes or TLBs.

A third version first sorts the indices before running the loop.
This ought to be closer to the streams benchmark's performance 
as it is streaming with holes in it.

### Discussion
Comparing these random access patterns to the streams benchmark for
a particular node could be interesting.

In serial, we don't have the problem of atomic updates in `histo`.
The parallel version in bale_classic as to use some form of atomicity 
to handle multiple threads concurrently update a particular entry in `counts[ ]`.

Running this simple version of `histo` on a node, with node level threads,
might reveal something about atomic updates to memory,
without the performance being dominated by the interconnection network.

### References
https://www.cs.virginia.edu/stream
