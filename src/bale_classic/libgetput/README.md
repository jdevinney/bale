# libgetput (A basic parallel library)

This library is meant to provide a very simple parallel programming layer. Namely it provides basic remote gets and puts, atomics, and a few collective functions (like barriers and reductions). Our history programming with UPC programming with this small set of operations led to a simple parallel programming model that could achieve good performance on machines of the past (Cray T3E for instance), but that suffered from poor performance on today's machines. For years, this made us nostalgic for computers that allowed for this style of programming without the performance hit. However, as our experience with bale has deepened, we now believe that parallel programmers need more than this. 

One feature of libgetput is that it acts as a wrapper for some of the simple functionality that is common between UPC and SHMEM. It allows most of bale to be easily compiled against UPC or SHMEM. 

As we mentioned above, as its name implies, libgetput provides basic remote get and put functions:
lgp_get_int64, lgp_put_int64 and lgp_getmem and lgp_putmem (which are
similar to SHMEM functions for single word gets and puts or more
general gets and puts of memory). One key distinction is that the
indexing in libgetput is UPC style indexing. That is we consider the
distributed array to be indexed in round-robin fashion: the first
element has affinity to PE 0, the next element has affinity to PE 1
and so on.

libgetput also supplies a variety of atomic functions, both fetching and non-fetching. Finally, libgetput provides some fundamental collectives: value-based reductions and barriers.