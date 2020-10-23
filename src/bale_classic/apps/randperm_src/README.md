# randperm

## Definition

The goal of randperm is to create, in parallel, a distributed array that holds a uniform random permutation of {0,...,N-1}. 

#### Serial algorithm

As a review, to get a random permutation in serial, a well known algorithm is Knuth's "shuffle" (also, called Fisher-Yates) algorithm.

```c
for(i=0; i<N; i++)
   perm[i] = i;
for(L=N-1; L>0; L--){
   r = lrand48() % (L+1);
   s = perm[L];
   perm[L] = perm[r];
   perm[r] = s;
}
```

## Discussion

The above serial algorithm doesn't obviously extend to
a parallel version. At the time we first wrote bale,
we didn't actually know any parallel algorithm for
random permutations.

We found a paper on the "dart throwing algorithm":

P.B.Gibbon, Y.Matias, and V.L.Ramachandran. Efficient low-contention Parallel algorithms. 
J. of Computer and System Sciences, 53:417-442, Dec 1992. 

Each PE is responsible for some slice of {0,...,N-1}. These items become
the "darts". A distributed "target" table, which is required to be at
least as large as N, but should be much larger in practice for the sake of
efficiency, is allocated with M entries and all PEs throw their darts at
random locations at the table. If a dart hits a location in the table that
has not yet been hit, that dart sticks and its index is recorded in the table
at the location. If a dart hits a location which already full, the dart must
be thrown again until it sticks. Once all darts are stuck, we have a random
permutation of the numbers of {0,...,N} by looking at the target array from
0,...,M-1 and recording the indices of stuck darts.

#### Why is it in bale?

This algorithm is fun, interesting, has an intuitive AGP version, and is amenable to aggregation, so we decided to code it up. 

After some time, it became clear to us that this algorithm was not optimal
in terms of remote communication. We implemented a superior version of the
algorithm in the alternates directory. Later we learned there are even
better (i.e. "communication-free") algorithms. See "Efficient Parallel
Random Sampling -- Vectorized, Cache-Efficient, and Online" by Sanders
et al. Rather than replace the dart throwing algorithm implementation,
we decided to keep it in the main line. For one thing, it relies on an
interesting pattern of communication: updates are made remotely and don't
always "succeed". Additionally, the process of evolving our understanding of
the randperm problem is exactly what a researcher goes through when attempting
to code up a novel algorithm. Making this process easier and more productive
is exactly the motivation of bale. We want to make it easier for parallel
programmers to tinker and experiment.

#### From the Book?

The AGP version of the main loop of randperm (where we rely on atomic compare
and swap) is quite nice looking. Also, in the AGP version, the redistribution
of the "darts" at the end once every dart has stuck somewhere, is also quite
simple, but does more communication than is necessary.
