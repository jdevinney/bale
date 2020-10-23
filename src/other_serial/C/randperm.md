## randperm

### Definition
We fill an array of `int64_t`'s with a flat random permutation.

### Algorithms
#### the Fisher-Yates algorithm
```
   fill the array, rand[ ], with indices 0 thru n-1. 
   for l=n-1, n-2, ... ,1
     swap rand[l] with a random entry in 0,1,...,l
```
By definition this picks "n balls from an urn without replacement".
This is a standard serial algorithm that is in fact a serial algorithm.
You have to process the entries from right to left one at a time.

#### the "dart throwing" algorithm
We pick a dart board (an array) that is bigger than the desired permutation,
say twice as big and fill the entries with -1. 
Then we randomly throw darts (numbers from 0 to `len-1`) 
at the dart board, re-throwing any dart that hits an entry that is already occupied (!= -1). 
Then we squeeze out the holes.

We picked the dartboard to be twice the size of the array 
so that even the last dart has a 50/50 chance of hitting an open entry.

### the sorting algorithm
We form an array of (index, key) pairs. Then we randomly fill the keys
and sort the array on the keys. Then we read the permutation from the indices.

NB. Repeated key are bad, but tolerated. 
They would be ok if ties were broken randomly or if doubles were real numbers. 

### Discussion
The Fisher-Yates algorithm must do one thing at a time, so it doesn't parallelize.

The dart throwing algorithm is here because it shadows the algorithm
in bale_classic.  It is in bale_classic because it is fun.
And because the AGP version is essentially the same as this serial version.

The sorting version seems like a reasonable parallel algorithm, but it 
not in bale_classic because its not that interesting nor fun.

### References
https://en.wikipedia.org/wiki/Fisher-Yates_shuffle
