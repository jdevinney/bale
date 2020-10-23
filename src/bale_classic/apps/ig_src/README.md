# indexgather

## Description

The indexgather app involves the asynchronous reading of random elements from a distributed array.

In SHMEM, indexgather looks like this:

```c
for(i = 0; i < N; i++)
   shmem_get(&target[i], &table[index[i]], sizeof(long), index[i] % NPES);
```

where table is a distributed array with M elements, index is a local array of indices mod M, 
and target is a local array where we record the results of the remote reads. Similar to [histogram](../histo_src/README.md), the intent of indexgather is that the PEs read asynchronously and without dependencies (allowing for aggregation).

## Discussion

Indexgather is another rather simple application. Just like histogram, it is important because it represents a pattern of communication that is frequently used in parallel applications. In terms of difficulty to aggregate though, it is a step up from histogram.

Below is indexgather written with exstack. Note that we need three distinct phases. 
  - In phase 1, we send out requests for the remote read. 
  - In phase 2, we process requests (pop them off in-buffers, 
 look up the requested index in the table, and push them back onto out-buffers). 
  -In phase 3, we receive our processed requests and record their values in our tgt array.

```c
while( exstack_proceed(ex, (i==l_num_req)) ) {
   i0 = i;  
   while(i < l_num_req) {
      l_indx = pckindx[i] >> 16;
      pe  = pckindx[i] & 0xffff;
      if(!exstack_push(ex, &l_indx, pe)) 
        break; 
      i++;
   }   

   exstack_exchange(ex);

   while(exstack_pop(ex, &idx , &fromth)) {
     idx  = ltable[idx];
     exstack_push(ex, &idx, fromth);
   }
   
   lgp_barrier();
   
   exstack_exchange(ex);
  
   for(j=i0; j<i; j++) {  // retrace the requests
      fromth = pckindx[j] & 0xffff;
      exstack_pop_thread(ex, &idx, (uint64_t)fromth);
      tgt[j] = idx;
   }
   lgp_barrier();
}
```

Because exstack is synchronous and we are gathering a single word, 
we can get reuse the same exstack for both the request and the response.
This is not possible with exstack2 or conveyors however since they are asynchronous. This makes indexgather (and any generalizations of it) especially interesting. It is extremely simple, yet frustratingly complicated to write code that aggregates the task.

