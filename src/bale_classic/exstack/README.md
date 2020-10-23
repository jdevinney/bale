# exstack

This library contains both exstack and exstack2.

**exstack**

exstack was originally written in 2011 and was our first attempt at an aggregation library. Exstack is synchronous and resembles Bulk Synchronous Programming. The main functions in its API (push and pop)
remain in its descendants (exstack2 and [conveyors](../convey/README.md)). In a typical exstack loop there are
three phases. 

 1. First, each PE pushes items onto its local out-buffers. 
    
 2. Once a PE sees that one of its out-buffers has become
    full, that PE goes to the "exchange" phase where it waits for all
    other PEs to join it. Once all PEs are in the exchange phase, all
    out-buffers are sent to their destination where they land in
    in-buffers. 

 3. PEs then enter into the pop phase where they pop items off
    of their in-buffers and do whatever computation is required. 

See [histogram](../apps/histo_src/README.md) for a simple example. While exstack is naive compared to conveyors, it still acheives very good performance on most of the bale apps.

**exstack2**

exstack2 was our attempt to improve both the performance and ease-of-use of exstack. In fact, bale grew out of our attempt to write exstack2. exstack2's main difference from exstack is the fact that exstack2 is asynchronous. PEs do not proceed in lock-step and there is no "exchange phase". The buffer sends happen automatically and individually as soon as a buffer is full. 

Both exstack and exstack2 require logic for what we call the "endgame". This refers to the fact that PEs must all continue to participate in any exstack or exstack2 loop until all PEs are "done". This is because, even if a PE is done pushing items to other PEs, it needs to stay ready to process other PEs pushes to it. The endgame of aggregation libraries (especially asynchronous versions) is an interesting topic in itself.

