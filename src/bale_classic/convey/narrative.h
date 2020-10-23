// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


/** \mainpage Conveyors

Conveyors are parallel objects that provide efficient, scalable mechanisms
for all-to-all communication.  The amount of data flowing between each pair
of processes is not fixed and may be small, or it may be large enough that
the data cannot all be stored in buffers simultaneously.  Each process
"pushes" items of data that it wishes to send, "pulls" items that it
receives, and repeatedly "advances" the conveyor to ensure forward progress.
These operations must be interleaved in a way that conforms with
\ref contract between the conveyor and its client processes.

By a "parallel object" we mean an object (in an abstract sense; this is C,
not C++) that spans multiple processes and some of whose operations
(methods) are collective in nature.  If the underlying model of parallel
programming supports a notion of "teams" or "communicators", then the set
of processes involved in a conveyor can be a subset of the processes in
the parallel program; otherwise all processes must participate.

Conveyors can be built upon several different programming models, including
SHMEM, MPI, and UPC.  The programming model does not affect the conveyor
API except by making certain types (subclasses) of conveyors unavailable.
For instance, when UPC is the communication substrate, the constructor for
synchronous "twohop" conveyors fails and returns \c NULL.


\section operations Operations

A conveyor, once created, has seven core operations: \b begin, \b push, \b pull,
\b unpull, \b advance, \b reset, and \b free.  Of these, \b begin, \b reset,
and \b free are collective, and calls to \b advance may be collective as well,
depending on the conveyor and its state.  The other operations are local.

The \b push, \b pull, and \b unpull operations can succeed or fail, and
their success or failure can be observed from their return values.  Which
operations a process has done, and their return values, determine the local
state of the conveyor.  The local state, in turn, determines which core
operations are legal.  An illegal operation always fails; it moves no data
and has no effect on the state of the conveyor. It does not crash, but it
may have side effects such as printing an error message.

\section states States and Transitions

<dl>
<dt> \e dormant </dt>
<dd> The legal operations are \b begin, \b reset (which does nothing), and
\b free.  The \b begin operation causes a transition to the \e working
state. </dd>
<dt> \e working </dt>
<dd> The legal operations are \b push, \b pull, \b unpull, and \b advance.
If \b advance is given the argument \c true, then the state transitions to
\e endgame, \e cleanup, or \e complete according to whether the return
value of \b advance is \e ok, \e near, or \e done.  Otherwise \b advance
returns \e ok and the state does not change. </dd>
<dt> \e endgame </dt>
<dd> The legal operations are \b pull, \b unpull, and \b advance (with
argument \c true).  If \b advance returns \e near or \e done, the state
changes to \e cleanup or \e complete respectively; otherwise it does not
change. </dd>
<dt> \e cleanup </dt>
<dd> The legal operations are \b pull, \b unpull, and \b advance (with
argument \c true).  If \b advance returns \e done, the state changes to
\e complete; otherwise it does not change. </dd>
<dt> \e complete </dt>
<dd> The only illegal operations are \b push and \b begin.  The \b pull
and \b unpull operations are legal but they fail, and \b advance returns
\e done.  The \b reset operation changes the state to \e dormant, and
\b free destroys the conveyor. </dd>
</dl>

A newly created conveyor is in the \e dormant state on every process.
The client can track the state by observing the return values from \b advance.

\section contract The Contract

To specify the contract between a conveyor and its client, we need the
notion of \b unpull and \b pull operations \e matching each other.  On a
given process, imagine matching up each successful \b pull with the most
recent preceding unmatched successful \b unpull, as if successful \b unpull
and \b pull operations were left and right parentheses, respectively.
(This relationship is not the most obvious one.  An unpull does not match
with the previous pull that it reverts; rather it matches with the
following pull that "repulls" the unpulled item.)

The client of a conveyor promises that each process will obey the following rules:
  - It will never call \b begin or \b reset or \b free when it is illegal,
    or in a non-collective way.
  - If the conveyor state is neither \e dormant nor \e complete, then the
    client will eventually attempt a \b pull at a time when there is no
    unmatched \b unpull.  (If this \b pull succeeds, it will be unmatched.)
  - If the conveyor state is neither \e dormant nor \e complete, then the
    client will eventually call \b advance with argument \c true.

Provided that the client upholds its part of the bargain, the conveyor promises:
  - If the conveyor is in state \e working, then repeated attempts to \b
    push an item, interleaved with \b advance operations, will eventually
    lead to success.
  - Every item that is successfully pushed is eventually delivered to
    exactly one unmatched \b pull operation on the desired destination
    process, prior the next \b reset or \b free.  Conversely, when a \b
    pull succeeds, the item it delivers is an item that was pushed to that
    PE after the most recent \b begin.
  - Items pushed on a particular process for delivery to a given process are
    delivered in the order that they are successfully pushed.  In other words,
    the channel between each pair of processes is first-in, first-out (FIFO).
  - If the immediately preceding operation on a process was a \b pull,
    then \b unpull succeeds.
  - Every item restored by a successful \b unpull is eventually delivered
    to exactly one subsequent matched \b pull, prior to the next \b reset
    or \b free.
  - If the conveyor is in state \e complete or \e cleanup on any process,
    then on no process is it in state \e working.  If the state is \e
    dormant on any process, then it is \e dormant on every process.
  - One of the following two properties holds:
     - Every call to \b advance is nonblocking.
     - If the conveyor is in state \e complete or \e cleanup on any process,
       then it is in state \e complete or \e cleanup on every process.  In this case
       \b advance is nonblocking.  Otherwise the conveyor is in state \e working or
       \e endgame on every process, and \b advance synchronizes all the processes.

The conveyor API does not currently promise anything about thread-safety.



\example histo.c
The \b histogram example illustrates one-way communication.  It uses a
single conveyor to accumulate counts in a distributed array.

\example gather.c
The \b indexgather example illustrates two-way communication.  It uses a
request conveyor and a reply conveyor simultaneously to gather values from
a distributed array.

\example bigather.c
A version of the \b indexgather example using a single biconveyor rather
than a pair of conveyors.

\example meld.c
The \b meld example illustrates a loop that combines data from two parallel
data structures. It relies on the strong progress guarantees of
biconveyors. If the loop were written with two pairs of conveyors rather
than two biconveyors, it would not work; it would be prone to deadlock.

*/
