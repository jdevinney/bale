Conveyors
=========

Conveyors are parallel objects that provide efficient, scalable mechanisms
for all-to-all communication.  The amount of data flowing between each pair
of processes is not fixed and may be small, or it may be large enough that
the data cannot all be stored in buffers simultaneously.  Each process
_pushes_ items of data that it wishes to send, _pulls_ items that it
receives, and repeatedly _advances_ the conveyor to ensure forward
progress.

This package implements conveyors in C.  By a "parallel object" we mean an
object (in an abstract sense) that spans multiple processes and some of
whose operations (methods) are collective in nature.  If the underlying
model of parallel programming supports a notion of a "current team" or a
"current communicator", then the set of processes involved in a conveyor
can be a subset of the processes in the parallel program; otherwise all
processes must participate.

Conveyors can be built upon several different programming models, including
SHMEM, MPI, and UPC.  The programming model does not affect the conveyor
API except by making certain types (subclasses) of conveyors unavailable.


Documentation
-------------

+ For installation instructions, see the INSTALL file.
+ For a full API reference, see `narrative.h` and `convey.h`.
  Better yet, build the package and view the Doxygen-created HTML.
  When using conveyors in a more complex way than illustrated by
  the sample code below, an understanding of the contract between
  the conveyor and its client is essential.
+ In the bale distribution, the paper `uconvey.pdf` (in the `docs`
  directory) offers a full discussion of the design of the conveyor
  package, the motivation behind it, and its performance characteristics.


Sample Code
-----------

To create a conveyor, call one of the `convey_new` constructors.  A
standard conveyor for fixed-sized items can be obtained thus:

    convey_t* c = convey_new(SIZE_MAX, 0, NULL, 0);

It should be reasonably efficient in terms of both memory usage and
communication throughput.

Once you have a conveyor, you can set up a "session" to transmit items of
a particular type.  A conveyor deals with blobs of bytes; it performs no
serialization or deserialization.  So its notion of a datatype is just a
size and an alignment.  If you have a datatype to transmit such as
`packet_t`, then in C11 you would start a conveyor session thus:

    convey_begin(c, sizeof(packet_t), alignof(packet_t));

assuming you had said `#include <stdalign.h>`.  If you are writing C99,
then you can safely write instead

    convey_begin(c, sizeof(packet_t), 0);

The `convey_begin` operation, like `convey_new`, is collective: all
processes must call it concurrently.

Now you need to write a loop to transmit data.


### Messages Without Replies

Suppose for simplicity that each process knows that it wants to send
`n` items, where `n` is a local variable.  Then you can write a conveyor
loop of the following form.

    packet_t p;
    int64_t i = 0, pe;
    while (convey_advance(c, i == n)) {
        for (; i < n; i++) {
            p = (packet_t) { ... };  // construct ith packet
            pe = ...;                // decide where to send it
            if (! convey_push(c, &p, &pe))
                break;
        }

        while (convey_pull(c, &p, &pe)) {
            // handle incoming packet p, sent by pe
        }
    }

The advance, push, and pull operations must be interleaved, and there
are many ways to do this, but the pattern above is common.  The
argument `i == n` tells `convey_advance` whether this process has
finished pushing its items.  Note that for conveyor operations, a
return value of 0 means "failure" or "stop", while a positive return
value means "success" or "continue".

Finally, the session must be terminated by performing a collective
operations: either `reset`, which preserves the conveyor for reuse in future
sessions, or `free`, which destroys it, or both.

    convey_reset(c);
    // ...and then sooner or later...
    convey_free(c);


### Queries With Replies

If some or all messages will trigger a response from the remote pe, then
you need to use two conveyors: one for queries and one for replies.
Multiplexing queries and replies onto a single conveyor is not compatible
with the defined conveyor behavior and could lead to deadlock.
Fortunately, conveyors are inexpensive to create and use.  (Alternatively,
you could use the experimental _biconveyor_ interface, which is new in
version 0.6.)

Assuming that you have created two conveyors `q` and `r`, and have defined
datatypes `query_t` and `reply_t`, then your query/reply loop may look
like this:

    convey_begin(q, sizeof(query_t), alignof(query_t));
    convey_begin(r, sizeof(reply_t), alignof(reply_t));

    query_t query;
    reply_t reply;
    int64_t i = 0, pe;

    while (convey_advance(r, !convey_advance(q, i == n))) {
        for (; i < n; i++) {
            query = (query_t) { ... };  // construct ith query
            pe = ...;                   // decide where to send it
            if (! convey_push(q, &query, &pe))
                break;
        }

        while (convey_pull(q, &query, &pe)) {
            reply = (reply_t) { ... };  // construct reply
            if (! convey_push(r, &reply, pe)) {
                convey_unpull(q);
                break;
            }
        }

        while (convey_pull(r, &reply, NULL)) {
            // handle the reply
        }
    }

    convey_reset(r);
    convey_reset(q);

In this situation, the `convey_unpull` operation is used to "put back" a
query whose reply cannot yet be sent.  Only the most recently pulled item
can be unpulled.  When pulling a reply, you can pass `NULL` as the third
argument of `convey_pull` to indicate that you do not care which process
sent the reply.  Finally, note that replies can arrive in any order.
Thus, in code of this type, the query and reply structures normally
include some sort of "token" or sequence number that is copied from
the query to the reply and that the receiver of the reply can use to
identify the corresponding query.  (Biconveyors remove the need for
such tokens; they deliver replies in the order of the corresponding
queries.)

### Alternative Idiom

Often it is desirable to separate the code that handles messages from the
code that sends messages, particularly if both operations are complex.
To achieve this separation, you must usually encapsulate the information
needed to handle messages in some structure, say of type `context_t`,
which may include the conveyor or conveyors involved.  You can then write
the following kind of code, illustrated here for messages without replies.

```
int absorb_messages(context_t* ctx, bool done) {
    packet_t p;
    int64_t from;
    while (convey_pull(ctx->conveyor, &p, &from)) {
        // handle incoming packet p
    }
    return convey_advance(ctx->conveyor, done);
}

    // to send messages:
    context_t _ctx = { ... };
    convey_begin(_ctx.conveyor, sizeof(packet_t), alignof(packet_t));

    for (int i = 0; i < n; i++) {
        p = (packet_t) { ... };  // construct ith packet
        pe = ...;                // decide where to send it
        while (! convey_push(_ctx.conveyor, &p, &pe))
            absorb_messages(&_ctx, false);
    }
    while (absorb_messages(&_ctx, true))
        ;

    convey_reset(_ctx.conveyor);
```


Additional Features
-------------------

Conveyors offer a few standard features not mentioned above, including
statistics gathering and some options for error handling.  In addition,
conveyors support three optional features that are important in some
applications:

+ **Variable-sized Items.** An _elastic_ conveyor that can handle arbitrary
  item sizes (up to a given bound) can be created with
  `convey_new_elastic`.  Items are then pushed using `convey_epush`
  and pulled using `convey_epull`.

+ **Steady Progress.** A standard conveyor does not flush its partially
  filled buffers until the applications signals that it has finished pushing
  items (by passing `true` to `convey_advance`). But a conveyor can be
  constructed with the option `convey_opt_PROGRESS`, in which case it
  does not require this signal to guarantee forward progress.

+ **Compression.** Some types of conveyors (but not elastic conveyors as
  yet) support an optional compression mechanism that squeezes out
  constant bits or bytes from transmitted items. An application can
  request this behavior using the constructor options `convey_opt_COMPRESS`.
  Whether it gains performance will depend on the amount of redundancy
  in the items and on the ratio of available compute cycles to
  communication bandwidth.
