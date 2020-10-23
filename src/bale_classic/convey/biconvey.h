// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


/** \file biconvey.h
 * The public API for biconveyors.
 *
 * A biconveyor is a kind of bidirectional conveyor. If a conveyor is
 * a way to implement active messages, then a biconveyor is a way to
 * implement active messages with replies.
 *
 * Client code can use a biconveyor to push queries to remote PEs and to
 * pull replies. Each query gives rise to exactly one reply, which returns
 * to the PE that pushed the query. Replies are pulled in the order that
 * the corresponding queries were pushed. Every biconveyor "session",
 * bounded by calls to \c biconvey_begin() and \c biconvey_reset(),
 * specifies a mechanism to compute the reply to each query. As with a
 * conveyor, pushes and pulls can fail, so they must be placed in a loop
 * with calls to the advance operation.
 *
 * There are no elastic biconveyors as yet; the queries and replies in
 * a session must have fixed size.
 *
 * A biconveyor is guaranteed to make progress as long as every
 * participating PE continues to make calls to \c biconvey_advance().
 * Under this condition, the reply to every pushed query will eventually
 * become available to be pulled. PEs need not perform either pushes or
 * pulls to ensure progress; thus all biconveyors are "steady".
 */

#ifndef BICONVEY_H
#define BICONVEY_H

#include "convey.h"

/** The opaque type of a biconveyor. */
typedef struct biconveyor biconvey_t;


/*** CONSTRUCTORS ***/

/** The generic constructor for biconveyors.
 *
 * This function is similar in all respects to \c convey_new() except that
 * it builds a biconveyor rather than a conveyor. Biconveyors are steady by
 * definition, so the \c convey_opt_PROGRESS option is ignored.
 *
 * \param[in] max_bytes   Desired limit on internal memory usage (per PE)
 * \param[in] n_local     Number of PEs in a local group (e.g., node or socket)
 * \param[in] alloc       Means of obtaining symmetric memory; can be \c NULL
 * \param[in] options     An OR of any set of values from \c enum #convey_option
 */
biconvey_t*
biconvey_new(size_t max_bytes, size_t n_local,
	     const convey_alc8r_t* alloc, uint64_t options);

/** Construct a synchronous biconveyor akin to \c exstack.
 *
 * Build a synchronous biconveyor that uses a single underlying simple conveyor
 * for both queries and replies. This is a collective operation; under mpp_utilV4,
 * the set of PEs involved is determined by the current communicator.
 *
 * \param[in] capacity    Capacity of each buffer, measured in bytes
 * \param[in] alloc       Means of obtaining symmetric memory; can be \c NULL
 * \param[in] a2a         Strategy (\c mpp_alltoall_t) for \c mpp_alltoallv_simple;
 * if \c NULL, the conveyor tries to use an underlying SHMEM or
 * MPI mechanism, and falls back to a default \c mpp_utilV4 strategy
 * \param[in] options     An OR of any set of values from \c enum #convey_option
 */
biconvey_t*
biconvey_new_simple(size_t capacity, const convey_alc8r_t* alloc,
		    const convey_mpp_a2a_t* a2a, uint64_t options);

/** Build an asynchronous biconveyor in which each message makes \a order hops.
 *
 * This is a collective operation; under mpp_utilV4, the set of PEs
 * involved is determined by the current communicator. It builds two tensor
 * conveyors whose configuration is determined by the given parameters, and
 * combines these to form a biconveyor.
 *
 * \param[in] capacity    Size of each internal buffer, in bytes
 * \param[in] order       Number of porters per conveyor (hops each item takes),
 *                        1 <= order <= 3
 * \param[in] n_local     Number of PEs in a local group (ignored if order = 1)
 * \param[in] n_buffers   Internal buffer multiplicity, must be a power of two;
 *                        values 1, 2, 4 are reasonable
 * \param[in] alloc       Means of obtaining symmetric memory; can be \c NULL
 * \param[in] options     Supports all options except \c convey_opt_SCATTER; see
 *                        \c convey_opt_PROGRESS; see enum #convey_option
 */
biconvey_t*
biconvey_new_tensor(size_t capacity, int order, size_t n_local, size_t n_buffers,
		    const convey_alc8r_t* alloc, uint64_t options);


/*** BICONVEYOR METHODS ***/

/** Prepare a biconveyor for pushes and pulls.
 *
 * This function sets the size of queries pushed by \c biconvey_push()
 * to \a query_bytes and of replies pulled by \c biconvey_pull() to
 * \a reply_bytes.  It also supplies the function \a answer that the
 * biconveyor will call to obtain the reply to each query, along with
 * an arbitrary \a context that the biconveyor will pass to \a answer.
 *
 * The \a answer function must allow its \a query and \a reply to have
 * arbitrary alignment.  In practice, this means that it must access
 * these arguments via \c memcpy().
 */
int biconvey_begin(biconvey_t* self, size_t query_bytes, size_t reply_bytes,
		   void (*answer)(const void* query, void* reply, void* context),
		   void* context);

/** Enqueue the given \a query for delivery to the given \a pe.
 *
 * This function is similar in all respects to the conveyor function
 * \c convey_push(). The \a query has the size \a query_bytes given in
 * the most recent call to \c biconvey_begin().
 */
int biconvey_push(biconvey_t* self, const void* query, int64_t pe);

/** Attempt to dequeue a reply.
 *
 * This function is similar to \c convey_pull() except that it does not
 * provide the caller with the identity of the PE that sent the reply. The
 * \a reply has the size \a reply_bytes given in the most recent call to
 * \c biconvey_begin(). Replies are pulled in the order of their
 * corresponding pushes.
 */
int biconvey_pull(biconvey_t* self, void* reply);

/** Make forward progress and check for completion.
 *
 * Externally, this function is similar in all respects to the conveyor
 * function \a convey_advance(). Internally, it does quite a bit more
 * work. In particular, it accepts incoming queries, calls the \a answer
 * function, transmits the resulting replies, and if necessary reorders
 * incoming replies.
 */
int biconvey_advance(biconvey_t* self, bool done);

/** Restore the biconveyor to a pristine state.
 *
 * This is a collective operation that works like \c convey_reset().
 */
int biconvey_reset(biconvey_t* self);

/** Destroy the biconveyor and release all memory that it allocated.
 *
 * This is a collective operation that works like \c convey_free().
 */
int biconvey_free(biconvey_t* self);


#endif
