// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


/** \file convey.h
 * The public API for conveyors.
 */

#ifndef CONVEY_H
#define CONVEY_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include "convey_alc8r.h"
#include "convey_codec.h"

/** Non-fatal return values from conveyor methods. */
enum convey_status {
  convey_DONE = 0,         ///< returned by convey_advance()
  convey_FAIL = 0,         ///< returned by convey_push(), convey_pull(), etc.
  convey_OK   = 1,         ///< returned by convey_push(), convey_advance(), etc.
  convey_NEAR = 2,         ///< returned by convey_advance()
};

/** Optional features of conveyors. */
enum convey_feature {
  convey_ELASTIC = 0x01,   ///< supports variable-sized items (epush/epull)
  convey_STEADY = 0x02,    ///< delivers all pushed items without requiring termination
  convey_THRIFTY = 0x04,   ///< supports buffer compression for remote transfers
};

/** Statistics tracked by conveyors. */
enum convey_statistic {
  convey_BEGINS = 0,
  convey_PUSHES,
  convey_PULLS,
  convey_UNPULLS,
  convey_ADVANCES,
  convey_COMMS,         // transfers of buffers
  convey_SYNCS,         // barrier, fence, quiet, etc.
  convey_BYTES,
  convey_CUMULATIVE,    // add to any of the above
};

/** The opaque type of a conveyor. */
typedef struct conveyor convey_t;

/** The structure returned by an elastic pull. */
typedef struct {
  size_t bytes;         ///< number of bytes in the pulled item 
  void* data;           ///< aligned pointer to the pulled item
  int64_t from;         ///< index of the PE that pushed it
} convey_item_t;

/** Version string, such as "0.5.0". */
extern const char convey_version[];


/*** GENERIC CONSTRUCTORS ***/

/** Options for conveyor constructors. */
enum convey_option {
  /// Turn off some error checks (e.g., that the \a pe argument of
  /// convey_push() is in range) to gain speed.
  convey_opt_RECKLESS = 0x01,
  /// Allocate and deallocate large buffers in \c convey_begin() and \c
  /// convey_reset(), respectively, rather than in the constructor and \c
  /// convey_free().
  convey_opt_DYNAMIC = 0x02,
  /// Do not report severe errors by printing a message; just return the
  /// error codes.
  convey_opt_QUIET = 0x04,
  /// Optimize for the case in which successive pushes tend to target
  /// random PEs, rather than dwelling on the same PE for a while.
  convey_opt_SCATTER = 0x08,
  /// Obtain a "steady" conveyor: one that ensures that all pushed items
  /// are delivered, provided only that every PE continues to call \c
  /// convey_pull() and \c convey_advance(); otherwise delivery is only
  /// guaranteed after all PEs have declared that they are done pushing.
  convey_opt_PROGRESS = 0x10,
  /// Print warning messages if the conveyor detects any unexpected
  /// condition that may reduce performance.
  convey_opt_ALERT = 0x20,
  /// Turn on data compression, if possible, for non-local hops. Standard
  /// compression methods will be used unless \c convey_set_codec() is
  /// called.
  convey_opt_COMPRESS = 0x40,
};

/** Expert options for conveyor constructors. */
enum convey_expert_option {
  /// Turn off optimizations that alter the format of internal buffers.
  /// This option is provided for conveyor testing and maintenance.
  convey_opt_STANDARD = 0x100000,
  /// Avoid the use of nonblocking SHMEM/UPC puts. This option is provided
  /// for testing and automatic tuning.
  convey_opt_BLOCKING = 0x200000,
};

/** The maximum alignment that a conveyor can guarantee. */
#define CONVEY_MAX_ALIGN (UINT64_C(64))

/** Strategies for \c mpp_alltoall_simple(); unlikely to be useful
 * 
 * Define \c convey_mpp_a2a_t as a synonym for \c mpp_alltoall_t so we don't
 * have to include mpp_utilV4.h (in UPC mode, for instance).
 */
typedef struct mpp_alltoall convey_mpp_a2a_t;

/** The generic constructor for ordinary conveyors.
 *
 * Build a conveyor, trying to maximize performance subject to the given
 * limit on the number of bytes per PE used internally by the conveyor.
 * If no available conveyor is compatible with this limit, then build
 * a conveyor that uses as little memory as possible. Setting \a max_bytes
 * to \c SIZE_MAX removes the limit; this value is generally recommended,
 * as efficient conveyors tend to use little memory anyway.
 *
 * This function uses simple heuristics, which may not be reliable, to
 * choose a conveyor and its parameters. If the environment variable
 * \c CONVEY_BUFFER_SIZE is set to a positive integer, then the conveyor
 * will prefer to transfer data in buffers of approximately that size.
 *
 * If \a n_local is zero, then the constructor will use a default value.
 * The value is obtained from the environment variable \c CONVEY_NLOCAL if
 * possible. Otherwise, the conveyor package tries to determine the number
 * of processes per node. If it fails, or if the resulting number does not
 * divide the number of processes, then \a n_local defaults to 1.
 *
 * \param[in] max_bytes   Desired limit on internal memory usage (per PE)
 * \param[in] n_local     Number of PEs in a local group (e.g., node or socket)
 * \param[in] alloc       Means of obtaining symmetric memory; can be \c NULL
 * \param[in] options     An OR of any set of values from \c enum #convey_option
 */

convey_t*
convey_new(size_t max_bytes, size_t n_local,
           const convey_alc8r_t* alloc, uint64_t options);

/** The generic constructor for elastic conveyors.
 *
 * Build a conveyor that supports variable-length items, up to a limit of
 * \a item_bound bytes per item. In other respects this function is the same
 * as \c convey_new().
 *
 * \param[in] item_bound  Maximum size of each item in bytes
 * \param[in] max_bytes   Desired limit on internal memory usage (per PE)
 * \param[in] n_local     Number of PEs in a local group (e.g., node or socket)
 * \param[in] alloc       Means of obtaining symmetric memory; can be \c NULL
 * \param[in] options     An OR of any set of values from \c enum #convey_option
 */

convey_t*
convey_new_elastic(size_t item_bound, size_t max_bytes, size_t n_local,
                   const convey_alc8r_t* alloc, uint64_t options);
                    


/*** CONVEYOR METHODS ***/

/** Prepare a conveyor for pushes and pulls.
 *
 * This function sets the item size for subsequent \c convey_push() and
 * \c convey_pull() operations to \a item_bytes, and it requests a certain
 * alignment for pointers supplied by \c convey_apull() and
 * \c convey_epull(). If \a align is not valid (it is not a power of 2, or
 * does not divide \a item_byts, or is larger than \c CONVEY_MAX_ALIGN),
 * then the required alignment is set to the largest power of 2, no greater
 * than \c CONVEY_MAX_ALIGN, that divides \a item_bytes. If the conveyor
 * was created with a non-NULL allocator, then the conveyor's alignment
 * guarantee is only as strong as the allocator's alignment guarantee.
 *
 * The state must be \e dormant on every process, and it changes to \e working.
 * This is a collective operation. It may check that the conveyor was correctly
 * built, e.g., that the same parameters were supplied by every process. The
 * return value is \c convey_OK on success, negative on failure.
 */
int convey_begin(convey_t* c, size_t item_bytes, size_t align);

/** Enqueue the given \a item for delivery to the given \a pe.
 *
 * This operation is legal only if the local state of \a c is \e working.
 * The return value is \c convey_OK on success, or \c convey_FAIL if the
 * item could not be pushed. A negative return value means that the
 * function was used incorrectly, e.g., the \a pe argument was out of
 * range, or the conveyor is not in the \e working state. If the return
 * value is \c convey_OK, then the item has been copied into the conveyor's
 * internal storage.
 */
int convey_push(convey_t* c, const void* item, int64_t pe);

/** Enqueue the given \a item of size \a bytes for delivery to the given \a pe.
 *
 * If the conveyor is not elastic (see enum #convey_feature), or if \a bytes
 * is greater than the conveyor's maximum item size, then then this function
 * returns a negative value. It is otherwise the same as \c convey_push().
 * Items of size zero (\a bytes == 0) are permitted, and will be delivered.
 */
int convey_epush(convey_t* c, size_t bytes, const void* item, int64_t pe);

/** Attempt to dequeue an item, optionally learning the PE it was sent \a from.
 *
 * This operation is always legal unless the local state is \e dormant. If
 * it succeeds, then it copies the incoming item into \a *item, writes
 * the index of the source PE into \a *from unless \a from is \c NULL,
 * and returns \c convey_OK. If it fails, it does not change \a *item or
 * \a *from, and it returns \c convey_FAIL. A negative return value
 * indicates that the operation is illegal or that an internal error occurred.
 *
 * The return value of \c convey_FAIL from \c convey_pull does not
 * necessarily imply that no items are available to pull. The only correct
 * way to ascertain that is to observe a return value of \c convey_DONE
 * from \c convey_advance.
 *
 * In an elastic conveyor, \c convey_pull() will fail unless the next item
 * to be pulled has the size used by \c convey_push().
 */
int convey_pull(convey_t* c, void* item, int64_t* from);

/** Attempt to dequeue an item, obtaining a properly aligned pointer to it.
 *
 * This function is very similar to \c convey_pull(), but instead of copying
 * the incoming item, it returns a pointer to it, which may be faster. The
 * pointer is aligned as determined by the most recent call to
 * \c convey_begin(). The pointer remains valid until the next successful pull
 * or unpull operation on the conveyor; after that, the content of the item
 * may change. If \c convey_apull() fails for any reason, then it returns
 * \c NULL and does not change \a *from.
 */
void* convey_apull(convey_t* c, int64_t* from);

/** Attempt to dequeue an item of arbitrary size.
 *
 * This function is similar to \c convey_apull() except that it returns a
 * status code and, if successful, writes information about the pulled item
 * into \a *result. Specifically, on success, \a result->bytes contains
 * the item size, \a result->data points to the item, and \a result->from
 * holds index of the source PE. If the conveyor is not elastic (see enum
 * #convey_feature), then the operation fails with a negative return value.
 */
int convey_epull(convey_t* c, convey_item_t* result);

/** Attempt to restore a previously pulled item to the conveyor.
 *
 * The item in question is the most recently pulled item that has not been
 * unpulled. If the unpull succeeds, then the item can be pulled again
 * later. This operation is legal whenever the local state is not \e
 * dormant, but it is only guaranteed to succeed if the last operation on
 * \a c was a successful \c convey_pull(), \c convey_apull(), or \c
 * convey_epull(). On success, the return value is \c convey_OK; on
 * failure, the return value is \c convey_FAIL.
 */
int convey_unpull(convey_t* c);

/** Make forward progress and check for completion.
 *
 * If \a done is \c true, declare that the caller will no longer push any
 * items (prior to the next \c convey_reset). In any case, ask the conveyor
 * to make progress. The return value is \c convey_DONE if the local state
 * has become \e complete, meaning that there are no items to pull and no
 * more items will be delivered (until after the next \c convey_reset). The
 * return value is \c convey_NEAR if the local state has become \e cleanup,
 * meaning that no more items will be delivered but items may remain to be
 * pulled. Otherwise the return value is \c convey_OK, or a negative value
 * if an error occurs -- for instance, the conveyor detects that the caller
 * has violated the contract. If the local state of \a c is not \e working
 * when \c convey_advance is called, then the \a done argument must be \c
 * true.
 */
int convey_advance(convey_t* c, bool done);

/** Restore the conveyor to a pristine state, aside from its statistics.
 *
 * This is a collective operation. The \c convey_reset operation is legal
 * if and only if the local state of \a c is \e complete. If the operation
 * is legal on every process and is called collectively, then it succeeds
 * and changes the local state to \e dormant. The return value is \c
 * convey_OK on success, negative on failure.
 */
int convey_reset(convey_t* c);

/** Destroy the conveyor and release all memory that it allocated.
 *
 * This is a collective operation. The \c convey_free operation is legal if
 * and only if the local state of \a c is \e complete or \e dormant. If the
 * operation is legal on every process and is called collectively, then it
 * succeeds. The return value is \c convey_OK on success, negative on
 * failure. Calling \c convey_free on a \c NULL conveyor always succeeds.
*/
int convey_free(convey_t* c);

/** Set up the methods used to compress and decompress buffers.
* 
* This is a collective operation that can only be performed when the state
* is \e dormant (that is, before the first \c convey_begin, or between
* \c convey_reset and \c convey_begin). If \a codec is \c NULL, then
* compression is disabled. The return value is \c convey_OK on success,
* \c convey_FAIL if the conveyor does not support compression, and a negative
* value if the state is wrong or some other severe error occurs.
*/
int convey_set_codec(convey_t* c, const convey_codec_t* codec, void* arg);

/** Obtain the size, in bytes, of the items that this conveyor handles.
 *
 * If the conveyor is \e dormant, then this function returns 0.
 */
size_t convey_item_size(convey_t* c);

/** Obtain a bitmask that describes the optional features supported by this conveyor.
 *
 * Such features (see enum #convey_feature) cannot weaken the promises made by the
 * conveyor as part of the contract, but they can strengthen them.
 */
uint64_t convey_features(convey_t* c);

/** Obtain one of a fixed set of statistics that conveyors may track.
 *
 * Some statistics are cleared by \c convey_reset; some accumulate over the
 * lifetime of the conveyor. A return value of -1 means that the conveyor does
 * not track that statistic. If the state of the conveyor is not \e dormant,
 * then the statistic may not be up to date.
 */
int64_t convey_statistic(convey_t* c, int which);

/** Convert an error code to an immutable string.
 *
 * Given a conveyor and a negative value returned by one of that conveyor's
 * methods, return a pointer to an immutable string that describes the error.
 */
const char* convey_error_string(convey_t* c, int error);


/*** SPECIFIC CONSTRUCTORS ***/

/** Construct a basic synchronous conveyor akin to \c exstack.
 *
 * Build a synchronous conveyor based on \c shmem_alltoallv, \c MPI_Alltoallv,
 * or \c mpp_alltoallv_simple. This is a collective operation; under mpp_utilV4,
 * the set of PEs involved is determined by the current communicator.
 *
 * \param[in] capacity    Capacity of each buffer, measured in bytes
 * \param[in] alloc       Means of obtaining symmetric memory; can be NULL
 * \param[in] a2a         Strategy (\c mpp_alltoall_t) for \c mpp_alltoallv_simple;
 * if \c NULL, the conveyor tries to use an underlying SHMEM or
 * MPI mechanism, and falls back to a default \c mpp_utilV4 strategy
 * \param[in] options     An OR of any set of values from enum #convey_option
 *
 * Conveyor Features: Perhaps \c convey_STEADY, depending on \a options.
 *
 * Error Conditions: If one of these errors occurs, then the constructor
 * returns \c NULL and may print an error message.
 *  - Capacity is zero
 *  - Allocator is \c NULL, communicator is not \c MPP_COMM_WORLD, and
 *    \c mpp_alloc() does not support teams
 *  - Unable to allocate memory
 *  - Strategy is \c NULL, library is MPI, and \a capacity * \c PROCS
 *    is too large (2 GiB or more)
 */

convey_t*
convey_new_simple(size_t capacity, const convey_alc8r_t* alloc,
                  const convey_mpp_a2a_t* a2a, uint64_t options);


/** Build a synchronous conveyor based on \c shmem_team_alltoallv or \c MPI_Alltoallv.
 *
 * \deprecated Twohop conveyors, even when they exist, are generally
 * dominated by tensor conveyors. They are not currently chosen by generic
 * constructors and will be removed in a future release.
 * 
 * This is a collective operation; under \c mpp_utilV4,
 * the set of PEs involved is determined by the current
 * communicator. The PEs are arranged in a matrix with \a row_procs PEs per
 * row, and the PE numbers within a row are consecutive. Each message
 * takes two hops: first within its source row and then within its
 * destination column. Twohop conveyors require either MPI or Cray SHMEM.
 *
 * \param[in] capacity    Minimum size in bytes for internal buffers
 * \param[in] row_procs   Number of processes per row (must divide \c PROCS)
 * \param[in] alloc       Means of obtaining symmetric memory; can be \c NULL
 * \param[in] options     Supports all options except \c convey_opt_SCATTER and
 *                        \c convey_opt_COMPRESS; see enum #convey_option
 *
 * Conveyor Features: Perhaps \c convey_STEADY, depending on \a options.
 */

convey_t*
convey_new_twohop(size_t capacity, size_t row_procs,
                  const convey_alc8r_t* alloc, uint64_t options)
#if defined(__GNUC__) && !defined(__UPC__)
__attribute__((deprecated))
#endif
;


/** Build an asynchronous conveyor in which each item makes \a order hops.
 *
 * This is a collective operation; under \c mpp_utilV4, the set of PEs
 * involved is determined by the current communicator. Tensor conveyors
 * are optimized for Cray SHMEM but can work atop SHMEM, UPC, or MPI.
 *
 * Internally, tensor conveyors use subsidiary objects called \e porters
 * that are akin to \c exstack2. Individual porters are not designed to
 * scale beyond a few hundred PEs, but by using multiple porters, the
 * overall conveyor can scale to hundreds of thousands of PEs. Each item
 * makes a certain number of hops called the \a order, where \a order is 1,
 * 2, or 3. If \a order is 2 or 3, then the first hop and (if \a order is 3)
 * last hop take place within a group of \a n_local consecutive PEs, which
 * should correspond to a hardware unit with fast communication, such as a
 * node or a socket. In these cases, \a n_local must divide the number of
 * participating PEs.
 *
 * \param[in] capacity    Size of each internal buffer, in bytes
 * \param[in] order       Number of porters (hops each item takes), 1 <= order <= 3
 * \param[in] n_local     Number of PEs in a local group (ignored if order = 1)
 * \param[in] n_buffers   Internal buffer multiplicity, must be a power of two;
 *                        values 1, 2, 4 are reasonable
 * \param[in] alloc       Means of obtaining symmetric memory; can be \c NULL
 * \param[in] options     Supports all options except \c convey_opt_SCATTER; see
 *                        enum #convey_option
 *
 * Conveyor Features: Perhaps \c convey_STEADY and/or \c convey_THRIFTY,
 * depending on the requested \a options.
 */

convey_t*
convey_new_tensor(size_t capacity, int order, size_t n_local, size_t n_buffers,
                  const convey_alc8r_t* alloc, uint64_t options);


/** Build an asynchronous elastic conveyor in which each item makes \a order hops.
 *
 * This is a collective operation; under \c mpp_utilV4, the set of PEs
 * involved is determined by the current communicator.
 *
 * Except for its first two arguments, this function is similar to \c
 * convey_new_tensor(). The \a monster_bytes parameter can be greater than
 * \a buffer_bytes. In this case the conveyor uses a separate mechanism, with
 * its own large buffers, for sending items of size greater than \a buffer_bytes.
 * Delivery of such large items may be less efficient, but they behave in
 * all other respects like small items.
 *
 * \param[in] buffer_bytes  Approximate size of each internal buffer
 * \param[in] monster_bytes Maximum item size for \c convey_epush()
 * \param[in] order         Number of porters (hops each item takes), 1 <= order <= 3
 * \param[in] n_local       Number of PEs in a local group (ignored if order = 1)
 * \param[in] n_buffers     Internal buffer multiplicity, must be a power of two;
 *                          values 1, 2, 4 are reasonable
 * \param[in] alloc         Means of obtaining symmetric memory; can be \c NULL
 * \param[in] options       Supports all options except \c convey_opt_SCATTER and
 *                          \c convey_opt_COMPRESS; see enum #convey_option
 *
 * Conveyor Features: Always \c convey_ELASTIC and sometimes \c convey_STEADY,
 * depending on the \a options.
 */

convey_t*
convey_new_etensor(size_t buffer_bytes, size_t monster_bytes,
                   int order, size_t n_local, size_t n_buffers,
                   const convey_alc8r_t* alloc, uint64_t options);


#endif
