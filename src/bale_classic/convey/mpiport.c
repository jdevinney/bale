// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of the conveyor package. For license information,
// see the LICENSE file in the top level directory of the distribution.


#include <stdlib.h>
#include "porter_impl.h"
#include "private.h"

// Keep track of porter IDs in order to have unique MPI tags
static uint64_t used_ids = 0;

// For simplicity, we assume 2 send buffers per link for now.
typedef struct mpi_porter {
  porter_t porter;
  // MPI-specific stuff
  MPI_Comm comm;
  MPI_Request* send_req;
  MPI_Request* recv_req;
  int porter_id;        // tag for message matching
  int mpi_error;        // keep track of first serious error
  // State and statistics
  int n_pending;        // count of buffers waiting to be taken
  int i_pending;        // which of these to consider next
  int* pending;         // sources of pending buffers (first n_pending valid)
  buffer_t* taken;      // buffer borrowed most recently
  int n_active;         // counts open incoming channels
  bool* finished;       // records closed incoming channels
} mpi_porter_t;


/*** Memory Functions ***/

static void*
symmetric_alloc(void* alc8r, size_t size, const char* tag, uint64_t value)
{
  void* ptr = NULL;
  int err = posix_memalign(&ptr, CONVEY_MAX_ALIGN, size);
  return err ? NULL : ptr;
}

static void
symmetric_free(void* alc8r, void* ptr)
{
  free(ptr);
}

static const convey_alc8r_t porter_alc8r = {
  .alc8r = NULL, .grab = &symmetric_alloc, .free = &symmetric_free,
};

static inline buffer_t*
mpiport_inbuf(porter_t* self, int pe)
{
  return (buffer_t*) (self->recv_buffers + pe * self->buffer_stride);
}


/*** Internal Methods ***/

static void
mpiport_test_code(mpi_porter_t* mpip, int rc)
{
  if (rc != MPI_SUCCESS && mpip->mpi_error == MPI_SUCCESS) {
    mpip->mpi_error = rc;
    // Report the first error
    char string[MPI_MAX_ERROR_STRING+1];
    int length;
    MPI_Error_string(rc, string, &length);
    string[length] = 0;
    mprint(MY_PROC, 0, "MPI ERROR: %s\n", string);
  }
}

// Here 0 <= source < n = mpip->porter.n_ranks
static void
mpiport_try_recv(mpi_porter_t* mpip, int source)
{
  porter_t* self = &mpip->porter;
  int rc = MPI_Irecv(mpiport_inbuf(self, source), self->buffer_bytes,
                     MPI_BYTE, self->relative[source],
                     mpip->porter_id, mpip->comm, &mpip->recv_req[source]);
  mpiport_test_code(mpip, rc);
}

static bool
mpiport_ready(porter_t* self, int dest, uint64_t emitted)
{
  return self->channels[dest].delivered >= emitted;
}

static bool
mpiport_send(porter_t* self, int dest, uint64_t level,
             size_t n_bytes, buffer_t* buffer, uint64_t signal)
{
  // Use buffer->start as a completion flag
  buffer->start = signal & 1;
  if (n_bytes == 0)
    n_bytes = sizeof(buffer_t);

  DEBUG_PRINT("sending %zu bytes%s to %d\n",
              n_bytes, (signal & 1) ? "!" : "", dest);
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  int rc = MPI_Isend(buffer, n_bytes, MPI_BYTE, self->relative[dest],
                     mpip->porter_id, mpip->comm, &mpip->send_req[dest]);
  mpiport_test_code(mpip, rc);

  self->send_count++;
  self->byte_count += n_bytes;
  return false;  // delivery is asynchronous
}

static bool
mpiport_progress(porter_t* self, int dest)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;

  if (dest >= 0) {
    // Check whether an outstanding send has completed
    MPI_Request* req = &mpip->send_req[dest];
    if (*req != MPI_REQUEST_NULL) {
      MPI_Status _status;
      int ready = 0;
      int rc = MPI_Test(req, &ready, &_status);
      mpiport_test_code(mpip, rc);
      if (!ready)
        return false;
      porter_record_delivery(self, dest, self->channels[dest].delivered + 1);
    }
    return true;
  }

  // Try to complete all outstanding sends
  int n_ranks = self->n_ranks;
  MPI_Status status[n_ranks];
  int index[n_ranks];
  int outcount;
  // FIXME: use MPI_STATUSES_IGNORE
  int rc = MPI_Testsome(n_ranks, mpip->send_req, &outcount, index, status);
  mpiport_test_code(mpip, rc);
  if (outcount != MPI_UNDEFINED)
    for (int i = 0; i < outcount; i++)
      porter_record_delivery(self, index[i],
                             self->channels[index[i]].delivered + 1);

  // Check whether all sends are complete
  bool done = true;
  for (int i = 0; done && i < n_ranks; i++) {
    channel_t* channel = &self->channels[i];
    done = (channel->delivered == channel->emitted);
  }
  return done;
}


/*** Public Methods ***/

static bool
mpiport_setup(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  size_t n = self->n_ranks;
  bool ok = !self->dynamic || porter_grab_buffers(self);

  if (ok) {
    int n_active = 0;
    for (int i = 0; i < n; i++) {
      bool active = (self->relative[i] >= 0);
      n_active += active;
      mpip->send_req[i] = MPI_REQUEST_NULL;
      mpip->recv_req[i] = MPI_REQUEST_NULL;
      mpip->finished[i] = !active;
    }

    // Current state
    mpip->mpi_error = MPI_SUCCESS;
    mpip->n_pending = 0;
    mpip->i_pending = 0;
    mpip->taken = NULL;
    mpip->n_active = n_active;
  }

  mpp_barrier(1);
  // Initiate nonblocking receive operations
  if (ok && mpip->n_active)
    for (int i = 0; i < n; i++) {
      if (! mpip->finished[i])
        mpiport_try_recv(mpip, i);
    }
  return ok;
}

static buffer_t*
mpiport_borrow(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;

  if (mpip->i_pending == mpip->n_pending) {
    int outcount;
    int rc = MPI_Testsome(self->n_ranks, mpip->recv_req, &outcount,
                          mpip->pending, MPI_STATUSES_IGNORE);
    mpiport_test_code(mpip, rc);

    // Check for closure messages here, not on returned buffers
    int j = 0;
    if (outcount != MPI_UNDEFINED)
      for (int i = 0; i < outcount; i++) {
        int source = mpip->pending[i];
        buffer_t* buffer = mpiport_inbuf(self, source);
        if (buffer->start) {
          mpip->finished[source] = true;
          mpip->n_active--;
          buffer->start = 0;
        }
        if (buffer->limit > 0)
          mpip->pending[j++] = source;
      }

    mpip->n_pending = j;
    mpip->i_pending = 0;
  }

  if (mpip->i_pending < mpip->n_pending) {
    int source = mpip->pending[mpip->i_pending];
    buffer_t* taken = mpiport_inbuf(self, source);
    taken->source = source;
    mpip->taken = taken;
    if (self->compress && taken->n_items > 0)
      porter_decompress(self, taken, source);
    return taken;
  }

  if (mpip->n_active == 0)
    self->drained = true;
  return NULL;
}

static void
mpiport_return(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  buffer_t* taken = mpip->taken;
  mpip->taken = NULL;

  if (taken->start >= taken->limit) {
    int i = mpip->i_pending++;
    int source = mpip->pending[i];
    if (! mpip->finished[source])
      mpiport_try_recv(mpip, source);
  }
}

static void
mpiport_reset(porter_t* self)
{
  mpp_barrier(1);
  if (self->dynamic)
    porter_free_buffers(self);
}

static void
mpiport_demolish(porter_t* self)
{
  mpi_porter_t* mpip = (mpi_porter_t*) self;
  if (mpip->porter_id >= 0)
    used_ids &= ~(UINT64_C(1) << mpip->porter_id);
  porter_free_buffers(self);
  free(mpip->finished);
  free(mpip->pending);
  free(mpip->recv_req);
  free(mpip->send_req);
}

static const porter_methods_t mpi_porter_methods = {
  .setup = &mpiport_setup,
  .borrow = &mpiport_borrow,
  .turnin = &mpiport_return,
  .reset = &mpiport_reset,
  .demolish = &mpiport_demolish,
  .ready = &mpiport_ready,
  .send = &mpiport_send,
  .progress = &mpiport_progress,
};


/*** Constructor and Destructor ***/

porter_t*
porter_new(int n, int32_t relative[n], int my_rank,
           size_t tag_bytes, size_t capacity, size_t multiplicity,
           const convey_alc8r_t* alloc, uint64_t options, int opcode)
{
  if (n <= 0 || multiplicity == 0 || (multiplicity & (multiplicity - 1)))
    return NULL;
  const bool dynamic = options & convey_opt_DYNAMIC;
  const bool steady = options & convey_opt_PROGRESS;
  const bool local = options & porter_opt_LOCAL;
  const bool compress = options & convey_opt_COMPRESS;

  mpi_porter_t* mpip = malloc(sizeof(mpi_porter_t));
  if (mpip == NULL)
    return NULL;
  porter_t* porter = &mpip->porter;

  if (alloc == NULL)
    alloc = &porter_alc8r;

  const size_t align_1 = CONVEY_MAX_ALIGN - 1;
  const size_t buffer_stride = (sizeof(buffer_t) + capacity + align_1) & ~align_1;
  *mpip = (mpi_porter_t) { .porter = {
      ._class_ = &mpi_porter_methods,
      .tag_bytes = tag_bytes, .buffer_stride = buffer_stride,
      .n_ranks = n, .my_rank = my_rank, .abundance = 1,  // FIXME
      .relative = relative, .alloc = *alloc, .opcode = opcode,
      .dynamic = dynamic, .inmult = false, },
    .porter_id = -1,
  };
  mpip->comm = mpp_comm_mpi(MPP_COMM_CURR);

  // Local allocations
  porter->send_areas = malloc(n * sizeof(area_t));
  porter->all_sent = malloc(n * sizeof(bool));
  porter->channels = malloc(n * sizeof(channel_t));
  if (steady)
    porter->waiting = malloc(n * sizeof(uint8_t));
  mpip->send_req = malloc(n * sizeof(MPI_Request));
  mpip->recv_req = malloc(n * sizeof(MPI_Request));
  mpip->pending = malloc(n * sizeof(int));
  mpip->finished = malloc(n * sizeof(bool));
  bool ok = (porter->send_areas && porter->all_sent && porter->channels &&
             (!steady || porter->waiting) && mpip->send_req && mpip->recv_req &&
             mpip->pending && mpip->finished);
  if (ok && !dynamic)
    ok = porter_grab_buffers(porter);
  if (ok && compress && !local)
    ok = porter_make_codata(porter);

  ok = mpp_and_long(ok);
  if (!ok) {
    // could report an out-of-memory error here
    porter_demolish(porter);
    return NULL;
  }

  // Agree on an unused ID for this porter
  uint64_t avail = mpp_and_long(~used_ids);
  if (avail == 0) {
    // could report an out-of-ids error here
    porter_demolish(porter);
    return NULL;
  }
  int id = _trailz(avail);
  used_ids |= UINT64_C(1) << id;
  mpip->porter_id = id;

  // Standardize the relative[] array
  int32_t n_procs = PROCS;
  for (int i = 0; i < n; i++)
    if (relative[i] < 0 || relative[i] >= n_procs)
      relative[i] = -1;

  // porter_setup() will erase everything
  return porter;
}
