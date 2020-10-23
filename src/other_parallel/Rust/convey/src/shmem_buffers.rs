//! Conveyor module for handling a shmem buffers
//!
//! Copyright (c) 2020, Institute for Defense Analyses
//! 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//!
//! All rights reserved.
//!
//! This file is part of Convey, a conveyor library for rust.  For
//! licence information see the LICENSE file in the top level dirctory
//! of the distribution.
use crate::Convey;
use shmem::atomic::Atomic;
use shmem::object::Object;
use std::cell::RefCell;

/// Shmem buffers handle communication to other ranks
pub struct ShmemBuffers {
    /// my rank
    pub my_rank: usize,
    /// the number of ranks
    pub num_ranks: usize,
    /// Recieve area, sized by number or ranks
    pub receive_buf: Object<u8>,
    /// Send area, also sized by number of ranks so we can use non-blocking sends
    pub send_buf: Object<u8>,
    /// valid flag for receive, on receiver, set by sender, cleared by receiver
    /// also communicates phase information
    pub receive_valid: Object<i64>,
    /// The goal we are gettng
    pub receive_goal: Object<u32>,
    /// busy flag for receive, on sender, set by receiver, cleared by sender
    pub receive_busy: Object<u64>,
}

struct ShmemBufferPool {
    items: RefCell<Vec<ShmemBuffers>>,
}

impl ShmemBufferPool {
    fn new() -> Self {
        Self {
            items: RefCell::new(Vec::new()),
        }
    }
    fn get(&self, convey: &Convey) -> ShmemBufferGuard {
        let item = match self.items.borrow_mut().pop() {
            Some(item) => item,
            None => {
                let num_buf_elements = convey.num_ranks * convey.send_max;
                let send_buf = convey
                    .shmem
                    .new_object(num_buf_elements)
                    .expect("shmem::new_object failed");
                let receive_buf = convey
                    .shmem
                    .new_object(num_buf_elements)
                    .expect("shmem::new_object failed");
                let receive_valid = convey
                    .shmem
                    .new_object(convey.num_ranks)
                    .expect("shmem::new_object failed");
                let receive_busy = convey
                    .shmem
                    .new_object(convey.num_ranks)
                    .expect("shmem::new_object failed");
                let receive_goal = convey
                    .shmem
                    .new_object(convey.num_ranks)
                    .expect("shmem::new_object failed");
                ShmemBuffers {
                    my_rank: convey.my_rank,
                    num_ranks: convey.num_ranks,
                    receive_buf,
                    send_buf,
                    receive_valid,
                    receive_busy,
                    receive_goal,
                }
            }
        };
        ShmemBufferGuard { inner: Some(item) }
    }
}

/// The actual shmem buffer
#[derive(Debug)]
pub struct ShmemBufferGuard {
    inner: Option<ShmemBuffers>,
}

impl Drop for ShmemBufferGuard {
    fn drop(&mut self) {
        let item: ShmemBuffers = self.inner.take().unwrap();
        //item.reset();
        REUSE_BUF.with(|rb| rb.items.borrow_mut().push(item));
    }
}

impl std::ops::Deref for ShmemBufferGuard {
    type Target = ShmemBuffers;
    fn deref(&self) -> &Self::Target {
        self.inner.as_ref().unwrap()
    }
}

impl std::ops::DerefMut for ShmemBufferGuard {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.inner.as_mut().unwrap()
    }
}

thread_local! {
    static REUSE_BUF:ShmemBufferPool = ShmemBufferPool::new();
}

// we break out shmem buffers so we can reuse them.  They are
// independent of the type of ConveySession
impl ShmemBuffers {
    /// generate a new ShmemBuffer, maybe reusing one which was just freed
    pub fn new(convey: &Convey) -> ShmemBufferGuard {
        convey.barrier();
        REUSE_BUF.with(|rb| rb.get(convey))
    }
    /// is a receiver available for to_rank from from_rank
    pub fn receive_is_available(&self, to_rank: usize, from_rank: usize) -> bool {
        if from_rank == self.my_rank {
            self.receive_busy.local_part()[to_rank] == 0
        } else {
            self.receive_busy.atomic_fetch(to_rank, from_rank).unwrap() == 0
        }
    }

    /// set receiver available for to_rank from from_rank
    pub fn set_receive_available(&self, to_rank: usize, from_rank: usize, val: bool) -> () {
        //println!("set receive avail {} {} {}", to_rank, from_rank, val);
        let busy = match val {
            true => 0,
            false => 1,
        };
        self.receive_busy
            .atomic_set(to_rank, busy, from_rank)
            .unwrap();
    }

    /// is a receiver valid for to_rank from from_rank
    pub fn receive_is_valid(&self, to_rank: usize, from_rank: usize) -> i64 {
        if to_rank == self.my_rank {
            self.receive_valid.local_part()[from_rank]
        } else {
            self.receive_valid.atomic_fetch(from_rank, to_rank).unwrap()
        }
    }

    /// set receiver valid for to_rank from from_rank
    pub fn set_receive_valid(&self, to_rank: usize, from_rank: usize, val: i64) -> () {
        self.receive_valid
            .atomic_set(from_rank, val, to_rank)
            .unwrap();
    }

    /// set goal for to_rank from from_rank
    pub fn set_my_goal_on_rank(&self, their_rank: usize, my_rank: usize, val: u32) -> () {
        self.receive_goal
            .atomic_set(my_rank, val, their_rank)
            .unwrap();
    }
    /// get goal for to_rank from from_rank
    pub fn get_their_goal_on_my_rank(&self, my_rank: usize, their_rank: usize) -> u32 {
        if my_rank == self.my_rank {
            self.receive_goal.local_part()[their_rank]
        } else {
            panic!("never get remote goal");
        }
    }
}

impl std::fmt::Debug for ShmemBuffers {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ShmemBuffers").finish()
    }
}
