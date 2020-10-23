#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Library for conveyors in rust
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Convey, a conveyor library for rust.  For
/// licence information see the LICENSE file in the top level dirctory
/// of the distribution.
use crate::session::ConveySession;
use serde::de::DeserializeOwned;
use serde::ser::Serialize;
use shmem::Shmem;
use std::cell::{Cell, RefCell};
use std::env;
use std::rc::Rc;

pub mod collect;
pub mod session;
pub mod shmem_buffers;

/// A Convey instance sets up communication and allocates buffers
pub struct Convey {
    num_ranks: usize,
    my_rank: usize,
    /// the maximum number of u64s to send
    pub send_max: usize,
    /// the minimum number of u64s to send efficiently
    pub send_min_eff: usize,
    shmem: Shmem,
    debug: Cell<bool>,
}

impl std::fmt::Debug for Convey {
    //  TODO: need to figure out how to debug buffers used, valid, etc
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ConveySession")
            .field("num_ranks", &self.num_ranks)
            .field("my_rank", &self.my_rank)
            .field("send_max", &self.send_max)
            .field("send_min_eff", &self.send_min_eff)
            //.field("shmem", &self._shmem)
            .finish()
    }
}

impl Convey {
    /// Initialize the conveyor.
    pub fn new() -> Result<Self, &'static str> {
        let shmem = Shmem::new().expect("shmem init failed");
        let num_ranks = shmem.n_pes();
        let my_rank = shmem.my_pe();
        let send_max = 8192;
        let send_min_eff = 1024;
        let debug = Cell::new(env::var("CONVEY_DEBUG").is_ok());

        Ok(Convey {
            num_ranks,
            my_rank,
            send_max,
            send_min_eff,
            shmem,
            debug,
        })
    }

    /// Return my rank in 0..num_ranks-1
    pub fn my_rank(&self) -> usize {
        self.my_rank
    }

    /// total number of ranks
    pub fn num_ranks(&self) -> usize {
        self.num_ranks
    }

    /// set or clear the debug flag
    pub fn debug(&self, val: bool) {
        self.debug.set(val);
    }

    /// read the debug flag
    pub fn is_debug(&self) -> bool {
        self.debug.get()
    }

    /// Convenience function: size on my rank
    pub fn per_my_rank(&self, n: usize) -> usize {
        (n + self.num_ranks - self.my_rank - 1) / self.num_ranks
    }

    /// Convenience function: size on my rank
    pub fn offset_rank(&self, n: usize) -> (usize, usize) {
        (n / self.num_ranks, n % self.num_ranks)
    }

    /// Begin a session.
    pub fn begin<'a, T: Copy + Serialize + DeserializeOwned>(
        pull_fn: impl FnMut(T, usize) + 'a,
    ) -> ConveySession<'a, T> {
        let session = ConveySession::new_with_pull(pull_fn);
        session
    }

    /// create a session without a pull_fn
    pub fn session<'a, T: Copy + Serialize + DeserializeOwned>() -> ConveySession<'a, T> {
        let session = ConveySession::new();
        session
    }

    /// execute a simple function
    pub fn simple<T, I>(map: I, pull_fn: impl FnMut(T, usize)) -> ()
    where
        T: Copy + Serialize + DeserializeOwned,
        I: Iterator<Item = (T, usize)>,
    {
        Convey::session().pull_fn(pull_fn).push_iter(map).finish();
    }

    /// execute a simple function with return
    pub fn simple_return<T1, T2, I>(
        map: I,
        mut pull_fn: impl FnMut(T1) -> T2,
        mut ret_fn: impl FnMut(usize, T2),
    ) -> ()
    where
        T1: Copy + Serialize + DeserializeOwned,
        T2: Copy + Serialize + DeserializeOwned,
        I: Iterator<Item = ((usize, T1), usize)>,
    {
        let mut session2 = Convey::begin(|item: (usize, T2), _from_rank: usize| {
            ret_fn(item.0, item.1);
        });
        session2.simple_done = false;
        {
            Convey::session()
                .pull_fn(|item: (usize, T1), from_rank: usize| {
                    let ret = pull_fn(item.1);
                    session2.push((item.0, ret), from_rank);
                })
                .push_iter(map)
                .finish();
        }
        session2.finish();
    }

    /// call to shmem barrier
    ///  TODO: Maybe rethink this, only real purpose is to do io/calls
    pub fn barrier(&self) {
        self.shmem.barrier();
    }
}

/// An RcVec is a vector which is wrapped in Rc<RefCell<>> This is
///  very useful in conveyors where information needs to be
///  transferred between a push loop and a pull loop
#[derive(Debug)]
pub struct RcVec<T> {
    vec: Rc<RefCell<Vec<T>>>,
}

impl<T> RcVec<T>
where
    T: Copy,
{
    /// Create a new RcVec
    pub fn new() -> Self {
        RcVec {
            vec: Rc::new(RefCell::new(Vec::new())),
        }
    }
    /// Push a value
    pub fn push(&self, val: T) {
        self.vec.borrow_mut().push(val);
    }
    /// Pop a value
    pub fn pop(&self) -> Option<T> {
        self.vec.borrow_mut().pop()
    }
    /// Get it's length
    pub fn len(&self) -> usize {
        self.vec.borrow().len()
    }
    /// get a value
    pub fn get_at(&self, index: usize) -> T {
        self.vec.borrow()[index]
    }
    /// set a value
    pub fn set_at(&self, index: usize, val: T) {
        self.vec.borrow_mut()[index] = val;
    }
}

/// A module which provides for testing convey and uses of convey.
/// This is necessary because our underlying shmem implementation does
/// not support multithreading, while the testing infrastructure uses
/// multithreading.  Therefore, we must cause only one test which
/// calls convey::new() to exist at a time.
pub mod testing_support;
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn it_works() {
        let mutex = super::testing_support::TestingMutex::new();
        let mut local: Vec<i64> = vec![0; 100];
        let index: Vec<usize> = vec![3, 15, 5, 33];
        {
            // This needs to be in a new block so that the session will be dropped
            //  before we access the outer data after the conveyor

            // local is moved into the session because it is accessed in the calback
            //   * Should be able to pass in a session which is the response
            //      converyor (for ig, etc), but still figuring this out
            let mut session = Convey::begin(|item: usize, _from_rank| {
                local[item] += 1;
            });

            for item in index {
                let rank = item % mutex.convey.num_ranks;
                let lindex = item / mutex.convey.num_ranks;
                session.push(lindex, rank);
            }

            session.finish();
        }
        assert_eq!(local[0], 0);
        assert_eq!(local[3], 1);
        assert_eq!(local[15], 1);
        assert_eq!(local[5], 1);
        assert_eq!(local[33], 1);
    }

    #[test]
    fn test_simple() {
        let mutex = super::testing_support::TestingMutex::new();
        // Our variables
        let mut local: Vec<i64> = vec![0; 100];
        let index: Vec<usize> = vec![3, 15, 5, 33];

        Convey::simple(
            index.iter().map(|x| mutex.convey.offset_rank(*x)),
            |item: usize, _from_rank| {
                local[item] += 1;
            },
        );

        assert_eq!(local[0], 0);
        assert_eq!(local[3], 1);
        assert_eq!(local[15], 1);
        assert_eq!(local[5], 1);
        assert_eq!(local[33], 1);
    }
}
