//! Conveyor module for handling a session
//!
//! Copyright (c) 2020, Institute for Defense Analyses
//! 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//!
//! All rights reserved.
//!
//! This file is part of Convey, a conveyor library for rust.  For
//! licence information see the LICENSE file in the top level dirctory
//! of the distribution.
//!
//! Conveyor module for handling a session
//! A session has a fixed data type and contains the buffering needed
//! for both sorting and recieving
//!
//! This version does shmem communication directly, bypassing most of Convey,
//! which results in less risk of deadlock, but uses more buffering

use crate::shmem_buffers::{ShmemBufferGuard, ShmemBuffers};
use crate::Convey;
use bincode::{deserialize, serialize};
use serde::de::DeserializeOwned;
use serde::ser::Serialize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

/// A ConveySession lasts while the type we are sending is fixed
///     Initiated with Convey::begin()
///     Ended when dropped
///     n.b. this needs to be a separate struct from Convey so we can set its type
pub struct ConveySession<'a, T> {
    /// so others can access
    convey: Convey,
    uid: usize,
    pull_fn: Option<Box<dyn FnMut(T, usize) + 'a>>,
    /// done means done, can't push after this
    pub simple_done: bool,
    /// always send partial work in advance()
    pub always_send: bool,
    done_push: bool,
    current_phase: u32,
    goal_phase: u32,
    phase_sent: Vec<u32>,
    current_phase_received: Vec<u32>,
    /// type level sorting queues, one per destination
    /// Ordering is fine as we push to back and iterate off front of vector
    to_push: Vec<Vec<T>>,
    send_max: usize,
    send_min_eff: usize,
    shmem: ShmemBufferGuard,
    last_report: Instant,
    overflows: u64,
}

impl<T> std::fmt::Debug for ConveySession<'_, T> {
    //  TODO: need to figure out how to debug both the pull_fn and valid/avail receive
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut push_len: Vec<usize> = Vec::new();
        for item in &self.to_push {
            push_len.push(item.len());
        }
        f.debug_struct("ConveySession")
            .field("uid", &self.uid)
            //.field("convey", &self.convey) // Only interesting thing in rank
            .field("rank", &self.convey.my_rank)
            .field("done_push", &self.done_push)
            .field("simple_done", &self.simple_done)
            .field("always_send", &self.always_send)
            .field("current_phase", &self.current_phase)
            .field("goal_phase", &self.goal_phase)
            .field("phase_sent", &self.phase_sent)
            .field("current_phase_received", &self.current_phase_received)
            //.field("shmem", &self.shmem)
            .field("to_push.len()", &push_len)
            .field("overflows", &self.overflows)
            .field("send_max", &self.send_max)
            .field("send_min_eff", &self.send_min_eff)
            .finish()
    }
}

impl<'a, T> ConveySession<'a, T>
where
    T: Copy + Serialize + DeserializeOwned,
{
    /// A new session
    ///
    /// Lifetime parameter assures convey struct lasts as long as we do
    pub fn new() -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(1);
        let mut to_push: Vec<Vec<T>> = Vec::new();
        let convey = Convey::new().expect("couldn't get convey instance");
        for _i in 0..convey.num_ranks {
            to_push.push(Vec::new());
        }
        let len = std::mem::size_of::<T>();
        let xfer_len = len + 8;
        let shmem = ShmemBuffers::new(&convey);
        let send_max = (convey.send_max - xfer_len) / len;
        let send_min_eff = (convey.send_min_eff - xfer_len) / len;
        ConveySession {
            pull_fn: None,
            uid: COUNTER.fetch_add(1, Ordering::Relaxed),
            simple_done: true,
            always_send: false,
            done_push: false,
            current_phase: 0,
            goal_phase: 2,
            phase_sent: vec![0; convey.num_ranks],
            current_phase_received: vec![0; convey.num_ranks],
            to_push: to_push,
            convey: convey,
            shmem: shmem,
            send_max: send_max,
            send_min_eff: send_min_eff,
            last_report: Instant::now(),
            overflows: 0,
        }
    }

    ///  New that provides a pull function within it
    pub fn new_with_pull(pull_fn: impl FnMut(T, usize) + 'a) -> Self {
        let mut ret = ConveySession::new();
        ret.pull_fn = Some(Box::new(pull_fn));
        ret
    }

    /// pass through debug function
    pub fn debug(&self, val: bool) {
        self.convey.debug(val);
    }

    /// Builder function to set always send, default is to defer small sends
    pub fn always_send(&mut self) -> &mut Self {
        self.always_send = true;
        self
    }

    /// Builder function to defer small sends, the default
    pub fn defer_small_sends(&mut self) -> &mut Self {
        self.always_send = false;
        self
    }

    /// Builder function to set simple_done flag, the default
    pub fn simple_done(&mut self) -> &mut Self {
        self.simple_done = true;
        self
    }

    /// Builder function to allow done to be withdrawn, default is simple_done
    pub fn allow_done_withdrawl(&mut self) -> &mut Self {
        self.simple_done = false;
        self
    }

    /// Builder function to add a pull function
    pub fn pull_fn(&mut self, pull_fn: impl FnMut(T, usize) + 'a) -> &mut Self {
        self.pull_fn = Some(Box::new(pull_fn));
        self
    }

    /// Builder function to do push iterators
    pub fn push_iter<I>(&mut self, map: I) -> &mut Self
    where
        I: Iterator<Item = (T, usize)>,
    {
        assert!(self.pull_fn.is_some());
        for item in map {
            self.push(item.0, item.1);
        }
        self
    }

    /// Builder termination function to terminate session
    pub fn finish(&mut self) -> () {
        self.last_report = Instant::now();
        //println!("enter finish {:?}", self);
        loop {
            if !self.advance(true) {
                break;
            }
            let d = self.last_report.elapsed();
            // deadlock? report after 5 seconds in finish
            if d.as_millis() > 5000 {
                println!("in finish {:?}", self);
                self.last_report = Instant::now();
            }
        }
        //println!("exit finish rank {:?}", self);
    }

    /// Push an item to a destination
    pub fn push(&mut self, item: T, dest: usize) -> () {
        /*
        println!(
            "uid: {} bytes on entry to push {:?}",
            self.uid,
            &self.receive_buf.local_part()[0..40],
        );
        */
        // If we push() and were done before, we need to seek phase + 2 now
        if self.done_push {
            if self.simple_done {
                panic!("Call to ConveySession::push() after done was indicated, configure to allow withdrawl of done");
            }
            self.done_push = false;
            self.goal_phase = self.current_phase + 2;
        }
        self.to_push[dest].push(item);
        //println!("push {} len {}", dest, self.to_push[dest].len());
        if self.to_push[dest].len() >= self.send_max {
            loop {
                self.advance(false);
                // don't let ourselves get rediculusly out of balance
                if self.to_push[dest].len() < self.send_max * 4 {
                    break;
                }
            }
        }
    }

    /// Advance the conveyor, return false when complete
    pub fn advance(&mut self, done: bool) -> bool {
        if self.convey.is_debug() {
            println!("enter advance {} {:?}", done, self);
        }
        // If we are done and are not in done_push state, increment our phase and set done_push state
        if done && !self.done_push {
            self.done_push = done;
            if self.simple_done {
                // if we are done using the simple definition, we are really done, set current_phase to 2
                self.current_phase = 2;
            }
        }
        // If we are not done and were done before, we need to seek phase + 2 now
        if !done && self.done_push {
            if self.simple_done {
                panic!("Call to ConveySession::advance(done:false) after done was indicated, configure to allow withdrawl of done");
            }
            self.goal_phase = self.current_phase + 2;
            self.done_push = done;
        }
        self.maybe_push();
        self.try_pull();

        if !done {
            return true;
        }
        // First check if we can advance current_phase
        // TODO: Make this a fold()
        for phase in &self.phase_sent {
            if *phase != self.current_phase {
                return true;
            }
        }

        for phase in &self.current_phase_received {
            if *phase < self.current_phase {
                return true;
            }
        }

        for q in &self.to_push {
            if q.len() > 0 {
                return true;
            }
        }

        if self.simple_done {
            // If we have sent everyone our data and phase and we have
            // received from everyone their phase in a "done" state, we are
            // done, exit now
            if self.convey.is_debug() {
                println!("exit advance(simple done){} {:?}", done, self);
            }
            return false;
        }
        // Once we have received and sent all of our current phase, we
        // are ready to advance, but only as far as our goal
        if self.current_phase < self.goal_phase {
            //println!("Advance current_phase {:?}", self);
            self.current_phase += 1;
        }

        for phase in &self.phase_sent {
            if *phase != self.goal_phase {
                return true;
            }
        }
        for phase in &self.current_phase_received {
            if *phase != self.goal_phase {
                return true;
            }
        }
        if self.convey.is_debug() {
            println!("exit advance {} {:?}", done, self);
        }
        return false;
    }

    /// Decide if we want to do communication
    ///  TODO: better decision about when to push
    pub fn maybe_push(&mut self) -> () {
        let my_rank = self.convey.my_rank;
        let done = self.done_push;
        let mut to_set_valid: Vec<usize> = Vec::new();
        //println!("enter maybe_push {:?}", self);

        // first check to see if we need to increase our goal
        if !self.simple_done {
            for to_rank in 0..self.convey.num_ranks {
                self.goal_phase = self
                    .goal_phase
                    .max(self.shmem.get_their_goal_on_my_rank(my_rank, to_rank));
            }
        }

        for to_rank in 0..self.convey.num_ranks {
            // see if they need to increase their goal
            if !self.simple_done
                && (self.shmem.get_their_goal_on_my_rank(my_rank, to_rank) < self.goal_phase)
            {
                self.shmem
                    .set_my_goal_on_rank(to_rank, my_rank, self.goal_phase);
            }
            let to_push = &self.to_push[to_rank];
            if ((done && (to_push.len() > 0 || self.phase_sent[to_rank] != self.current_phase))
                || to_push.len()
                    > if self.always_send {
                        0
                    } else {
                        self.send_min_eff
                    })
                && self.shmem.receive_is_available(to_rank, my_rank)
            {
                self.shmem.set_receive_available(to_rank, my_rank, false);

                let s = self.shmem.send_buf.local_part();
                let srange = if to_push.len() > self.send_max {
                    let v = serialize(&to_push[0..self.send_max]).unwrap();
                    assert!(v.len() <= self.convey.send_max);
                    let srange =
                        self.convey.send_max * to_rank..(self.convey.send_max * to_rank) + v.len();
                    s[srange.clone()].copy_from_slice(&v);
                    self.to_push[to_rank] = to_push[self.send_max..].to_vec();
                    self.overflows += 1;
                    srange
                } else {
                    let v = serialize(to_push).unwrap();
                    assert!(v.len() <= self.convey.send_max);
                    let srange =
                        self.convey.send_max * to_rank..(self.convey.send_max * to_rank) + v.len();
                    s[srange.clone()].copy_from_slice(&v);
                    self.phase_sent[to_rank] = self.current_phase;
                    self.to_push[to_rank] = Vec::new();
                    srange
                };

                if self.convey.is_debug() {
                    println!(
                        "send uid {} from {} to {} range {:?} offset {} bytes {:?} remaining_len {} ",
                        self.uid,
                        my_rank,
                        to_rank,
                        srange,
                        my_rank * self.convey.send_max,
                        &self.shmem.send_buf.local_part()[srange.start..srange.start + 20],
                        self.to_push[to_rank].len(),
                    );
                }

                self.shmem
                    .receive_buf
                    .put_with_offset_and_range(
                        my_rank * self.convey.send_max,
                        &self.shmem.send_buf,
                        srange,
                        to_rank,
                    )
                    .expect("put failure");
                to_set_valid.push(to_rank);
            }
        }
        self.convey.shmem.fence();

        for to_rank in &to_set_valid {
            let val = if done && self.current_phase == self.phase_sent[*to_rank] {
                -((self.current_phase + 1) as i64)
            } else {
                (self.current_phase + 1) as i64
            };
            self.shmem.set_receive_valid(*to_rank, my_rank, val);
        }
        self.convey.shmem.quiet();
    }

    /// Try to pull an item off the queue
    pub fn try_pull(&mut self) -> () {
        let to_rank = self.convey.my_rank;
        //println!("enter try_pull {:?}", self);
        for from_rank in 0..self.convey.num_ranks {
            let val = self.shmem.receive_is_valid(to_rank, from_rank);
            if val != 0 {
                let phase = if val < 0 {
                    if self.simple_done {
                        // if we are in simple_done and we have reached "done", we are done, skip phase to 2
                        2
                    } else {
                        ((-val) - 1) as u32
                    }
                } else {
                    (val - 1) as u32
                };
                self.current_phase_received[from_rank] = phase;
                self.shmem.set_receive_valid(to_rank, from_rank, 0);
                let range =
                    (from_rank * self.convey.send_max)..((from_rank + 1) * self.convey.send_max);

                if self.convey.is_debug() {
                    println!(
                        "receive uid {} from {} to {} len {} range {:?} bytes {:?}",
                        self.uid,
                        from_rank,
                        to_rank,
                        0,
                        range,
                        &self.shmem.receive_buf.local_part()[range.start..range.start + 20],
                    );
                }

                let v: Vec<T> = deserialize(&self.shmem.receive_buf.local_part()[range])
                    .expect("deserializtion failed?");
                self.shmem.set_receive_available(to_rank, from_rank, true);

                for item in v {
                    (self.pull_fn.as_deref_mut().unwrap())(item, from_rank)
                }
            }
        }
        self.convey.shmem.quiet();
    }

    /// Convenience functions, transfer to convey impl, size on my rank for size
    pub fn per_my_rank(&self, n: usize) -> usize {
        self.convey.per_my_rank(n)
    }

    /// Convenience function: transfer to convey impl, ofset and rank of an index
    pub fn offset_rank(&self, n: usize) -> (usize, usize) {
        self.convey.offset_rank(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_builder() {
        let mutex = crate::testing_support::TestingMutex::new();
        let mut local: Vec<i64> = vec![0; 100];
        let index: Vec<usize> = vec![3, 15, 5, 33];
        Convey::session()
            .pull_fn(|item: usize, _from_rank| {
                local[item] += 1;
            })
            .push_iter(index.iter().map(|x| mutex.convey.offset_rank(*x)))
            .finish();

        assert_eq!(local[0], 0);
        assert_eq!(local[3], 1);
        assert_eq!(local[15], 1);
        assert_eq!(local[5], 1);
        assert_eq!(local[33], 1);
    }
}
