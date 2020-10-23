//! Permutation Library, a part of the Sparse Matrix Library
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For license information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use convey_hpc::collect::CollectValues;
use convey_hpc::Convey;
use rand::distributions::{Distribution, Uniform};
use std::cell::{Ref, RefCell, RefMut};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::rc::Rc;

/// A permutation needs a convey and the permutation itself
pub struct Perm {
    convey: Convey,
    perm: Vec<usize>,
}

impl std::fmt::Debug for Perm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let my_rank = self.convey.my_rank();
        f.debug_struct("Perm")
            .field("rank", &my_rank)
            .field("len", &self.perm.len())
            .finish()
    }
}

impl Perm {
    /// new takes a vec to start with
    pub fn new(v: Vec<usize>) -> Self {
        let convey = Convey::new().unwrap();
        Perm {
            perm: v,
            convey: convey,
        }
    }

    /// Create an identity permutation
    pub fn identity(n: usize) -> Self {
        let convey = Convey::new().unwrap();
        let mut perm = vec![0; convey.per_my_rank(n)];
        let my_rank = convey.my_rank();
        let num_ranks = convey.num_ranks();

        for i in 0..perm.len() {
            perm[i] = i * num_ranks + my_rank;
        }
        Perm { perm, convey }
    }

    /// Get the length of this permutation
    pub fn len(&self) -> usize {
        self.perm.len()
    }

    /// get the entry from a permutation
    pub fn entry(&self, index: usize) -> usize {
        self.perm[index]
    }

    /// get a reference to the permutation itself
    pub fn perm(&mut self) -> &mut Vec<usize> {
        &mut self.perm
    }

    /// predicate: is this a real permutation
    pub fn is_perm(&self) -> bool {
        let mut flags: Vec<bool> = vec![false; self.perm.len()];
        //println!("is_perm{:?}", self.perm);
        Convey::simple(
            self.perm.iter().map(|val| self.convey.offset_rank(*val)),
            |offset: usize, _from_rank: usize| {
                flags[offset] = true;
            },
        );
        let flags_unset = flags
            .iter()
            .fold(0, |acc, x| if *x { acc } else { acc + 1 });
        let total_flags_unset = flags_unset.reduce_sum();
        if total_flags_unset > 0 {
            println!("not a perm, {} unset", total_flags_unset);
        }
        total_flags_unset == 0
    }

    /// Create a random permutation
    pub fn random(n: usize, seed: i64) -> Self {
        Perm::random_darts_resubmit_local(n, seed)
    }

    /// Create a permutation from a dart-throwing vector
    fn from_target(target: Vec<i64>, n: usize) -> Self {
        let convey = Convey::new().unwrap();
        let mut perm = vec![0; convey.per_my_rank(n)];
        let count = target
            .iter()
            .fold(0_u64, |acc, x| if *x != -1_i64 { acc + 1 } else { acc });
        let offset = count.reduce_prior_sum() as usize;
        let mut pos = 0;
        let num_ranks = convey.num_ranks();
        let mut rank = offset % num_ranks;
        {
            let mut session = Convey::begin(|item: i64, _from_rank| {
                perm[pos] = item as usize;
                pos += 1;
            });

            for item in &target {
                if *item >= 0 {
                    session.push(*item, rank);
                    rank += 1;
                    if rank >= num_ranks {
                        rank = 0;
                    }
                }
            }
            session.finish();
        }
        Perm::new(perm)
    }

    /// writes the first and last part of an perm to the specified file
    /// # Arguments
    /// * a       the array
    /// * maxdisp the number of entries that are written, 0 means everything,
    ///            otherwise write the first and last maxdisp/2 entries
    /// * filename the filename to written to
    pub fn dump(&self, maxdisp: usize, filename: &str) -> Result<(), std::io::Error> {
        let path = Path::new(&filename);
        let mut file = File::create(path)?;

        let safe_disp = if maxdisp <= self.perm.len() && maxdisp > 0 {
            maxdisp / 2
        } else {
            self.perm.len() / 2
        };

        for entry in &self.perm[0..safe_disp] {
            write!(file, "{}\n", entry)?;
        }
        for entry in &self.perm[self.perm.len() - safe_disp..self.perm.len()] {
            write!(file, "{}\n", entry)?;
        }
        Ok(())
    }

    /// version of random_darts which take darts which are invalid (already used)
    /// and resubmit them on the local thread to throw again
    pub fn random_darts_resubmit_local(n: usize, _seed: i64) -> Self {
        let identity = Perm::identity(n);
        let convey = identity.convey;
        let queue = Rc::new(RefCell::new(identity.perm));
        let m = n * 2;
        let mut target = vec![-1_i64; convey.per_my_rank(m)];

        // initialize the random number generator
        let mut rng = rand::thread_rng();
        let die = Uniform::from(0..m as usize);

        {
            let mut session = Convey::begin(|item: (usize, usize), _from_rank| {
                let idx = item.0;
                if target[idx] == -1_i64 {
                    target[idx] = item.1 as i64;
                } else {
                    queue.borrow_mut().push(item.1);
                }
            });
            session.simple_done = false;

            while session.advance(get_len(queue.borrow()) == 0) {
                let mut batch = get_batch(queue.borrow_mut(), 200);
                //println!("pushing {}/{} items", batch.len(), get_len(queue.borrow()));
                while let Some(item) = batch.pop() {
                    let (offset, rank) = convey.offset_rank(die.sample(&mut rng));
                    //println!("pushing i{} o{} r{} ", item, offset, rank);
                    session.push((offset, item), rank);
                }
            }
        }

        Perm::from_target(target, n)
    }

    /// version of random_darts which take darts which are invalid (already used)
    /// and return them to the sender, who in turn will resend them
    pub fn random_darts_return_rejects(n: usize, _seed: i64) -> Self {
        let identity = Perm::identity(n);
        let convey = identity.convey;
        let queue = Rc::new(RefCell::new(identity.perm));
        let m = n * 2;
        let mut target = vec![-1_i64; convey.per_my_rank(m)];

        // initialize the random number generator
        let mut rng = rand::thread_rng();
        let die = Uniform::from(0..m as usize);
        {
            // In a separate block so session2 life ends before target use needed
            let session2 = Rc::new(RefCell::new(Convey::begin(|item: usize, _from_rank| {
                queue.borrow_mut().push(item);
            })));
            session2.borrow_mut().simple_done = false;
            {
                // In a separate block so session1 life ends before target use needed
                let mut session1 = Convey::begin(|item: (usize, usize), from_rank| {
                    if target[item.0] == -1_i64 {
                        target[item.0] = item.1 as i64;
                    } else {
                        session2.borrow_mut().push(item.1, from_rank);
                    }
                });

                session1.simple_done = false;
                loop {
                    let s1continue = session1.advance(get_len(queue.borrow()) == 0);
                    if !s1continue {
                        if !session2.borrow_mut().advance(true) {
                            break;
                        }
                    }
                    let mut batch = get_batch(queue.borrow_mut(), 20);
                    //println!("pushing {}/{} items", batch.len(), get_len(queue.borrow()));
                    while let Some(item) = batch.pop() {
                        let (offset, rank) = convey.offset_rank(die.sample(&mut rng));
                        session1.push((offset, item), rank);
                    }
                }
            }
        }
        Perm::from_target(target, n)
    }
}

/// helper function: get a batch of work from a queue
///   These help us avoid borrow checker problems.
fn get_batch(mut queue: RefMut<'_, Vec<usize>>, limit: usize) -> Vec<usize> {
    let mut ret: Vec<usize> = Vec::new();
    for _i in 0..queue.len().min(limit) {
        ret.push(queue.pop().unwrap());
    }
    ret
}

/// helper function: get length from a queue
///   These help us avoid borrow checker problems.
fn get_len(queue: Ref<'_, Vec<usize>>) -> usize {
    queue.len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use convey_hpc::testing_support::TestingMutex;
    #[test]
    fn perm() {
        let mutex = TestingMutex::new();
        let perm = Perm::identity(100);
        println!("identity perm: {:?}", perm.perm);
        assert_eq!(perm.perm.len(), mutex.convey.per_my_rank(100));
        drop(mutex);
    }
    #[test]
    fn is_perm() {
        let mutex = TestingMutex::new();
        let perm = Perm::identity(100);
        println!("identity perm: {:?}", perm.perm);
        assert_eq!(perm.is_perm(), true);
        drop(mutex);
    }
    #[test]
    fn is_perm_fail() {
        let mutex = TestingMutex::new();
        let mut perm = Perm::identity(100);
        perm.perm[0] += 1;
        println!("identity perm: {:?}", perm.perm);
        assert_eq!(perm.is_perm(), false);
        drop(mutex);
    }
    #[test]
    fn random_darts_resubmit_local() {
        let mutex = TestingMutex::new();
        let perm = Perm::random_darts_resubmit_local(2000, 0);
        //println!("random perm(resubmit): {:?}", perm.perm);
        assert_eq!(perm.is_perm(), true);
        assert_eq!(perm.perm.len(), mutex.convey.per_my_rank(2000));
        drop(mutex);
    }
    #[test]
    fn random_darts_return_rejects() {
        let mutex = TestingMutex::new();
        let perm = Perm::random_darts_return_rejects(100, 0);
        println!("random perm(rejects): {:?}", perm.perm);
        assert_eq!(perm.is_perm(), true);
        assert_eq!(perm.perm.len(), mutex.convey.per_my_rank(100));
        drop(mutex);
    }
}
