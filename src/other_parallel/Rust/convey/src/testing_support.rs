use crate::Convey;
use lazy_static::lazy_static;
use std::sync::{Mutex, MutexGuard};

lazy_static! {
    static ref LOCK: Mutex<i32> = Mutex::new(0);
}

/// A strut to hold our open mutex and open convey instance
#[derive(Debug)]
pub struct TestingMutex<'a> {
    /// Tests can use this convey instance
    pub convey: Convey,
    // Test should *not* use this, so keep it private
    _data: MutexGuard<'a, i32>,
}

impl<'a> TestingMutex<'a> {
    /// Create a new TestingMutex instance
    pub fn new() -> TestingMutex<'a> {
        // it is important to get the convey structure first due to
        // requirement that we always have one shmem instance open.
        // There could actually be a race condition here, so maybe we
        // need to have a slight delay or synchronization somewhere
        let convey = Convey::new().unwrap();
        let data = LOCK.lock().unwrap();
        TestingMutex {
            convey: convey,
            _data: data,
        }
    }
}
