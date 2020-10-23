//! module which implements shmem operations not related to data
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Convey, a conveyor library for rust.  For
/// licence information see the LICENSE file in the top level dirctory
/// of the distribution.
use crate::error::Error::NewAfterDrop;
use crate::object::{GlobalObject, Object};
use crate::Shmem;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;

/// To improve efficiency, we pre-allocate some work arrays on a per-thread basis
///  TODO: reimplement using Rc approach, eliminates need for lifetime
#[derive(Debug)]
pub struct ShmemThreadLocal {
    /// sync work array passed to collectives
    pub sync_obj: Object<i64>,
    /// the size of the work object TODO: Remove this, not needed
    pub work_size: usize,
    /// the work array passed to shmem collectives
    pub work_i64: Object<i64>,
}

impl ShmemThreadLocal {
    /// Initialize the thread local storage
    pub fn new() -> Self {
        let sync_size = shmem_sys::SHMEM_SYNC_SIZE as usize;
        let work_size = shmem_sys::SHMEM_REDUCE_MIN_WRKDATA_SIZE as usize;
        ShmemThreadLocal {
            sync_obj: Object::<i64>::new(sync_size).unwrap(),
            work_size: work_size,
            work_i64: Object::<i64>::new(work_size).unwrap(),
        }
    }
}

use std::cell::RefCell;
thread_local! {
    /// Instantiate the thread local data
    pub static SHMEM_WORK: RefCell<ShmemThreadLocal> = RefCell::new(ShmemThreadLocal::new());
}

static SHMEM_CREATED: AtomicBool = AtomicBool::new(false);
use lazy_static::lazy_static;
// count of crate::new() so we only finalize after last drop
lazy_static! {
    static ref SHMEM_COUNTER: Mutex<i32> = Mutex::new(0);
}

impl Shmem {
    /// Create a new instance of the shmem object
    ///
    ///  Key thing here is that we can call shmem_init() only once but may have
    ///  multiple customers of Shmem in a program (even across threads), so we need to have
    ///  a lock-protected code and a counter of calls so we only call shmem_finalize() on the
    ///  last drop().
    ///
    ///  This also leads to an error if you call new(), drop() and then new() again as we will
    ///  have already called shmem_finalize() and we can't call shmem_init() again.  Sadness.
    pub fn new() -> crate::Result<Self> {
        let mut num = SHMEM_COUNTER.lock().unwrap();
        let first =
            SHMEM_CREATED.compare_exchange(false, true, Ordering::Relaxed, Ordering::Relaxed);
        if first.is_ok() {
            // SAFETY: Call to shmem library with no arguments
            unsafe {
                shmem_sys::shmem_init();
            }
        } else {
            if *num == 0 {
                return Err(NewAfterDrop);
            }
        }
        *num += 1;
        Ok(Shmem {})
    }
    /// the number of active Shmem instances
    pub fn active_counter() -> usize {
        *SHMEM_COUNTER.lock().unwrap() as usize
    }
    /// wraper for my_pe()
    pub fn my_pe(&self) -> usize {
        Self::my_pe_raw()
    }
    /// wraper for my_pe(), callable without reference to an instance
    pub fn my_pe_raw() -> usize {
        // SAFETY: Call to shmem library with no arguments
        unsafe { shmem_sys::shmem_my_pe() as usize }
    }
    /// wraper for n_pes()
    pub fn n_pes(&self) -> usize {
        Self::n_pes_raw()
    }
    /// wraper for n_pes(), callable without reference to an instance
    pub fn n_pes_raw() -> usize {
        // SAFETY: Call to shmem library with no arguments
        unsafe { shmem_sys::shmem_n_pes() as usize }
    }
    /// create a new object
    pub fn new_object<T: Copy>(&self, elements: usize) -> crate::Result<Object<T>> {
        Object::new(elements)
    }
    /// create a new globalobject
    pub fn new_global_object<T: Copy>(
        &self,
        elements: usize,
        stride: usize,
    ) -> crate::Result<GlobalObject<T>> {
        GlobalObject::new(elements, stride)
    }
    /// wrapper for shmem_barrier_all()
    pub fn barrier(&self) {
        // SAFETY: Call to shmem library with no arguments
        unsafe { shmem_sys::shmem_barrier_all() };
    }
    /// wrapper for shmem_fence()
    pub fn fence(&self) {
        // SAFETY: Call to shmem library with no arguments
        unsafe { shmem_sys::shmem_fence() };
    }
    /// wrapper for shmem_quiet()
    pub fn quiet(&self) {
        // SAFETY: Call to shmem library with no arguments
        unsafe { shmem_sys::shmem_quiet() };
    }
}

impl Drop for Shmem {
    fn drop(&mut self) {
        // Only call finalize on last drop()
        let mut num = SHMEM_COUNTER.lock().unwrap();
        *num -= 1;
        if *num == 0 {
            // SAFETY: Call to shmem library with no arguments
            unsafe {
                shmem_sys::shmem_finalize();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atomic::Atomic;
    // Need to use only one test because tests are done in parallel threads
    // and we can only call shmem_init() and shmem_finalize() once.  This
    // implementation handles the case of multiple calls to crate::new() but
    // the problem occurs when one set of tests finishes completely before
    // others start.
    // we do test this error is raised by a final allocation after the
    // others are raised
    #[test]
    fn one_and_only_shmem_test() {
        {
            let shmem1 = Shmem::new().expect("initialization failed 1");
            {
                let shmem2 = Shmem::new().expect("initialization failed 2");
                // first make sure we can call new more than once
                assert_eq!(shmem1.my_pe(), 0);
                assert_eq!(shmem2.my_pe(), 0);
            }
            let a = shmem1.new_object::<i64>(10).expect("allocation error");
            let lp = a.local_part();
            lp[2] = 40;
            a.atomic_add(2, 2, 0).expect("atomic add error");
            shmem1.barrier();
            assert_eq!(lp[2], 42);
        }
        // final test, make sure we get new after drop error
        let shmem3 = Shmem::new();
        match shmem3 {
            Ok(_) => panic!("expected NewAfterDrop error"),
            Err(NewAfterDrop) => (),
            _ => panic!("Wrong error, should be NewAfterDrop"),
        }
    }
}
