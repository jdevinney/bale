//! This package implements the rust interface to shmem-sys, which in turn
//!  is a wrapper for OpenShmem implmentations, currently version 1.4

#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

use shmem::Shmem;

pub mod collect;
pub mod error;
pub mod object;

use std::rc::Rc;

/// Generic result type for this library
pub type Result<T> = std::result::Result<T, crate::error::Error>;

/// Our main struct, keeps track of global stuff we need
#[derive(Debug)]
pub struct Pshmem {
    /// The index of the current rank
    pub rank: usize,
    /// The total number of ranks
    pub size: usize,
    shmem: Shmem,
}

/// To improve efficiency, we pre-allocate some work arrays on a per-thread basis
///  We use Rc block so no need to use lifetimes
#[derive(Debug)]
pub struct PshmemThreadLocal {
    /// i64 thread local work block
    work_i64: Rc<Object<i64>>,
    /// u64 thread local work block
    work_u64: Rc<Object<u64>>,
    /// i32 thread local work block
    work_i32: Rc<Object<i32>>,
    /// u32 thread local work block
    work_u32: Rc<Object<u32>>,
    /// isize thread local work block
    work_isize: Rc<Object<isize>>,
    /// isize thread local work block
    work_usize: Rc<Object<usize>>,
}

impl PshmemThreadLocal {
    /// initialize thread local storage
    pub fn new(work_size: usize) -> Self {
        PshmemThreadLocal {
            work_i64: Rc::new(Object::<i64>::new(work_size).unwrap()),
            work_u64: Rc::new(Object::<u64>::new(work_size).unwrap()),
            work_i32: Rc::new(Object::<i32>::new(work_size).unwrap()),
            work_u32: Rc::new(Object::<u32>::new(work_size).unwrap()),
            work_isize: Rc::new(Object::<isize>::new(work_size).unwrap()),
            work_usize: Rc::new(Object::<usize>::new(work_size).unwrap()),
        }
    }
}

use std::cell::RefCell;
thread_local! {
    /// The actual thread local instance
    pub static PSHMEM_WORK: RefCell<PshmemThreadLocal> = RefCell::new(PshmemThreadLocal::new(100));
}

use crate::object::{GlobalObject, Object};

// eventually condition this with #[cfg(shmem)], for example
impl Pshmem {
    /// Create a new Pshmem instance
    pub fn new() -> Self {
        let shmem = Shmem::new().expect("shmem init failed");
        let rank = shmem.my_pe();
        let size = shmem.n_pes();
        Pshmem { rank, size, shmem }
    }
    /// We call the index in Pshmem space to be it's rank
    pub fn rank(&self) -> usize {
        self.rank
    }
    /// Number of Pshmem spaces in our program
    pub fn size(&self) -> usize {
        self.size
    }
    /// Barrier: all PEs must enter before any leaves
    ///  TODO: Figure out what this means on a per thread basis
    pub fn barrier(&self) {
        self.shmem.barrier();
    }
    /// Create a new object
    pub fn new_object<T: Copy>(&self, elements: usize) -> crate::Result<Object<T>> {
        Object::new(elements)
    }
    /// Create a new global object
    pub fn new_global_object<T: Copy>(
        &self,
        elements: usize,
        blocking: usize,
    ) -> crate::Result<GlobalObject<T>> {
        GlobalObject::new(elements, blocking)
    }
}
