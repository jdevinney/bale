//! Module to handle shmem collectives, this is largely a bunch of macros to
//! instantiate all of the variants of atomics and types
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Convey, a conveyor library for rust.  For
/// licence information see the LICENSE file in the top level dirctory
/// of the distribution.
use crate::error::Error::BoundsExceeded;
use crate::object::Object;
use crate::shmem::SHMEM_WORK;
use crate::Shmem;

/// This trait presents collective operations of shmem to Objects
///   TODO: not all are implemented yet
///   TODO: only impl for i64 so far
pub trait Collect {
    /// Our associated type, will be filled in by impl
    type T;
    /// Compute the minumum over all PEs
    fn min_to_all(&self, offset: usize, source: &[Self::T]) -> crate::Result<()>;
    /// Compute the maximum over all PEs
    fn max_to_all(&self, offset: usize, source: &[Self::T]) -> crate::Result<()>;
    /// Compute the sum over all PEs
    fn sum_to_all(&self, offset: usize, source: &[Self::T]) -> crate::Result<()>;
}

// Convenience macros to create definitions for a given type
macro_rules! at_decl {
    ($fn: ident, $shfn: ident, $shtype: ident) => {
        fn $fn(&self, offset: usize, source: &[Self::T]) -> crate::Result<()> {
            if offset + source.len() > self.local_part().len() {
                return Err(BoundsExceeded);
            }
            let dest = self.void_ptr_with_offset(offset).unwrap();
            SHMEM_WORK.with(|sw| {
                if sw.borrow().work_size < (source.len() + 1) / 2 {
                    todo!(
                        "Need to implement collectives of size {} > {}",
                        sw.borrow().work_size,
                        (source.len() + 1) / 2
                    );
                }
                unsafe {
                    Ok(shmem_sys::$shfn(
                        dest as *mut $shtype,
                        source.as_ptr() as *const $shtype,
                        source.len() as shmem_sys::nreduce_t,
                        0,
                        1,
                        Shmem::n_pes_raw() as i32,
                        sw.borrow().work_i64.local_part().as_ptr() as *mut $shtype,
                        sw.borrow().sync_obj.local_part().as_ptr() as *mut i64,
                    ))
                }
            })
        }
    };
}

impl Collect for Object<i64> {
    type T = i64;
    at_decl!(min_to_all, shmem_longlong_min_to_all, i64);
    at_decl!(max_to_all, shmem_longlong_min_to_all, i64);
    at_decl!(sum_to_all, shmem_longlong_min_to_all, i64);
}

impl Collect for Object<u64> {
    type T = u64;
    at_decl!(min_to_all, shmem_longlong_min_to_all, i64);
    at_decl!(max_to_all, shmem_longlong_min_to_all, i64);
    at_decl!(sum_to_all, shmem_longlong_min_to_all, i64);
}

impl Collect for Object<usize> {
    type T = usize;
    at_decl!(min_to_all, shmem_longlong_min_to_all, i64);
    at_decl!(max_to_all, shmem_longlong_min_to_all, i64);
    at_decl!(sum_to_all, shmem_longlong_min_to_all, i64);
}

impl Collect for Object<i32> {
    type T = i32;
    at_decl!(min_to_all, shmem_int_min_to_all, i32);
    at_decl!(max_to_all, shmem_int_min_to_all, i32);
    at_decl!(sum_to_all, shmem_int_min_to_all, i32);
}

impl Collect for Object<u32> {
    type T = u32;
    at_decl!(min_to_all, shmem_int_min_to_all, i32);
    at_decl!(max_to_all, shmem_int_min_to_all, i32);
    at_decl!(sum_to_all, shmem_int_min_to_all, i32);
}
