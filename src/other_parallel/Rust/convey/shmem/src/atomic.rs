//! Module to handle shmem atomics, this is largely a bunch of macros to
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
use crate::error::Error::InvalidPE;
use crate::object::{GlobalObject, Object};
use crate::Shmem;
use std::os::raw::c_int;

/// This trait presents atomic operations of shmem to Objects
///   TODO: not all are implemented yet
pub trait Atomic {
    /// Our associated type, will be filled in by impl
    type T;
    /// Fetch a remote element, return value
    fn atomic_fetch(&self, offset: usize, pe: usize) -> crate::Result<Self::T>;
    /// Set a remote element to a value
    fn atomic_set(&self, offset: usize, value: Self::T, pe: usize) -> crate::Result<()>;
    /// Add a value to a remote element
    fn atomic_add(&self, offset: usize, value: Self::T, pe: usize) -> crate::Result<()>;
    /// Add a value to a remote element and return old value
    fn atomic_fetch_add(&self, offset: usize, value: Self::T, pe: usize) -> crate::Result<Self::T>;
    /// Increment remote element
    fn atomic_inc(&self, offset: usize, pe: usize) -> crate::Result<()>;
    /// Increment remote element and return old value
    fn atomic_fetch_inc(&self, offset: usize, pe: usize) -> crate::Result<Self::T>;
}

// Convenience macros to create definitions for a given type
macro_rules! at_declv {
    ($fn: ident, $shfn: ident, $ret: ty) => {
        fn $fn(&self, offset: usize, value: Self::T, pe: usize) -> crate::Result<$ret> {
            if pe >= Shmem::n_pes_raw() {
                return Err(InvalidPE);
            }
            let dest = self.void_ptr_with_offset(offset)?;
            //SAFETY: call to shmem with checked arguments
            unsafe { Ok(shmem_sys::$shfn(dest as *mut Self::T, value, pe as c_int)) }
        }
    };
}
macro_rules! at_decl {
    ($fn: ident, $shfn: ident, $ret: ty) => {
        fn $fn(&self, offset: usize, pe: usize) -> crate::Result<$ret> {
            if pe >= Shmem::n_pes_raw() {
                return Err(InvalidPE);
            }
            let dest = self.void_ptr_with_offset(offset)?;
            //SAFETY: call to shmem with checked arguments
            unsafe { Ok(shmem_sys::$shfn(dest as *mut Self::T, pe as c_int) as $ret) }
        }
    };
}

impl Atomic for Object<i64> {
    type T = i64;
    at_decl!(atomic_fetch, shmem_int64_atomic_fetch, i64);
    at_declv!(atomic_set, shmem_int64_atomic_set, ());
    at_declv!(atomic_add, shmem_int64_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_int64_atomic_fetch_add, i64);
    at_decl!(atomic_inc, shmem_int64_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_int64_atomic_fetch_inc, i64);
}

impl Atomic for Object<u64> {
    type T = u64;
    at_decl!(atomic_fetch, shmem_uint64_atomic_fetch, u64);
    at_declv!(atomic_set, shmem_uint64_atomic_set, ());
    at_declv!(atomic_add, shmem_uint64_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_uint64_atomic_fetch_add, u64);
    at_decl!(atomic_inc, shmem_uint64_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_uint64_atomic_fetch_inc, u64);
}

impl Atomic for Object<i32> {
    type T = i32;
    at_decl!(atomic_fetch, shmem_int32_atomic_fetch, i32);
    at_declv!(atomic_set, shmem_int32_atomic_set, ());
    at_declv!(atomic_add, shmem_int32_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_int32_atomic_fetch_add, i32);
    at_decl!(atomic_inc, shmem_int32_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_int32_atomic_fetch_inc, i32);
}

impl Atomic for Object<u32> {
    type T = u32;
    at_decl!(atomic_fetch, shmem_uint32_atomic_fetch, u32);
    at_declv!(atomic_set, shmem_uint32_atomic_set, ());
    at_declv!(atomic_add, shmem_uint32_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_uint32_atomic_fetch_add, u32);
    at_decl!(atomic_inc, shmem_uint32_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_uint32_atomic_fetch_inc, u32);
}

impl Atomic for Object<isize> {
    type T = isize;
    at_decl!(atomic_fetch, shmem_ptrdiff_atomic_fetch, isize);
    at_declv!(atomic_set, shmem_ptrdiff_atomic_set, ());
    at_declv!(atomic_add, shmem_ptrdiff_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_ptrdiff_atomic_fetch_add, isize);
    at_decl!(atomic_inc, shmem_ptrdiff_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_ptrdiff_atomic_fetch_inc, isize);
}

//bindgen is currently making usize be u64. This will unfortuntately
//require a cast to u64 for the result returns
impl Atomic for Object<usize> {
    type T = u64;
    at_decl!(atomic_fetch, shmem_size_atomic_fetch, u64);
    at_declv!(atomic_set, shmem_size_atomic_set, ());
    at_declv!(atomic_add, shmem_size_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_size_atomic_fetch_add, u64);
    at_decl!(atomic_inc, shmem_size_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_size_atomic_fetch_inc, u64);
}

/// This trait presents atomic operations of shmem to GlobalObjects
///   TODO: not all are implemented yet
pub trait GlobalAtomic {
    /// Our associated type, will be filled in by impl
    type T;
    /// Fetch a remote element, return value
    fn atomic_fetch(&self, offset: usize) -> crate::Result<Self::T>;
    /// Set a remote element to a value
    fn atomic_set(&self, offset: usize, value: Self::T) -> crate::Result<()>;
    /// Add a value to a remote element
    fn atomic_add(&self, offset: usize, value: Self::T) -> crate::Result<()>;
    /// Add a value to a remote element and return old value
    fn atomic_fetch_add(&self, offset: usize, value: Self::T) -> crate::Result<Self::T>;
    /// Increment remote element
    fn atomic_inc(&self, offset: usize) -> crate::Result<()>;
    /// Increment remote element and return old value
    fn atomic_fetch_inc(&self, offset: usize) -> crate::Result<Self::T>;
}

// Convenience macros to create definitions for a given type
macro_rules! at_declv {
    ($fn: ident, $shfn: ident, $ret: ty) => {
        fn $fn(&self, offset: usize, value: Self::T) -> crate::Result<$ret> {
            let (pe, loffset) = self.pe_loffset(offset);
            let dest = self.void_ptr_with_offset(loffset)?;
            //SAFETY: call to shmem with checked arguments
            unsafe { Ok(shmem_sys::$shfn(dest as *mut Self::T, value, pe as c_int)) }
        }
    };
}
macro_rules! at_decl {
    ($fn: ident, $shfn: ident, $ret: ty) => {
        fn $fn(&self, offset: usize) -> crate::Result<$ret> {
            let (pe, loffset) = self.pe_loffset(offset);
            let dest = self.void_ptr_with_offset(loffset)?;
            //SAFETY: call to shmem with checked arguments
            unsafe { Ok(shmem_sys::$shfn(dest as *mut Self::T, pe as c_int)) }
        }
    };
}

impl GlobalAtomic for GlobalObject<i64> {
    type T = i64;
    at_decl!(atomic_fetch, shmem_int64_atomic_fetch, i64);
    at_declv!(atomic_set, shmem_int64_atomic_set, ());
    at_declv!(atomic_add, shmem_int64_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_int64_atomic_fetch_add, i64);
    at_decl!(atomic_inc, shmem_int64_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_int64_atomic_fetch_inc, i64);
}

impl GlobalAtomic for GlobalObject<u64> {
    type T = u64;
    at_decl!(atomic_fetch, shmem_uint64_atomic_fetch, u64);
    at_declv!(atomic_set, shmem_uint64_atomic_set, ());
    at_declv!(atomic_add, shmem_uint64_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_uint64_atomic_fetch_add, u64);
    at_decl!(atomic_inc, shmem_uint64_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_uint64_atomic_fetch_inc, u64);
}

impl GlobalAtomic for GlobalObject<i32> {
    type T = i32;
    at_decl!(atomic_fetch, shmem_int32_atomic_fetch, i32);
    at_declv!(atomic_set, shmem_int32_atomic_set, ());
    at_declv!(atomic_add, shmem_int32_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_int32_atomic_fetch_add, i32);
    at_decl!(atomic_inc, shmem_int32_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_int32_atomic_fetch_inc, i32);
}

impl GlobalAtomic for GlobalObject<u32> {
    type T = u32;
    at_decl!(atomic_fetch, shmem_uint32_atomic_fetch, u32);
    at_declv!(atomic_set, shmem_uint32_atomic_set, ());
    at_declv!(atomic_add, shmem_uint32_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_uint32_atomic_fetch_add, u32);
    at_decl!(atomic_inc, shmem_uint32_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_uint32_atomic_fetch_inc, u32);
}

impl GlobalAtomic for GlobalObject<isize> {
    type T = isize;
    at_decl!(atomic_fetch, shmem_ptrdiff_atomic_fetch, isize);
    at_declv!(atomic_set, shmem_ptrdiff_atomic_set, ());
    at_declv!(atomic_add, shmem_ptrdiff_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_ptrdiff_atomic_fetch_add, isize);
    at_decl!(atomic_inc, shmem_ptrdiff_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_ptrdiff_atomic_fetch_inc, isize);
}

impl GlobalAtomic for GlobalObject<usize> {
    type T = u64;
    at_decl!(atomic_fetch, shmem_size_atomic_fetch, u64);
    at_declv!(atomic_set, shmem_size_atomic_set, ());
    at_declv!(atomic_add, shmem_size_atomic_add, ());
    at_declv!(atomic_fetch_add, shmem_size_atomic_fetch_add, u64);
    at_decl!(atomic_inc, shmem_size_atomic_inc, ());
    at_decl!(atomic_fetch_inc, shmem_size_atomic_fetch_inc, u64);
}
