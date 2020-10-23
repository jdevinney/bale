//! Module to implement objects and global objects
//!  This is the heart of the shmem presentation to rust.
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Convey, a conveyor library for rust.  For
/// licence information see the LICENSE file in the top level dirctory
/// of the distribution.
use crate::error::Error::{BoundsExceeded, InvalidPE, UnsupportedBlocking};
use crate::Shmem;
use std::ffi::c_void;
use std::ops::Range;
use std::os::raw::c_int;

/// A shared object, indexed by local offset and pe
///  never implement Clone or Copy for this as the underlying memory is allocated via shmalloc()
#[derive(Debug)]
pub struct Object<T> {
    /// The local part of the underlying shmem array as a rust slice
    local_ptr: *mut T,
    local_elements: usize,
    is_subset: bool,
}

impl<T> Object<T>
where
    T: Copy,
{
    /// get a new shmem memory block as specified and build a slice representation of it
    pub fn new(local_elements: usize) -> crate::Result<Self> {
        let bytes = local_elements * std::mem::size_of::<T>();
        // SAFETY:
        //  - size of allocation is computed above
        //  - shmalloc returns uninitialized memory, so we initialize it
        //  - local_ptr is raw, so it can only be used in unsafe, see
        //    other SAFETY notes in this module
        let ptr = unsafe {
            let ptr = shmem_sys::shmalloc(bytes as u64) as *mut u8;
            let sl = std::slice::from_raw_parts_mut(ptr, bytes);
            sl.iter_mut().for_each(|x| *x = 0);
            shmem_sys::shmem_barrier_all();
            ptr as *mut T
        };

        //println!("new:{:?}, elements:{}, len:{}", ptr, local_elements, bytes);
        Ok(Object {
            local_ptr: ptr,
            local_elements: local_elements,
            is_subset: false,
        })
    }

    /// Make an object out of a subset of an object's local elements, changing type
    pub fn subset<T1: Copy>(&self, range: Range<usize>) -> crate::Result<Object<T1>> {
        let from_len = std::mem::size_of::<T>();
        let to_len = std::mem::size_of::<T1>();
        if range.end > self.local_elements || range.start >= range.end {
            return Err(BoundsExceeded);
        }
        // SAFETY: this is copying a range of a previous safe allocation, and the range is checked
        let local_ptr = unsafe { self.local_ptr.add(range.start) as *mut T1 };
        let local_elements = ((range.end - range.start) * from_len) / to_len;
        /*
            println!(
                "subset old:{:?}, new:{:?}, oldelements:{}, new_elements:{}, from_len:{}, to_len:{}",
                self.local_ptr, local_ptr, self.local_elements, local_elements, from_len, to_len
            );
        */
        Ok(Object {
            local_ptr,
            local_elements,
            is_subset: true,
        })
    }

    /// Make a mutable pointer to the local part of the struture
    pub fn local_part(&self) -> &mut [T] {
        // SAFETY: creates a correctly sized slice
        // Potential UNSAFETY: The shmem model allows multiple mutators of memory
        //   borrow checking.  Argument might want to be &mut self so we enforce
        //   borrow checking within the rust portion of a program, even if shmem
        //   is inherently unsafe due to remote puts and gets going on behind
        //   the scenes
        unsafe { std::slice::from_raw_parts_mut(self.local_ptr, self.local_elements) }
    }

    /// private (to crate) helper function to safely compute mutable pointer to send to shmem routine
    pub(crate) fn void_ptr_with_offset(&self, offset: usize) -> crate::Result<*mut c_void> {
        if offset > self.local_elements {
            return Err(BoundsExceeded);
        }
        // SAFETY: bounds are checked, result can only be used in an unsafe block
        //         Look to use cases for SAFETY messages
        Ok(unsafe { self.local_ptr.add(offset) as *mut c_void })
    }

    /// shmem get operations, with variants for passing offset and ranges
    /// get a slice of memory from a remote pe
    pub fn get(&self, source: &Object<T>, pe: usize) -> crate::Result<()> {
        self.get_with_offset_and_range(0, source, 0..source.local_elements, pe)
    }
    /// get variant when a destination offset is specifed
    pub fn get_with_offset(
        &self,
        offset: usize,
        source: &Object<T>,
        pe: usize,
    ) -> crate::Result<()> {
        self.get_with_offset_and_range(offset, source, 0..source.local_elements, pe)
    }

    /// get variant when a destination offset and source range is specifed
    pub fn get_with_offset_and_range(
        &self,
        offset: usize,
        source: &Object<T>,
        range: Range<usize>,
        pe: usize,
    ) -> crate::Result<()> {
        // bounds checking on dest base on offset and source.len()
        let num = range.end - range.start;
        if num <= 0 {
            return Ok(()); // nothing to do
        }
        if offset + num > self.local_elements {
            return Err(BoundsExceeded);
        }
        if range.end > source.local_elements {
            return Err(BoundsExceeded);
        }
        if pe >= Shmem::n_pes_raw() {
            return Err(InvalidPE);
        }
        let dest = self.void_ptr_with_offset(offset)?;
        let source = source.void_ptr_with_offset(range.start)?;
        // SAFETY: call to shmem function with checked arguments
        unsafe {
            shmem_sys::shmem_getmem(dest, source, num as u64, pe as c_int);
            Ok(())
        }
    }

    /// get variant where source is a GlobalObject
    pub fn get_global_with_offset_and_range(
        &self,
        offset: usize,
        source: &GlobalObject<T>,
        range: Range<usize>,
        pe: usize,
    ) -> crate::Result<()> {
        // bounds checking on dest base on offset and source.len()
        let num = range.end - range.start;
        if num <= 0 {
            return Ok(()); // nothing to do
        }
        if offset + num > self.local_elements {
            return Err(BoundsExceeded);
        }
        if range.end > source.local_elements {
            return Err(BoundsExceeded);
        }
        if pe >= Shmem::n_pes_raw() {
            return Err(InvalidPE);
        }
        let dest = self.void_ptr_with_offset(offset)?;
        let source = source.void_ptr_with_offset(range.start)?;
        // SAFETY: call to shmem function with checked arguments
        unsafe {
            shmem_sys::shmem_getmem(dest, source, num as u64, pe as c_int);
            Ok(())
        }
    }

    /// Present shmem put operations
    /// put a slice of memory to a remote pe for any type
    pub fn put(&self, source: &Object<T>, pe: usize) -> crate::Result<()> {
        self.put_with_offset_and_range(0, source, 0..source.local_elements, pe)
    }
    /// put variant with destination offset
    pub fn put_with_offset(
        &self,
        offset: usize,
        source: &Object<T>,
        pe: usize,
    ) -> crate::Result<()> {
        self.put_with_offset_and_range(offset, source, 0..source.local_elements, pe)
    }
    /// put variant with destination offset and source range
    pub fn put_with_offset_and_range(
        &self,
        offset: usize,
        source: &Object<T>,
        range: Range<usize>,
        pe: usize,
    ) -> crate::Result<()> {
        // bounds checking on dest base on offset and source.len()
        let num = range.end - range.start;
        if num <= 0 {
            return Ok(()); // nothing to do
        }
        if offset + num > self.local_elements {
            return Err(BoundsExceeded);
        }
        if range.end > source.local_elements {
            return Err(BoundsExceeded);
        }
        if pe >= Shmem::n_pes_raw() {
            return Err(InvalidPE);
        }
        let dest = self.void_ptr_with_offset(offset)?;
        let source = source.void_ptr_with_offset(range.start)?;
        // SAFETY: call to shmem function with checked arguments
        unsafe {
            shmem_sys::shmem_putmem(dest, source, num as u64, pe as c_int);
            Ok(())
        }
    }

    /// Do a shmem_put with non blocking.  Must call fence() before valid
    pub fn put_nbi_with_offset_and_range(
        &self,
        offset: usize,
        source: &Object<T>,
        range: Range<usize>,
        pe: usize,
    ) -> crate::Result<()> {
        // bounds checking on dest base on offset and source.len()
        let num = range.end - range.start;
        if num <= 0 {
            return Ok(()); // nothing to do
        }
        if offset + num > self.local_elements {
            return Err(BoundsExceeded);
        }
        if range.end > source.local_elements {
            return Err(BoundsExceeded);
        }
        if pe >= Shmem::n_pes_raw() {
            return Err(InvalidPE);
        }
        let dest = self.void_ptr_with_offset(offset)?;
        let source = source.void_ptr_with_offset(range.start)?;
        // SAFETY: call to shmem function with checked arguments
        unsafe {
            shmem_sys::shmem_putmem_nbi(dest, source, num as u64, pe as c_int);
            Ok(())
        }
    }
}

impl<T> Drop for Object<T> {
    fn drop(&mut self) {
        if self.is_subset {
            // don't free the subset objects, they will be freed by the full object
            // TODO: how safe is this?
            return;
        }
        if Shmem::active_counter() > 0 {
            // SAFETY: call to shmem function with argument unchanged from initial shmem alloc
            unsafe { shmem_sys::shfree(self.local_ptr as *mut c_void) }
        } else {
            // Can't panic here, as static state is freed after shmem_finalize().
            // Don't know how to fix this
            //panic!("Dropped object after finalize");
        }
    }
}

/// A shared global object, indexed by a global index
///  never implement Clone or Copy for this as the underlying memory is allocated via shmalloc()
#[derive(Debug)]
pub struct GlobalObject<T> {
    /// The local part of the underlying shmem array as a rust slice
    local_ptr: *mut T,
    local_elements: usize,
    /// The total elements in the user's request
    ///   n.b. this could be less than n_pes()*localpart.len()
    total_elements: usize,
    /// The type of blocking used
    ///  0 means pe index increments slowest ("block")
    ///  1 means pe index increments fastest ("cyclic")
    ///  TODO: >1 means block cyclic
    blocking: usize,
}

impl<T> GlobalObject<T>
where
    T: Copy,
{
    /// get a new shmem memory block as specified and build a slice representation of it
    ///  n.b., must round up to even number of element per PE
    pub fn new(total_elements: usize, blocking: usize) -> crate::Result<Self> {
        let local_elements = (total_elements + Shmem::n_pes_raw() - 1) / Shmem::n_pes_raw();
        let bytes = local_elements * std::mem::size_of::<T>();
        // SAFETY:
        //  - size of allocation is computed above
        //  - shmalloc returns uninitialized memory, so we initialize it
        //  - local_ptr is raw, so it can only be used in unsafe, see
        //    other SAFETY notes in this module
        let ptr = unsafe {
            let ptr = shmem_sys::shmalloc(bytes as u64) as *mut u8;
            let sl = std::slice::from_raw_parts_mut(ptr, bytes);
            sl.iter_mut().for_each(|x| *x = 0);
            shmem_sys::shmem_barrier_all();
            ptr as *mut T
        };
        if blocking > 1 {
            return Err(UnsupportedBlocking);
        }
        Ok(GlobalObject {
            local_ptr: ptr,
            local_elements,
            total_elements,
            blocking,
        })
    }

    /// Make a mutable pointer to the local part of the struture
    /// This is safe because the underlying memory is owned by shmalloc()
    pub fn local_part(&self) -> &mut [T] {
        // SAFETY: creates a correctly sized slice
        // Potential UNSAFETY: The shmem model allows multiple mutators of memory
        //   borrow checking.  Argument might want to be &mut self so we enforce
        //   borrow checking within the rust portion of a program, even if shmem
        //   is inherently unsafe due to remote puts and gets going on behind
        //   the scenes
        unsafe { std::slice::from_raw_parts_mut(self.local_ptr, self.local_elements) }
    }

    /// private to crate helper function to safely compute mutable pointer to send to shmem routine
    pub(crate) fn void_ptr_with_offset(&self, offset: usize) -> crate::Result<*mut c_void> {
        if offset > self.total_elements {
            return Err(BoundsExceeded);
        }
        // SAFETY: bounds are checked, result can only be used in an unsafe block
        //         Look to use cases for SAFETY messages
        Ok(unsafe { self.local_ptr.add(offset) as *mut c_void })
    }

    /// helper function to calculate local offsets and pes
    pub fn pe_loffset(&self, offset: usize) -> (usize, usize) {
        if self.blocking == 0 {
            (offset / self.local_elements, offset % self.local_elements)
        } else {
            (offset % Shmem::n_pes_raw(), offset / Shmem::n_pes_raw())
        }
    }

    /// Present shmem get operations to rust
    /// get a slice of memory from a remote pe for any type
    pub fn get(&self, source: &GlobalObject<T>) -> crate::Result<()> {
        Ok(self.get_with_range(source, 0..source.local_elements)?)
    }
    /// get variant where a source range is specifed
    pub fn get_with_range(
        &self,
        source: &GlobalObject<T>,
        range: Range<usize>,
    ) -> crate::Result<()> {
        if range.end - range.start > self.total_elements {
            return Err(BoundsExceeded);
        }
        let (pe, loffset) = source.pe_loffset(range.start);
        let src = source.void_ptr_with_offset(loffset)?;
        let dest = self.void_ptr_with_offset(0)?;
        // SAFETY: call to shmem function with checked arguments
        unsafe {
            shmem_sys::shmem_getmem(dest, src, (range.end - range.start) as u64, pe as c_int);
            Ok(())
        }
    }

    /// Present shmem put operations to rust
    /// put a slice of memory to a remote pe for any type
    pub fn put(&self, offset: usize, source: &GlobalObject<T>) -> crate::Result<()> {
        // bounds checking on dest base on offset and source.len()
        if offset + source.local_elements > self.total_elements {
            return Err(BoundsExceeded);
        }
        let (pe, loffset) = self.pe_loffset(offset);
        let dest = self.void_ptr_with_offset(loffset)?;
        // SAFETY: call to shmem function with checked arguments
        unsafe {
            shmem_sys::shmem_putmem(
                dest as *mut c_void,
                source.local_ptr as *const c_void,
                source.local_elements as u64,
                pe as c_int,
            );
            Ok(())
        }
    }
}

impl<T> Drop for GlobalObject<T> {
    fn drop(&mut self) {
        if Shmem::active_counter() > 0 {
            // SAFETY: call to shmem function with argument unchanged from initial shmem alloc
            unsafe { shmem_sys::shfree(self.local_ptr as *mut c_void) }
        } else {
            panic!("Dropped global object after finalize");
        }
    }
}
