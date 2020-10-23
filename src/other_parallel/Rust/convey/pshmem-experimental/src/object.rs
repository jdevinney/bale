//! Module to implement objects from underlying shmem crate
//! Eventually this will be cfg'd for other underlying implementation crates

use crate::error::Error::{BoundsExceeded, UnsupportedBlocking};
use crate::PSHMEM_WORK;
use shmem::Shmem;
use std::ops::Range;
use std::rc::Rc;

/// A shared object, indexed by local offset and pe
#[derive(Debug)]
pub struct Object<T> {
    /// the implementation object
    pub shmem_object: shmem::object::Object<T>,
}

impl<T> Object<T>
where
    T: Copy,
{
    /// create a new object
    pub fn new(local_elements: usize) -> crate::Result<Self> {
        let shmem_object = shmem::object::Object::new(local_elements)?;
        Ok(Object { shmem_object })
    }

    /// return a slice to the local part
    pub fn local_part(&self) -> &mut [T] {
        self.shmem_object.local_part()
    }

    /// get a slice of memory from a remote rank for any type
    pub fn get(&self, source: &Object<T>, source_rank: usize) -> crate::Result<()> {
        self.get_with_offset_and_range(0, source, 0..source.local_part().len(), source_rank)
    }

    /// get variant with destination offset
    pub fn get_with_offset(
        &self,
        offset: usize,
        source: &Object<T>,
        source_rank: usize,
    ) -> crate::Result<()> {
        self.get_with_offset_and_range(offset, source, 0..source.local_part().len(), source_rank)
    }

    /// get variant with destination offset and source range
    pub fn get_with_offset_and_range(
        &self,
        offset: usize,
        source: &Object<T>,
        range: Range<usize>,
        source_rank: usize,
    ) -> crate::Result<()> {
        let num = range.end - range.start;
        if num <= 0 {
            return Ok(()); // nothing to do
        }
        if offset + num > self.local_part().len() {
            return Err(BoundsExceeded);
        }
        if range.end > source.local_part().len() {
            return Err(BoundsExceeded);
        }
        if offset + source.local_part().len() > self.shmem_object.local_part().len() {
            return Err(BoundsExceeded);
        }
        Ok(self.shmem_object.get_with_offset_and_range(
            offset,
            &source.shmem_object,
            range,
            source_rank,
        )?)
    }

    /// get variant for global objects
    pub fn get_global_with_range(
        &self,
        source: &GlobalObject<T>,
        range: Range<usize>,
    ) -> crate::Result<()> {
        if source.blocking > 0 {
            // we cannot do a gather
            return Err(UnsupportedBlocking);
        }
        let num = range.end - range.start;
        if num <= 0 {
            return Ok(()); // nothing to do
        }
        if source.local_part().len() > self.shmem_object.local_part().len() {
            return Err(BoundsExceeded);
        }
        let (source_rank, source_loffset) = source.pe_loffset(range.start);
        Ok(self.shmem_object.get_global_with_offset_and_range(
            0,
            &source.shmem_global_object,
            source_loffset..source_loffset + num,
            source_rank,
        )?)
    }

    /// put a slice of memory to a remote rank for any type
    pub fn put(&self, dest: &Object<T>, dest_rank: usize) -> crate::Result<()> {
        self.put_with_offset_and_range(0, dest, 0..dest.local_part().len(), dest_rank)
    }

    /// put variant with source offset
    pub fn put_with_offset(
        &self,
        offset: usize,
        dest: &Object<T>,
        dest_rank: usize,
    ) -> crate::Result<()> {
        self.put_with_offset_and_range(offset, dest, 0..dest.local_part().len(), dest_rank)
    }

    /// put variant with source offset and destination rnage
    pub fn put_with_offset_and_range(
        &self,
        offset: usize,
        dest: &Object<T>,
        range: Range<usize>,
        dest_rank: usize,
    ) -> crate::Result<()> {
        let num = range.end - range.start;
        if num <= 0 {
            return Ok(()); // nothing to do
        }
        if offset + num > self.local_part().len() {
            return Err(BoundsExceeded);
        }
        if range.end > dest.local_part().len() {
            return Err(BoundsExceeded);
        }
        if offset + dest.local_part().len() > self.shmem_object.local_part().len() {
            return Err(BoundsExceeded);
        }
        Ok(self.shmem_object.put_with_offset_and_range(
            offset,
            &dest.shmem_object,
            range,
            dest_rank,
        )?)
    }
}

/// Fetchers get data from other ranks and return it as either values or Vec
pub trait Fetchers {
    /// Our associated type, will be filled in by impl
    type T;
    /// fetch a remote element and return as value
    fn fetch(&self, source_rank: usize) -> crate::Result<Self::T>;
    /// fetch variant with offset
    fn fetch_with_offset(&self, offset: usize, source_rank: usize) -> crate::Result<Self::T>;
    /// fetch a remote slice and return as Vec
    fn fetch_slice(&self, source_rank: usize) -> crate::Result<Vec<Self::T>>;
    /// fetch_slcie variant with range
    fn fetch_slice_with_range(
        &self,
        range: Range<usize>,
        source_rank: usize,
    ) -> crate::Result<Vec<Self::T>>;
}

// This is a bit of a pain, but we need to instantiate fetchers on a per-type basis
// because they each use a different work array.  The following macro helps us do
// this fairly tersely
macro_rules! instantiate_fetchers {
    ($ty: ty, $work: ident) => {
        impl Fetchers for Object<$ty> {
            type T = $ty;
            fn fetch(&self, source_rank: usize) -> crate::Result<Self::T> {
                self.fetch_with_offset(0, source_rank)
            }
            fn fetch_with_offset(
                &self,
                offset: usize,
                source_rank: usize,
            ) -> crate::Result<Self::T> {
                let slice = PSHMEM_WORK.with({ |sw| Rc::clone(&sw.borrow().$work) });
                slice.get_with_offset(offset, self, source_rank)?;
                Ok(slice.local_part()[0])
            }

            fn fetch_slice(&self, source_rank: usize) -> crate::Result<Vec<Self::T>> {
                self.fetch_slice_with_range(0..self.local_part().len(), source_rank)
            }

            fn fetch_slice_with_range(
                &self,
                range: Range<usize>,
                source_rank: usize,
            ) -> crate::Result<Vec<Self::T>> {
                let slice = PSHMEM_WORK.with({ |sw| Rc::clone(&sw.borrow().$work) });
                // try to avoid creating/destroying new object
                if slice.local_part().len() >= range.end - range.start {
                    slice.get_with_offset_and_range(0, self, range, source_rank)?;
                    let mut ret: Vec<Self::T> = Vec::new();
                    ret.copy_from_slice(slice.local_part());
                    Ok(ret)
                } else {
                    let slice = Object::new(range.end - range.start)?;
                    slice.get_with_offset_and_range(0, self, range, source_rank)?;
                    let mut ret: Vec<Self::T> = Vec::new();
                    ret.copy_from_slice(slice.local_part());
                    Ok(ret)
                }
            }
        }
    };
}

instantiate_fetchers!(i64, work_i64);
instantiate_fetchers!(u64, work_u64);
instantiate_fetchers!(i32, work_i32);
instantiate_fetchers!(u32, work_u32);
instantiate_fetchers!(isize, work_isize);
instantiate_fetchers!(usize, work_usize);

/// A global shared object, indexed by global index
#[derive(Debug)]
pub struct GlobalObject<T> {
    /// Blocking strategy
    ///  0 means pe index increments slowest ("block")
    ///  1 means pe index increments fastest ("cyclic")
    ///  TODO: >1 means block cyclic
    blocking: usize,
    /// the implementation object
    shmem_global_object: shmem::object::GlobalObject<T>,
}

impl<T> GlobalObject<T>
where
    T: Copy,
{
    /// create a new global object
    pub fn new(total_elements: usize, blocking: usize) -> crate::Result<Self> {
        // Round to the next higher number of elements if not even
        let shmem_global_object = shmem::object::GlobalObject::new(total_elements, blocking)?;
        if blocking > 1 {
            return Err(UnsupportedBlocking);
        }
        Ok(GlobalObject {
            blocking,
            shmem_global_object,
        })
    }

    /// get local piece as slice
    pub fn local_part(&self) -> &mut [T] {
        self.shmem_global_object.local_part()
    }

    /// helper function to calculate local offsets and pes
    pub fn pe_loffset(&self, offset: usize) -> (usize, usize) {
        if self.blocking == 0 {
            (
                offset / self.local_part().len(),
                offset % self.local_part().len(),
            )
        } else {
            (offset % Shmem::n_pes_raw(), offset / Shmem::n_pes_raw())
        }
    }
}

/// Global Fetchers get data from other ranks and return it as either values or Vec
pub trait GlobalFetchers {
    /// Our associated type, will be filled in by impl
    type T;
    /// fetch a remote element and return as value
    fn fetch(&self) -> crate::Result<Self::T>;
    /// fetch variant with offset
    fn fetch_with_offset(&self, offset: usize) -> crate::Result<Self::T>;
    /// fetch a remote slice and return as Vec
    fn fetch_slice(&self) -> crate::Result<Vec<Self::T>>;
    /// fetch_slcie variant with range
    fn fetch_slice_with_range(&self, range: Range<usize>) -> crate::Result<Vec<Self::T>>;
}

// This is a bit of a pain, but we need to instantiate fetchers on a per-type basis
// because they each use a different work array.  The following macro helps us do
// this fairly tersely
macro_rules! instantiate_global_fetchers {
    ($ty: ty, $work: ident) => {
        impl GlobalFetchers for GlobalObject<$ty> {
            type T = $ty;
            fn fetch(&self) -> crate::Result<Self::T> {
                self.fetch_with_offset(0)
            }
            fn fetch_with_offset(&self, offset: usize) -> crate::Result<Self::T> {
                //let slice = GlobalObject::new(1, 0)?;
                let slice = PSHMEM_WORK.with({ |sw| Rc::clone(&sw.borrow().$work) });
                slice.get_global_with_range(self, offset..offset + 1)?;
                Ok(slice.local_part()[0])
            }

            fn fetch_slice(&self) -> crate::Result<Vec<Self::T>> {
                self.fetch_slice_with_range(0..self.local_part().len())
            }

            fn fetch_slice_with_range(&self, range: Range<usize>) -> crate::Result<Vec<Self::T>> {
                //let slice = GlobalObject::new(range.end - range.start, 0)?;
                let slice = PSHMEM_WORK.with({ |sw| Rc::clone(&sw.borrow().$work) });
                // try to avoid creating/destroying new object
                if slice.local_part().len() >= range.end - range.start {
                    slice.get_global_with_range(self, range)?;
                    let mut ret: Vec<Self::T> = Vec::new();
                    ret.copy_from_slice(slice.local_part());
                    Ok(ret)
                } else {
                    let slice = Object::new(range.end - range.start)?;
                    slice.get_global_with_range(self, range)?;
                    let mut ret: Vec<Self::T> = Vec::new();
                    ret.copy_from_slice(slice.local_part());
                    Ok(ret)
                }
            }
        }
    };
}

instantiate_global_fetchers!(i64, work_i64);
instantiate_global_fetchers!(u64, work_u64);
instantiate_global_fetchers!(i32, work_i32);
instantiate_global_fetchers!(u32, work_u32);
instantiate_global_fetchers!(isize, work_isize);
instantiate_global_fetchers!(usize, work_usize);
