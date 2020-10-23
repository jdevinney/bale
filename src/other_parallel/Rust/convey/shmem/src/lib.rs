#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]
//! This package implements the rust interface to shmem-sys, which in turn
//!  is a wrapper for OpenShmem implmentations, currently version 1.4
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Convey, a conveyor library for rust.  For
/// licence information see the LICENSE file in the top level dirctory
/// of the distribution.

/// Generic result type for this library
pub type Result<T> = std::result::Result<T, error::Error>;

/// Our instance struct, currently no per-instance state
#[derive(Debug)]
pub struct Shmem {}

pub mod atomic;
pub mod collect;
pub mod error;
pub mod object;
pub mod shmem;

// Tests are in shmem.rs.  Unfortunately cannot have tests in other modules.
//    (See test notes in shmem.rs)  Todo: this could be done
