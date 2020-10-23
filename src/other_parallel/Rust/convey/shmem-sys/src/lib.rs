// Some allows to prevent warnings
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]
//
// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of Convey, a conveyor library for rust.  For
// licence information see the LICENSE file in the top level dirctory
// of the distribution.

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

// bother: need a different type in various versions of shmem.h
#[cfg(cray)]
pub type nreduce_t = size_t;

#[cfg(not(cray))]
pub type nreduce_t = i32;
