//! Bale Serial Sparsemat error library
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For license information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use std::error;
use std::fmt;
use std::time::SystemTimeError;

#[derive(Debug)]
pub enum SparseMatError {
    Sparsemat(ParseMmError),
    Io(std::io::Error),
    Time(std::time::SystemTimeError),
    ParseInt(std::num::ParseIntError),
    ParseFloat(std::num::ParseFloatError),
    Re(regex::Error),
}

#[derive(Debug)]
pub struct ParseMmError {
    detail: String,
}

impl ParseMmError {
    pub fn new(msg: String) -> ParseMmError {
        ParseMmError { detail: msg }
    }
}

impl fmt::Display for ParseMmError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.detail)
    }
}

impl error::Error for ParseMmError {
    fn description(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for SparseMatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            SparseMatError::Sparsemat(ref e) => e.fmt(f),
            // This is a wrapper, so defer to the underlying types' implemenntation of `fmt`.
            SparseMatError::Io(ref e) => e.fmt(f),
            SparseMatError::Time(ref e) => e.fmt(f),
            SparseMatError::ParseInt(ref e) => e.fmt(f),
            SparseMatError::ParseFloat(ref e) => e.fmt(f),
            SparseMatError::Re(ref e) => e.fmt(f),
        }
    }
}

impl error::Error for SparseMatError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        match *self {
            SparseMatError::Sparsemat(ref e) => Some(e),
            // The cause is the underlying implementation error type. Is implicitly
            // cast to the trait object `&error::Error`. This works because the
            // underlying type already implements the `Error` trait.
            SparseMatError::Io(ref e) => Some(e),
            SparseMatError::Time(ref e) => Some(e),
            SparseMatError::ParseInt(ref e) => Some(e),
            SparseMatError::ParseFloat(ref e) => Some(e),
            SparseMatError::Re(ref e) => Some(e),
        }
    }
}

// Implement the conversion from `SystemTimeError` and `io::Error` to `SparseMatError`.
// This will be automatically called by `?` if a `SystemTimeError`
// needs to be converted into a `SparseMatError`.
impl From<SystemTimeError> for SparseMatError {
    fn from(err: SystemTimeError) -> SparseMatError {
        SparseMatError::Time(err)
    }
}

impl From<std::io::Error> for SparseMatError {
    fn from(err: std::io::Error) -> SparseMatError {
        SparseMatError::Io(err)
    }
}

impl From<std::num::ParseIntError> for SparseMatError {
    fn from(err: std::num::ParseIntError) -> SparseMatError {
        SparseMatError::ParseInt(err)
    }
}

impl From<std::num::ParseFloatError> for SparseMatError {
    fn from(err: std::num::ParseFloatError) -> SparseMatError {
        SparseMatError::ParseFloat(err)
    }
}

impl From<regex::Error> for SparseMatError {
    fn from(err: regex::Error) -> SparseMatError {
        SparseMatError::Re(err)
    }
}
