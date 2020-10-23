//! Module to handle errors that occur in shmem
//!  Derived from https://blog.burntsushi.net/rust-error-handling/

/// Enum to keep all the error types we and our depenencies use
#[derive(Debug)]
pub enum Error {
    /// An error that occured in a called function
    Io(std::io::Error),
    /// We cannot re-establish a shmem_init() after shmem_finalize() was called
    ///  This should only occur if a user of our package calls new() multiple times
    NewAfterDrop,
    /// Some bound on an object::Object or object::GlobalObject was exceeded
    BoundsExceeded,
    /// Tried to reference an invalid PE
    InvalidPE,
    /// We don't support arbitrary blocking (yet?)
    UnsupportedBlocking,
}

impl std::fmt::Display for Error {
    /// Format one of our errors for display to user
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Error::Io(ref err) => err.fmt(f),
            Error::NewAfterDrop => write!(f, "Shmem::new() after last Shem::drop()"),
            Error::BoundsExceeded => write!(f, "Shmem bounds exceeded on remote operation"),
            Error::InvalidPE => write!(f, "PE number exceeds num_pes()"),
            Error::UnsupportedBlocking => write!(f, "Unsupported blocking factor"),
        }
    }
}

impl From<std::io::Error> for Error {
    /// Pull the io error into our space
    fn from(err: std::io::Error) -> Error {
        Error::Io(err)
    }
}
