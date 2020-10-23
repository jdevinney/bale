//! Module to handle errors that occur in shmem
//!  Derived from https://blog.burntsushi.net/rust-error-handling/

/// Enum to keep all the error types we and our dependencies use
#[derive(Debug)]
pub enum Error {
    /// An error that occured in a called function by shmem
    Io(std::io::Error),
    /// An error that occured in a called function in shmem
    Shmem(shmem::error::Error),
    /// Some bound on an object::Object or object::GlobalObject was exceeded
    BoundsExceeded,
    /// We don't support arbitrary blocking (yet?)
    UnsupportedBlocking,
}

impl std::fmt::Display for Error {
    /// Format one of our errors for display to user
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match *self {
            Error::Io(ref err) => err.fmt(f),
            Error::Shmem(ref err) => err.fmt(f),
            Error::BoundsExceeded => write!(f, "Shmem bounds exceeded on remote operation"),
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

impl From<shmem::error::Error> for Error {
    /// Pull the shmem error into our space
    fn from(err: shmem::error::Error) -> Error {
        Error::Shmem(err)
    }
}
