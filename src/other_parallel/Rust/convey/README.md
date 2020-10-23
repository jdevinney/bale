# convey_private
Conveyor Implementation for Rust

This project is to explore implementing conveyors in the Rust
programming language.  This is inspired by many projects in the PGAS
research community including UPC, OpenShmem, Bale, and Conveyors.

**Note: This repository is currently private, please do not redistribute.**
Open an 'issue' to get access for others.

**Note: Expect rapid changes in code organization and interfaces, this is highly experimental.**

## Cray Build Instructions 

To build this package on Cray, see [Cray
Instructions](README_cray.md).  Ignore what follows here as it is
currently complicated.

## Non-Cray Build Instructions 

To build this package you will need an implementation of openshmem,
version 1.4.  If not, you will need to set the environment
variable SHMEM_PATH to the path where openshmem was installed.  The
default is /usr/local.  Also you will need a version of llvm in your
path, specifically the build needs to be able to run 'llvm-config --prefix'

To build this just say `cargo build --release --examples` and then you can run it with:

```
oshrun -n 4 target/release/examples/histo_convey
```

There are sub-crates for adapting to implementation technologies this is
build to top of.

`shmem` and `shmem-sys` implement an adaptation to the OpenShmem (1.4)
system.  `shmem-sys` is the wrapper for a 1.4-confomrant OpenShmem
library.  `shmem` provides the a subset of these features to `pshmem`.
See the [Shmem README.md](shmem/README.md) for more details, including
some examples of using this directly.

