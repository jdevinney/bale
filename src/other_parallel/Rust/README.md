# Bale Applications written in Rust

This is the implementation of the Bale applications, using the Rust
programming language.

## Cray Installation

If you are attempting to build on a Cray system, this is currently
tricky.  See the [Cray Instructions](README_cray.md) and ignore 
the rest of these instructions.

## Non-Cray Installation

It depends on the rust conveyor library [available
here](https://github.com/wwc559/convey_private). It is 
strongly advised that you do a *standalone build and test* of 
`convey_private` *first*, following the instructions in that library.
This can be a *bit tricky* as you probably need to specify where the
underlying parallel communication library is installed on
your system, via the `SHMEM_PATH` environment variable.

This build is setup to run correctly if you place `convey_private`
in the same directory where you placed bale.  Modify Cargo.toml,
delta_stepping/Cargo.toml, spmat/Cargo.toml, triangle/Cargo.toml, 
and toposort/Cargo.toml if you have placed it in a different location.

To build this and say `cargo build --release --workspace` and then you can 
run any of the bale apps with:

```
oshrun -n 4 target/release/delta_stepping
oshrun -n 4 target/release/histo_convey
oshrun -n 4 target/release/ig_convey
oshrun -n 4 target/release/permute_convey
oshrun -n 4 target/release/randperm_convey
oshrun -n 4 target/release/toposort
oshrun -n 4 target/release/triangle
```



