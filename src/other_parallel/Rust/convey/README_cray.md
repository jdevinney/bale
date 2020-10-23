# Installation instructions for Cray

The cray build is complicated, at least on the NERSC Cori system.  The
crux of the probem is that some Rust crates will not build with the
llvm module loaded while others will not build without it loaded.
This is unfortunate and we will look for solutions that work better
going forward.

1. Build with default modules loaded, shmem-sys will fail
```
    cargo build --release
```
2. source the fix-modules.sh script
```
    source scripts/fix-modules.sh
```	
3. rebuild, it should work this time.
```
    cargo build --release --workspace
    cargo build --release --examples
```
To run the application:

```
srun -n 4 target/release/examples/histo_convey
```

If you are building bale_private/src/other_parallel/Rust, the
instructions are similar, so just build there (on a fresh
login, without having loaded fix-modules.sh)


