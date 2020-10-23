# Cray Build Instructions

We depend on the rust conveyor library [available
here](https://github.com/wwc559/convey_private). Install this
first but DO NOT BUILD.

The Cray build is complicated, at least on the NERSC Cori system.
The crux of the probem is that some Rust crates will not build with the
llvm module loaded while others will not build without it loaded.
This is unfortunate and we will look for solutions that work
better going forward.

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
srun -n 4 target/release/delta_stepping
srun -n 4 target/release/histo_convey
srun -n 4 target/release/ig_convey
srun -n 4 target/release/permute_convey
srun -n 4 target/release/randperm_convey
srun -n 4 target/release/toposort
srun -n 4 target/release/triangle
```

If you are building convey_private, the instructions are similar, so
just build there (on a fresh login, without having loaded
fix-modules.sh)


