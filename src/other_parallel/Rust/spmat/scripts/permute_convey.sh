#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=1
#SBATCH --nodes=2
#SBATCH --tasks-per-node=16
#SBATCH --constraint=haswell

srun $HOME/Rust/pshmem_private/target/release/examples/permute_convey
