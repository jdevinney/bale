#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=1
#SBATCH --nodes=8
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell

srun $HOME/Rust/pshmem_private/target/release/toposort -n 20000
