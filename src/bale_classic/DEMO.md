# Table of Contents

* [Building on some common platforms](#Building-on-some-common-platforms)
  * Linux SMP OpenMPI/OSHMEM
  * Linux SMP Sandia OpenSHMEM
  * Linux SMP with GNU or CLANG UPC
  * Cray XC30
  * Mac OSX
* [Running Examples](#Running-Examples)



## Building on some common platforms

For this document, let BALEDIR = the directory containing this file.

## ... on a Linux SMP with OpenMPI/OSHMEM

You should have a reasonably modern C complier (gcc 8.3.0 seems to be working well for us) and 'oshcc' in your path. You also need OpenMPI 4.x.x or higher.

```bash
cd $BALEDIR
export PLATFORM=linux_oshmem
export CC=oshcc
./bootstrap.sh
./make_bale -s
```

This builds everything in $BALEDIR/build_linux_oshmem. Binaries appear in $BALEDIR/build_linux_oshmem/bin.

**Note**: Due to a bug in OpenMPI/OpenSHMEM, exstack2 and conveyor apps sometimes hang. We recommend you avoid running them under oshmem by running application with the implementation mask set to 3 (-M 3).

## ... on a Linux SMP with Sandia OpenSHMEM

You should have a reasonably modern C complier (gcc 8.3.0 seems to be working well for us) and 'oshcc' in your path. You also need SOS 1.4.5 or higher, preferrably built on XPMEM. See note below.

Otherwise, exactly the same instructions as above substituting "linux_sos" for PLATFORM.

**Note**: Due to a bug in SOS built using CMA transport layer, exstack2 and conveyor apps sometimes hang. We recommend you avoid running them under SOS by running application with the implementation mask set to 3 (-M 3).

## ... on a Linux SMP with GNU UPC (GUPC) or Clang UPC (CUPC)

You should have a reasonably modern C complier (gcc 8.3.0 seems to be working well for us) and gupc/upc in your path for GUPC or CUPC. You should also have GUPC 5.2.0.1 or Clang-UPC 3.9.1.1.

```bash
cd $BALEDIR
export PLATFORM=linux_gupc
export UPC=gupc (for GUPC)
export UPC=upc (for CUPC)
unset CC
./bootstrap.sh
./make_bale -u
```

This builds everything in $BALEDIR/build_linux_oshmem. Binaries appear in $BALEDIR/build_linux_oshmem/bin.

## ... on a Cray XC30

- If you want to use Cray UPC:
    - Use the PrgEnvCray module to access the Cray UPC compiler.
    - `export UPC="cc -hupc"`
- If you want to use Cray SHMEM: Use PrgEnvCray and load the cray-shmem module.
- If you want to use Cray OpenSHMEMX: Use PrgEnvgnu and load cray-openshmemx. 

```bash
cd $BALEDIR
export PLATFORM=xc30
./bootstrap.sh
./make_bale -u
```

## ... on Mac OS X

This is a bit tricky :)

1. install osx compiling environment, bring up the mac app store and install xcode
   After installation you need to accept the license:

   ```bash
   sudo xcodebuild -license
   ```

2. to get autotools, best to use brew (https://brew.sh):

   ```bash
   brew install autoconf
   brew install automake
   brew install libtool
   ```

3. to get clang-upc, which seems to run well, go to:
   https://clangupc.github.io/clang-upc/install.html
   and follow the directions. You will need to do

   ```bash
   brew install cmake
   ```

to get this to work.  On recent OS X versions (10.15, maybe 10.14) you need
    to pass another argument to cmake:

     `-DDEFAULT_SYSROOT:STRING="$(xcrun --show-sdk-path)"`

4. to use shmem, SOS openShmem seems best for now.  Go to:
   https://github.com/Sandia-OpenSHMEM/SOS/wiki/OFI-Build-Instructions
   and follow the directions.  For recent versions OS X (10.15, maybe 10.14) you
   need to set an environment variable before runnig anything will work.
   `export SHMEM_OFI_DOMAIN=lo0`
   performance is horrid because it only has a socket based implementation

5. when using the runit.sh script with clang, you need to say

   ```bash
   ./runit.sh -l $PWD/clang_upc_run.sh -c 2
   ```

   Note that -c should be small but not less than 2, otherwise some things hang

## Running Examples

Try running a simple test (remember to use -M 3 with OpenMPI/oshmem or SOS)
### with oshrun
```bash
oshrun -n 4 $BALEDIR/build_$PLATFORM/bin/histo
```
### with slurm
```bash
srun -n 4 $BALEDIR/build_$PLATFORM/bin/histo
```
### with gupc
```bash
$BALEDIR/build_$PLATFORM/bin/histo -n 4
```

You should see something like this...

```bash

***************************************************************
Bale Version 3.00 (UPC 201311): 2110-08-17.12:50
Running command on 4 PEs: ../build_ucs3_gupc/bin/histo
***************************************************************

num_updates_per_pe: 100000
table_size_per_pe: 1000
Standard options:
----------------------------------------------------
buf_cnt (buffer size)    (-b): 1024
seed                     (-s): 122222
cores_per_node           (-c): 0
Models Mask              (-M): 15

       AGI:    0.009 seconds     0.000 GB/s injection bandwidth
   Exstack:    0.003 seconds     0.000 GB/s injection bandwidth
  Exstack2:    0.004 seconds     0.000 GB/s injection bandwidth
  Conveyor:    0.012 seconds     0.000 GB/s injection bandwidth
```

## Running a suite of tests

The run_apps.py script allows you to launch a suite of tests of one or more bale apps. Run with the --help option to see the usage.

## Analyzing tests

All bale apps can output profiling data into a json file. Use the run_apps.py script to run multiple runs and collect the data into one consolidated json file.  We have included an example Jupyter notebook (plot_results.ipynb) that demonstrates reading the json file and plotting some results. 