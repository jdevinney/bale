#!/bin/sh
#
#
#  Copyright(C) 2018-2020, Institute for Defense Analyses
#  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
#
#  All rights reserved.
#  
#  This file is a part of Bale.  For license information see the
#  LICENSE file in the top level directory of the distribution.
#  

# this script runs all of the bale apps (C versions) with trivial parameters
# It exits if any application exits abnormally.
#
cores_per_node=""
quick_run=0

while getopts ":c:qM:" opt; do
    case $opt in
        q ) quick_run=1;;
        \? ) echo 'usage: runall [-q]'
             exit 1
    esac
done

# Let make see that they are all up to date
make 
#

if [ $quick_run == 1 ]; then
    # this makes the tests run a little quicker, if you want to run a longer
    # set of tests, use a second argument (can be anything!)
    options+=" -n 1000"
fi

echo; echo XXXXXXXXXXXX demo_spmat XXXXXXXXXXXXXXX
./demo_spmat


for app in histo ig topo randperm permute_matrix transpose_matrix triangles unionfind sssp
do
    # just run the command with -h
    echo; echo XXXXXXXXXXXX $app XXXXXXXXXXXXXXX
    echo;
    ./$app 
    echo;
done
