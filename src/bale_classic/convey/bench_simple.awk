# Copyright (c) 2020, Institute for Defense Analyses
# 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
#
# All rights reserved.
#
# This file is part of the conveyor package. For license information,
# see the LICENSE file in the top level directory of the distribution.

BEGIN {
    if (buf == 0) {
        print "set approximate buffer bytes with -vbuf=..."
        exit
    }
    if (ppn <= 0 || ppn % 2 == 1) {
        print "set procs/node to an even number with -v ppn=..."
        exit
    }
    if (max <= 0) {
        print "set maximum nodes with -v max=..."
        exit
    }

    bytes[0] = 8; bytes[1] = 16; bytes[2] = 32; bytes[3] = 128
    print "#!/bin/sh -e"
    print "export MPP_INIT_QUIET=1"
    nodes = 2; power = 2
    for (i = 0; nodes <= max; i++) {
        print "$NLAUNCH", ppn * nodes, "./alltoall? -- <<EOF"
        for (j = 0; j < 4; j++) {
            size = bytes[j]
# each process wants to send 256 MiB
            load = 268435456 / (nodes * ppn * size)
            cap = buf
# make sure buffers occupy at most 1GB
            while (cap * nodes * 32 * 2 > 1000 * 1000 * 1000) {
                cap = cap / 2
            }
            for (k = 0; k < 2; k++) {
                option = (k == 0) ? "" : "-x "
                printf("-c%d -n3 %ssimple %d %d\n", cap, option, load, size)
            }
        }
        print "EOF"
        if (i % 2 == 1) {
            power *= 2
            nodes = power
        } else {
            if (power >= 32) nodes = 45 * power / 32
            else if (power >= 8) nodes = 11 * power / 8
            else nodes = 3 * power / 2
        }
    }
    print "exit 0"
}
