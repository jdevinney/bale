# Copyright (c) 2020, Institute for Defense Analyses
# 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
#
# All rights reserved.
#
# This file is part of the conveyor package. For license information,
# see the LICENSE file in the top level directory of the distribution.

BEGIN {
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
            if (nodes <= 64) {
                printf("-b2 -n3 vector %d %d\n", load, size);
            }
            if (nodes <= 128) {
                printf("-b2 -n3 -t%d matrix %d %d\n", ppn / 2, load, size);
            }
            printf("-b2 -n3 -t%d matrix %d %d\n", ppn, load, size);
            if (nodes > ppn) {
                printf("-b2 -n3 -t%d tensor %d %d\n", ppn / 2, load, size);
            }
            if (nodes >= 1024) {
                printf("-b2 -n3 -t%d tensor %d %d\n", ppn, load, size);
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
