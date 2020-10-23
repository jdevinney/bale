# Copyright (c) 2020, Institute for Defense Analyses
# 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
#
# All rights reserved.
#
# This file is part of the conveyor package. For license information,
# see the LICENSE file in the top level directory of the distribution.


use strict;

open(CONFIG, "config.h") || die "cannot open config.h for reading\n";
my @config = <CONFIG>;
close(CONFIG);
my $enabled = grep /ENABLE_NONBLOCKING 1/, @config;
my $old_enabled = $enabled;

my @lines = <>;
my @blines = grep /^1/, @lines;
my @nblines = grep /^0/, @lines;
die "incomplete data for tuning\n" if $#blines != $#nblines;
die "insufficient data for tuning\n" if $#blines < 10;

my @bdata = map { (split)[2] } @blines;
my @nbdata = map { (split)[2] } @nblines;

if (! $enabled) {
    for (my $i = 0; $i < @bdata; $i++) {
        $bdata[$i] = ($bdata[$i] + $nbdata[$i]) / 2;
    }
}
my ($peak, $value) = &findPeak(\@bdata);
if ($enabled) {
    my ($npeak, $nvalue) = &findPeak(\@nbdata);
    if ($nvalue >= 0.95 * $value) {
        $peak = $npeak;
        $value = $nvalue;
        @blines = @nblines;
    }
    else {
        $enabled = 0;
    }
}

my $bufsiz = (split ' ',$blines[$peak])[1];
my $old_bufsiz = (join '', grep /CONVEY_BUFFER_SIZE [0-9]+\s*$/, @config);
$old_bufsiz =~ s/^.*CONVEY_BUFFER_SIZE ([0-9]+)\s$/\1/;

printf "Estimated bandwidth is %.1f MB/sec/PE\n", ($value);
print "Setting ENABLE_NONBLOCKING to $enabled (was $old_enabled)\n";
print "Setting CONVEY_BUFFER_SIZE to $bufsiz (was $old_bufsiz)\n";
print "Rerun 'make' to build these values into the library.\n";

my $config = join '', @config;
$config =~ s/ENABLE_NONBLOCKING [0-1]/ENABLE_NONBLOCKING $enabled/gm;
$config =~ s/CONVEY_BUFFER_SIZE [0-9]+\s*$/CONVEY_BUFFER_SIZE $bufsiz/gm;
open(CONFIG, ">config.h") || die "cannot open config.h for writing\n";
print CONFIG $config;
close(CONFIG);

exit;


sub findPeak {
    my ($array) = @_;
    my $peak = 0;
    my $value = 0;
    for (my $i = 1; $i < @$array - 1; $i++) {
        my $avg = ($$array[$i-1] + $$array[$i] + $$array[$i+1]) / 3;
        if ($avg > $value) {
            $value = $avg;
            $peak = $i;
        }
    }
    return ($peak, $value);
}
