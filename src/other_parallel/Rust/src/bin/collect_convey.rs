#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]
//! Example program to do collectives in conveyors
//
// Copyright (c) 2020, Institute for Defense Analyses
// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
//
// All rights reserved.
//
// This file is part of Bale.  For license information see the
// LICENSE file in the top level dirctory of the distribution.
//
use clap::{App, Arg};
use convey_hpc::collect::CollectValues;
use convey_hpc::Convey;
use std::time::Instant;

fn main() {
    let matches = App::new("collect")
        .version("0.1.0")
        .about("test of collectives")
        .arg(
            Arg::with_name("collectives")
                .short("c")
                .long("collectives")
                .takes_value(true)
                .help("the number of collective calls per rank"),
        )
        .arg(
            Arg::with_name("verbose")
                .short("-v")
                .long("verbose")
                .takes_value(false)
                .help("increase the amount of verbosity"),
        )
        .get_matches();

    let collectives: usize = matches
        .value_of("collectives")
        .unwrap_or("100")
        .parse()
        .expect("bad collectives arg");
    let verbose: u64 = matches.occurrences_of("verbose");

    let convey = Convey::new().expect("bad");
    do_collect_convey(&convey, collectives, verbose);
    do_collect_shmem(&convey, collectives, verbose);
}

fn do_collect_convey(convey: &Convey, collectives: usize, verbose: u64) {
    let value: usize = 42;
    let now = Instant::now();
    let my_rank = convey.my_rank();
    let num_ranks = convey.num_ranks();

    let mut total_result = 0;

    for _i in 0..collectives {
        total_result += value.reduce_sum();
    }
    let d = now.elapsed();
    if verbose > 0 || my_rank == 0 {
        println!(
            "convey pe{}/{}, {} reductions, {} reductions, {} msec",
            my_rank,
            num_ranks,
            collectives,
            total_result,
            d.as_millis(),
        );
    }
}

fn do_collect_shmem(convey: &Convey, collectives: usize, verbose: u64) {
    let _value: usize = 42;
    let my_rank = convey.my_rank();
    let num_ranks = convey.num_ranks();

    let now = Instant::now();

    let total_result = 0;

    for _i in 0..collectives {
        //total_result += Shmem::sum_to_all(value);
    }
    let d = now.elapsed();
    if verbose > 0 || my_rank == 0 {
        println!(
            "shmem pe{}/{}, {} reductions, {} reductions, {} msec",
            my_rank,
            num_ranks,
            collectives,
            total_result,
            d.as_millis(),
        );
    }
}
