#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]
//! Example program to exercise random permutatons in conveyors
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
use convey_hpc::Convey;
use spmat::perm::Perm;
use std::time::Instant;

fn main() {
    let matches = App::new("randperm")
        .version("0.1.0")
        .about("test of random permutation")
        .arg(
            Arg::with_name("numrows")
                .short("n")
                .long("numrows")
                .takes_value(true)
                .help("the number of rows per rank"),
        )
        .arg(
            Arg::with_name("verbose")
                .short("-v")
                .long("verbose")
                .takes_value(false)
                .help("increase the amount of verbosity"),
        )
        .arg(
            Arg::with_name("quiet")
                .short("-q")
                .long("quiet")
                .takes_value(false)
                .help("make quiet"),
        )
        .get_matches();

    let numrows_per_rank: usize = matches
        .value_of("numrows")
        .unwrap_or("1000000")
        .parse()
        .expect("bad numrows arg");
    let verbose: u64 = matches.occurrences_of("verbose");

    let convey = Convey::new().expect("failed conveyor initializtoin");
    let my_rank = convey.my_rank();
    let num_ranks = convey.num_ranks();
    let quiet = matches.is_present("quiet") || (verbose == 0 && my_rank > 0);

    let numrows = numrows_per_rank * num_ranks;

    if !quiet {
        println!("Running randperm on {} ranks", num_ranks);
        println!("Number of rows per rank     (-n) {}", numrows_per_rank);
    }

    let mut perm;
    for i in 0..2 {
        let now = Instant::now();
        if i == 0 {
            if !quiet {
                print!("   using random_darts_resubmit_local: ");
            }
            perm = Perm::random_darts_resubmit_local(numrows, 0);
        } else {
            if !quiet {
                print!("   using random_darts_return_rejects: ");
            }
            perm = Perm::random_darts_resubmit_local(numrows, 0);
        }
        let d = now.elapsed();
        if !quiet {
            println!(" {} milliseconds", d.as_millis());
        }
        if !perm.is_perm() {
            println!("ERROR: check_result failed");
        }
    }
}
