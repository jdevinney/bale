#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]
//! Program to do permutations in conveyors
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
use spmat::SparseMat;
use std::time::Instant;

fn main() {
    let matches = App::new("randperm")
        .version("0.1.0")
        .about("test of permute matrix")
        .arg(
            Arg::with_name("numrows")
                .short("n")
                .long("numrows")
                .takes_value(true)
                .help("the number of rows per rank"),
        )
        .arg(
            Arg::with_name("nz_per_row")
                .short("Z")
                .long("nz_per_row")
                .takes_value(true)
                .help("Number of non-zeros per row"),
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
        .unwrap_or("10000")
        .parse()
        .expect("bad numrows arg");
    let verbose: u64 = matches.occurrences_of("verbose");
    let mut erdos_renyi_prob: f64 = matches
        .value_of("er_prob")
        .unwrap_or("0.0")
        .parse()
        .expect("er_prob: not a float");
    let mut nz_per_row: f64 = matches
        .value_of("nz_per_row")
        .unwrap_or("35")
        .parse()
        .expect("nz_per_row: not a float");

    let convey = Convey::new().expect("conveyor initializtion failed");
    let my_rank = convey.my_rank();
    let num_ranks = convey.num_ranks();
    let quiet = matches.is_present("quiet") || (verbose == 0 && my_rank > 0);
    let numrows = numrows_per_rank * num_ranks;

    if erdos_renyi_prob == 0.0 {
        // use nz_per_row to get erdos_renyi_prob
        erdos_renyi_prob = (2.0 * nz_per_row) / (numrows - 1) as f64;
        if erdos_renyi_prob > 1.0 {
            erdos_renyi_prob = 1.0;
        }
    } else {
        // use erdos_renyi_prob to get nz_per_row
        nz_per_row = erdos_renyi_prob * numrows as f64
    }

    if !quiet {
        println!("Running randperm on {} ranks", num_ranks);
        println!("Number of rows per rank     (-n) {}", numrows_per_rank);
        println!("Avg # non-eros per row      (-Z) {}", nz_per_row);
        println!("Erdos-Renyi edge probabilty (-e) {}", erdos_renyi_prob);
    }

    let rp = Perm::random(numrows, 0);
    let cp = Perm::random(numrows, 0);

    let mat = SparseMat::gen_erdos_renyi_graph(numrows, erdos_renyi_prob, true, 2, 0);

    for i in 0..1 {
        let now = Instant::now();
        if i == 0 {
            if !quiet {
                print!("   using permute_matrix: ");
            }
            let _outmat = mat.permute(&rp, &cp);
        }
        let d = now.elapsed();
        if !quiet {
            println!(" {} milliseconds", d.as_millis());
        }
    }
}
