#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Bale Serial Toposort application
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For license information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use clap::{App, Arg};
use toposort::generate_toposort_input;
use toposort::TopoSort;

// Demo application that finds an upper triangular form for a matrix.
// That is, we are given a matrix that is a random row and column
// permutation of a an upper triangular matrix (with ones on the
// diagonal).  This algorithm finds a row and column permutation that
// would return it to an upper triangular form.

fn main() {
    let matches = App::new("TopoSort")
        .version("0.1.0")
        .about("Implements a test of TopoSort")
        .arg(
            Arg::with_name("numrows")
                .short("n")
                .long("numrows")
                .takes_value(true)
                .help("The number of rows in test matrix"),
        )
        .arg(
            Arg::with_name("er_prob")
                .short("e")
                .long("erdos_renyi_prob")
                .takes_value(true)
                .help("Probability of an edge in the erdos renyi graph"),
        )
        .arg(
            Arg::with_name("dump_files")
                .short("d")
                .long("dump_files")
                .takes_value(false)
                .help("Produce short dumps as the algorithm progresses"),
        )
        .arg(
            Arg::with_name("quiet")
                .short("q")
                .long("quiet")
                .takes_value(false)
                .help("produce less chatty output"),
        )
        .get_matches();

    // input args, just constants for now
    let numrows: usize = matches
        .value_of("numrows")
        .unwrap_or("20000")
        .parse()
        .expect("numrows: not an integer");
    let erdos_renyi_prob: f64 = matches
        .value_of("er_prob")
        .unwrap_or("0.1")
        .parse()
        .expect("er_prob: not an float");

    let seed = 12346;
    let quiet = matches.is_present("quiet");
    let dump_files = matches.is_present("dump_files");

    if !quiet {
        println!("creating input matrix for toposort");
    }

    let mat = generate_toposort_input(numrows, erdos_renyi_prob, seed, dump_files);

    if !quiet {
        println!("input matrix stats:");
        mat.stats();
    }

    if dump_files {
        mat.dump(20, "mat.out").expect("could not write mat.out");
    }

    mat.write_mm_file("topo_mat.mm")
        .expect("could not write topo_mat.mm");

    let tmat = mat.transpose();

    if dump_files {
        tmat.dump(20, "trans.out")
            .expect("could not write trans.out");
    }

    tmat.write_mm_file("topo_tmat.mm")
        .expect("could not write topo_tmat.mm");

    if !quiet {
        println!("Running toposort on mat (and tmat) ...");
    }

    let mut matret;
    for i in 0..2 {
        if i == 0 {
            if !quiet {
                print!("   using generic toposort: ");
            }
            matret = mat.toposort_queue(&tmat);
        } else {
            if !quiet {
                print!("   using loop toposort: ");
            }
            matret = mat.toposort_loop(&tmat);
        }
        if !mat.check_result(&matret, dump_files) {
            println!("ERROR: check_result failed");
        }
        if !quiet {
            println!(" {} seconds", matret.laptime)
        }
    }
}
