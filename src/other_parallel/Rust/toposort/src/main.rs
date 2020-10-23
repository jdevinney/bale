//! Toposort application
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
use convey_hpc::Convey;
use toposort::generate_toposort_input;
use toposort::TopoSort;

// Demo application that finds an upper triangular form for a matrix.
// That is, we are given a matrix that is a random row and column permutation
// of a an upper triangular matrix (with ones on the diagonal).
// This algorithm finds a row and column permutation that would return it
// to an upper triangular form.

fn main() {
    let matches = App::new("TopoSort")
        .version("0.1.0")
        .about("Implements a test of TopoSort")
        .arg(
            Arg::with_name("numrows_per_rank")
                .short("n")
                .long("numrows_per_rank")
                .takes_value(true)
                .help("The number of rows per rank"),
        )
        .arg(
            Arg::with_name("nz_per_row")
                .short("Z")
                .long("nz_per_row")
                .takes_value(true)
                .help("Number of non-zeros per row"),
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
    let numrows_per_rank: usize = matches
        .value_of("numrows_per_rank")
        .unwrap_or("20000")
        .parse()
        .expect("numrows_per_rank: not an integer");
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

    let seed = 12346;
    let convey = Convey::new().expect("conveyor initializtion failed");
    let my_rank = convey.my_rank();
    let num_ranks = convey.num_ranks();
    let quiet = matches.is_present("quiet") || my_rank > 0;
    let dump_files = matches.is_present("dump_files");

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
        println!("Running toposort on {} ranks", num_ranks);
        println!("Number of rows per rank     (-n) {}", numrows_per_rank);
        println!("Avg # non-eros per row      (-Z) {}", nz_per_row);
        println!("Erdos-Renyi edge probabilty (-e) {}", erdos_renyi_prob);

        println!("creating input matrix for toposort");
    }

    let mat = generate_toposort_input(numrows, erdos_renyi_prob, seed, dump_files);

    if !quiet {
        println!(
            "input matrix: {} rows, {} nonzeros",
            mat.numrows(),
            mat.nnz()
        );
    }

    if dump_files {
        mat.dump(20, "mat.out").expect("could not write mat.out");
        mat.write_mm_file("topo_mat.mm")
            .expect("could not write topo_mat.mm");
    }

    let tmat = mat.transpose();

    if dump_files {
        tmat.dump(20, "trans.out")
            .expect("could not write trans.out");
        tmat.write_mm_file("topo_tmat.mm")
            .expect("could not write topo_tmat.mm");
    }

    if !quiet {
        println!("Running toposort on mat (and tmat) ...");
    }

    let mut matret;
    for i in 0..2 {
        if i == 0 {
            if !quiet {
                print!("   using queue toposort: ");
            }
            matret = mat.toposort_queue(&tmat);
        } else {
            if !quiet {
                print!("   using queue2 toposort: ");
            }
            matret = mat.toposort_queue2(&tmat);
        }
        if !mat.check_result(&matret, dump_files) {
            println!("ERROR: check_result failed");
        }
        if !quiet {
            println!(" {} seconds", matret.laptime)
        }
    }
}
