#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Deltastepping applications
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For license information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use chrono::{DateTime, Local};
use clap::{App, Arg};
use convey_hpc::Convey;
use delta_stepping::DeltaStepping;
use spmat::SparseMat;

// Application that finds shortest path lengths from a single source
// in a directed graph, using the Meyer/Sanders delta-stepping
// algorithm.

// Find single-source shortest path lengths in a directed graph.
//
// First we generate the problem by making a random directed graph
// with edge costs c(v,w).
//
// Each vertex has a "tentative distance" during the algorithm. The
// source has tentative distance 0, and all other vertices initially
// have have tentative distance inf. We proceed by "relaxing" edges
// (in a clever order): relaxing edge (v,w) changes the tentative
// distance of w to min(tent(w), tent(v) + c(v,w)). Eventually each
// vertex's tent() distance becomes final, or "settled"; initially
// only the source is settled.
//
// Unsettled vertices with tent() < inf are kept in "buckets" by
// tent() value; bucket i contains vertices with tent() at least
// i*\Delta and less than (i+1)*\Delta, where \Delta is a parameter.
//
// The algorithm has three nested loops.
//
// The outer (serial) loop is over buckets; an iteration processes
// vertices in the lowest nonempty bucket until it is empty.
//
// The middle (serial) loop is over "phases"; a phase consists of
// removing all the vertices in the active bucket and relaxing all the
// "light" edges out of them (an edge is "light" if it has cost at
// most \Delta, "heavy" otherwise). The edge relaxations in a phase
// may cause vertices to enter the active bucket; phases continue
// until the bucket is empty.  At that point all the vertices that
// were in that bucket are settled.  Following the light-edge phases,
// one more phase relaxes all the heavy edges from vertices deleted
// from the active bucket; this cannot cause any vertices to go into
// that bucket.
//
// The inner (parallel) loop implements the edge relaxations in a
// single phase.  Those relaxations can be done in any order, provided
// they are done atomically.

fn main() {
    // figure out parallel environment
    let convey = Convey::new().expect("conveyor initializtion failed");
    let my_rank = convey.my_rank();
    let num_ranks = convey.num_ranks();

    // parse the command line arguments
    let matches = App::new("DeltaStepping")
        .version("0.1.0")
        .about("Implements a test of DeltaStepping")
        .arg(
            Arg::with_name("numrows")
                .short("n")
                .long("numrows")
                .takes_value(true)
                .help("The number of rows in test matrix"),
        )
        .arg(
            Arg::with_name("source_vtx")
                .short("s")
                .long("source_vtx")
                .takes_value(true)
                .help("The number of the source vertex"),
        )
        .arg(
            Arg::with_name("er_prob")
                .short("e")
                .long("erdos_renyi_prob")
                .takes_value(true)
                .help("Probability of an edge in the erdos renyi graph"),
        )
        .arg(
            Arg::with_name("input_file")
                .short("i")
                .long("input_file")
                .takes_value(true)
                .help("Matrix Market input file for graph with edge weights"),
        )
        .arg(
            Arg::with_name("forced_delta")
                .short("f")
                .long("forced_delta")
                .takes_value(true)
                .help("Bucket width for delta-stepping, override algorithm's choice"),
        )
        .arg(
            Arg::with_name("dump_files")
                .short("d")
                .long("dump_files")
                .takes_value(false)
                .help("Write the matrix to sssp.mm and the output distances to sssp.dst"),
        )
        .arg(
            Arg::with_name("quiet")
                .short("q")
                .long("quiet")
                .takes_value(false)
                .help("produce less chatty output"),
        )
        .arg(
            Arg::with_name("trace")
                .short("t")
                .long("trace")
                .takes_value(false)
                .help("write data structure trace file from each thread"),
        )
        .get_matches();

    // input args
    let numrows: usize = matches
        .value_of("numrows")
        .unwrap_or("10")
        .parse()
        .expect("numrows: not an integer");
    let source_vtx: usize = matches
        .value_of("source_vtx")
        .unwrap_or("0")
        .parse()
        .expect("source_vtx: not an integer");
    let erdos_renyi_prob: f64 = matches
        .value_of("er_prob")
        .unwrap_or("0.3")
        .parse()
        .expect("er_prob: not a float");
    let input_file: &str = matches.value_of("input_file").unwrap_or("NONE");
    let forced_delta: f64 = matches
        .value_of("forced_delta")
        .unwrap_or("0.0")
        .parse()
        .expect("forced_delta: not a float");
    let seed = 12346; // the random-number seed is actually never used
    let quiet = matches.is_present("quiet") || my_rank > 0;
    let trace = matches.is_present("trace");
    let dump_files = matches.is_present("dump_files");

    // done with options, now do it

    let mut mat: SparseMat;
    if matches.is_present("input_file") {
        mat = SparseMat::read_mm_file(input_file).expect("can't read MatrixMarket file");
    } else {
        let mode = 3; // mode 3 means nonsymmetric matrix, directed graph (not acyclic)
        mat = SparseMat::gen_erdos_renyi_graph(numrows, erdos_renyi_prob, false, mode, seed);
        mat.randomize_values();
    };

    if !quiet {
        let now: DateTime<Local> = Local::now();
        println!(
            "Running delta_stepping on {} from source_vtx {} using {} PEs at {}",
            if matches.is_present("input_file") {
                input_file
            } else {
                "random matrix"
            },
            source_vtx,
            num_ranks,
            now
        );
        println!("\nInput matrix stats:");
        mat.stats();
        println!("");
    }
    if dump_files {
        mat.write_mm_file("sssp.mm")
            .expect("could not write sssp.mm");
    }

    let matret = mat.delta_stepping(
        source_vtx,
        if forced_delta == 0.0 {
            None
        } else {
            Some(forced_delta)
        },
        quiet,
        trace,
    );

    if dump_files {
        matret.write_dst("sssp.dst").expect("results write error");
    }
    let checks = mat.check_result(&matret, input_file, quiet);

    if !quiet {
        println!("\nResult of check_result is {}", checks);
    }
}
