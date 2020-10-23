#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Bale Serial deltastepping applicaton
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
use delta_stepping::DeltaStepping;
use itertools::join;
use regex::Regex;
use sparsemat::SparseMat;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

// Application that finds shortest path lengths from a single source
// in a directed graph, using the Meyer/Sanders delta-stepping
// algorithm.  This serial Rust implementation is based on the serial
// C in bale, and is intended as a first step toward a Rust/Conveyors
// version.

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
// The inner (potentially parallel) loop implements the edge
// relaxations in a single phase.  Those relaxations can be done in
// any order, provided they are done atomically.

fn main() {
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
            Arg::with_name("source")
                .short("s")
                .long("source")
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

    // input args
    let numrows: usize = matches
        .value_of("numrows")
        .unwrap_or("10")
        .parse()
        .expect("numrows: not an integer");
    let source: usize = matches
        .value_of("source")
        .unwrap_or("2")
        .parse()
        .expect("source: not an integer");
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
    let quiet = matches.is_present("quiet");
    let dump_files = matches.is_present("dump_files");

    let mut mat: SparseMat;
    if matches.is_present("input_file") {
        mat = SparseMat::read_mm_file(input_file).expect("can't read MatrixMarket file");
    } else {
        mat = SparseMat::erdos_renyi_graph(numrows, erdos_renyi_prob, false, seed);
        mat.randomize_values();
    };

    if !quiet {
        println!("input matrix stats:");
        mat.stats();
    }

    if dump_files {
        mat.dump(20, "mat.out").expect("could not write mat.out");
    }

    mat.write_mm_file("sssp_mat.mm")
        .expect("could not write sssp_mat.mm");

    if !quiet {
        let now: DateTime<Local> = Local::now();
        println!(
            "Running delta_stepping on {} from source {} at {} ...",
            if matches.is_present("input_file") {
                input_file
            } else {
                "random matrix"
            },
            source,
            now
        );
    }

    let matret = mat.delta_stepping(
        source,
        if forced_delta == 0.0 {
            None
        } else {
            Some(forced_delta)
        },
    );
    if dump_files {
        // matret.dump(0, "results.out").expect("results write error");
        matret.dump_wts("results.wts").expect("results write error");
    }

    // hack to check against an output file if it's there
    // this would have been like 6 lines in Python ... just sayin' ...
    if matches.is_present("input_file") {
        let re = Regex::new(r"\.").unwrap();
        let mut tokens: Vec<&str> = re.split(input_file).collect();
        if let Some(_) = tokens.pop() {
            tokens.push("wts");
        }
        let check_file = join(&tokens, ".");
        if let Ok(fp) = File::open(&check_file) {
            let reader = BufReader::new(fp);
            let mut check_dist: Vec<f64> = vec![];
            for line in reader.lines() {
                let d = line
                    .expect("can't read check file")
                    .parse::<f64>()
                    .expect("can't read check file");
                check_dist.push(d);
            }
            let check_nv = check_dist[0] as usize;
            if (check_nv != mat.numrows()) || (check_nv != check_dist.len() - 1) {
                println!(
                    "check file problem: mat.numrows = {}, check says nv = {}, check has {} distances",
                    mat.numrows(), check_nv, check_dist.len()-1
                );
            } else {
                let mut diff: f64 = 0.0;
                let mut csum: f64 = 0.0;
                for v in 0..check_nv {
                    if check_dist[v + 1].is_finite() || matret.distance[v].is_finite() {
                        diff += (check_dist[v + 1] - matret.distance[v]).powi(2);
                    }
                    if check_dist[v + 1].is_finite() {
                        csum += check_dist[v + 1].powi(2);
                    }
                }
                if diff <= f64::EPSILON.sqrt() * csum {
                    println!("\nCORRECT! relative diff = {}", diff / csum);
                } else {
                    println!("\nDISAGREE! relative diff = {}", diff / csum);
                }
            }
        } else {
            println!("No ground truth file '{}' to check against", check_file);
        }
    } else {
        println!("No ground truth file to check against");
    }

    if !mat.check_result(&matret, dump_files) {
        println!("ERROR: check_result failed");
    }
    if !quiet {
        println!(" {} seconds", matret.laptime)
    }
}
