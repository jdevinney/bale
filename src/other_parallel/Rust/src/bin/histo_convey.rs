#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]
//! Program to do histogram in conveyors
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
use rand::distributions::{Distribution, Uniform};
use std::time::Instant;

fn main() {
    let matches = App::new("histo")
        .version("0.1.0")
        .about("test of histogram")
        .arg(
            Arg::with_name("buckets")
                .short("b")
                .long("buckets")
                .takes_value(true)
                .help("then number of buckets per rank"),
        )
        .arg(
            Arg::with_name("updates")
                .short("u")
                .long("updates")
                .takes_value(true)
                .help("the number of updates per rank"),
        )
        .arg(
            Arg::with_name("ranks_per_node")
                .short("rpn")
                .long("ranks_per_node")
                .takes_value(true)
                .help("The number of ranks per node, used to scale injection bandwidth"),
        )
        .arg(
            Arg::with_name("verbose")
                .short("-v")
                .long("verbose")
                .takes_value(false)
                .help("increase the amount of verbosity"),
        )
        .get_matches();

    let buckets: usize = matches
        .value_of("buckets")
        .unwrap_or("4000")
        .parse()
        .expect("bad buckets arg");
    let updates: u64 = matches
        .value_of("updates")
        .unwrap_or("40000")
        .parse()
        .expect("bad updates arg");
    let mut ranks_per_node: usize = matches
        .value_of("ranks_per_node")
        .unwrap_or("32")
        .parse()
        .expect("bad ranks_per_node arg");
    let verbose: u64 = matches.occurrences_of("verbose");

    let convey = Convey::new().expect("conveyor initializtion failed");

    let num_ranks = convey.num_ranks();
    if ranks_per_node > num_ranks {
        ranks_per_node = num_ranks;
    }

    do_histo_convey(&convey, buckets, updates, ranks_per_node, verbose);
    do_histo_convey_simple(&convey, buckets, updates, ranks_per_node, verbose);
}

/// the main work function using flexible conveyors
fn do_histo_convey(
    convey: &Convey,
    buckets: usize,
    updates: u64,
    ranks_per_node: usize,
    verbose: u64,
) {
    let me = convey.my_rank();
    let num = convey.num_ranks();
    if verbose > 0 {
        println!("Hello, world from rank {} of {}!", me, num);
    }

    // a random number generator, uniform over the total number of
    // buckets
    let mut rng = rand::thread_rng();
    let die = Uniform::from(0..buckets * num as usize);

    // create our data structure before the block, so we can access later
    let mut local: Vec<i64> = vec![0; buckets.max(5)];
    let mut total_updates: u64 = 0;

    let now = Instant::now();
    {
        // Always put the session in a new block, as you will
        // not be able to able to local after conveyor is done
        let mut session = Convey::begin(|item: usize, _from_rank| {
            //if from_rank != me {
            //   println!("to {} from {} received item {}", me, from_rank, item);
            //}
            local[item] += 1;
            total_updates += 1;
        });

        //println!("{:?}", session);
        // do the updates
        for _i in 0..updates {
            let index = die.sample(&mut rng);
            let rank = index % num;
            let offset = index / num;
            //println!("to {} from {}  sent item {}", rank, me, offset);
            session.push(offset, rank);
        }
        session.finish();
    }
    let d = now.elapsed();
    if verbose > 0 || me == 0 {
        println!(
            "pe{}, {} updates, {} msec, {} mb/node-sec first 5 buckets {:?}",
            me,
            total_updates,
            d.as_millis(),
            (total_updates * 8 * ranks_per_node as u64) / d.as_micros() as u64,
            &local[0..5]
        );
    }
}

/// the main work function using simple conveyors
fn do_histo_convey_simple(
    convey: &Convey,
    buckets: usize,
    updates: u64,
    ranks_per_node: usize,
    verbose: u64,
) {
    let me = convey.my_rank();
    let num = convey.num_ranks();
    if verbose > 0 {
        println!("Hello, world from rank {} of {}!", me, num);
    }

    let mut rng = rand::thread_rng();
    let die = Uniform::from(0..buckets * num as usize);

    // create our data structure before the block, so we can access later
    let mut local: Vec<i64> = vec![0; buckets.max(5)];
    let mut total_updates: u64 = 0;

    let now = Instant::now();
    Convey::simple(
        (0..updates).map(|_x| convey.offset_rank(die.sample(&mut rng))),
        |item: usize, _from_rank| {
            local[item] += 1;
            total_updates += 1;
        },
    );

    let d = now.elapsed();
    if verbose > 0 || me == 0 {
        println!(
            "pe{}, {} updates, {} msec, {} mb/node-sec first 5 buckets {:?}",
            me,
            total_updates,
            d.as_millis(),
            (total_updates * 8 * ranks_per_node as u64) / d.as_micros() as u64,
            &local[0..5]
        );
    }
}
