//! Triangle counting application
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
use convey_hpc::collect::CollectValues;
use convey_hpc::Convey;
use triangle::Triangle;

fn main() {
    let matches = App::new("Triangle")
        .version("0.1.0")
        .about("Implements a test of Triangle Counting")
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
            Arg::with_name("gen_kron")
                .short("K")
                .long("gen_kron")
                .takes_value(true)
                .help("Kroniker product graph size"),
        )
        .arg(
            Arg::with_name("alg")
                .short("a")
                .long("alg")
                .takes_value(true)
                .help("algorithm variant, 0 is L & L * U, 1 is L & U * L"),
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
    let alg: u16 = matches
        .value_of("alg")
        .unwrap_or("0")
        .parse()
        .expect("alg: not an integer");
    let gen_kron = matches.is_present("gen_kron");
    let gen_kron_str = matches.value_of("gen_kron").unwrap_or("");

    let seed = 12346;
    let convey = Convey::new().expect("Conveyor initialzation failed");
    let my_rank = convey.my_rank();
    let num_ranks = convey.num_ranks();
    let quiet = matches.is_present("quiet") || my_rank > 0;
    let _dump_files = matches.is_present("dump_files");

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
        println!("Running triangle on {} ranks", num_ranks);
        println!("Number of rows per rank     (-n) {}", numrows_per_rank);
        println!("Algorithm                   (-a) {}", alg);
        if gen_kron {
            println!("Avg # non-zeros per row     (-Z) {}", nz_per_row);
            println!("Erdos-Renyi edge probabilty (-e) {}", erdos_renyi_prob);
        } else {
            println!("Generating Kronicker Prod   (-K) {}", gen_kron_str);
        }
        println!("creating input matrix for triangle");
    }

    let (mat, correct_answer) = if gen_kron {
        let mut iter = gen_kron_str.split(" ");
        let mode: u16 = iter.next().unwrap_or("0").parse().expect("bad mode arg");
        let mut args: Vec<u16> = Vec::new();
        for item in iter {
            args.push(item.parse().expect("bad kron value"));
        }
        // Build up the kronecker mode
        let ans = match mode {
            0 => 0.0,
            1 => {
                let mut ans = 1.0;
                for item in &args {
                    ans *= 3 as f64 * *item as f64 + 1.0;
                }
                ans *= 1.0 / 6.0;
                let mut x = 1.0;
                for item in &args {
                    x *= *item as f64 + 1.0;
                }
                ans - 0.5 * x + 1.0 / 3.0
            }
            2 => {
                (1.0 / 6.0) * 4.0f64.powf(args.len() as f64) - 2.0f64.powf((args.len() - 1) as f64)
                    + 1.0 / 3.0
            }
            _ => panic!("unsupported kronecker mode"),
        };

        let half = args.len() / 2;
        (
            triangle::generate_kronecker_graph(&args[..half], &args[half..], mode),
            ans.round(),
        )
    } else {
        (
            spmat::SparseMat::gen_erdos_renyi_graph(numrows, erdos_renyi_prob, false, 1, seed),
            0.0,
        )
    };
    if !quiet {
        println!(
            "input matrix: {} rows, {} nonzeros, expected count: {}",
            mat.numrows(),
            mat.nnz(),
            correct_answer
        );
    }

    let tmat = mat.transpose();
    let upper = if alg == 1 { Some(&tmat) } else { None };

    if !quiet {
        println!("Running triangle on mat (and tmat) ...");
    }

    let mut triret;
    for i in 0..1 {
        if i == 0 {
            if !quiet {
                print!("   using triangle_push: ");
            }
            triret = mat.triangle_push(upper);
        } else {
            if !quiet {
                print!("   using triangle_pull: ");
            }
            triret = mat.triangle_pull(upper);
        }
        let tot_shref = triret.sh_refs.reduce_sum();
        let tot_tri_cnt = triret.tri_cnt.reduce_sum();
        if !quiet {
            println!(
                " {:.4} seconds, {} Shrefs, {} Triangles",
                triret.laptime, tot_shref, tot_tri_cnt
            )
        }
    }
}
