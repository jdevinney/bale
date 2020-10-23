use clap::{App, Arg};
use rand::distributions::{Distribution, Uniform};
use shmem::atomic::{Atomic, GlobalAtomic};
use shmem::Shmem;

fn main() {
    let matches = App::new("histo")
        .version("0.1.0")
        .about("test of histogram")
        .arg(
            Arg::with_name("buckets")
                .short("b")
                .long("buckets")
                .takes_value(true)
                .help("then number of buckets"),
        )
        .arg(
            Arg::with_name("updates")
                .short("u")
                .long("updates")
                .takes_value(true)
                .help("the number of updates"),
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

    let shmem = Shmem::new().expect("shmem initializtion failed");
    do_histo_shmem(&shmem, buckets, updates);
    do_histo_shmem_global(&shmem, buckets, updates);
}

fn do_histo_shmem(shmem: &Shmem, buckets: usize, updates: u64) {
    let me = shmem.my_pe();
    let num = shmem.n_pes();
    //println!("Hello, world from PE {} of {}!", me, num);

    let mut rng = rand::thread_rng();
    let die = Uniform::from(0..buckets as usize);

    let histo = shmem
        .new_object::<i64>(buckets / num)
        .expect("allocation error");

    for _i in 0..updates {
        let index = die.sample(&mut rng);
        let pe = index % num;
        let offset = index / num;
        histo.atomic_inc(offset, pe).expect("failed atomic inc");
    }
    shmem.barrier();
    println!("pe{}, first bucket {}", me, histo.local_part()[0]);
}

fn do_histo_shmem_global(shmem: &Shmem, buckets: usize, updates: u64) {
    let me = shmem.my_pe();
    //let num = shmem.n_pes();
    //println!("Hello, world from PE {} of {}!", me, num);

    let mut rng = rand::thread_rng();
    let die = Uniform::from(0..buckets as usize);

    let histo = shmem
        .new_global_object::<i64>(buckets, 0)
        .expect("allocation error");

    for _i in 0..updates {
        let index = die.sample(&mut rng);
        histo.atomic_inc(index).expect("failed atomic inc");
    }
    shmem.barrier();
    println!("pe{}, first bucket {}", me, histo.local_part()[0]);
}
