//! Program to do histogram in conveyors
use clap::{App, Arg};
use convey_hpc::Convey;
use rand::distributions::{Distribution, Uniform};

fn main() {
    let matches = App::new("ig")
        .version("0.1.0")
        .about("test of index gather")
        .arg(
            Arg::with_name("requests")
                .short("r")
                .long("requests")
                .takes_value(true)
                .help("then number of requests"),
        )
        .arg(
            Arg::with_name("table_entries")
                .short("t")
                .long("table_entries")
                .takes_value(true)
                .help("the size of the table"),
        )
        .get_matches();

    let requests: u64 = matches
        .value_of("requests")
        .unwrap_or("4000")
        .parse()
        .expect("bad requests arg");
    let table_entries: usize = matches
        .value_of("table_entries")
        .unwrap_or("40000")
        .parse()
        .expect("bad table_size arg");

    let convey = Convey::new().expect("shmem initializtion failed");
    do_ig_convey(&convey, table_entries, requests);
    do_ig_convey_simple(&convey, table_entries, requests);
}

fn do_ig_convey(convey: &Convey, table_entries: usize, requests: u64) {
    let me = convey.my_rank();
    let num = convey.num_ranks();
    println!("Hello, world from rank {} of {}!", me, num);

    let mut rng = rand::thread_rng();
    let die = Uniform::from(0..table_entries);

    // create our data structure before the block, so we can access later
    let ltable: Vec<i64> = vec![42 + me as i64; (table_entries + num - 1) / num];
    let mut tgt: Vec<i64> = vec![0; ((requests as usize + num - 1) / num).max(10)];
    let mut total_requests: i64 = 0;
    let mut total_returns: i64 = 0;

    {
        // Always put the session in a new block, as you will not be
        // able to able to local content after conveyor is done
        //  Could also use an explicit drop
        let mut session2 = Convey::begin(|item: (usize, i64), _from_rank| {
            tgt[item.0] = item.1;
            total_returns += 1;
        });

        //println!("{:?}", session2);
        {
            let mut session1 = Convey::begin(|item: (usize, usize), from_rank| {
                session2.push((item.0, ltable[item.1]), from_rank);
                total_requests += 1;
            });

            //println!("{:?}", session1);
            // do the updates
            for i in 0..requests / num as u64 {
                let index = die.sample(&mut rng);
                let rank = index % num;
                let offset = index / num;
                session1.push((i as usize, offset), rank);
            }
            session1.finish();
        }
        session2.finish();
    }
    println!(
        "pe{}, {} requests, {} returns first 10 tgt {:?}",
        me,
        total_requests,
        total_returns,
        &tgt[0..10]
    );
}

fn do_ig_convey_simple(convey: &Convey, table_entries: usize, requests: u64) {
    let me = convey.my_rank();
    let num = convey.num_ranks();
    println!("Hello, world from rank {} of {}!", me, num);

    let mut rng = rand::thread_rng();
    let die = Uniform::from(0..table_entries);

    // create our data structure before the block, so we can access later
    let ltable: Vec<i64> = vec![42 + me as i64; (table_entries + num - 1) / num];
    let mut tgt: Vec<i64> = vec![0; ((requests as usize + num - 1) / num).max(10)];
    let mut total_requests: i64 = 0;
    let mut total_returns: i64 = 0;

    Convey::simple_return(
        (0..(requests / num as u64)).map(|x| {
            let (offset, rank) = convey.offset_rank(die.sample(&mut rng));
            ((x as usize, offset), rank)
        }),
        |item: usize| {
            total_requests += 1;
            ltable[item]
        },
        |index: usize, val: i64| {
            tgt[index] = val;
            total_returns += 1;
        },
    );

    println!(
        "pe{}, {} requests, {} returns first 10 tgt {:?}",
        me,
        total_requests,
        total_returns,
        &tgt[0..10]
    );
}
