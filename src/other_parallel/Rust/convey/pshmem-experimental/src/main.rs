use pshmem::object::Fetchers;
use pshmem::Pshmem;

fn main() {
    let pshmem = Pshmem::new();
    let me = pshmem.rank();
    let size = pshmem.size();

    // A shared object with 1 usize element per rank
    let x = pshmem.new_object::<usize>(1).expect("allocation error");
    // the local slice of this data
    let lp = x.local_part();

    lp[0] = pshmem.rank() + 42;

    pshmem.barrier();
    println!(
        "Hello Pshmem World from rank {} of {}, a nice number is {}",
        me,
        pshmem.size(),
        x.fetch((me + 1) % size).expect("fetch failed")
    );
}
