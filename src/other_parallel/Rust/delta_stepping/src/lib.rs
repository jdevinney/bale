#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Deltastepping
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
use convey_hpc::collect::CollectValues;
use convey_hpc::Convey;
use itertools::join;
use regex::Regex;
use serde::{Deserialize, Serialize};
use spmat::wall_seconds;
use spmat::SparseMat;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::BufRead;
use std::io::BufReader;
use std::io::{Error, Write};
use std::ops::Range;
use std::path::Path;

/// A helper function for dumping only part of a big data structure.
/// This should really go somewhere else than the delta_stepper lib.
pub fn display_ranges(max_disp: usize, num_items: usize) -> Vec<Range<usize>> {
    let mut ranges: Vec<Range<usize>> = Vec::new();
    if max_disp <= num_items && max_disp > 0 {
        ranges.push(0..max_disp / 2);
        ranges.push(num_items - max_disp / 2..num_items);
    } else {
        ranges.push(0..num_items);
    }
    ranges
}

/// Output structure for single-source shortest path
#[derive(Debug, Clone)]
pub struct SsspInfo {
    /// source size
    pub source: usize,
    /// time to complete
    pub laptime: f64,
    /// distance vector
    pub distance: Vec<f64>,
}

impl SsspInfo {
    /// Dump output distances to a file
    pub fn dump(&self, max_disp: usize, filename: &str) -> Result<(), Error> {
        let path = Path::new(&filename);
        let mut file = OpenOptions::new().write(true).create(true).open(path)?;
        writeln!(
            file,
            "=========================================================="
        )?;
        let now: DateTime<Local> = Local::now();
        writeln!(file, "Final Distances at {}", now)?;

        write!(file, "vtx: dist\n")?;
        for r in display_ranges(max_disp, self.distance.len()) {
            for v in r {
                write!(file, "{}: {}\n", v, self.distance[v])?;
            }
        }
        Ok(())
    }

    /// Write output distances to a file, collecting them all on rank 0
    pub fn write_dst(&self, filename: &str) -> Result<(), Error> {
        let convey = Convey::new().expect("conveyor initializtion failed");
        let my_rank = convey.my_rank();
        let num_ranks = convey.num_ranks();
        let num_vtxs = self.distance.len().reduce_sum();
        let mut all_distances: Vec<f64> = vec![0.0; if my_rank == 0 { num_vtxs } else { 0 }];
        #[derive(Copy, Clone, Debug, Serialize, Deserialize)]
        struct Item {
            vtx: usize,
            dst: f64,
        }
        {
            let mut session = Convey::begin(|item: Item, _from_rank| {
                all_distances[item.vtx] = item.dst;
            });
            for (i, d) in self.distance.iter().enumerate() {
                session.push(
                    Item {
                        vtx: num_ranks * i + my_rank, // should be session.global_index(i)
                        dst: *d,
                    },
                    0,
                );
            }
            session.finish();
        }
        if my_rank == 0 {
            let path = Path::new(&filename);
            let mut file = OpenOptions::new().write(true).create(true).open(path)?;
            writeln!(file, "{}", num_vtxs)?;
            for v in all_distances {
                writeln!(file, "{}", v)?;
            }
        }
        Ok(())
    }
}

/// A potential edge relaxation to be examined
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
struct Request {
    w_g: usize, // head of edge being relaxed, global vertex index
    dist: f64,  // new distance from source to w_g using that edge
                //  v_g: usize,  // could include tail of edge (v_g,w_g) in request to build shortest path tree
}

/// A struct and methods for all the data structures in delta stepping
/// Each bucket is a circular doubly-linked list, linked by prev_elt and next_elt,
/// indexed by an "elt" that is either a vertex in the bucket, or a "bucket header".
/// A vertex's "elt" number is its vertex number, locally indexed on this rank.
/// The header of bucket i (with reuse) is "elt" number nvtxs_this_rank + i.

struct BucketSearcher<'a> {
    graph: &'a SparseMat,           // the graph being searched
    delta: f64,                     // width of a bucket
    num_buckets: usize,             // number of actual buckets, taking into account reuse
    tentative_dist: Vec<f64>,       // current tentative distance from source to this vtx
    prev_elt: Vec<usize>, // back link in list of vertices in each bucket (including bucket header)
    next_elt: Vec<usize>, // forward link in list of vertices in each bucket (including header)
    activated: Vec<bool>, // has this vtx ever been activated?
    vtx_bucket: Vec<Option<usize>>, // what bucket if any is this vtx in?
    bucket_size: Vec<usize>, // number of vertices (from this rank) in this bucket
    bucket_header: Vec<usize>, // which elt is this bucket's header? (just to make code clearer)

                          // parallel notes:
                          // Vertices belong to PEs in round-robin order.
                          // Each PE has its own copy of every bucket, but only containing vertices it owns.
                          //     (Say nv >= 100 * num_buckets * THREADS. Then storage for buckets is relatively small.)
                          // There is no shared memory; every array is either local or split and locally indexed.
                          // Arrays bucket_size and bucket_header are local to each PE.
                          // Arrays indexed only by vtx are split among PEs, with local indexing, but are not shared:
                          //     tentative_dist, activated, bucket.
                          // Arrays prev_elt and next_elt are also split, linking together only the vtxs on the local PE.
                          //     The bucket-header nodes at the end of prev_elt and next_elt have copies on each PE.
                          // The parallelism all happens in relax_requests, where a request to relax an edge with head w_g
                          //     gets conveyed to the PE that owns w_g.
}

impl<'a> BucketSearcher<'a> {
    /// Create a bucket structure for a weighted graph
    fn new(graph: &SparseMat, delta: f64) -> BucketSearcher<'_> {
        let nvtxs_this_rank = graph.numrows_this_rank();

        let mut my_max_edge_len: f64 = 0.0;
        if let Some(edge_len) = &graph.value {
            for c in edge_len {
                if *c > my_max_edge_len {
                    my_max_edge_len = *c;
                }
            }
        } else {
            panic!("Graph must have edge weights (values)");
        }
        let max_edge_len = my_max_edge_len.reduce_max();
        // upper bound on number of buckets we will ever need at the same time
        let num_buckets = (max_edge_len / delta).ceil() as usize + 1;
        // tentative distances all start out infinite, including source
        let tentative_dist = vec![f64::INFINITY; nvtxs_this_rank];
        // circular linked lists have room for the bucket headers,
        // and initially every list is empty (every element points to itself).
        let prev_elt: Vec<usize> = (0..nvtxs_this_rank + num_buckets).collect();
        let next_elt: Vec<usize> = (0..nvtxs_this_rank + num_buckets).collect();
        // initially no vtx has ever been activated.
        let activated = vec![false; nvtxs_this_rank];
        // initially every vtx is in no bucket.
        let vtx_bucket = vec![Option::<usize>::None; nvtxs_this_rank];
        // initially every bucket contains no vertices.
        let bucket_size = vec![0; num_buckets];
        // immutable bucket_header is just for clearer code
        let bucket_header: Vec<usize> = (nvtxs_this_rank..nvtxs_this_rank + num_buckets).collect();

        graph.barrier();

        BucketSearcher {
            graph,
            delta,
            num_buckets,
            tentative_dist,
            prev_elt,
            next_elt,
            activated,
            vtx_bucket,
            bucket_size,
            bucket_header,
        }
    }

    /// Append bucket structure state to the file trace.pe#.out (which will get big)
    fn trace(&self, max_disp: usize, title: &str, nums: Vec<usize>) -> Result<(), Error> {
        let filename = format!("trace.{}.out", self.graph.my_rank());
        let path = Path::new(&filename);
        let mut file = OpenOptions::new().append(true).create(true).open(path)?;
        let nvtxs_this_rank = self.graph.numrows_this_rank();
        let now: DateTime<Local> = Local::now();
        writeln!(
            file,
            "=========================================================="
        )?;
        write!(
            file,
            "BucketSearcher: rank {} of {}: {}",
            self.graph.my_rank(),
            self.graph.num_ranks(),
            title
        )?;
        for n in nums {
            write!(file, " {}", n)?;
        }
        writeln!(file, " at {}\n", now)?;
        writeln!(
            file,
            "nvtxs_this_rank={}, num_buckets={}, delta={}",
            nvtxs_this_rank, self.num_buckets, self.delta
        )?;
        writeln!(
            file,
            "elt: global_vtx prev_elt next_elt vtx_bucket activated tentative_dist"
        )?;
        for r in display_ranges(max_disp, nvtxs_this_rank) {
            for v in r {
                writeln!(
                    file,
                    "{}: {} {} {} {} {} {}",
                    v,
                    self.graph.global_index(v),
                    self.prev_elt[v],
                    self.next_elt[v],
                    if let Some(b) = self.vtx_bucket[v] {
                        b.to_string()
                    } else {
                        "N".to_string()
                    },
                    self.activated[v],
                    self.tentative_dist[v]
                )?;
            }
        }
        for v in nvtxs_this_rank..nvtxs_this_rank + self.num_buckets {
            writeln!(file, "{}: {} {}", v, self.prev_elt[v], self.next_elt[v])?;
        }
        writeln!(file, "bucket (bucket_size on this rank): elt elt ...")?;
        for r in display_ranges(max_disp, self.num_buckets) {
            for b in r {
                write!(file, "{} ({}):", b, self.bucket_size[b])?;
                let mut v = self.next_elt[self.bucket_header[b]];
                while v != self.bucket_header[b] {
                    write!(file, " {}", v)?;
                    v = self.next_elt[v];
                }
                write!(file, "\n")?;
            }
        }
        writeln!(file, " ")?;
        Ok(())
    }

    /// total size of a bucket over all ranks
    fn global_bucket_size(&self, bucket: usize) -> usize {
        let ret = self.bucket_size[bucket].reduce_sum();
        ret
    }

    /// find the next nonempty bucket after start_bucket, % num_buckets, if any
    fn next_nonempty_bucket(&self, start_bucket: usize) -> Option<usize> {
        let mut my_steps = 1;
        while my_steps < self.num_buckets
            && self.bucket_size[(start_bucket + my_steps) % self.num_buckets] == 0
        {
            my_steps += 1;
        }
        let steps = my_steps.reduce_min();
        if steps < self.num_buckets {
            Some((start_bucket + steps) % self.num_buckets)
        } else {
            None
        }
    }

    /// what bucket does a vtx with this tentative distance go in?
    fn home_bucket(&self, dist: f64) -> usize {
        ((dist / self.delta).floor() as usize) % self.num_buckets
    }

    /// remove a vtx from the bucket it's in. (harmless if not in a bucket). w is a local vtx index.
    fn remove_from_bucket(&mut self, w: usize) {
        if let Some(b) = self.vtx_bucket[w] {
            self.prev_elt[self.next_elt[w]] = self.prev_elt[w];
            self.next_elt[self.prev_elt[w]] = self.next_elt[w];
            self.prev_elt[w] = w;
            self.next_elt[w] = w;
            self.vtx_bucket[w] = None;
            self.bucket_size[b] -= 1;
        }
    }

    /// insert a vtx into a bucket it's not in. (vtx must not be in a bucket) w is a local vtx index.
    fn place_in_bucket(&mut self, w: usize, new_bucket: usize) {
        assert!(self.vtx_bucket[w] == None);
        assert!(self.prev_elt[w] == w);
        assert!(self.next_elt[w] == w);
        self.prev_elt[w] = self.bucket_header[new_bucket];
        self.next_elt[w] = self.next_elt[self.bucket_header[new_bucket]];
        self.prev_elt[self.next_elt[w]] = w;
        self.next_elt[self.prev_elt[w]] = w;
        self.vtx_bucket[w] = Some(new_bucket);
        self.bucket_size[new_bucket] += 1;
    }

    /// make a list of all relaxation requests from light edges with tails in bucket b
    fn find_light_requests(&self, b: usize) -> Vec<Request> {
        let mut requests: Vec<Request> = Vec::new();
        if let Some(edge_len) = &self.graph.value {
            let mut v = self.next_elt[self.bucket_header[b]];
            while v != self.bucket_header[b] {
                for adj in self.graph.offset[v]..self.graph.offset[v + 1] {
                    let vw_len = edge_len[adj];
                    if vw_len <= self.delta {
                        // light edge
                        requests.push(Request {
                            w_g: self.graph.nonzero[adj], // nonzero[] contains global indices
                            dist: self.tentative_dist[v] + vw_len,
                        });
                    }
                }
                v = self.next_elt[v];
            }
        } else {
            panic!("Graph must have edge weights (values)");
        }
        requests
    }

    /// make a list of all relaxation requests from heavy edges with tails on vtx_list
    fn find_heavy_requests(&self, vtx_list: Vec<usize>) -> Vec<Request> {
        let mut requests: Vec<Request> = Vec::new();
        if let Some(edge_len) = &self.graph.value {
            for v in vtx_list {
                for adj in self.graph.offset[v]..self.graph.offset[v + 1] {
                    let vw_len = edge_len[adj];
                    if vw_len > self.delta {
                        // heavy edge
                        requests.push(Request {
                            w_g: self.graph.nonzero[adj], // nonzero[] contains global indices
                            dist: self.tentative_dist[v] + vw_len,
                        });
                    }
                }
            }
        } else {
            panic!("Graph must have edge weights (values)");
        }
        requests
    }

    /// return vertices in bucket that have not ever been activated (removed from an active bucket)
    fn newly_active_vertices(&self, b: usize) -> Vec<usize> {
        let mut new_vtxs: Vec<usize> = Vec::new();
        let mut v = self.next_elt[self.bucket_header[b]];
        while v != self.bucket_header[b] {
            if !self.activated[v] {
                new_vtxs.push(v);
            }
            v = self.next_elt[v];
        }
        new_vtxs
    }

    /// remove all vertices from the bucket and mark them activated
    fn empty_bucket(&mut self, b: usize) {
        let header = self.bucket_header[b];
        let mut v = self.next_elt[header];
        while v != header {
            let w = self.next_elt[v];
            self.next_elt[v] = v;
            self.prev_elt[v] = v;
            self.vtx_bucket[v] = None;
            self.activated[v] = true;
            v = w;
        }
        self.next_elt[header] = header;
        self.prev_elt[header] = header;
        self.bucket_size[b] = 0;
    }

    /// relax all the requests from this phase, in parallel
    fn relax_requests(&mut self, requests: Vec<Request>) {
        // convey the request r=(w_g,d) to the PE that owns vtx w_g, and have it call relax()
        let mut session = Convey::begin(|item: Request, _from_rank| {
            self.relax(item);
        });
        for r in requests {
            let rank = session.offset_rank(r.w_g).1;
            session.push(r, rank);
        }
        session.finish();
    }

    /// relax an incoming edge to vtx r.w_g with new source distance r.dist,
    ///     and rebucket r.w_g if necessary.
    /// this functions like a user-defined atomic: since it is called by r.w_g's PE,
    ///     there is no race on tentative_dist[r.w_g] .
    fn relax(&mut self, r: Request) {
        let w = self.graph.local_index(r.w_g); // panics if argument is not on this rank
        if r.dist < self.tentative_dist[w] {
            let new_bucket = self.home_bucket(r.dist);
            if let Some(old_bucket) = self.vtx_bucket[w] {
                // w was in a bucket,
                if old_bucket != new_bucket {
                    // move w from old to new bucket
                    self.remove_from_bucket(w);
                    self.place_in_bucket(w, new_bucket);
                }
            } else {
                // w was not in a bucket, put it in new bucket
                self.place_in_bucket(w, new_bucket);
            }
            self.tentative_dist[w] = r.dist;
        }
    }
}

/// Trait to implement deltastepping extensions to SparseMat
pub trait DeltaStepping {
    /// Deltastepping function
    fn delta_stepping(
        &self,
        source: usize,
        forced_delta: Option<f64>,
        quiet: bool,
        trace: bool,
    ) -> SsspInfo;
    /// check results
    fn check_result(&self, info: &SsspInfo, input_file: &str, quiet: bool) -> bool;
    /// compare vectors, treating infinity carefully
    fn sq_rel_diff(&self, a: &Vec<f64>, b: &Vec<f64>) -> f64;
}

impl DeltaStepping for SparseMat {
    /// This implements parallel delta stepping.
    /// # Arguments: source vertex, optional overriding delta value
    fn delta_stepping(
        &self,
        source: usize,
        forced_delta: Option<f64>,
        quiet: bool,
        trace: bool,
    ) -> SsspInfo {
        assert!(self.numrows() == self.numcols());
        assert!(source < self.numrows());

        let t1 = wall_seconds();

        // choose a value for delta, the bucket width
        let delta;
        if let Some(d) = forced_delta {
            delta = d;
        } else {
            let (_my_mindeg, my_maxdeg, _my_sumdeg) =
                self.rowcounts().fold((self.numcols(), 0, 0), |acc, x| {
                    (acc.0.min(x), acc.1.max(x), acc.2 + x)
                });
            let maxdeg = my_maxdeg.reduce_max();
            delta = 1.0 / (maxdeg as f64);
        }
        if !quiet {
            println!(
                "delta_stepping: nvtxs = {}, nedges = {}, delta = {}",
                self.numrows(),
                self.nnz(),
                delta
            );
        }

        // initialize buckets, activated flags, etc.
        let mut searcher = BucketSearcher::new(&self, delta);

        // use relax to set tentative_dist[source] to 0, which also puts it in bucket 0
        if self.offset_rank(source).1 == self.my_rank() {
            searcher.relax(Request {
                w_g: source,
                dist: 0.0,
            });
        }

        if trace {
            searcher
                .trace(20, "after relax source", vec![source])
                .expect("bucket trace write failed");
        }

        // outer loop: for each nonempty bucket in order ...
        let mut outer = 0;
        let mut active_bucket_if_any = Some(0);
        while let Some(active_bucket) = active_bucket_if_any {
            if !quiet {
                println!(
                    "\nouter loop iteration {}: active_bucket = {}",
                    outer, active_bucket
                );
            }

            let mut removed: Vec<usize> = Vec::new(); // vertices removed from active bucket, R in paper

            // middle loop: while active bucket (B[i] in paper) is not empty ...
            let mut phase = 0;
            while searcher.global_bucket_size(active_bucket) > 0 {
                if !quiet {
                    println!(
                        "middle loop iteration {}: active_bucket {} has {} vtxs on rank 0",
                        phase, active_bucket, searcher.bucket_size[active_bucket]
                    );
                }

                // find light edges with tails in active bucket;
                // empty active bucket, keeping a set "removed" of unique vtxs removed from this bucket;
                // relax light edges, which may put some removed and other vtxs into the active bucket.

                let requests = searcher.find_light_requests(active_bucket);
                removed.append(&mut searcher.newly_active_vertices(active_bucket));
                searcher.empty_bucket(active_bucket);
                searcher.relax_requests(requests);

                if trace {
                    searcher
                        .trace(20, "end of middle iter", vec![outer, phase])
                        .expect("bucket trace write failed");
                }
                phase += 1;
            } // end of middle looop

            // relax heavy edges with tails in removed set, which cannot add vtxs to active bucket
            let requests = searcher.find_heavy_requests(removed);
            searcher.relax_requests(requests);

            if trace {
                searcher
                    .trace(20, "end of outer iter", vec![outer])
                    .expect("bucket trace write failed");
            }
            outer += 1;

            active_bucket_if_any = searcher.next_nonempty_bucket(active_bucket);
        } // end of outer loop

        if !quiet {
            println!("\nDid {} iterations of outer loop.", outer);
        }

        // return the info struct, which will now own the distance array
        SsspInfo {
            distance: searcher.tentative_dist,
            source: source,
            laptime: wall_seconds() - t1,
        }
    }

    /// Square of relative 2-norm difference between two distributed vectors, with inf-inf treated as 0.
    /// This should be standard somewhere, but I can't find it.
    fn sq_rel_diff(&self, a: &Vec<f64>, b: &Vec<f64>) -> f64 {
        assert!(a.len() == b.len());
        let mut l_diff: f64 = 0.0;
        let mut l_csum: f64 = 0.0;
        for v in 0..a.len() {
            if a[v].is_finite() || b[v].is_finite() {
                l_diff += (a[v] - b[v]).powi(2);
            }
            if a[v].is_finite() {
                l_csum += a[v].powi(2);
            }
        }
        let diff = l_diff.reduce_sum();
        let csum = l_csum.reduce_sum();

        if diff == 0.0 {
            0.0
        } else {
            diff / csum
        }
    }

    /// check the result of delta stepping
    ///
    /// # Arguments
    /// * info data from the run to check
    /// * input file name
    /// * quiet flag
    fn check_result(&self, info: &SsspInfo, input_file: &str, quiet: bool) -> bool {
        if !quiet {
            println!(
                "\ncheck_result: source is {}, rank 0 time is {}",
                info.source, info.laptime
            );
        }
        let mut l_unreachable = 0;
        let mut l_max_dist: f64 = 0.0;
        let mut l_sum_dist: f64 = 0.0;
        for v in 0..self.numrows_this_rank() {
            if info.distance[v].is_finite() {
                l_max_dist = f64::max(l_max_dist, info.distance[v]);
                l_sum_dist += info.distance[v];
            } else {
                l_unreachable += 1;
            }
        }
        let unreachable = l_unreachable.reduce_sum();
        let max_dist: f64 = l_max_dist.reduce_max();
        let sum_dist = l_sum_dist.reduce_sum();
        if !quiet {
            println!(
                "unreachable vertices: {}; max finite distance: {}; avg finite distance: {}",
                unreachable,
                max_dist,
                sum_dist / (self.numrows() as f64 - unreachable as f64)
            );
        }
        // check against ground truth file if it's there
        if input_file == "NONE" {
            if !quiet {
                println!("No ground truth file to check against");
            }
        } else {
            let re = Regex::new(r"\.").unwrap();
            let mut tokens: Vec<&str> = re.split(input_file).collect();
            if let Some(_) = tokens.pop() {
                tokens.push("dst");
            }
            let check_file = join(&tokens, ".");
            if let Ok(fp) = File::open(&check_file) {
                let reader = BufReader::new(fp);
                let mut check_dst: Vec<f64> = vec![];
                for (index, line) in reader.lines().enumerate() {
                    let d = line
                        .expect("can't read check file")
                        .parse::<f64>()
                        .expect("can't read check file");
                    if index > 0 && self.my_rank() == self.offset_rank(index - 1).1 {
                        check_dst.push(d);
                    }
                }
                let diff = self.sq_rel_diff(&check_dst, &info.distance);
                if diff <= f64::EPSILON.sqrt() {
                    if !quiet {
                        println!("\nCORRECT! squared relative diff = {}", diff);
                    }
                } else {
                    if !quiet {
                        println!("\nDISAGREE! squared relative diff = {}", diff);
                    }
                    return false;
                }
            } else {
                println!("No ground truth file '{}' to check against", check_file);
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::SparseMat;
    use crate::DeltaStepping;
    use convey_hpc::testing_support::TestingMutex;

    /*
        #[test]
        #[should_panic]
        fn delta_stepping_rand1() {
            let mutex = TestingMutex::new();
            let numrows = 100;
            let erdos_renyi_prob = 0.05;
            let unit_diag = false;
            let mode = 3;
            let seed = 1;
            let mat = SparseMat::gen_erdos_renyi_graph(numrows, erdos_renyi_prob, unit_diag, mode, seed);
            let source = 0;
            let forced_delta = None;
            let quiet = true;
            let trace = false;
            let _matret = mat.delta_stepping(source, forced_delta, quiet, trace);
            drop(mutex);
        }
    */
    #[test]
    fn delta_stepping_rand2() {
        let mutex = TestingMutex::new();
        let numrows = 100;
        let erdos_renyi_prob = 0.05;
        let unit_diag = false;
        let mode = 3;
        let seed = 1;
        let mut mat =
            SparseMat::gen_erdos_renyi_graph(numrows, erdos_renyi_prob, unit_diag, mode, seed);
        mat.randomize_values();
        let source = 0;
        let forced_delta = None;
        let quiet = true;
        let trace = false;
        let matret_a = mat.delta_stepping(source, forced_delta, quiet, trace);
        for delta in vec![0.01, 0.1, 1.0, 100.0] {
            let matret_b = mat.delta_stepping(source, Some(delta), quiet, trace);
            let diff = mat.sq_rel_diff(&matret_a.distance, &matret_b.distance);
            assert!(diff <= f64::EPSILON.sqrt());
        }
        drop(mutex);
    }
    #[test]
    fn delta_stepping_sparse() {
        let mutex = TestingMutex::new();
        let input_file = "../../../../example_matrices/sparse100.mm";
        let mat = SparseMat::read_mm_file(input_file).expect("can't read input file");
        let source = 2;
        let forced_delta = None;
        let quiet = true;
        let trace = false;
        let matret_a = mat.delta_stepping(source, forced_delta, quiet, trace);
        assert!(mat.check_result(&matret_a, input_file, quiet));
        for delta in vec![0.01, 0.1, 1.0, 100.0] {
            let matret_b = mat.delta_stepping(source, Some(delta), quiet, trace);
            let diff = mat.sq_rel_diff(&matret_a.distance, &matret_b.distance);
            assert!(diff <= f64::EPSILON.sqrt());
        }
        drop(mutex);
    }
}
