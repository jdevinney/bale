#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Bale Serial deltastepping library
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
use sparsemat::wall_seconds;
use sparsemat::SparseMat;
use std::fs::OpenOptions;
use std::io::Error;
use std::io::Write;
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
    /// the resulgint distance vector
    pub distance: Vec<f64>,
    /// how big the source was
    pub source: usize,
    /// time to complete
    pub laptime: f64,
}

impl SsspInfo {
    /// Dump output distances to a file
    pub fn dump(&self, max_disp: usize, filename: &str) -> Result<(), Error> {
        //2 outs but perm has 1??
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
                // should use zip 0-0
                write!(file, "{}: {}\n", v, self.distance[v])?;
            }
        }
        Ok(())
    }

    /// Dump output distances to a file foo.wts 
    pub fn dump_wts(&self, filename: &str) -> Result<(), Error> {
        let path = Path::new(&filename);
        let mut file = OpenOptions::new().write(true).create(true).open(path)?;
        writeln!(file, "{}", self.distance.len())?;
        for v in &self.distance {
            writeln!(file, "{}", v)?;
        }
        Ok(())
    }
}

// A potential edge relaxation to be examined
struct Request {
    w: usize, // head of edge being relaxed
    dist: f64, // new distance from source to w using that edge
              //  v: usize,    // could include tail of edge (v,w) in request to build shortest path tree
}

// A struct and methods for all the data structures in delta stepping
// Each bucket is a circular doubly-linked list, linked by prev_elt and next_elt,
// indexed by an "elt" that is either a vertex in the bucket, or a "bucket header".
// A vertex's "elt" number is its vertex number.
// The header of bucket i (with reuse) is "elt" number num_buckets + i.

struct BucketSearcher<'a> {
    graph: &'a SparseMat,           // the graph being searched
    delta: f64,                     // width of a bucket
    num_buckets: usize,             // number of actual buckets, taking into account reuse
    tentative_dist: Vec<f64>,       // current tentative distance from source to this vtx
    prev_elt: Vec<usize>, // back link in list of vertices in each bucket (including bucket header)
    next_elt: Vec<usize>, // forward link in list of vertices in each bucket (including header)
    activated: Vec<bool>, // has this vtx ever been activated?
    vtx_bucket: Vec<Option<usize>>, // what bucket if any is this vtx in?
    bucket_size: Vec<usize>, // number of vertices in this bucket
    bucket_header: Vec<usize>, // which elt is this bucket's header? (just to make code clearer)

                          // iterator notes:
                          // BucketSearcher should also have an iterator that updates the active bucket, skipping to the
                          // next nonempty bucket (% num_buckets) and stopping when all buckets are empty.
                          // (But I'm not going to put active_bucket in the struct.)
                          // Also there should be an iterator that yields the vertices in a bucket,
                          // and keeps working when its current vertex is removed from the bucket.

                          // parallel notes:
                          // Vertices belong to PEs in round-robin order.
                          // Each PE has its own copy of every bucket, but only containing vertices it owns.
                          // (Say nv >= 100 * num_buckets * THREADS. Then storage for buckets is relatively small.)
                          // Arrays bucket_size and bucket_header are local to each PE.
                          // Arrays indexed only by vtx are split among PEs (but don't have to be shared, I think):
                          //     tentative_dist, activated, bucket.
                          // Arrays prev_elt and next_elt are also split, linking together only the vtxs on the local PE.
                          // The bucket-header nodes at the end of prev_elt and next_elt have copies on each PE.
                          // The parallelism all happens in relax_requests, where a request to relax an edge with head w
                          // gets conveyed to the PE that owns w.
}

// An iterator for BucketSearcher that finds the next nonempty bucket
impl<'a> IntoIterator for &'a BucketSearcher<'a> {
    type Item = usize;
    type IntoIter = ActiveBucketIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        ActiveBucketIterator {
            searcher: self,
            index: 0,
        }
    }
}

struct ActiveBucketIterator<'a> {
    searcher: &'a BucketSearcher<'a>,
    index: usize,
}

impl<'a> Iterator for ActiveBucketIterator<'a> {
    type Item = usize;
    fn next(&mut self) -> Option<usize> {
        let start_index = self.index;
        while self.searcher.bucket_size[self.index] == 0 {
            self.index = (self.index + 1) % self.searcher.num_buckets;
            if self.index == start_index {
                return None;
            }
        }
        Some(self.index)
    }
}

impl<'a> BucketSearcher<'a> {
    // Create a bucket structure for a weighted graph
    fn new(graph: &SparseMat, delta: f64) -> BucketSearcher<'_> {
        // I guess this uses lifetime elision
        let nv = graph.numrows();
        let mut max_edge_len: f64 = 0.0;
        if let Some(edge_len) = &graph.value {
            // 0-0 yuck! and max_edge_len = edge_len.max(); fails
            for c in edge_len {
                if *c > max_edge_len {
                    // 0-0 the compiler wants *c here, which I don't understand
                    max_edge_len = *c;
                }
            }
        } else {
            panic!("Graph must have edge weights (values)");
        }
        println!("max_edge_len = {}", max_edge_len);
        // upper bound on number of buckets we will ever need at the same time
        let num_buckets = (max_edge_len / delta).ceil() as usize + 1;
        // tentative distances all start out infinite, including source
        let tentative_dist = vec![f64::INFINITY; nv];
        // circular linked lists have room for the bucket headers,
        // and initially every list is empty (every element points to itself).
        let prev_elt: Vec<usize> = (0..nv + num_buckets).collect();
        let next_elt: Vec<usize> = (0..nv + num_buckets).collect();
        // initially no vtx has ever been activated.
        let activated = vec![false; nv];
        // initially every vtx is in no bucket.
        let vtx_bucket = vec![Option::<usize>::None; nv];
        // initially every bucket contains no vertices.
        let bucket_size = vec![0; num_buckets];
        // immutable bucket_header is just for clearer code
        let bucket_header: Vec<usize> = (nv..nv + num_buckets).collect();

        //barrier here

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

    // Dump bucket structure state to a file
    fn dump(
        &self,
        max_disp: usize,
        filename: &str,
        title: &str,
        nums: Vec<usize>,
    ) -> Result<(), Error> {
        let path = Path::new(&filename);
        let mut file = OpenOptions::new().append(true).create(true).open(path)?;
        let nv = self.graph.numrows();
        let now: DateTime<Local> = Local::now();
        writeln!(
            file,
            "=========================================================="
        )?;
        write!(file, "BucketSearcher: {}", title)?;
        for n in nums {
            write!(file, " {}", n)?;
        }
        writeln!(file, " at {}\n", now)?;
        writeln!(
            file,
            "nv={}, num_buckets={}, delta={}",
            nv, self.num_buckets, self.delta
        )?;
        writeln!(
            file,
            "elt: prev_elt next_elt vtx_bucket activated tentative_dist"
        )?;
        for r in display_ranges(max_disp, nv) {
            for v in r {
                // should be firstlist(r).zip(secondlist(r)) 0-0
                writeln!(
                    file,
                    "{}: {} {} {} {} {}",
                    v,
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
        for v in nv..nv + self.num_buckets {
            // should use zip 0-0
            writeln!(file, "{}: {} {}", v, self.prev_elt[v], self.next_elt[v])?;
        }
        writeln!(file, "bucket (bucket_size): elt elt ...")?;
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

    // find the next nonempty bucket after start_bucket, % num_buckets, if any
    fn next_nonempty_bucket(&self, start_bucket: usize) -> Option<usize> {
        let mut new_bucket = (start_bucket + 1) % self.num_buckets;
        while self.bucket_size[new_bucket] == 0 {
            new_bucket = (new_bucket + 1) % self.num_buckets;
            if new_bucket == start_bucket {
                return None;
            }
        }
        Some(new_bucket)
    }

    // what bucket does a vtx with this tentative distance go in?
    fn home_bucket(&self, dist: f64) -> usize {
        ((dist / self.delta).floor() as usize) % self.num_buckets
    }

    // remove a vtx from the bucket it's in. (harmless if not in a bucket)
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

    // insert a vtx into a bucket it's not in. (vtx must not be in a bucket)
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

    // make a list of all relaxation requests from light edges with tails in bucket
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
                            w: self.graph.nonzero[adj],
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

    // make a list of all relaxation requests from heavy edges with tails on vtx_list
    fn find_heavy_requests(&self, vtx_list: Vec<usize>) -> Vec<Request> {
        let mut requests: Vec<Request> = Vec::new();
        if let Some(edge_len) = &self.graph.value {
            for v in vtx_list {
                for adj in self.graph.offset[v]..self.graph.offset[v + 1] {
                    let vw_len = edge_len[adj];
                    if vw_len > self.delta {
                        // heavy edge
                        requests.push(Request {
                            w: self.graph.nonzero[adj],
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

    // return vertices in bucket that have not been activated(removed from an active bucket) before
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

    // remove all vertices from the bucket and mark them activated
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

    // relax all the requests from this phase (could be parallel)
    fn relax_requests(&mut self, requests: Vec<Request>) {
        // convey the request r=(w,d) to the PE that owns vtx w here, and have it call relax
        for r in requests {
            self.relax(r);
        }
        // barrier here (maybe also barrier at beginning of relax_requests?)
    }

    // relax an incoming edge to vtx r.w with new source distance r.dist, and rebucket r.w if necessary
    // this will be called by r.w's PE, so there is no race on tent[r.w]
    fn relax(&mut self, r: Request) {
        if r.dist < self.tentative_dist[r.w] {
            let new_bucket = self.home_bucket(r.dist);
            if let Some(old_bucket) = self.vtx_bucket[r.w] {
                // r.w was in a bucket,
                if old_bucket != new_bucket {
                    // move r.w from old to new bucket
                    self.remove_from_bucket(r.w);
                    self.place_in_bucket(r.w, new_bucket);
                }
            } else {
                // r.w was not in a bucket, put it in new bucket
                self.place_in_bucket(r.w, new_bucket);
            }
            self.tentative_dist[r.w] = r.dist;
        }
    }
}

/// Trait to extend sparsemat for delta_stepping
pub trait DeltaStepping {
    /// the implementation
    fn delta_stepping(&self, source: usize, forced_delta: Option<f64>) -> SsspInfo;
    /// check function
    fn check_result(&self, info: &SsspInfo, dump_files: bool) -> bool;
}

impl DeltaStepping for SparseMat {
    /// This implements the sequential AGI variant of delta stepping.
    /// # Argument: source vertex
    fn delta_stepping(&self, source: usize, forced_delta: Option<f64>) -> SsspInfo {
        assert!(self.numrows() == self.numcols());
        assert!(source < self.numrows());

        let t1 = wall_seconds().expect("wall second error");

        let (_mindeg, maxdeg, _sumdeg) = self.rowcounts().fold((self.numcols(), 0, 0), |acc, x| {
            (acc.0.min(x), acc.1.max(x), acc.2 + x)
        });

        // choose a value for delta, the bucket width
        let delta;
        if let Some(d) = forced_delta {
            delta = d;
        } else {
            delta = 1.0 / (maxdeg as f64);
        }
        println!(
            "delta_stepping: nvtxs = {}, nedges = {}, delta = {}",
            self.numrows(),
            self.offset[self.numrows()],
            delta
        );

        // initialize buckets, activated flags, etc.
        let mut searcher = BucketSearcher::new(&self, delta);

        // use relax to set tent(source) to 0, which also puts it in bucket 0
        searcher.relax(Request {
            w: source,
            dist: 0.0,
        });

        searcher
            .dump(20, "trace.out", "after relax source", vec![source])
            .expect("bucket dump failed");

        // outer loop: for each nonempty bucket in order ...
        let mut outer = 0;
        let mut active_bucket_if_any = Some(0);
        while let Some(active_bucket) = active_bucket_if_any {
            println!(
                "\nouter loop iteration {}: active_bucket = {}",
                outer, active_bucket
            );

            let mut removed: Vec<usize> = Vec::new(); // vertices removed from active bucket, R in paper

            // middle loop: while active bucket (B[i] in paper) is not empty ...
            let mut phase = 0;
            while searcher.bucket_size[active_bucket] > 0 {
                println!(
                    "middle loop iteration {}: active_bucket has {} vtxs",
                    phase, searcher.bucket_size[active_bucket]
                );

                // find light edges with tails in active bucket;
                // empty active bucket, keeping a set "removed" of unique vtxs removed from this bucket;
                // relax light edges, which may put some removed and other vtxs into the active bucket.

                let requests = searcher.find_light_requests(active_bucket);
                removed.append(&mut searcher.newly_active_vertices(active_bucket));
                searcher.empty_bucket(active_bucket);
                searcher.relax_requests(requests);

                searcher
                    .dump(20, "trace.out", "end of middle iter", vec![outer, phase])
                    .expect("bucket dump failed");
                phase += 1;
            } // end of middle looop

            // relax heavy edges with tails in removed set, which cannot add vtxs to active bucket
            let requests = searcher.find_heavy_requests(removed);
            searcher.relax_requests(requests);

            searcher
                .dump(20, "trace.out", "end of outer iter", vec![outer])
                .expect("bucket dump failed");
            outer += 1;

            active_bucket_if_any = searcher.next_nonempty_bucket(active_bucket);
        } // end of outer loop

        println!("\nDid {} iterations of outer loop.", outer);

        // return the info struct, which will now own the distance array
        SsspInfo {
            distance: searcher.tentative_dist,
            source: source,
            laptime: wall_seconds().expect("wall second error") - t1,
        }
    }

    /// check the result of delta stepping
    ///
    /// # Arguments
    /// * info data from the run to check
    /// * dump_files debugging flag
    fn check_result(&self, info: &SsspInfo, dump_files: bool) -> bool {
        println!(
            "\ncheck_result: source is {}, dump_files is {}",
            info.source, dump_files
        );
        if dump_files {
            info.dump(20, "dist.out").expect("info dump error");
        }
        let mut unreachable = 0;
        let mut max_dist: f64 = 0.0;
        let mut sum_dist: f64 = 0.0;
        for v in 0..self.numrows() {
            if info.distance[v].is_finite() {
                max_dist = f64::max(max_dist, info.distance[v]);
                sum_dist += info.distance[v];
            } else {
                unreachable += 1;
            }
        }
        println!(
            "unreachable vertices: {}; max finite distance: {}; avg finite distance: {}",
            unreachable,
            max_dist,
            sum_dist / ((self.numrows() - unreachable) as f64)
        );
        true
    }
}
