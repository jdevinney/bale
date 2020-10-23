#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Toposort library
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For license information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use convey_hpc::collect::CollectValues;
use convey_hpc::{Convey, RcVec};
use serde::{Deserialize, Serialize};
use spmat::perm::Perm;
use spmat::wall_seconds;
use spmat::SparseMat;

/// Generate a matrix that is the a random permutation of a sparse uppper triangular matrix.
///
/// # Arguments
/// * numrows the number of rows (and columns) in the produced matrix
/// * erdos_renyi_prob the probability that an entry in the matrix is non-zero
/// * seed the seed for random number generator that determines the original
///   matrix and the permutations
/// * dump_files is a debugging flag
///
/// Make the upper triangular matrix is the upper half of an Erdös–Renyi random graph
/// We force the diagonal entry and pick the other edges with probability erdos_renyi_prob
/// then randomly permute the rows and the columns.  
/// The toposort algorithm takes this matrix and finds one of the possibly many row and
/// column permutations that would bring the matrix back to an upper triangular form.
pub fn generate_toposort_input(
    numrows: usize,
    edge_prob: f64,
    seed: i64,
    dump_files: bool,
) -> SparseMat {
    let numcols = numrows;
    let tmat = SparseMat::gen_erdos_renyi_graph(numrows, edge_prob, true, 2, seed);

    if dump_files {
        tmat.dump(20, "orig_tri.out")
            .expect("could not write orig_tri.out");
    }

    assert_eq!(tmat.is_upper_triangular(true), true);

    // get random row and column permutations
    let rperminv = Perm::random(numcols, 1234);
    let cperminv = Perm::random(numcols, 5678);

    assert_eq!(cperminv.is_perm(), true);
    assert_eq!(rperminv.is_perm(), true);

    if dump_files {
        rperminv
            .dump(20, "rperm.out")
            .expect("could not write rperm.out");
        rperminv
            .dump(20, "cperm.out")
            .expect("could not write cperm.out");
    }

    let mat = tmat.permute(&rperminv, &cperminv);

    if dump_files {
        mat.dump(20, "perm.out").expect("could not write perm.out");
    }
    mat
}

/// Information returned from toposort, the two permutations and time
#[derive(Debug)]
pub struct TopoInfo {
    rperm: Perm,
    cperm: Perm,
    /// time to complete
    pub laptime: f64,
}

impl TopoInfo {
    fn new(len: usize) -> TopoInfo {
        let rperm = Perm::identity(len);
        let cperm = Perm::identity(len);
        TopoInfo {
            rperm,
            cperm,
            laptime: 0.0,
        }
    }
}

/// a trait for Toposort extensions to SparseMat
pub trait TopoSort {
    /// first version of toposort algorithm
    fn toposort_queue(&self, tmat: &SparseMat) -> TopoInfo;
    /// second version of toposort algorithm
    fn toposort_queue2(&self, tmat: &SparseMat) -> TopoInfo;
    /// third version of toposort algorithm
    fn toposort_loop(&self, tmat: &SparseMat) -> TopoInfo;
    /// checker
    fn check_result(&self, info: &TopoInfo, dump_files: bool) -> bool;
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
enum TopoSend {
    ToColq {
        offset: usize,
        lev: usize,
    },
    MarkRow {
        row: usize,
        global_col: usize,
        lev: usize,
    },
}

impl TopoSort for SparseMat {
    /// This routine implements the agi variant of toposort
    /// # Arguments
    /// * tmat the transpose of mat
    fn toposort_queue(&self, tmat: &SparseMat) -> TopoInfo {
        let nr = self.per_my_rank(self.numrows());
        let nc = self.per_my_rank(self.numcols());
        let mut ret = TopoInfo::new(self.numrows());

        // row queue just has row
        let row_queue: RcVec<usize> = RcVec::new();
        // column queue carries level
        let col_queue: RcVec<(usize, usize)> = RcVec::new();

        // sum of all the nonzeros in a row
        let nz_sum: RcVec<usize> = RcVec::new();
        // count of all the nonzeros in a row
        let mut nz_cnt: Vec<usize> = vec![0; nr];

        let level: RcVec<usize> = RcVec::new();
        let mut matched_col: Vec<usize> = vec![0; nr];

        let mut num_levels = 0;

        for row in 0..nr {
            nz_sum.push(0);
            level.push(0);
            nz_cnt[row] = self.offset[row + 1] - self.offset[row];
            if nz_cnt[row] == 1 {
                row_queue.push(row);
            }
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                let t = nz_sum.get_at(row);
                nz_sum.set_at(row, t + *nz);
            }
        }

        let check = nz_cnt.iter().fold(0, |acc, x| acc + x);
        assert_eq!(check, self.nonzero.len());

        let t1 = wall_seconds();

        let mut r_and_c_done = 0;
        {
            let mut session = Convey::begin(|item: TopoSend, _from_rank| {
                match item {
                    TopoSend::ToColq { offset, lev } => {
                        // we have revieved a transposed row, push it on the
                        // column queue with the associated level
                        //println!("off {} lev {}", offset, lev);
                        col_queue.push((offset, lev));
                    }
                    TopoSend::MarkRow {
                        row,
                        global_col,
                        lev,
                    } => {
                        // we have recieved a row entry from the
                        // transpose_matrix update the row_sum for the
                        // row in the non_transpose matrix to indicate
                        // we have seen it
                        //println!("row {} gcol {} lev {}", row, global_col, level);
                        let t = nz_sum.get_at(row);
                        nz_sum.set_at(row, t - global_col);
                        nz_cnt[row] -= 1;
                        // if we are on a new level, update the level count
                        if lev >= level.get_at(row) {
                            level.set_at(row, lev + 1);
                            if lev + 1 > num_levels {
                                num_levels = lev + 1;
                            }
                        }
                        // if we have received all the elements for
                        // this row but one, we are now eligible to go
                        // on the row queue
                        if nz_cnt[row] == 1 {
                            row_queue.push(row);
                        }
                    }
                }
            });

            session.always_send = true;

            while session.advance(r_and_c_done == nr + nc) {
                // for each degree 1 row (which is from the row queue)
                // we send the column for that row to the rank which
                // holds the transpose of the column so it can add it
                // to the column queue for that rank
                while let Some(row) = row_queue.pop() {
                    let transpose_row = nz_sum.get_at(row);
                    let lev = level.get_at(row);
                    matched_col[row] = transpose_row;
                    let (offset, rank) = self.offset_rank(transpose_row);
                    session.push(TopoSend::ToColq { offset, lev }, rank);
                    r_and_c_done += 1;
                }

                // for each pivot in the column queue, we look up the
                // related nonzeros in the transpos matrix and, for
                // each of them, send it to the rank which holds that
                // row in the non-transpose matrix so we can subtract
                // the elements off the row_sum
                while let Some((col, lev)) = col_queue.pop() {
                    for nz in &tmat.nonzero[tmat.offset[col]..tmat.offset[col + 1]] {
                        let global_row = *nz;
                        let global_col = col * self.num_ranks() + self.my_rank();
                        let (row, rank) = self.offset_rank(global_row);
                        session.push(
                            TopoSend::MarkRow {
                                row,
                                global_col,
                                lev,
                            },
                            rank,
                        );
                    }
                    r_and_c_done += 1;
                }
                if false {
                    println!(
                        "rank {} tsq1b rq {} cq {} {} {}",
                        self.my_rank(),
                        row_queue.len(),
                        col_queue.len(),
                        r_and_c_done,
                        nr + nc
                    );
                    if wall_seconds() - t1 > 10.0 {
                        session.debug(true);
                    }
                }
            }
        }

        //println!("tsq1");
        num_levels += 1;
        num_levels = num_levels.reduce_max();

        let mut level_sizes: Vec<usize> = vec![0; num_levels];
        let mut level_start: Vec<usize> = vec![0; num_levels];

        for row in 0..nr {
            level_sizes[level.get_at(row)] += 1;
        }

        let mut total = 0;
        for i in 0..num_levels {
            let size = level_sizes[i];
            level_start[i] = total + size.reduce_prior_sum();
            level_sizes[i] = size.reduce_sum();
            //println!("start[{}] {}, sizes {}", i, level_start[i], level_sizes[i]);
            total += level_sizes[i];
        }

        // we need a copy of rp to avoid borrow problems in conveyor that follows
        let mut rp = Perm::identity(self.numrows());
        for i in 0..nr {
            let lev = level.get_at(i);
            assert!(lev <= num_levels);
            rp.perm()[i] = (self.numrows() - 1) - level_start[lev];
            level_start[lev] += 1;
        }
        //println!("tsq2");

        Convey::session()
            .pull_fn(|item: (usize, usize), _from_rank| {
                ret.cperm.perm()[item.0] = item.1;
            })
            .push_iter((0..nr).map(|i| {
                let pos = rp.entry(i);
                let (col, rank) = self.offset_rank(matched_col[i]);
                ((col, pos), rank)
            }))
            .finish();
        ret.rperm = rp;
        ret.laptime = wall_seconds() - t1;
        //println!("tsq3");
        ret
    }

    /// This routine implements the agi variant of toposort
    /// # Arguments
    /// * tmat the transpose of mat
    fn toposort_queue2(&self, tmat: &SparseMat) -> TopoInfo {
        let nr = self.per_my_rank(self.numrows());
        let nc = self.per_my_rank(self.numcols());
        let mut ret = TopoInfo::new(self.numrows());

        // row queue just has row
        let row_queue: RcVec<usize> = RcVec::new();
        // column queue carries level
        let col_queue: RcVec<(usize, usize)> = RcVec::new();

        // sum of all the nonzeros in a row
        let nz_sum: RcVec<usize> = RcVec::new();
        // count of all the nonzeros in a row
        let mut nz_cnt: Vec<usize> = vec![0; nr];

        let level: RcVec<usize> = RcVec::new();
        let mut matched_col: Vec<usize> = vec![0; nr];

        let mut num_levels = 0;

        for row in 0..nr {
            nz_sum.push(0);
            level.push(0);
            nz_cnt[row] = self.offset[row + 1] - self.offset[row];
            if nz_cnt[row] == 1 {
                row_queue.push(row);
            }
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                let t = nz_sum.get_at(row);
                nz_sum.set_at(row, t + *nz);
            }
        }

        let t1 = wall_seconds();

        let mut r_and_c_done = 0;
        {
            let mut to_colq = Convey::begin(|item: (usize, usize), _from_rank| {
                // we have revieved a transposed row, push it on the
                // column queue with the associated level
                //println!("off {} lev {}", offset, lev);
                col_queue.push((item.0, item.1));
            });

            to_colq.always_send = true;

            {
                let mut mark_row = Convey::begin(|item: (usize, usize, usize), _from_rank| {
                    let row = item.0;
                    let global_col = item.1;
                    let lev = item.2;
                    // we have recieved a row entry from the
                    // transpose_matrix update the row_sum for the
                    // row in the non_transpose matrix to indicate
                    // we have seen it
                    //println!("row {} gcol {} lev {}", row, global_col, level);
                    let t = nz_sum.get_at(row);
                    nz_sum.set_at(row, t - global_col);
                    nz_cnt[row] -= 1;
                    // if we are on a new level, update the level count
                    if lev >= level.get_at(row) {
                        level.set_at(row, lev + 1);
                        if lev + 1 > num_levels {
                            num_levels = lev + 1;
                        }
                    }
                    // if we have received all the elements for
                    // this row but one, we are now eligible to go
                    // on the row queue
                    if nz_cnt[row] == 1 {
                        row_queue.push(row);
                    }
                });

                mark_row.always_send = true;
                loop {
                    let a1 = to_colq.advance(r_and_c_done == nr + nc);
                    let a2 = mark_row.advance(r_and_c_done == nr + nc);
                    if !a1 && !a2 {
                        break;
                    }

                    // for each degree 1 row (which is from the row queue)
                    // we send the column for that row to the rank which
                    // holds the transpose of the column so it can add it
                    // to the column queue for that rank
                    while let Some(row) = row_queue.pop() {
                        let transpose_row = nz_sum.get_at(row);
                        let lev = level.get_at(row);
                        matched_col[row] = transpose_row;
                        let (offset, rank) = self.offset_rank(transpose_row);
                        to_colq.push((offset, lev), rank);
                        r_and_c_done += 1;
                    }
                    //println!("tsq1a {}", r_and_c_done);

                    // for each pivot in the column queue, we look up the
                    // related nonzeros in the transpos matrix and, for
                    // each of them, send it to the rank which holds that
                    // row in the non-transpose matrix so we can subtract
                    // the elements off the row_sum
                    while let Some((col, lev)) = col_queue.pop() {
                        for nz in &tmat.nonzero[tmat.offset[col]..tmat.offset[col + 1]] {
                            let global_row = *nz;
                            let global_col = col * self.num_ranks() + self.my_rank();
                            let (row, rank) = self.offset_rank(global_row);
                            mark_row.push((row, global_col, lev), rank);
                        }
                        r_and_c_done += 1;
                    }
                    //println!("tsq1b {}", r_and_c_done);
                }
            }
        }

        num_levels += 1;
        num_levels = num_levels.reduce_max();

        let mut level_sizes: Vec<usize> = vec![0; num_levels];
        let mut level_start: Vec<usize> = vec![0; num_levels];

        for row in 0..nr {
            level_sizes[level.get_at(row)] += 1;
        }

        let mut total = 0;
        for i in 0..num_levels {
            let size = level_sizes[i];
            level_start[i] = total + size.reduce_prior_sum();
            level_sizes[i] = size.reduce_sum();
            //println!("start[{}] {}, sizes {}", i, level_start[i], level_sizes[i]);
            total += level_sizes[i];
        }

        // we need a copy of rp to avoid borrow problems in conveyor that follows
        let mut rp = Perm::identity(self.numrows());
        for i in 0..nr {
            let lev = level.get_at(i);
            assert!(lev <= num_levels);
            rp.perm()[i] = (self.numrows() - 1) - level_start[lev];
            level_start[lev] += 1;
        }

        Convey::session()
            .pull_fn(|item: (usize, usize), _from_rank| {
                ret.cperm.perm()[item.0] = item.1;
            })
            .push_iter((0..nr).map(|i| {
                let pos = rp.entry(i);
                let (col, rank) = self.offset_rank(matched_col[i]);
                ((col, pos), rank)
            }))
            .finish();
        ret.rperm = rp;
        ret.laptime = wall_seconds() - t1;
        ret
    }

    /// This routine implements the agi variant of toposort
    /// # Arguments
    /// * tmat the transpose of mat
    fn toposort_loop(&self, _tmat: &SparseMat) -> TopoInfo {
        todo!();
    }

    /// check the result toposort
    ///
    ///  check that the permutations are in fact permutations and the check that applying
    ///  them to the original matrix yields an upper triangular matrix
    /// # Arguments
    /// * info data from the run to check
    /// * dump_files debugging flag
    fn check_result(&self, info: &TopoInfo, dump_files: bool) -> bool {
        if !info.cperm.is_perm() {
            println!("cperm not a permutation");
            return false;
        }

        if !info.rperm.is_perm() {
            println!("rperm not a permutation");
            return false;
        }

        let mat2 = self.permute(&info.rperm, &info.cperm);
        if !mat2.is_upper_triangular(false) {
            println!("permutation is not upper triangular");
            return false;
        }
        if dump_files {
            mat2.dump(20, "mat2.out").expect("io error");
        }
        true
    }
}
