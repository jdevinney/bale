#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Bale Serial Toposort library
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For license information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use sparsemat::wall_seconds;
use sparsemat::Perm;
use sparsemat::SparseMat;

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
    let tmat = SparseMat::erdos_renyi_tri(numrows, edge_prob, false, true, seed);

    if dump_files {
        tmat.dump(20, "orig_tri.out")
            .expect("could not write orig_tri.out");
    }

    assert_eq!(tmat.is_upper_triangular(), true);

    // get random row and column permutations
    let rperminv = Perm::random(numcols, 1234);
    let cperminv = Perm::random(numcols, 5678);
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

/// return values used by main
#[derive(Debug, Clone)]
pub struct TopoInfo {
    rperm: Perm,
    cperm: Perm,
    /// the time to complete
    pub laptime: f64,
}

impl TopoInfo {
    fn new(len: usize) -> TopoInfo {
        let rperm = Perm::new(len);
        let cperm = Perm::new(len);
        TopoInfo {
            rperm,
            cperm,
            laptime: 0.0,
        }
    }
}

/// a trait to extend SparseMat for toposort
pub trait TopoSort {
    /// queue implementation
    fn toposort_queue(&self, tmat: &SparseMat) -> TopoInfo;
    /// loop implemetation
    fn toposort_loop(&self, tmat: &SparseMat) -> TopoInfo;
    /// check results
    fn check_result(&self, info: &TopoInfo, dump_files: bool) -> bool;
}

impl TopoSort for SparseMat {
    /// This routine implements the agi variant of toposort
    /// # Arguments
    /// * tmat the transpose of mat
    fn toposort_queue(&self, tmat: &SparseMat) -> TopoInfo {
        let nr = self.numrows();
        let nc = self.numcols();
        let mut ret = TopoInfo::new(nr);

        let mut queue: Vec<usize> = vec![0; nr];
        let mut rowtrck: Vec<i64> = vec![0; nr];
        let mut start = 0;
        let mut end = 0;

        for row in 0..nr {
            rowtrck[row] = 0;
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                rowtrck[row] += (1 << 32) + *nz as i64;
            }
            if rowtrck[row] >> 32 == 1 {
                queue[end] = row;
                end += 1;
            }
        }

        // we a pick a row with a single nonzero = col.
        // setting rperm[pos] = row and cprem[pos] = col
        // in a sense 'moves' the row and col to the bottom
        // right corner of matrix.
        // Next, we cross out that row and col by decrementing
        //  the rowcnt for any row that contains that col
        // repeat
        //
        let t1 = wall_seconds().expect("wall second error");

        let mut n_pivots = 0;
        while start < end {
            let row = queue[start];
            start += 1;
            let col = (rowtrck[row] & 0xFFFF) as usize; // see cool trick

            //println!("row {} col {} queue {:?} rowtrck {:x?}", row, col, queue, rowtrck );
            ret.rperm.perm()[row] = nr - 1 - n_pivots;
            ret.cperm.perm()[col] = nc - 1 - n_pivots;
            n_pivots += 1;

            // look at this column (tmat's row) to find all the rows that hit it
            for nz in &tmat.nonzero[tmat.offset[col]..tmat.offset[col + 1]] {
                let trow = *nz;
                assert!(trow < self.numrows());
                rowtrck[trow] -= (1 << 32) + col as i64;
                if (rowtrck[trow] >> 32) == 1 {
                    queue[end] = trow;
                    end += 1;
                }
            }
        }

        ret.laptime = wall_seconds().expect("wall second error") - t1;

        assert!(n_pivots == nr);
        ret
    }

    /// This routine implements the agi variant of toposort
    /// # Arguments
    /// * tmat the transpose of mat
    fn toposort_loop(&self, tmat: &SparseMat) -> TopoInfo {
        let nr = self.numrows();
        let nc = self.numcols();
        let mut ret = TopoInfo::new(nr);

        let mut rowtrck: Vec<i64> = vec![0; nr];

        for row in 0..nr {
            rowtrck[row] = 0;
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                rowtrck[row] += (1 << 32) + *nz as i64;
            }
        }

        // we a pick a row with a single nonzero = col.
        // setting rperm[pos] = row and cprem[pos] = col
        // in a sense 'moves' the row and col to the bottom
        // right corner of matrix.
        // Next, we cross out that row and col by decrementing
        //  the rowcnt for any row that contains that col
        // repeat
        //
        let t1 = wall_seconds().expect("wall second error");

        let mut n_pivots = 0;
        while n_pivots < nr {
            for row in 0..nr {
                if (rowtrck[row] >> 32) == 1 {
                    let col = (rowtrck[row] & 0xFFFF) as usize; // see cool trick
                    ret.rperm.perm()[row] = nr - 1 - n_pivots;
                    ret.cperm.perm()[col] = nc - 1 - n_pivots;
                    n_pivots += 1;
                    rowtrck[row] = 0;

                    // look at this column (tmat's row) to find all the rows that hit it
                    for nz in &tmat.nonzero[tmat.offset[col]..tmat.offset[col + 1]] {
                        let trow = *nz;
                        assert!((trow) < self.numrows());
                        rowtrck[trow] -= (1 << 32) + col as i64;
                    }
                }
            }
        }

        ret.laptime = wall_seconds().expect("wall second error") - t1;

        assert!(n_pivots == nr);
        ret
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
        if !mat2.is_upper_triangular() {
            println!("permutation is not upper triangular");
            return false;
        }
        if dump_files {
            mat2.dump(20, "mat2.out").expect("couldn't write mat2.out");
        }
        true
    }
}
