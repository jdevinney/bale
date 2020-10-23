#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Triangle library
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For license information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use convey_hpc::Convey;
use spmat::wall_seconds;
use spmat::SparseMat;

/// Generate a distributed graph that is the product of a collection
/// of star graphs. This is done in two stages. In the first stage,
/// the list of stars (parameterized by an integer m, K_{1,m}) is
/// split in half and each half-list forms a local adjacency matrix
/// for the Kronecker product of the stars (matrices B and C). Then
/// the two local matrices are combined to form a distributed
/// adjacency matrix for the Kronecker product of B and C.
pub fn generate_kronecker_graph(b_spec: &[u16], c_spec: &[u16], mode: u16) -> SparseMat {
    let bmat = SparseMat::gen_local_mat_from_stars(b_spec, mode);
    let cmat = SparseMat::gen_local_mat_from_stars(c_spec, mode);

    if bmat.my_rank() == 0 {
        println!(
            "Generating Mode {} Kronecker Product graph (A = B X C) with parameters: {:?} X {:?}",
            mode, b_spec, c_spec
        );
        println!("B has {} rows/cols and {} nnz", bmat.numrows(), bmat.nnz());
        println!("C has {} rows/cols and {} nnz", cmat.numrows(), cmat.nnz());
    }
    bmat.kron_prod_dist(&cmat, true)
}

/// Information that comes from trangle counting
#[derive(Debug)]
pub struct TriInfo {
    /// number of triangles
    pub tri_cnt: u64,
    /// number of remote ops to find
    pub sh_refs: u64,
    /// time to find
    pub laptime: f64,
}

impl TriInfo {
    fn new() -> TriInfo {
        TriInfo {
            tri_cnt: 0,
            sh_refs: 0,
            laptime: 0.0,
        }
    }
}

/// A trait to extend SparseMat for Triangle Counting
pub trait Triangle {
    /// push variant of triangle counting
    fn triangle_push(&self, upper: Option<&SparseMat>) -> TriInfo;
    /// pull variant of triangle counting
    fn triangle_pull(&self, upper: Option<&SparseMat>) -> TriInfo;
    /// generate a star graph
    fn gen_star(m: u16, mode: u16) -> SparseMat;
    /// generate a local matrix from star graph
    fn gen_local_mat_from_stars(spec: &[u16], mode: u16) -> SparseMat;
    /// generate a local kron graph
    fn kron_prod(&self, other: &SparseMat) -> SparseMat;
    /// generate a dist kron graph
    fn kron_prod_dist(&self, other: &SparseMat, mode: bool) -> SparseMat;
}

impl Triangle for SparseMat {
    /// This routine implements the push variant of triangle counting
    /// # Arguments
    /// * tmat the transpose of mat
    fn triangle_push(&self, upper: Option<&SparseMat>) -> TriInfo {
        let mut ret = TriInfo::new();
        let nr = self.numrows_this_rank();

        let t1 = wall_seconds();
        let mut success_count = 0;
        let mut num_pushed = 0;
        {
            let mut session = Convey::begin(|item: (usize, usize), _from_rank| {
                let w = item.0;
                let vj = item.1;
                let mat = match upper {
                    Some(m) => m,
                    None => self,
                };
                for nz in mat.offset[vj]..mat.offset[vj + 1] {
                    if w == mat.nonzero[nz] {
                        success_count += 1;
                        break;
                    }
                    if w < mat.nonzero[nz] {
                        // we have passed our slot
                        break;
                    }
                }
            });
            for row in 0..nr {
                for nz in self.offset[row]..self.offset[row + 1] {
                    let global_col = self.nonzero[nz];
                    let (vj, rank) = self.offset_rank(global_col);

                    let (mat, earlyout) = match upper {
                        Some(m) => (m, false),
                        None => (self, true),
                    };
                    for nz2 in mat.offset[row]..mat.offset[row + 1] {
                        let w = mat.nonzero[nz2];
                        if earlyout && w > global_col {
                            break;
                        }
                        session.push((w, vj), rank);
                        num_pushed += 1;
                    }
                }
            }
            session.finish();
        }
        ret.laptime = wall_seconds() - t1;
        ret.sh_refs = num_pushed;
        ret.tri_cnt = success_count;
        ret
    }

    fn triangle_pull(&self, _upper: Option<&SparseMat>) -> TriInfo {
        todo!();
    }
    fn gen_star(m: u16, mode: u16) -> SparseMat {
        let n = m as usize;
        let mut star = SparseMat::new_local(n + 1, n + 1, 2 * n + (if mode > 0 { 1 } else { 0 }));
        let mut pos = 0;

        if mode == 1 {
            star.nonzero[pos] = 0;
            pos += 1;
        }

        for i in 0..n {
            star.nonzero[pos] = i + 1;
            pos += 1;
        }
        star.offset[1] = pos;

        for i in 1..(n + 1) {
            star.nonzero[pos] = 0;
            pos += 1;
            star.offset[i + 1] = pos;
        }
        if mode == 2 {
            star.nonzero[pos] = n;
            pos += 1;
            star.offset[n + 1] = pos;
        }

        star
    }

    fn gen_local_mat_from_stars(spec: &[u16], mode: u16) -> SparseMat {
        let veclen = spec.len();
        if veclen == 1 {
            SparseMat::gen_star(spec[0], mode)
        } else {
            let mut mats: Vec<SparseMat> = Vec::new();

            for i in 0..veclen {
                mats.push(SparseMat::gen_star(spec[i], mode));
            }

            if veclen == 2 {
                mats[0].kron_prod(&mats[1])
            } else {
                mats.push(mats[0].kron_prod(&mats[1]));
                for i in 0..veclen - 3 {
                    mats.push(mats[2 + i].kron_prod(&mats[veclen + i]));
                }
                mats[2 + veclen - 3].kron_prod(&mats[veclen + veclen - 3])
            }
        }
    }

    fn kron_prod(&self, other: &SparseMat) -> SparseMat {
        let numrows = self.numrows() * other.numrows();
        let nnz = self.nnz() * other.nnz();
        if false {
            println!(
                "kpl r:{} n:{} x r:{} n:{} -> r:{} n:{}",
                self.numrows(),
                self.nnz(),
                other.numrows(),
                other.nnz(),
                numrows,
                nnz
            );
        }
        let mut ret = SparseMat::new_local(numrows, numrows, nnz);

        let mut offset: Vec<usize> = Vec::new();

        offset.push(0);
        for row in 0..self.numrows() {
            let d1 = self.offset[row + 1] - self.offset[row];
            for i in 0..other.numrows() {
                let d2 = other.offset[i + 1] - other.offset[i];
                offset.push(d1 * d2);
            }
        }
        for row in 0..numrows {
            offset[row + 1] = offset[row + 1] + offset[row];
            ret.offset[row + 1] = offset[row + 1];
        }

        for row in 0..self.numrows() {
            for j in self.offset[row]..self.offset[row + 1] {
                let col = self.nonzero[j];
                for otherrow in 0..other.numrows() {
                    for k in other.offset[otherrow]..other.offset[otherrow + 1] {
                        let index = row * other.numrows() + otherrow;
                        ret.nonzero[offset[index]] = col * other.numcols() + other.nonzero[k];
                        offset[index] += 1;
                    }
                }
            }
        }
        ret
    }
    /// generates a global matrix from two local or global matrices
    fn kron_prod_dist(&self, other: &SparseMat, lower: bool) -> SparseMat {
        let numrows = self.numrows_this_rank() * other.numrows_this_rank();
        let convey = Convey::new().expect("fixme");
        let my_rank = convey.my_rank();
        let num_ranks = convey.num_ranks();

        let mut nnz_this_rank = 0;

        for l_row in 0..convey.per_my_rank(numrows) {
            let g_row = (l_row * num_ranks) + my_rank;
            let row_b = g_row / other.numrows_this_rank();
            let row_c = g_row % other.numrows_this_rank();
            nnz_this_rank += (self.offset[row_b + 1] - self.offset[row_b])
                * (other.offset[row_c + 1] - other.offset[row_c]);
            if false {
                println!(
                    "nnz {} g {} b {} {} c {} {}",
                    nnz_this_rank,
                    g_row,
                    row_b,
                    self.offset[row_b + 1] - self.offset[row_b],
                    row_c,
                    other.offset[row_c + 1] - other.offset[row_c]
                );
            }
        }

        let mut ret = SparseMat::new(numrows, numrows, nnz_this_rank);

        let mut pos = 0;

        ret.offset[0] = 0;

        for l_row in 0..convey.per_my_rank(numrows) {
            let g_row = (l_row * num_ranks) + my_rank;
            let row_b = g_row / other.numrows_this_rank();
            let row_c = g_row % other.numrows_this_rank();
            for j in self.offset[row_b]..self.offset[row_b + 1] {
                let col_b = self.nonzero[j];
                for k in other.offset[row_c]..other.offset[row_c + 1] {
                    let col = col_b * other.numcols() + other.nonzero[k];
                    if !lower || col < g_row {
                        ret.nonzero[pos] = col;
                        pos += 1;
                    } else {
                        break;
                    }
                }
            }
            ret.offset[l_row + 1] = pos;
        }

        // the initial allocation of nonzeros was an upper bound, if
        // lower was set it could be shorter, adjust here
        ret.truncate_nonzeros(pos);
        ret
    }
}
