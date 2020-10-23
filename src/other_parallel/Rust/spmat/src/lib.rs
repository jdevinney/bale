#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Sparse Matrix Library
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For license information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use crate::perm::Perm;
use convey_hpc::collect::CollectValues;
use convey_hpc::Convey;
use rand::Rng;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Write;
use std::io::{Error, ErrorKind};
use std::rc::Rc;

use std::path::Path;

use std::time::{SystemTime, UNIX_EPOCH};
/// A routine to give access to the wall clock timer on most UNIX-like systems.
///    Uses rust's SystemTime.
pub fn wall_seconds() -> f64 {
    let n = SystemTime::now().duration_since(UNIX_EPOCH).unwrap();
    (n.as_secs() as f64) + (n.as_micros() as f64) * 1.0e-6
}

/// An element of a sparse matrix as a triple
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
struct Entry {
    row: usize,
    col: usize,
    val: f64,
}

/// A sparse matrix.  Data is stored as a num_rows_this_rank size
/// vector of offsets in the nonzero vector, which has column numbers
/// for each row.
pub struct SparseMat {
    numrows: usize, // the total number of rows in the matrix
    numrows_this_rank: usize,
    numcols: usize, // the nonzeros have values between 0 and numcols
    nnz: usize,     // total number of nonzeros in the matrix
    nnz_this_rank: usize,
    /// the row offsets into the array of nonzeros, size is nrows+1,
    /// offsets[nrows] is nnz
    pub offset: Vec<usize>,
    /// nonzero list
    pub nonzero: Vec<usize>, // the global array of nonzeros
    /// nonzero values, if present
    pub value: Option<Vec<f64>>, // the global array of nonzero values, optional
    convey: Option<Convey>, // a conveyor for our use
}

impl std::fmt::Debug for SparseMat {
    // 0-0 should add values
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let debug_rows = self.numrows_this_rank.min(10);
        let mut rows: Vec<(usize, Vec<usize>)> = Vec::new();
        for i in 0..debug_rows {
            let start = self.offset[i];
            let end = self.offset[i + 1];
            let debug_cols = (end - start).min(10);
            let mut row: Vec<usize> = Vec::new();
            for j in start..start + debug_cols {
                row.push(self.nonzero[j]);
            }
            rows.push((end - start, row));
        }
        f.debug_struct("SparseMat")
            .field("rank", &self.my_rank())
            .field("numrows", &self.numrows)
            .field("numcols", &self.numcols)
            .field("nnz", &self.nnz)
            .field("nnz.len()", &self.nonzero.len())
            .field("offset.len()", &self.offset.len())
            .field("first_part", &rows)
            .finish()
    }
}

impl SparseMat {
    /// new creates a distributed sparse matrix across all ranks using
    /// number of rows, number of columns, and number of rows on
    /// current rank
    pub fn new(numrows: usize, numcols: usize, nnz_this_rank: usize) -> Self {
        let convey = Convey::new().expect("fixme");
        let total_nnz = nnz_this_rank.reduce_sum();
        SparseMat {
            numrows: numrows,
            numrows_this_rank: convey.per_my_rank(numrows),
            numcols: numcols,
            nnz: total_nnz,
            nnz_this_rank: nnz_this_rank,
            offset: vec![0; convey.per_my_rank(numrows) + 1],
            nonzero: vec![0; nnz_this_rank],
            value: None,
            convey: Some(convey),
        }
    }

    /// new_with_values creates a distributed sparse matrix intialized with zero values
    pub fn new_with_values(numrows: usize, numcols: usize, nnz_this_rank: usize) -> Self {
        let convey = Convey::new().expect("fixme");
        let total_nnz = nnz_this_rank.reduce_sum();
        SparseMat {
            numrows: numrows,
            numrows_this_rank: convey.per_my_rank(numrows),
            numcols: numcols,
            nnz: total_nnz,
            nnz_this_rank: nnz_this_rank,
            offset: vec![0; convey.per_my_rank(numrows) + 1],
            nonzero: vec![0; nnz_this_rank],
            value: Some(vec![0.0_f64; nnz_this_rank]),
            convey: Some(convey),
        }
    }

    /// new_local creates a local sparse matrix on each rank using
    /// number of rows, number of columns, and number of rows on
    /// current rank
    pub fn new_local(numrows: usize, numcols: usize, nnz_this_rank: usize) -> Self {
        SparseMat {
            numrows: numrows,
            numrows_this_rank: numrows,
            numcols: numcols,
            nnz: nnz_this_rank,
            nnz_this_rank: nnz_this_rank,
            offset: vec![0; numrows + 1],
            nonzero: vec![0; nnz_this_rank],
            value: None,
            convey: None,
        }
    }

    /// new_local creates a local sparse matrix, with values on each rank using
    /// number of rows, number of columns, and number of rows on
    /// current rank
    pub fn new_local_with_values(numrows: usize, numcols: usize, nnz_this_rank: usize) -> Self {
        SparseMat {
            numrows: numrows,
            numrows_this_rank: numrows,
            numcols: numcols,
            nnz: nnz_this_rank,
            nnz_this_rank: nnz_this_rank,
            offset: vec![0; numrows + 1],
            nonzero: vec![0; nnz_this_rank],
            value: Some(vec![0.0_f64; nnz_this_rank]),
            convey: None,
        }
    }

    /// collect all of a distributed sparse matrix to a local matrix
    /// on rank 0.  must be called collectively by all ranks, returns
    /// empty matrix of same dims on all ranks but 0.
    pub fn to_local(&self) -> SparseMat {
        todo!();
    }

    /// distribute a local sparse matrix on rank 0 to a distributed
    /// sparse matrix.  must be called collectively by all ranks,
    /// though only rank 0 supplies the input matrix.
    pub fn to_distributed(&self) -> SparseMat {
        todo!();
    }

    /// randomize the values of a sparse matrix
    pub fn randomize_values(&mut self) {
        let mut value = Vec::new();
        let mut rng = rand::thread_rng();
        for _ in 0..self.nnz_this_rank {
            value.push(rng.gen::<f64>());
        }
        self.value = Some(value);
    }

    /// an iterator over row counts
    pub fn rowcounts<'a>(&'a self) -> impl Iterator<Item = usize> + 'a {
        self.offset[0..self.numrows_this_rank]
            .iter()
            .zip(&self.offset[1..self.numrows_this_rank + 1])
            .map(|(a, b)| b - a)
    }

    /// convenience function to return my rank
    pub fn my_rank(&self) -> usize {
        if self.convey.is_some() {
            self.convey.as_ref().unwrap().my_rank()
        } else {
            0
        }
    }

    /// Convenience function returning num_ranks
    pub fn num_ranks(&self) -> usize {
        if self.convey.is_some() {
            self.convey.as_ref().unwrap().num_ranks()
        } else {
            1
        }
    }

    /// Convenience function returning the number of elements on
    /// my_rank for a total number of elements
    pub fn per_my_rank(&self, n: usize) -> usize {
        if let Some(convey) = &self.convey {
            convey.per_my_rank(n)
        } else {
            n
        }
    }

    /// Convenience function returning the offset and rank in a
    /// distributed sparsemat for a given index
    pub fn offset_rank(&self, n: usize) -> (usize, usize) {
        if let Some(convey) = &self.convey {
            convey.offset_rank(n)
        } else {
            (n, 0)
        }
    }

    /// global row index to local row index on this rank
    pub fn local_index(&self, n: usize) -> usize {
        let (offset, rank) = self.offset_rank(n);
        if rank != self.my_rank() {
            panic!("attempt to get local index for row not on this rank");
        }
        offset
    }

    /// local row index on this rank to global row index
    pub fn global_index(&self, n: usize) -> usize {
        self.num_ranks() * n + self.my_rank()
    }

    /// inspector: num_rows()
    pub fn numrows(&self) -> usize {
        self.numrows
    }

    /// inspector: num_rows_this_rank()
    pub fn numrows_this_rank(&self) -> usize {
        self.numrows_this_rank
    }

    /// inspector: num_cols()
    pub fn numcols(&self) -> usize {
        self.numcols
    }

    /// inspector: nnz()
    pub fn nnz(&self) -> usize {
        self.nnz
    }

    /// barrier that works for local SparseMat too
    pub fn barrier(&self) {
        if let Some(convey) = &self.convey {
            convey.barrier()
        }
    }

    /// prints some stats of a sparse matrix
    pub fn stats(&self) -> () {
        let val_str = if let Some(_) = self.value {
            "(with values)"
        } else {
            "(pattern only)"
        };
        if let Some(_) = self.convey {
            if self.my_rank() != 0 {
                return ();
            }
            println!("    distributed matrix {}", val_str);
            println!("    num_ranks = {}", self.num_ranks());
            println!("    numrows.0 = {}", self.numrows_this_rank);
            println!("    nnz.0     = {}", self.nnz_this_rank);
        } else {
            println!("    local matrix {} on rank {}", val_str, self.my_rank());
        }
        println!("    numrows   = {}", self.numrows);
        println!("    numcols   = {}", self.numcols);
        println!("    nnz       = {}", self.nnz);

        // compute min, max, sum all at once for efficiency, only once thru itereator
        let (mindeg, maxdeg, sumdeg) = self.rowcounts().fold((self.numcols, 0, 0), |acc, x| {
            (acc.0.min(x), acc.1.max(x), acc.2 + x)
        });

        let avgdeg = sumdeg as f64 / self.numrows as f64;

        println!(
            "    min, avg, max degree (rank 0) = {}, {}, {}",
            mindeg, avgdeg, maxdeg
        );
    }

    /// dumps a sparse matrix to a file in a ASCII format
    /// # Arguments
    /// * maxrows the number of rows that are written, 0 means everything,
    ///           otherwise write the first and last maxrows/2 rows
    /// * filename the filename to written to
    pub fn dump(&self, maxrows: usize, filename: &str) -> Result<(), Error> {
        let path = Path::new(&filename);
        let mut file = File::create(path)?;

        let use_maxrow: bool = if maxrows < self.numrows && maxrows > 0 {
            true
        } else {
            false
        };
        let start_row = if use_maxrow {
            self.numrows - maxrows / 2
        } else {
            self.numrows
        };
        let stop_row = if use_maxrow {
            maxrows / 2
        } else {
            self.numrows
        };

        writeln!(file, "\n--------- offsets:")?;
        for off in &self.offset[0..stop_row + 1] {
            write!(file, "{} ", off)?;
        }
        if use_maxrow {
            write!(file, "... ")?;
            for off in &self.offset[start_row..self.numrows] {
                write!(file, "{}", off)?;
            }
        }

        writeln!(
            file,
            "{}",
            match &self.value {
                None => "\n--------- row col:",
                Some(_) => "\n--------- row col val:",
            },
        )?;

        for i in 0..stop_row {
            for k in self.offset[i]..self.offset[i + 1] {
                if let Some(value) = &self.value {
                    writeln!(file, "{} {} {}", i, self.nonzero[k], value[k])?;
                } else {
                    writeln!(file, "{} {}", i, self.nonzero[k])?;
                }
            }
        }
        if use_maxrow {
            write!(file, ".\n.\n.\n")?;
            for i in start_row..self.numrows {
                for k in self.offset[i]..self.offset[i + 1] {
                    if let Some(value) = &self.value {
                        writeln!(file, "{} {} {}", i, self.nonzero[k], value[k])?;
                    } else {
                        writeln!(file, "{} {}", i, self.nonzero[k])?;
                    }
                }
            }
        }
        Ok(())
    }

    /// read a sparse matrix from a file in a MatrixMarket ASCII format
    /// # Arguments
    /// * filename the file to be read
    pub fn read_mm_file(filename: &str) -> Result<SparseMat, Error> {
        let path = Path::new(&filename);
        let fp = File::open(path)?;
        let mut reader = BufReader::new(fp);
        SparseMat::read_mm(&mut reader)
    }

    /// read a sparse matrix from an arbitrary reader
    pub fn read_mm<R>(reader: &mut R) -> Result<SparseMat, Error>
    where
        R: BufRead,
    {
        let line = reader.lines().next().unwrap_or(Ok("".to_string()))?;
        let re1 = Regex::new(r"^%%MatrixMarket *matrix *coordinate *pattern*").expect("bad regex");
        let re2 = Regex::new(r"^%%MatrixMarket *matrix *coordinate *real*").expect("bad regex");
        if !re1.is_match(&line) && !re2.is_match(&line) {
            return Err(Error::new(
                ErrorKind::Other,
                format!("invalid MatrixMarket header {}", line),
            ));
        }
        let has_values = re2.is_match(&line);

        let line = reader.lines().next().unwrap_or(Ok("".to_string()))?;
        let re = Regex::new(r"^(\d*)\s*(\d*)\s*(\d*)\s*").expect("bad regex");
        let cleaned = re.replace_all(&line, "$1 $2 $3");
        let inputs: Vec<String> = cleaned.split(" ").map(|x| x.to_string()).collect();
        let nr: usize = inputs[0].parse().expect("bad row dimension");
        let nc: usize = inputs[1].parse().expect("bad column dimension");
        let nnz: usize = inputs[2].parse().expect("bad nnz");
        if nr == 0 || nc == 0 {
            panic!(format!("bad matrix sizes {} {} {}", nr, nc, nnz));
        }
        let mut matrix = SparseMat::new(nr, nc, 0); // we don't know nnz_this_rank yet
        let mut elts: Vec<Entry> = vec![];
        for line in reader.lines() {
            let line = line.expect("read error");
            let inputs: Vec<&str> = line.split_whitespace().collect();
            let row: usize = inputs[0].parse::<usize>().expect("bad row number") - 1;
            let col: usize = inputs[1].parse::<usize>().expect("bad column number") - 1;
            let val: f64 = if has_values {
                inputs[2].parse().expect("bad element value")
            } else {
                1.0
            };
            if row >= nr || col >= nc {
                panic!(format!("bad matrix nonzero coordinates {} {}", row, col));
            }
            if matrix.offset_rank(row).1 == matrix.my_rank() {
                elts.push(Entry {
                    row: matrix.local_index(row),
                    col: col,
                    val: val,
                });
            }
        }
        matrix.nnz_this_rank = elts.len();
        matrix.nnz = matrix.nnz_this_rank.reduce_sum();
        if matrix.nnz != nnz {
            panic!(format!(
                "incorrect number of nonzeros {} != {}",
                matrix.nnz, nnz
            ));
        }

        elts.sort_by(|a, b| {
            if a.row == b.row {
                a.col.cmp(&b.col)
            } else {
                a.row.cmp(&b.row)
            }
        });

        let mut values: Vec<f64> = vec![];
        let mut row: usize = 0;
        for (i, elt) in elts.iter().enumerate() {
            // first adjust so we are on correct row
            while row < elt.row {
                row += 1;
                matrix.offset[row] = i;
            }
            assert!(elt.row == row);
            matrix.nonzero.push(elt.col);
            if has_values {
                values.push(elt.val);
            }
        }
        for i in row + 1..matrix.numrows_this_rank + 1 {
            matrix.offset[i] = matrix.nnz_this_rank;
        }
        if has_values {
            matrix.value = Some(values);
        }
        Ok(matrix)
    }

    /// write a sparse matrix to a file in a MatrixMarket ASCII format
    /// # Arguments
    /// * filename the filename to written to
    pub fn write_mm_file(&self, filename: &str) -> Result<(), Error> {
        // use a new conveyor, not the one in the SparseMat, to avoid borrow issues
        let my_rank = self.my_rank();
        let file = if my_rank > 0 {
            Rc::new(RefCell::new(None))
        } else {
            let path = Path::new(&filename);
            Rc::new(RefCell::new(Some(File::create(path)?)))
        };
        if let Some(mut f) = file.borrow().as_ref() {
            if let Some(_) = &self.value {
                writeln!(f, "%%MatrixMarket matrix coordinate real general")
                    .expect("can't write .mm file");
            } else {
                writeln!(f, "%%MatrixMarket matrix coordinate pattern general")
                    .expect("can't write .mm file");
            }
            writeln!(f, "{} {} {}", self.numrows, self.numcols, self.nnz)
                .expect("can't write .mm file");
        };
        {
            let mut session = Convey::begin(|entry: Entry, _from_rank| {
                // only rank 0 will ever get asked to do this
                if let Some(mut f) = file.borrow().as_ref() {
                    if let Some(_) = &self.value {
                        writeln!(f, "{} {} {}", entry.row, entry.col, entry.val)
                            .expect("can't write .mm file");
                    } else {
                        writeln!(f, "{} {}", entry.row, entry.col).expect("can't write .mm file");
                    }
                }
            });
            for i in 0..self.numrows_this_rank {
                let i_g = self.global_index(i);
                for k in self.offset[i]..self.offset[i + 1] {
                    session.push(
                        Entry {
                            row: i_g + 1,
                            col: self.nonzero[k] + 1,
                            val: if let Some(value) = &self.value {
                                value[k]
                            } else {
                                0.0
                            },
                        },
                        0, // always push to rank 0, which writes the file
                    );
                }
            }
            session.finish();
        }
        Ok(())
    }

    /// perdicate: is the SparseMat lower triangular
    pub fn is_lower_triangular(&self, unit_diagonal: bool) -> bool {
        let mut lower_cnt = 0;
        let mut diag_missing_cnt = 0;
        for i in 0..self.numrows_this_rank {
            let global_row = i * self.num_ranks() + self.my_rank();
            let mut pivot = false;
            for nz in &self.nonzero[self.offset[i]..self.offset[i + 1]] {
                if *nz < global_row {
                    lower_cnt += 1;
                } else if *nz == global_row {
                    pivot = true;
                }
            }
            if !pivot {
                diag_missing_cnt += 1;
            }
        }
        let total_lower = lower_cnt.reduce_sum();
        let total_diag_missing = if unit_diagonal {
            diag_missing_cnt.reduce_sum()
        } else {
            0
        };
        total_lower != 0 || total_diag_missing != 0
    }

    /// perdicate: is the SparseMat upper triangular
    pub fn is_upper_triangular(&self, unit_diagonal: bool) -> bool {
        let mut lower_cnt = 0;
        let mut diag_missing_cnt = 0;
        for i in 0..self.numrows_this_rank {
            let global_row = i * self.num_ranks() + self.my_rank();
            let mut pivot = false;
            for nz in &self.nonzero[self.offset[i]..self.offset[i + 1]] {
                if *nz > global_row {
                    lower_cnt += 1;
                } else if *nz == global_row {
                    pivot = true;
                }
            }
            if !pivot {
                diag_missing_cnt += 1;
            }
        }
        let total_lower = lower_cnt.reduce_sum();
        let total_diag_missing = if unit_diagonal {
            diag_missing_cnt.reduce_sum()
        } else {
            0
        };
        total_lower != 0 || total_diag_missing != 0
    }

    /// createa a new SparseMat which is the permutation of rows and
    /// columns as defined by given permutations
    pub fn permute(&self, rperminv: &Perm, cperminv: &Perm) -> Self {
        if let Some(_) = &self.value {
            todo!()
        }
        let mut rowcounts = vec![0_usize; self.numrows_this_rank];
        assert_eq!(rperminv.len(), self.numrows_this_rank);
        assert_eq!(cperminv.len(), self.numrows_this_rank);
        assert!(rperminv.is_perm());
        assert!(cperminv.is_perm());
        let mut nnz = 0_usize;
        Convey::simple(
            (0..rperminv.len()).map(|idx| {
                let (offset, rank) = self.offset_rank(rperminv.entry(idx));
                let cnt = self.offset[idx + 1] - self.offset[idx];
                ((offset, cnt), rank)
            }),
            |item: (usize, usize), _from_rank| {
                rowcounts[item.0] = item.1;
                nnz += item.1;
            },
        );

        assert_eq!(self.nnz, nnz.reduce_sum());
        let mut permuted = SparseMat::new(self.numrows, self.numcols, nnz);
        permuted.offset[0] = 0;

        // step 2: set up new offsets
        for i in 1..permuted.numrows_this_rank + 1 {
            permuted.offset[i] = permuted.offset[i - 1] + rowcounts[i - 1];
        }
        //println!("{:?}{:?}", permuted, permuted.offset);
        assert_eq!(
            permuted.offset[permuted.numrows_this_rank],
            permuted.nnz_this_rank
        );
        // step 3: distrbute nonzeros
        let mut wrkoff = vec![0_usize; self.numrows_this_rank];
        {
            let mut session = Convey::begin(|item: (usize, usize), _from_rank| {
                let index = permuted.offset[item.0] + wrkoff[item.0];
                permuted.nonzero[index] = item.1;
                wrkoff[item.0] += 1;
            });
            let mut row = 0;
            for i in 0..self.nnz_this_rank {
                while i == self.offset[row + 1] {
                    row += 1;
                }
                let (offset, rank) = self.offset_rank(rperminv.entry(row));
                session.push((offset, self.nonzero[i]), rank);
            }
            session.finish();
        }

        // step 4: do column permutation (essentailly indexgather)
        {
            let my_nnz = permuted.nnz_this_rank;
            let cloned_nonzero = permuted.nonzero.clone();
            Convey::simple_return(
                (0..my_nnz).map(|i| {
                    let (offset, rank) = self.offset_rank(cloned_nonzero[i]);
                    ((i, offset), rank)
                }),
                |item: usize| cperminv.entry(item),
                |index: usize, val: usize| {
                    permuted.nonzero[index] = val;
                },
            );
        }
        permuted
    }

    /// Create a new SparseMat which is the transpose of the given SparseMat
    pub fn transpose(&self) -> Self {
        if let Some(_) = &self.value {
            todo!()
        }
        let mut colcnt = vec![0_usize; self.per_my_rank(self.numrows)];
        let mut nnz = 0_usize;

        // distributed calculation of column counts, with resulting local nonzeros
        Convey::simple(
            self.nonzero.iter().map(|nz| self.offset_rank(*nz)),
            |item: usize, _from_rank| {
                colcnt[item] += 1;
                nnz += 1;
            },
        );

        let mut trans = SparseMat::new(self.numcols, self.numrows, nnz);
        if self.nnz != trans.nnz {
            println!("self: {:?} colcnt: {:?}", self, colcnt);
        }
        assert_eq!(self.nnz, trans.nnz);

        trans.offset[0] = 0;
        for idx in 1..self.offset.len() {
            trans.offset[idx] = trans.offset[idx - 1] + colcnt[idx - 1];
        }

        let mut wrkoff = vec![0_usize; trans.numrows];
        {
            let mut session = Convey::begin(|item: (usize, usize), _from_rank| {
                let index = trans.offset[item.0] + wrkoff[item.0];
                trans.nonzero[index] = item.1;
                wrkoff[item.0] += 1;
            });
            let mut row = 0;
            for i in 0..self.nonzero.len() {
                while i == self.offset[row + 1] {
                    row += 1;
                }
                let (offset, rank) = self.offset_rank(self.nonzero[i]);
                let tcol = row * self.num_ranks() + self.my_rank();
                session.push((offset, tcol), rank);
            }
            session.finish();
        }
        trans
    }

    /// Predicate: are two SparseMats equal?
    pub fn compare(&self, other: &SparseMat) -> bool {
        if (self.value == None && other.value != None)
            || (self.value != None && other.value == None)
        {
            println!("presence of values differ {:?} {:?}", self, other);
            return false;
        }
        if (self.numcols != other.numcols)
            || (self.numrows != other.numrows)
            || (self.numrows_this_rank != other.numrows_this_rank)
            || (self.nnz != other.nnz)
            || (self.nnz_this_rank != other.nnz_this_rank)
        {
            println!("counts differ {:?} {:?}", self, other);
            return false;
        }
        if self.offset != other.offset {
            println!("offsets differ {:?} {:?}", self, other);
            return false;
        }
        // need to sort nonzeros before compare 0-0 jg: why?
        for row in 0..self.numrows_this_rank {
            // we already know self & other offsets are the same
            let range = self.offset[row]..self.offset[row + 1];
            if let Some(sval) = &self.value {
                if let Some(oval) = &other.value {
                    let mut self_sorted: Vec<_> = self.nonzero[range.clone()]
                        .iter()
                        .zip(&sval[range.clone()])
                        .collect();
                    let mut other_sorted: Vec<_> = other.nonzero[range.clone()]
                        .iter()
                        .zip(&oval[range.clone()])
                        .collect();
                    self_sorted.sort_by(|a, b| b.0.cmp(&a.0));
                    other_sorted.sort_by(|a, b| b.0.cmp(&a.0));
                    if self_sorted != other_sorted {
                        println!("nonzeros differ {:?} {:?}", self, other);
                        return false;
                    }
                }
            } else {
                let mut self_sorted: Vec<_> = self.nonzero[range.clone()].iter().collect();
                let mut other_sorted: Vec<_> = other.nonzero[range.clone()].iter().collect();
                self_sorted.sort();
                other_sorted.sort();
                if self_sorted != other_sorted {
                    println!("nonzeros differ {:?} {:?}", self, other);
                    return false;
                }
            }
        }
        true
    }

    /// Create a new SparseMat which is the sum of 2 given
    pub fn add(&self, other: &SparseMat) -> Self {
        if let Some(_) = &self.value {
            todo!()
        }
        // need to check for errors
        let mut sum = SparseMat::new(
            self.numrows,
            self.numcols,
            self.nnz_this_rank + other.nnz_this_rank,
        );
        let mut nnz = 0;
        for row in 0..self.numrows_this_rank {
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                sum.nonzero[nnz] = *nz;
                nnz += 1;
            }
            for nz in &other.nonzero[other.offset[row]..other.offset[row + 1]] {
                sum.nonzero[nnz] = *nz;
                nnz += 1;
            }
            sum.offset[row + 1] = nnz;
        }
        sum
    }

    /// truncate the nonzero array and update as appropraite
    pub fn truncate_nonzeros(&mut self, last: usize) {
        self.nonzero.truncate(last);
        self.nnz_this_rank = last;
        self.nnz = last.reduce_sum();
    }

    /// Generate an erdos_renyi random sparse matrix
    pub fn gen_erdos_renyi_graph(
        num_vert: usize,
        prob: f64,
        unit_diag: bool,
        mode: u8,
        seed: i64,
    ) -> Self {
        if mode == 0 {
            // symmetric matrix, undirected graph
            let upper = SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, false, seed);
            let lower = upper.transpose();
            upper.add(&lower)
        } else if mode == 1 {
            // lower triangular matrix
            SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, true, seed)
        } else if mode == 2 {
            // upper triangular matrix
            SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, false, seed)
        } else {
            // nonsymmetric matrix, directed graph
            let upper = SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, false, seed);
            let lower = SparseMat::erdos_renyi_tri(num_vert, prob, unit_diag, true, seed + 1);
            upper.add(&lower)
        }
    }

    /// Worker function: generate an erdos_renyi triangular random sparse matrix
    fn erdos_renyi_tri(
        num_vert: usize,
        prob: f64,
        unit_diag: bool,
        lower: bool,
        _seed: i64,
    ) -> Self {
        let mut tri = SparseMat::new(num_vert, num_vert, 0);
        // let mut tri = SparseMat::new_local(num_vert, num_vert, 0); // jg: make local mtx for test
        let l_max = (std::u32::MAX as f64).ln();
        let d = (1.0 - prob).ln();
        let numrows = num_vert;
        let mut row = tri.my_rank();
        let mut col = if lower { 0 } else { tri.my_rank() + 1 };
        let num_ranks = tri.num_ranks();

        // if we are upper triangular and unit_diag, we need to set
        // the first row
        if unit_diag && !lower {
            tri.nonzero.push(tri.my_rank());
        }
        while row < num_vert {
            let mut r = rand::random::<u32>();
            if r == std::u32::MAX {
                r -= 1;
            }

            col += (1.0 + ((((std::u32::MAX - r) as f64).ln() - l_max) / d).floor()) as usize;
            while (col >= (if lower { row } else { numrows })) && row < numrows {
                if lower {
                    if unit_diag {
                        tri.nonzero.push(row);
                    }
                    col = col - row;
                    row += num_ranks;
                    tri.offset[row / num_ranks] = tri.nonzero.len();
                } else {
                    row += num_ranks;
                    col = row + 1 + col - numrows;
                    tri.offset[row / num_ranks] = tri.nonzero.len();
                    if (row < numrows) && unit_diag {
                        tri.nonzero.push(row);
                    }
                }
            }
            if row < numrows {
                tri.nonzero.push(col);
            }
        }
        tri.offset[row / num_ranks] = tri.nonzero.len();
        tri.nnz_this_rank = tri.nonzero.len();
        tri.nnz = tri.nnz_this_rank.reduce_sum();
        tri
    }
}

/// Bring in the permutation library
pub mod perm;

#[cfg(test)]
mod tests {
    use super::SparseMat;
    use crate::Perm;
    use convey_hpc::testing_support::TestingMutex;
    #[test]
    fn rand_mat_tri_upper() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, false, false, 0);
        assert_eq!(mat.is_lower_triangular(false), false);
        assert_eq!(mat.is_upper_triangular(false), true);
        drop(mutex);
    }
    #[test]
    fn rand_mat_tri_lower() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, false, true, 0);
        assert_eq!(mat.is_lower_triangular(false), true);
        assert_eq!(mat.is_upper_triangular(false), false);
        drop(mutex);
    }
    #[test]
    fn rand_mat_tri_upper_unit() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, true, false, 0);
        assert_eq!(mat.is_lower_triangular(true), false);
        assert_eq!(mat.is_upper_triangular(true), true);
        drop(mutex);
    }
    #[test]
    fn rand_mat_tri_lower_unit() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, true, true, 0);
        assert_eq!(mat.is_lower_triangular(true), true);
        assert_eq!(mat.is_upper_triangular(true), false);
        drop(mutex);
    }
    #[test]
    fn rand_mat_perm1() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(1000, 0.05, true, false, 0);
        let rperm = Perm::random(1000, 0);
        let cperm = Perm::random(1000, 0);
        let _permuted = mat.permute(&rperm, &cperm);
        assert_eq!(mat.is_lower_triangular(true), false);
        assert_eq!(mat.is_upper_triangular(true), true);
        //assert_eq!(permuted.is_lower_triangular(true), false);
        //assert_eq!(permuted.is_upper_triangular(true), false);
        drop(mutex);
    }
    /*
        #[test]
        #[should_panic]
        fn rand_mat_perm2() {
            let mutex = TestingMutex::new();
            let mut mat = SparseMat::erdos_renyi_tri(1000, 0.05, true, false, 0);
            mat.randomize_values();
            let rperm = Perm::random(1000, 0);
            let cperm = Perm::random(1000, 0);
            let _permuted = mat.permute(&rperm, &cperm);
            assert_eq!(mat.my_rank(), mutex.convey.my_rank);
        }
    */
    #[test]
    fn transpose1() {
        let mutex = TestingMutex::new();
        let mat = SparseMat::erdos_renyi_tri(1000, 0.05, true, false, 0);
        let tmat = mat.transpose();
        let ttmat = tmat.transpose();
        println!("mat:{:?}", mat);
        println!("tmat:{:?}", tmat);
        println!("ttmat:{:?}", ttmat);
        assert_eq!(mat.compare(&tmat), false);
        assert_eq!(mat.compare(&ttmat), true);
        drop(mutex);
    }
    /*
        #[test]
        #[should_panic]
        fn transpose2() {
            let mutex = TestingMutex::new();
            let mut mat = SparseMat::erdos_renyi_tri(1000, 0.05, true, false, 0);
            mat.randomize_values();
            let _ = mat.transpose();
            assert_eq!(mat.my_rank(), mutex.convey.my_rank);
        }
    */
    #[test]
    fn write_read_mm1() {
        let mutex = TestingMutex::new();
        assert!(SparseMat::read_mm_file("not_there.mm").is_err());
        drop(mutex);
    }
    #[test]
    fn write_read_mm2() {
        let mutex = TestingMutex::new();
        let mat_a = SparseMat::erdos_renyi_tri(100, 0.05, true, false, 0);
        let f1 = format!(
            "{}/test_write_read_mm2.mm",
            std::env::temp_dir().to_str().unwrap()
        );
        mat_a.write_mm_file(&f1).expect("failed write");
        let mat_b = SparseMat::read_mm_file(&f1).expect("failed read");
        assert!(mat_b.value == None);
        assert_eq!(mat_a.compare(&mat_b), true);
        drop(mutex);
    }
    #[test]
    fn write_read_mm3() {
        let mutex = TestingMutex::new();
        let mut mat_a = SparseMat::erdos_renyi_tri(100, 0.05, true, false, 0);
        let f1 = format!(
            "{}/test_write_read_mm3.mm",
            std::env::temp_dir().to_str().unwrap()
        );
        mat_a.randomize_values();
        mat_a.write_mm_file(&f1).expect("failed write");

        let mat_b = SparseMat::read_mm_file(&f1).expect("failed read");
        assert!(mat_b.value != None);
        assert_eq!(mat_a.compare(&mat_b), true);
        drop(mutex);
    }
}
