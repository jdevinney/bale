#![warn(
    missing_docs,
    future_incompatible,
    missing_debug_implementations,
    rust_2018_idioms
)]

//! Bale Serial Sparsemat library
///
/// Copyright (c) 2020, Institute for Defense Analyses
/// 4850 Mark Center Drive, Alexandria, VA 22311-1882; 703-845-2500
///
/// All rights reserved.
///
/// This file is part of Bale.  For license information see the
/// LICENSE file in the top level dirctory of the distribution.
///
use rand::Rng;
use regex::Regex;
use std::fs::File;
use std::io::{BufRead, Write};
use std::io::{BufReader, BufWriter};
use std::path::Path;

/// Implement our own error handling extention for bad formats of matrix market
mod err;
use err::SparseMatError::Sparsemat;
use err::{ParseMmError, SparseMatError};

/// Generic result return in this library
type Result<T> = std::result::Result<T, SparseMatError>;

/// A structure to hold a sparse matrix
#[derive(Debug, Clone)]
pub struct SparseMat {
    numrows: usize, // the total number of rows in the matrix
    numcols: usize, // the nonzeros have values between 0 and numcols
    nnz: usize,     // total number of nonzeros in the matrix
    /// the row offsets into the arrays nonzeros and values, size is nrows+1,
    /// offset[nrows] is nnz
    pub offset: Vec<usize>,
    /// The nonzero column indicies
    pub nonzero: Vec<usize>,
    /// The values, if present
    pub value: Option<Vec<f64>>,
}

/// struct to hold the state used to iterate across the row of a sparse matrix.
#[derive(Debug)]
pub struct NextNz {
    mat: SparseMat,
    start: usize,
    idx: usize,
    stop: usize,
    col: usize,
    val: Option<f64>,
}

use std::time::{SystemTime, UNIX_EPOCH};
/// A routine to give access to the wall clock timer on most UNIX-like systems.
///    Uses rust's SystemTime.
pub fn wall_seconds() -> Result<f64> {
    let n = SystemTime::now().duration_since(UNIX_EPOCH)?;
    Ok((n.as_secs() as f64) + (n.as_micros() as f64) * 1.0e-6)
}

/// A permutation
#[derive(Debug, Clone)]
pub struct Perm {
    perm: Vec<usize>,
}

use rand::seq::SliceRandom;
use rand::thread_rng;

impl Perm {
    /// a new permutation of some size
    pub fn new(n: usize) -> Perm {
        let perm: Vec<usize> = (0..n).collect();
        Perm { perm }
    }

    /// Get the length of this permutation
    pub fn len(&self) -> usize {
        self.perm.len()
    }

    /// get the entry from a permutation
    pub fn entry(&self, index: usize) -> usize {
        self.perm[index]
    }

    /// get a reference to the permutation itself
    pub fn perm(&mut self) -> &mut Vec<usize> {
        &mut self.perm
    }

    /// create an int64_t array which holds a uniform random permutation
    /// # Arguments
    /// * n the length of the global array
    /// * _seed seed for the random number generator (currently unused)
    ///
    /// Rust suffle implements the standard serial algorithm, known at least as
    ///     Fisher-Yates or Knuth shuffle, to generate the uniform permutation.
    pub fn random(n: usize, _seed: i64) -> Perm {
        let mut ret = Perm::new(n);
        ret.perm.shuffle(&mut thread_rng());
        ret
    }

    /// create a permutation which will reverse the given perm
    pub fn inverse(&self) -> Perm {
        let mut ret = Perm::new(self.perm.len());
        ret.perm = ret
            .perm
            .iter()
            .map(|&x| self.perm.iter().position(|&y| y == x).expect("not found"))
            .collect();
        ret
    }

    /// checks that an array is in fact a permutation
    /// # Arguments
    /// * perm pointer to the array
    ///
    ///Every element in the flag array will equal 1 iff perm is a permutation.
    pub fn is_perm(&self) -> bool {
        let len: usize = self.perm.len();
        let mut flag: Vec<usize> = vec![0; len];

        for entry in &self.perm {
            if *entry >= len {
                return false;
            }
            flag[*entry as usize] += 1;
        }

        for entry in flag {
            if entry != 1 {
                return false;
            }
        }
        true
    }

    /// writes the first and last part of an perm to the specified file
    /// # Arguments
    /// * a       the array
    /// * maxdisp the number of entries that are written, 0 means everything,
    ///            otherwise write the first and last maxdisp/2 entries
    /// * filename the filename to written to
    pub fn dump(&self, maxdisp: usize, filename: &str) -> Result<()> {
        let path = Path::new(&filename);
        let mut file = File::create(path)?;

        let safe_disp = if maxdisp <= self.perm.len() && maxdisp > 0 {
            maxdisp / 2
        } else {
            self.perm.len() / 2
        };

        for entry in &self.perm[0..safe_disp] {
            write!(file, "{}\n", entry)?;
        }
        for entry in &self.perm[self.perm.len() - safe_disp..self.perm.len()] {
            write!(file, "{}\n", entry)?;
        }
        Ok(())
    }
}

impl SparseMat {
    /// a new sparse matrix without values
    pub fn new(numrows: usize, numcols: usize, nnz: usize) -> SparseMat {
        let offset: Vec<usize> = vec![0; numrows + 1];
        let nonzero: Vec<usize> = vec![0; nnz];
        SparseMat {
            numrows,
            numcols,
            nnz,
            offset,
            nonzero,
            value: None,
        }
    }

    /// a new sparse matrix with values
    pub fn new_with_values(numrows: usize, numcols: usize, nnz: usize) -> SparseMat {
        let offset: Vec<usize> = vec![0; numrows + 1];
        let nonzero: Vec<usize> = vec![0; nnz];
        let value: Option<Vec<f64>> = Some(vec![0.0; nnz]);
        SparseMat {
            numrows,
            numcols,
            nnz,
            offset,
            nonzero,
            value,
        }
    }

    /// inspector fn for numrows
    pub fn numrows(&self) -> usize {
        self.numrows
    }

    /// inspector fn for numcols
    pub fn numcols(&self) -> usize {
        self.numcols
    }

    /// give a sparse matrix random (uniform [0,1]) values
    pub fn randomize_values(&mut self) {
        let mut value = Vec::new();
        let mut rng = rand::thread_rng();
        for _ in 0..self.nnz {
            value.push(rng.gen::<f64>());
        }
        self.value = Some(value);
    }

    /// returns an iterator over row counts, very useful in this library
    pub fn rowcounts<'a>(&'a self) -> impl Iterator<Item = usize> + 'a {
        self.offset[0..self.numrows]
            .iter()
            .zip(&self.offset[1..self.numrows + 1])
            .map(|(a, b)| b - a)
    }

    /// dumps a sparse matrix to a file in a ASCII format
    /// # Arguments
    /// * maxrows the number of rows that are written, 0 means everything,
    ///           otherwise write the first and last maxrows/2 rows
    /// * filename the filename to written to
    pub fn dump(&self, maxrows: usize, filename: &str) -> Result<()> {
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
            for off in &self.offset[start_row..self.numrows + 1] {
                write!(file, "{} ", off)?;
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

    /// prints some stats of a sparse matrix
    pub fn stats(&self) -> () {
        println!("    numrows  = {}", self.numrows);
        println!("    numcols  = {}", self.numcols);
        println!("    nnz      = {}", self.nnz);
        if let Some(_) = self.value {
            println!("    matrix with values");
        } else {
            println!("    matrix pattern only");
        }

        // compute min, max, sum all at once for efficiency, only once thru itereator
        let (mindeg, maxdeg, sumdeg) = self.rowcounts().fold((self.numcols, 0, 0), |acc, x| {
            (acc.0.min(x), acc.1.max(x), acc.2 + x)
        });

        let avgdeg = sumdeg as f64 / self.numrows as f64;

        println!(
            "    min, avg, max degree = {}, {}, {}",
            mindeg, avgdeg, maxdeg
        );
    }

    /// writes a sparse matrix to a file in a MatrixMarket ASCII formats
    /// # Arguments
    /// * filename the filename to written to
    pub fn write_mm_file(&self, filename: &str) -> Result<()> {
        let path = Path::new(&filename);
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        self.write_mm(&mut writer)
    }

    /// write mm to any writer
    pub fn write_mm<W>(&self, writer: &mut W) -> Result<()>
    where
        W: Write,
    {
        if let Some(value) = &self.value {
            writeln!(writer, "%%MatrixMarket matrix coordinate real general")?;
            writeln!(writer, "{} {} {}", self.numrows, self.numcols, self.nnz)?;
            for i in 0..self.numrows {
                for k in self.offset[i]..self.offset[i + 1] {
                    writeln!(writer, "{} {} {}", i + 1, self.nonzero[k] + 1, value[k])?;
                }
            }
        } else {
            writeln!(writer, "%%MatrixMarket matrix coordinate pattern general")?;
            writeln!(writer, "{} {} {}", self.numrows, self.numcols, self.nnz)?;
            for i in 0..self.numrows {
                for nz in &self.nonzero[self.offset[i]..self.offset[i + 1]] {
                    writeln!(writer, "{} {}", i + 1, nz + 1)?;
                }
            }
        }
        Ok(())
    }

    /// read a sparse matrix from a file in a MatrixMarket ASCII format
    /// # Arguments
    /// * filename the file to be read
    pub fn read_mm_file(filename: &str) -> Result<SparseMat> {
        let path = Path::new(&filename);
        let fp = File::open(path)?;
        let mut reader = BufReader::new(fp);
        SparseMat::read_mm(&mut reader)
    }

    /// read from any BufRead
    pub fn read_mm<R>(reader: &mut R) -> Result<SparseMat>
    where
        R: BufRead,
    {
        let line = reader.lines().next().unwrap_or(Ok("".to_string()))?;
        let re1 = Regex::new(r"^%%MatrixMarket *matrix *coordinate *pattern*")?;
        let re2 = Regex::new(r"^%%MatrixMarket *matrix *coordinate *real*")?;
        if !re1.is_match(&line) && !re2.is_match(&line) {
            return Err(Sparsemat(ParseMmError::new(format!(
                "invalid header {}",
                line
            ))));
        }
        let has_values = re2.is_match(&line);

        let line = reader.lines().next().unwrap_or(Ok("".to_string()))?;
        let re = Regex::new(r"^(\d*)\s*(\d*)\s*(\d*)\s*")?;
        let cleaned = re.replace_all(&line, "$1 $2 $3");
        let inputs: Vec<String> = cleaned.split(" ").map(|x| x.to_string()).collect();
        let nr: usize = inputs[0].parse()?;
        let nc: usize = inputs[1].parse()?;
        let nnz: usize = inputs[2].parse()?;
        if nr == 0 || nc == 0 || nnz == 0 {
            return Err(Sparsemat(ParseMmError::new(format!(
                "bad matrix sizes {} {} {}",
                nr, nc, nnz
            ))));
        }

        #[derive(Debug)]
        struct Elt {
            row: usize,
            col: usize,
            val: f64,
        }

        let mut elts: Vec<Elt> = vec![];
        for line in reader.lines() {
            let line = line.expect("read error");
            let inputs: Vec<&str> = line.split_whitespace().collect();
            let row: usize = inputs[0].parse()?;
            let col: usize = inputs[1].parse()?;
            let val: f64 = if has_values { inputs[2].parse()? } else { 1.0 };
            if row == 0 || row > nr || col == 0 || col > nc {
                return Err(Sparsemat(ParseMmError::new(format!(
                    "bad matrix nonzero coordinates {} {}",
                    row, col
                ))));
            }
            elts.push(Elt {
                row: row - 1,
                col: col - 1,
                val: val,
            }); // MatrixMarket format is 1-up, not 0-up
        }
        if elts.len() != nnz {
            return Err(Sparsemat(ParseMmError::new(format!(
                "incorrect number of nonzeros {} != {}",
                elts.len(),
                nnz
            ))));
        }

        elts.sort_by(|a, b| {
            if a.row == b.row {
                a.col.cmp(&b.col)
            } else {
                a.row.cmp(&b.row)
            }
        });

        let mut ret = SparseMat::new(nr, nc, nnz);
        let mut values: Vec<f64> = vec![];
        let mut row: usize = 0;
        for (i, elt) in elts.iter().enumerate() {
            // first adjust so we are on correct row
            while row < elt.row {
                row += 1;
                ret.offset[row] = i;
            }
            assert!(elt.row == row);
            ret.nonzero[i] = elt.col;
            if has_values {
                values.push(elt.val);
            }
        }
        for i in row + 1..ret.numrows + 1 {
            ret.offset[i] = ret.nnz;
        }
        if has_values {
            ret.value = Some(values);
        }
        Ok(ret)
    }

    /// apply row and column permutations to a sparse matrix
    /// # Arguments
    /// * rperminv pointer to the global array holding the inverse of the row permutation
    /// *cperminv pointer to the global array holding the inverse of the column permutation
    ///     rperminv[i] = j means that row i of A goes to row j in matrix Ap
    ///     cperminv[i] = j means that col i of A goes to col j in matrix Ap
    pub fn permute(&self, rperminv: &Perm, cperminv: &Perm) -> SparseMat {
        let mut ap = SparseMat::new(self.numrows, self.numcols, 0);
        let rperm = rperminv.inverse();

        if let Some(_) = &self.value {
            todo!()
        }

        // fill in permuted rows with permuted columns
        for i in 0..ap.numrows {
            let row = rperm.perm[i];
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                ap.nonzero.push(cperminv.perm[*nz]);
            }
            ap.offset[i + 1] = ap.nonzero.len();
        }
        ap.nnz = ap.nonzero.len();
        ap
    }

    /// produce the transpose of a sparse matrix
    pub fn transpose(&self) -> SparseMat {
        let mut colcounts: Vec<usize> = vec![0; self.numcols + 1];
        // histogram the column counts of A into colcounts
        for nz in &self.nonzero[0..self.nnz] {
            colcounts[*nz] += 1;
        }

        if let Some(_) = &self.value {
            todo!()
        }

        let mut at = SparseMat::new(self.numcols, self.numrows, self.nnz);

        // Build new offsets, convert colcounts to cumulative
        for i in 0..at.numrows {
            at.offset[i + 1] = at.offset[i] + colcounts[i];
            colcounts[i] = at.offset[i]
        }
        colcounts[at.numrows] = at.nnz;

        // redistribute the nonzeros
        //This is still tidy, right?
        for row in 0..self.numrows {
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                at.nonzero[colcounts[*nz]] = row;
                colcounts[*nz] += 1;
            }
        }
        at
    }

    /// checks that the sparse matrix is (strictly, i.e. zero on diagonal) lower triangluar
    pub fn is_lower_triangular(&self) -> bool {
        for row in 0..self.numrows {
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                if *nz > row {
                    return false;
                }
            }
        }
        true
    }

    /// checks that the sparse matrix is (strictly) upper triangluar
    pub fn is_upper_triangular(&self) -> bool {
        for row in 0..self.numrows {
            for nz in &self.nonzero[self.offset[row]..self.offset[row + 1]] {
                if *nz < row {
                    return false;
                }
            }
        }
        true
    }

    /// sort the non-zeros in each row of a sparse matrix
    pub fn sort_nonzeros(&mut self) -> () {
        if let Some(_) = &self.value {
            todo!()
        }
        for row in 0..self.numrows {
            self.nonzero[self.offset[row]..self.offset[row + 1]].sort_by(|a, b| a.cmp(b));
        }
    }

    /// compare the structs that hold two sparse matrices
    /// # Arguments
    /// * rmat pointer to the right sparse matrix
    pub fn compare(&self, rmat: &SparseMat) -> bool {
        if self.value == None && rmat.value != None {
            println!("(self.value = None)  != (rmat.value = Some(_))",);
        } else if self.value != None && rmat.value == None {
            println!("(self.value = Some(_))  != (rmat.value = None)",);
        } else if self.numrows != rmat.numrows {
            println!(
                "(self.numrows = {})  != (rmat.numrows = {})",
                self.numrows, rmat.numrows
            );
        } else if self.numrows != rmat.numrows {
            println!(
                "(self.numrows = {})  != (rmat.numrows = {})",
                self.numrows, rmat.numrows
            );
        } else if self.numcols != rmat.numcols {
            println!(
                "(self.numcols = {})  != (rmat.numcols = {})",
                self.numcols, rmat.numcols
            );
        } else if self.nnz != rmat.nnz {
            println!("(self.nnz = {})  != (rmat.nnz = {})", self.nnz, rmat.nnz);
        } else if self.nnz != rmat.nnz {
            println!("(self.lnnz = {}) != (rmat.lnnz = {})", self.nnz, rmat.nnz);
        } else if self.offset[0] != 0 || rmat.offset[0] != 0 || (self.offset[0] != rmat.offset[0]) {
            println!(
                "(self.offset[0] = {})  != (rmat.offset[0] = {})",
                self.offset[0], rmat.offset[0]
            );
        } else {
            for row in 0..self.numrows {
                if self.offset[row + 1] != rmat.offset[row + 1] {
                    println!(
                        "(self.offset[{}] = {})  != (rmat.offset[{}] = {})",
                        row + 1,
                        self.offset[row + 1],
                        row + 1,
                        rmat.offset[row + 1]
                    );
                    return false;
                }
            }

            for i in 0..self.nnz {
                if self.nonzero[i] != rmat.nonzero[i] {
                    println!(
                        "(self.nonzero[{}] = {})  != (rmat.nonzero[{}] = {})",
                        i, self.nonzero[i], i, rmat.nonzero[i]
                    );
                    return false;
                }
            }

            if let Some(sval) = &self.value {
                if let Some(rval) = &rmat.value {
                    for i in 0..self.nnz {
                        if (sval[i] - rval[i]).abs() > 10.0 * f64::EPSILON * sval[i].abs() {
                            println!(
                                "(self.value[{}] = {})  != (rmat.value[{}] = {})",
                                i, sval[i], i, rval[i]
                            );
                            return false;
                        }
                    }
                }
            }
            // Got to end, all good!
            return true;
        }
        // Something went wrong in one of the if clauses above, return false
        false
    }

    // TODO: this whole part for generating random matrices needs to be reworked.
    // - It really can only make square matrices (cause we only cared about graphs)
    ///  Generates the upper or lower half of the adjacency matrix for an Erdos-Renyi random graph.
    ///   This subroutine uses ALG1 from the paper "Efficient Generation of Large Random Networks"
    ///   by Batageli and Brandes appearing in Physical Review 2005. Instead of flipping a coin for each potential edge
    ///   this algorithm generates a sequence of "gaps" between 1s in the upper or lower triangular portion of the
    ///   adjancecy matrix using a geometric random variable.
    ///
    /// # Arguments
    /// * numrows The total number of vertices in the graph (rows in the matrix).
    /// * p The probability that each non-loop edge is present.
    /// * mode  enum ER_TRIANGLE:
    ///            ER_TRI_L, ER_TRI_U, ER_TRI_LWD, ER_TRI_UWD for
    ///             lower triangle, upper triangle, lower with diagonal, upper with diagonal, resp.
    /// * seed A random seed.
    pub fn erdos_renyi_tri(
        numrows: usize,
        p: f64,
        lower: bool,
        diag: bool,
        _seed: i64,
    ) -> SparseMat {
        let l_max = (std::u32::MAX as f64).ln();
        let d = (1.0 - p).ln();

        let mut mat = SparseMat::new(numrows, numrows, 0);
        let mut row = 0;
        let mut col = if lower { 0 } else { 1 };

        if diag && !lower {
            mat.nonzero.push(0);
        }
        while row < numrows {
            let mut r = rand::random::<u32>();
            if r == std::u32::MAX {
                r -= 1;
            }

            col += (1.0 + ((((std::u32::MAX - r) as f64).ln() - l_max) / d).floor()) as usize;

            while (col >= (if lower { row } else { numrows })) && row < numrows {
                if lower {
                    if diag {
                        mat.nonzero.push(row);
                    }
                    col = col - row;
                    row += 1;
                    mat.offset[row] = mat.nonzero.len();
                } else {
                    row += 1;
                    col = row + 1 + col - numrows;
                    mat.offset[row] = mat.nonzero.len();
                    if (row < numrows) && diag {
                        mat.nonzero.push(row);
                    }
                }
            }
            if row < numrows {
                mat.nonzero.push(col);
            }
        }
        mat.offset[row] = mat.nonzero.len();
        mat.nnz = mat.nonzero.len();
        mat
    }

    /// Generates the upper or lower half of the adjacency matrix for an Erdos-Renyi random graph.
    ///  This is here mostly for illustration cause it is so slow.
    ///  It flips a coin for each possible edge
    /// # Arguments
    /// * numrows The total number of vertices in the graph (rows in the matrix).
    /// * p The probability that each non-loop edge is present.
    /// * mode  enum ER_TRIANGLE:
    ///             ER_TRI_L, ER_TRI_U, ER_TRI_LWD, ER_TRI_UWD for
    ///             lower triangle, upper triangle, lower with diagonal, upper with diagonal, resp.
    /// * seed A random seed.
    pub fn naive_erdos_renyi_tri(
        numrows: usize,
        p: f64,
        lower: bool,
        diag: bool,
        _seed: f64,
    ) -> SparseMat {
        let p32: u32 = (p * std::u32::MAX as f64) as u32;
        let mut mat = SparseMat::new(numrows, numrows, 0);
        let numcols = numrows;

        for row in 0..numrows {
            if lower {
                for col in 0..row {
                    if rand::random::<u32>() < p32 {
                        mat.nonzero.push(col);
                    }
                }
                if diag {
                    mat.nonzero.push(row);
                }
            } else {
                if diag {
                    mat.nonzero.push(row);
                }
                for col in row + 1..numcols {
                    if rand::random::<u32>() < p32 {
                        mat.nonzero.push(col);
                    }
                }
            }
            mat.offset[row + 1] = mat.nonzero.len();
        }
        mat.nnz = mat.nonzero.len();
        mat
    }

    /// Generates an Erdos-Renyi random graph
    ///  Puts together appropriate triangluar matrices created by erdos_renyi_tri
    /// # Arguments
    /// * numrows The total number of vertices in the graph (rows in the matrix).
    /// * p The probability that each non-loop edge is present.
    /// * mode  enum ER_GRAPH:
    ///              ER_GRAPH_SIMPLE, ER_GRAPH_DIRECT for a simple (undirected) graph, simple directed graph, resp
    /// * seed A random seed.
    pub fn erdos_renyi_graph(numrows: usize, p: f64, simple: bool, seed: i64) -> SparseMat {
        let u = SparseMat::erdos_renyi_tri(numrows, p, false, false, seed);
        let l = if simple {
            u.transpose()
        } else {
            SparseMat::erdos_renyi_tri(numrows, p, true, false, seed)
        };
        let mut mat = SparseMat::new(numrows, numrows, 0);

        for row in 0..numrows {
            for nz in &l.nonzero[l.offset[row]..l.offset[row + 1]] {
                mat.nonzero.push(*nz);
            }
            for nz in &u.nonzero[u.offset[row]..u.offset[row + 1]] {
                mat.nonzero.push(*nz);
            }
            mat.offset[row + 1] = mat.nonzero.len();
        }
        mat.nnz = mat.nonzero.len();
        mat
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_perm1() {
        let perm = Perm::new(4);
        assert_eq!(true, perm.is_perm());
    }

    #[test]
    fn is_perm2() {
        let mut perm = Perm::new(4);
        perm.perm[1] = 3;
        assert_eq!(false, perm.is_perm());
    }

    #[test]
    fn is_perm3() {
        let mut perm = Perm::new(4);
        perm.perm[1] = 4;
        assert_eq!(false, perm.is_perm());
    }

    #[test]
    fn rand_perm1() {
        let perm = Perm::random(100, 0);
        assert_eq!(true, perm.is_perm());
    }

    #[test]
    fn read_mm1() {
        assert!(SparseMat::read_mm_file("not_there.mm").is_err());
    }

    #[test]
    fn read_mm2() {
        let mat = SparseMat::read_mm_file("good/10_ring.mm");
        match mat {
            Ok(_m) => (),
            Err(e) => panic!("{}", e),
        }
    }

    #[test]
    fn read_mm3() {
        let mat = SparseMat::read_mm_file("good/10_ring_values.mm");
        match mat {
            Ok(_m) => (),
            Err(e) => panic!("{}", e),
        }
    }

    #[test]
    fn write_mm1() {
        let f1 = format!("{}/10_ring1.mm", std::env::temp_dir().to_str().unwrap());
        let mat = SparseMat::read_mm_file("good/10_ring.mm").expect("failed read");
        mat.write_mm_file(&f1).expect("failed write");
        let mat1 = SparseMat::read_mm_file(&f1).expect("failed read");
        assert_eq!(mat.compare(&mat1), true);
    }

    #[test]
    fn write_mm2() {
        let f1 = format!(
            "{}/10_ring1_values.mm",
            std::env::temp_dir().to_str().unwrap()
        );
        let mat = SparseMat::read_mm_file("good/10_ring_values.mm").expect("failed read");
        mat.write_mm_file(&f1).expect("failed write");
        let mat1 = SparseMat::read_mm_file(&f1).expect("failed read");
        assert_eq!(mat.compare(&mat1), true);
    }

    #[test]
    fn dump1() {
        let f1 = format!("{}/10_ring.dump", std::env::temp_dir().to_str().unwrap());
        let mat = SparseMat::read_mm_file("good/10_ring.mm").expect("failed read");
        mat.dump(mat.numrows, &f1).expect("failed dump");
        // How to test this is right?
    }

    #[test]
    fn compare1() {
        let mat = SparseMat::read_mm_file("good/10_ring.mm").expect("failed read");
        assert_eq!(mat.compare(&mat), true);
        let mat1 = SparseMat::read_mm_file("good/10_ring_values.mm").expect("failed read");
        assert_eq!(mat.compare(&mat1), false);
    }

    #[test]
    fn triangular1() {
        let mat = SparseMat::read_mm_file("good/10_ring.mm").expect("failed read");
        assert_eq!(mat.is_lower_triangular(), true);
        assert_eq!(mat.is_upper_triangular(), false);
    }

    #[test]
    fn transpose1() {
        let mat = SparseMat::read_mm_file("good/10_ring.mm").expect("failed read");
        let tmat = mat.transpose();
        let ttmat = tmat.transpose();
        println!("mat:{:?}", mat);
        println!("tmat:{:?}", tmat);
        println!("ttmat:{:?}", ttmat);
        assert_eq!(mat.compare(&tmat), false);
        assert_eq!(mat.compare(&ttmat), true);
    }

    #[test]
    #[should_panic]
    fn transpose2() {
        let mat = SparseMat::read_mm_file("good/10_ring_values.mm").expect("failed read");
        let _ = mat.transpose();
    }

    #[test]
    fn perm1() {
        let mat = SparseMat::read_mm_file("good/10_ring.mm").expect("failed read");
        let rperm = Perm::random(10, 0);
        let irperm = rperm.inverse();
        let cperm = Perm::random(10, 0);
        let icperm = cperm.inverse();
        let pmat = mat.permute(&rperm, &cperm);
        let ppmat = pmat.permute(&irperm, &icperm);
        println!("rperm:{:?}", rperm);
        println!("irperm:{:?}", irperm);
        println!("cperm:{:?}", cperm);
        println!("icperm:{:?}", icperm);
        println!("mat:{:?}", mat);
        println!("pmat:{:?}", pmat);
        println!("ppmat:{:?}", ppmat);
        assert_eq!(pmat.compare(&mat), false);
        assert_eq!(ppmat.compare(&mat), true);
    }

    #[test]
    #[should_panic]
    fn perm2() {
        let mat = SparseMat::read_mm_file("good/10_ring_values.mm").expect("failed read");
        let rperm = Perm::random(10, 0);
        let cperm = Perm::random(10, 0);
        let _ = mat.permute(&rperm, &cperm);
    }

    #[test]
    fn rand_mat_tri1() {
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, true, false, 0);
        assert_eq!(mat.is_lower_triangular(), true);
        assert_eq!(mat.is_upper_triangular(), false);
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, true, true, 0);
        assert_eq!(mat.is_lower_triangular(), true);
        assert_eq!(mat.is_upper_triangular(), false);
    }

    #[test]
    fn rand_mat_tri2() {
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, false, false, 0);
        assert_eq!(mat.is_lower_triangular(), false);
        assert_eq!(mat.is_upper_triangular(), true);
        let mat = SparseMat::erdos_renyi_tri(100, 0.05, false, true, 0);
        assert_eq!(mat.is_lower_triangular(), false);
        assert_eq!(mat.is_upper_triangular(), true);
    }

    #[test]
    fn rand_mat_graph1() {
        let mat = SparseMat::erdos_renyi_graph(100, 0.05, false, 0);
        assert_eq!(mat.is_lower_triangular(), false);
        assert_eq!(mat.is_upper_triangular(), false);
        let mat = SparseMat::erdos_renyi_graph(100, 0.05, true, 0);
        assert_eq!(mat.is_lower_triangular(), false);
        assert_eq!(mat.is_upper_triangular(), false);
    }
}
