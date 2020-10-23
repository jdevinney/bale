/******************************************************************
//
//
//  Copyright(C) 2020, Institute for Defense Analyses
//  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
// 
//
//  All rights reserved.
//  
//   This file is a part of Bale.  For license information see the
//   LICENSE file in the top level directory of the distribution.
//  
// 
 *****************************************************************/ 
use Random;
use CyclicDist;
use BlockDist;
use PrivateDist;
use Barriers;
use Time;
use spmat;

//
// This is mostly a direct copy of topo.chpl with triangle counting loops
// replacing the toposort algorithm.

const numTasksPerLocale = if dataParTasksPerLocale > 0 then dataParTasksPerLocale
                                                       else here.maxTaskPar;

const numTasks = numLocales * numTasksPerLocale;

proc triangles_pretty(L: Sparsemat) {
  var th_per_locale = numTasksPerLocale;
  var rp = newCyclicArr(0..#L.nr, int);
  var cp = newCyclicArr(0..#L.nc, int);
  var pos: atomic int;
  pos.write(L.nr-1);
  
  // triangle main
  // L is the lower triangular matrix submatrix of the adjacency matrix
  // really computing +/+/( L.&(L*U)) 
  var totcnt = newCyclicArr(0..#numLocales, int);
  coforall l in Locales do on l{
    var cnt: [0..#numTasksPerLocale] int;
    coforall t in 0..#th_per_locale{
      var my_rank = l.id*th_per_locale + t;
      for row in L.my_task_rows(my_rank){    
        for col in L.row_nz(row){
          // count the intersection of the columns in row(row) and row(col) of L
          for wc in L.row_nz(col) {
            for wr in L.row_nz(row){
              if (wc == wr) {
                //writeln("tri: ",row, col, wc);
                cnt[t] += 1;
              }
              if (wr > wc ) {break;}
            }
          }
        }
      }
    }
    totcnt[l.id] = (+ reduce cnt);
  }
  var count = (+ reduce totcnt);
  return(count);
}

proc triangles(L: Sparsemat) {
  var th_per_locale = numTasksPerLocale;
  var rp = newCyclicArr(0..#L.nr, int);
  var cp = newCyclicArr(0..#L.nc, int);
  var pos: atomic int;
  pos.write(L.nr-1);
  
  // triangle main
  // L is the lower triangular matrix submatrix of the adjacency matrix
  // really computing +/+/( L.&(L*U))

  var totcnt = newCyclicArr(0..#numLocales, int);
  coforall l in Locales do on l{
    var cnt: [0..#numTasksPerLocale] int;
    coforall t in 0..#th_per_locale{
      var my_rank = l.id*th_per_locale + t;
      for row in L.my_task_rows(my_rank){
        for col in L.row_nz(row){
          // count the intersection of the columns in row(row) and row(col) of L
          var start = L.row_start(row);
          for wc in L.row_nz(col) {
            for i in start..(L.row_end(row)-1) {
              var wr = L.nonzero[i];
              if (wc == wr) {
                //writeln("row= ", row, " col= ", col, " w= ", wc);
                cnt[t] += 1;
                start = i+1;
                break;
              }
              if (wr > wc ) {
                start = i;
                break;
              }
            }
          }
        }
      }
    }
    totcnt[l.id] = (+ reduce cnt);
  }
  var count = (+ reduce totcnt);
  return(count);
}

/*****************************************************************/
/*                           MAIN                                */
/*****************************************************************/

config const N = 12;
config const p = 0.1;
config const seed = 1038;
config const detail = 3;

var t: Timer;


var A = generate_erdos_renyi_half(N, p, seed, false);
//A.print(detail);

t.start();
var numtris = triangles(A);
//var numtris2 = triangles_pretty(A);
//assert(numtris == numtris2);

t.stop();
writeln("Found ", numtris," triangles ");

writeln("Triangle counting time ", t.elapsed()," secs");t.clear();


