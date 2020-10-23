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
use PeekPoke;
use Time;
use spmat;

// How data is laid out
// nonzero array is BLOCK distributed
// ROWS are cyclic distributed among locales
// offset array is cyclic distributed
// rows are BLOCK distributed amongs tasks (by convention)


// QUESTIONS, thoughts, suggestions:
// 1. Review iterator story... in transpose and permute_matrix I could use an iterator that iterates over all nonzeros properly (done)
// 2. Learn more about returning a tuple (toposort) (done)
// 3. it would be nice if you could demote atomics so that they were no longer atomic and could be treated normally
// 4. Is there a way to do a global (+ reduce x) in the middle of a SPMD like loop?
// 5. I think it would be nice to have a better way to get task locale subdomains

const numTasksPerLocale = if dataParTasksPerLocale > 0 then dataParTasksPerLocale
                                                       else here.maxTaskPar;
const numTasks = numLocales * numTasksPerLocale;




proc toposort(A: Sparsemat) {
  var th_per_locale = numTasksPerLocale;
  var rp = newCyclicArr(0..#A.nr, int);
  var cp = newCyclicArr(0..#A.nc, int);
  var pos: atomic int;
  pos.write(A.nr-1);
  
  var At = A.transpose();

  var D = newCyclicDom(0..#A.nr);
  var rowstat: [D] atomic int;

  // initialize rowstat array (for each row it tracks rowcnt and sum of live col indices)
  forall row in A.rows(){
    var sum = 0;
    var cnt = 0;
    for col in A.row_nz(row){
      sum += col;
      cnt += 1;
    }
    rowstat[row].add(sum << 32 | cnt);
  }

  // toposort main
  coforall l in Locales do on l{
      coforall t in 0..#th_per_locale{
        var my_rank = l.id*th_per_locale + t;
        var rows_left = A.num_rows_per_task[my_rank];
        while(rows_left > 0) {
          for row in A.my_task_rows(my_rank){
            var info = rowstat[row].read();
            if((info & 0x0000FFFF) == 1){
              var col = info >> 32;
                
              // reserve place in rp and cp              
              var newpos = pos.fetchAdd(-1);
              rp[newpos] = row;
              cp[newpos] = col;
                
              // update stats for the rows in this column
              for row2 in At.row_nz(col){
                rowstat[row2].sub(info);
              }
              rows_left -= 1;
            }
          }
          chpl_task_yield();
        }
      }    
    }

    /* var b = new Barrier(numLocales*th_per_locale); */
  
    /* // populate the queue with degree one rows */
    /* coforall l in Locales do on l{ */
    /*     coforall t in 0..#th_per_locale{         */
    /*       var my_rank = l.id*th_per_locale + t; */
        
    /*       for row in A.my_task_rows(my_rank){ */
    /*         var sum = 0; */
    /*         var cnt = 0; */
    /*         for col in A.row_nz(row){ */
    /*           sum += col; */
    /*           cnt += 1; */
    /*         } */
    /*         rowstat[row].add(sum << 32 | cnt); */
    /*       } */

    /*       b.barrier(); */

    /*       var rows_left = A.num_rows_per_task[my_rank]; */
      
    /*       while(rows_left > 0) { */
    /*         for row in A.my_task_rows(my_rank){ */
    /*           var info = rowstat[row].read(); */
    /*           if((info & 0x0000FFFF) == 1){ */
    /*             var col = info >> 32; */
              
    /*             // reserve place in rp and cp               */
    /*             var newpos = pos.fetchAdd(-1); */
    /*             rp[newpos] = row; */
    /*             cp[newpos] = col; */
              
    /*             // update stats for the rows in this column */
    /*             for row2 in At.row_nz(col){ */
    /*               rowstat[row2].sub(info);  */
    /*             } */
    /*             rows_left -= 1; */
    /*           } */
    /*         } */
    /*         chpl_task_yield(); */
    /*       } */
    /*     } */
    /*   } */
  
  return(rp, cp);
}

proc toposort2(A: Sparsemat){

  var th_per_locale = numTasksPerLocale;
  var rp = newCyclicArr(0..#A.nr, int);
  var cp = newCyclicArr(0..#A.nc, int);
  var pos: atomic int;
  pos.write(A.nr-1);
  
  var At = A.transpose();

  var D = newCyclicDom(0..#A.nr);
  var rowstat: [D] atomic int;

  // initialize rowstat array (for each row it tracks rowcnt and sum of live col indices)
  forall row in A.rows(){
    var sum = 0;
    for col in A.row_nz(row){
      sum += (A.nc + col);
    }
    rowstat[row].write(sum);
  }

  // toposort main
  coforall l in Locales do on l{
      coforall t in 0..#th_per_locale{
        var my_rank = l.id*th_per_locale + t;
        var rows_left = A.num_rows_per_task[my_rank];
        var queue: [0..#rows_left] int;
        var nc = A.nc;
        var done_val = nc*10;
        while(rows_left > 0) {
          
          var start = 0, end = 0;
          var d1 = 0;
          for row in A.my_task_rows(my_rank){            
            var info = rowstat[row].read();
            if(info < 2*nc){
              d1+=1;
              queue[end] = row;
              end += 1;
            }
          }

          var newpos = pos.fetchSub(d1);
          rows_left -= d1;
          
          while(start < end){
            var row = queue[start];
            start+=1;
            var info = rowstat[row].read();
            var col = info - nc;
            rp[newpos] = row;
            cp[newpos] = col;
            newpos -= 1;
            for row2 in At.row_nz(col){
              rowstat[row2].sub(info);
            }
            rowstat[row].write(done_val);//so we don't look at this row again
          }
          chpl_task_yield();
        }
      }    
    }
  return(rp, cp);
}

config const detail = 0;
config const timers = true;

var t: Timer;
proc printTimer(section) {
  if !timers then return;
  t.stop();
  writeln(section, " time ", t.elapsed() ," secs");
  t.clear();
  t.start();
}
/*****************************************************************/
/*                           MAIN                                */
/*****************************************************************/

config const nPerTask = 10;
config const N = nPerTask * numTasks;
//config const p = 0.1;
config const Z = 30.0;
config const seed = 1038;

var p: real;
p = (2.0*Z)/(N-1);
if(p > 1.0){p = 1.0;}

t.start();
/* generate random graph (it will be lower-triagular)*/
var At = generate_erdos_renyi_half(N, p, seed, true);

/* transpose At to get upper-triangular matrix */
var A = At.transpose();
A.print(detail);
printTimer("Matrix generation");

/* get random row and columns permutations */
var rowperm = rand_permp(N, seed*2);
var colperm = rand_permp(N, seed*3);
if(detail > 1){
  writeln("rowperm = ", rowperm);
  writeln("colperm = ", colperm);
}
printTimer("Permutation generation");

/* permute the matrix */
var Ap = permute_matrix(A, rowperm, colperm);
Ap.print(detail);
printTimer("Permute Matrix");

/* run toposort! */
var (rp, cp) = toposort(Ap);
if(detail > 1){
  writeln("topo rowperm = ",rp);
  writeln("topo colperm = ",cp);
}
printTimer("Toposort");

/* verify the result */
var App = permute_matrix(Ap, rp, cp);
App.print(detail);
printTimer("Permute Matrix(2)");

if(App.is_unit_upper_triangular()){
  writeln("Success!");
}else{
  writeln("Fail!");
 }
//for nz in App{
//writeln("this");
//}
