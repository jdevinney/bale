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
use CyclicDist, BlockDist;
use Random;
use Assert;
use Time;
config const N=10; // number of updates
config const M=10; // size of table

// allocate main table and array of random ints
const Mspace = {0..M-1};
const D = Mspace dmapped Cyclic(startIdx=Mspace.low);
var A: [D] atomic int;
const Nspace = {0..(N*numLocales - 1)};
const D2 = Nspace dmapped Block(Nspace);
var rindex: [D2] int;
var t: Timer;
t.start();

/* set up loop: populate the rindex array with random numbers mod M */
fillRandom(rindex, 208); // the 208 is a seed
forall r in rindex{
  r = mod(r, M);
}

t.stop();
writeln("Set up time: ", t.elapsed());
t.clear();
t.start();

/* In this code, we present 3 ways to write histogram in Chapel */

/* first, using a simple iterator over the array */
forall r in rindex{
  A[r].add(1); //atomic add
}

t.stop();
writeln("Loop 1: ", t.elapsed());
t.clear();
t.start();

/* We can write the main loop in a more node-centric way though. */
coforall loc in Locales do on loc do{
    forall i in rindex.localSubdomain(){
      A[rindex[i]].add(1);
    }
  }

t.stop();
writeln("Loop 2: ", t.elapsed());
t.clear();
t.start();

/* Most economical of all, we can also write the main loop in this “vector” way: */
A[rindex].add(-2);

t.stop();
writeln("Loop 3: ", t.elapsed());
t.start();

/* make sure all the updates happened correctly */
forall r in rindex{
  assert(A[r].read() == 0);
}