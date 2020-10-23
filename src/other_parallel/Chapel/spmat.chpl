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

// How data is laid out
// nonzero array is BLOCK distributed
// ROWS are cyclic distributed among locales
// offset array is cyclic distributed
// rows are BLOCK distributed amongs tasks (by convention)



const numTasksPerLocale = if dataParTasksPerLocale > 0 then dataParTasksPerLocale
                                                       else here.maxTaskPar;
const numTasks = numLocales * numTasksPerLocale;

iter task_local_inds(arr, my_task_rank, num_tasks){
  for i in arr.localSubdomain(){
    if(((i/numLocales) % num_tasks) == my_task_rank){
      yield(i);
    }
  }
}

// this function is agnostic of layout (BLOCK or CYCLIC)
proc my_num_things(nthings:int, my_rank:int, total_ranks:int){
  return((nthings + total_ranks - my_rank - 1)/total_ranks);
}

proc my_first_thing(nthings:int, my_rank:int, total_ranks:int, layout:int){
  if(layout==0){//BLOCK
    return(my_rank*(nthings/total_ranks) + min(my_rank, mod(nthings,total_ranks)));
  }else{//CYCLIC
    return(my_rank);
  }
}

/* The sparse matrix class : inspired by sparsemat_t in bale_classic. */

class Sparsemat{
  // you could try to make this a record if performance suffers
  var nr: int;
  var nc: int;
  var nnz: int;
  var off_dom;
  var nnz_dom;  
  var offset : [off_dom] int; 
  var nonzero : [nnz_dom] int;
  var first_row_per_task;
  var num_rows_per_task;
  var row_stride_per_task;
  
  /******************************************************************/
  proc init(nr:int, nc:int){
    var th_per_locale = numTasksPerLocale;
    var ntasks = th_per_locale*numLocales;
    this.nr  = nr;
    this.nc = nc;
    this.off_dom = {0..0};
    this.nnz_dom = {0..0};
    this.first_row_per_task  = newCyclicArr(0..#ntasks, int);
    this.num_rows_per_task   = newCyclicArr(0..#ntasks, int);
    this.row_stride_per_task = newCyclicArr(0..#ntasks, int);
    coforall l in Locales do on l{
        var l_nr = my_num_things(nr, l.id, numLocales);
        coforall t in 0..#th_per_locale{
          var my_rank = l.id*th_per_locale + t;
          this.first_row_per_task[my_rank] = l.id + my_first_thing(l_nr, t, th_per_locale, 0)*numLocales;
          this.num_rows_per_task[my_rank] = my_num_things(l_nr, t, th_per_locale);
          this.row_stride_per_task[my_rank] = numLocales;
        }
      }
  }
  /******************************************************************/
  proc init(nr:int, nc:int, nnz_per_locale:[]int){
    var max_nnz_per_locale = (max reduce nnz_per_locale);   
    var th_per_locale = numTasksPerLocale;
    var ntasks = th_per_locale*numLocales;
    this.nr  = nr;
    this.nc = nc;
    this.nnz = (+ reduce nnz_per_locale);
    this.off_dom = newCyclicDom(0..#(nr + numLocales));
    this.nnz_dom = newBlockDom(0..#(max_nnz_per_locale*numLocales));
    this.first_row_per_task  = newCyclicArr(0..#ntasks, int);
    this.num_rows_per_task   = newCyclicArr(0..#ntasks, int);
    this.row_stride_per_task = newCyclicArr(0..#ntasks, int);
    coforall l in Locales do on l{
        var l_nr = my_num_things(nr, l.id, numLocales);
        var local_off_dom = off_dom.localSubdomain();
        coforall t in 0..#th_per_locale{
          var my_rank = l.id*th_per_locale + t;
          this.first_row_per_task[my_rank] = local_off_dom.alignedLow + my_first_thing(l_nr, t, th_per_locale, 0)*local_off_dom.stride;
          this.num_rows_per_task[my_rank] = my_num_things(l_nr, t, th_per_locale);
          this.row_stride_per_task[my_rank] = local_off_dom.stride;
        }
      }
  }
  /******************************************************************/
  proc getrowcounts(): []int{
    var rc = newCyclicArr(0..#this.nr, int);
    forall row in this.rows(){
      rc[row] = this.row_start(row + numLocales) - this.row_start(row);
    }
    return(rc);
  }
  /******************************************************************/
  proc get_offsets_from_rowcounts(rowcounts){
    //nonzero is BLOCK distributed
    coforall l in Locales do on l{
        ref lrowcounts = rowcounts.localSlice[rowcounts.localSubdomain()];
        const local_offset_inds = this.offset.localSubdomain();
        ref loffset = this.offset.localSlice[local_offset_inds((l.id + numLocales)..)];
        loffset = (+ scan lrowcounts);
	this.offset.localSlice[local_offset_inds] += this.nonzero.localSubdomain().alignedLow; //add block offset
      }    
  }
  /******************************************************************/
  proc transpose(){
    var th_per_locale = numTasksPerLocale;
    
    // get column counts
    var D = newCyclicDom(0..#this.nc);
    var colcounts: [D] atomic int;
    var cc: [D] int;

    forall (_, col) in this{
      colcounts[col].add(1);
    }

    // copy colcounts into non-atomic array
    forall i in 0..#colcounts.size{
      cc[i] = colcounts[i].peek();
    }

    // allocate transpose matrix
    var lnnz: [0..#numLocales] int;
    coforall l in Locales do on l{
        lnnz[l.id] = (+ reduce cc.localSlice[cc.localSubdomain()]);
      }
    var At = new owned Sparsemat(this.nc, this.nr, lnnz);
    At.get_offsets_from_rowcounts(cc);

    // get an temp. offset array that is atomic
    D = newCyclicDom(0..#(At.nr + numLocales));
    var tmpoffset: [D] atomic int;
    forall i in 0..#At.offset.size{
      tmpoffset[i].write(At.offset[i]);
    }

    // populate the nonzero array of At
    forall (row, col) in this{
      var pos = tmpoffset[col].fetchAdd(1);
      At.nonzero[pos] = row;
    }

    return(At);
  }

  /******************************************************************/
  /* determine if a matrix is unit-upper-triangular */
  proc is_unit_upper_triangular(): bool{
    var pivots = 0;
    forall (row, col) in this with (+ reduce pivots){
      if(col < row){pivots += this.nr+1;}
      if(col == row){pivots += 1;}
    }
    if(pivots != this.nr){return(false);}
    return(true);
  }
  
  /******************************************************************/
  proc row_start(row: int){
    return(this.offset[row]);
  }
  /******************************************************************/
  proc row_end(row: int){
    return(this.offset[row + numLocales]); // this depends on the rows being CYCLIC dist
  }
  /******************************************************************/
  /*                           ITERATORS                            */
  /******************************************************************/
  iter row_nz(row: int) ref{ //serial iterator over nonzeros in a row
    var rs = this.row_start(row);
    var re = this.row_start(row + numLocales);// this depends on the rows being CYCLIC dist
    if(re < rs){
      writeln("Hey!!!", rs, " ", re);
    }
    for i in (rs..(re - 1)){ // this depends on nonzero being BLOCK distributed
      ref nz = this.nonzero[i];
      yield(nz);
    }
  }

  /******************************************************************/
  iter my_task_rows(my_rank:int){
    for rowind in (0..#this.num_rows_per_task[my_rank]){
      yield(rowind*this.row_stride_per_task[my_rank] + this.first_row_per_task[my_rank]);
    }
  }
  /******************************************************************/

  iter these(){ //serial class iterator
    writeln("Using serial CLASS iterator!");
    for row in (0..#this.nr){
      for col in this.row_nz(row){
        yield((row,col));
      }
    }
  }
  
  /******************************************************************/
  iter rows(){ //serial row iterator
    writeln("Using serial row iterator!");
    for row in (0..#this.nr){
      yield(row);
    }
  }

  /******************************************************************/
  iter nzs() ref{ //serial nonzero (column indices) iterator
    writeln("Using serial nonzeros iterator");
    for r in this.rows(){
      for c in this.row_nz(r){        
        yield(c);
      }
    }
  }
  /******************************************************************/
  iter rows(param tag: iterKind) where tag == iterKind.standalone{
    //parallel row iterator
    var th_per_locale = numTasksPerLocale;
    coforall l in Locales do on l{
        coforall t in (0..#th_per_locale){
          var my_rank = l.id*th_per_locale + t;
          for row in this.my_task_rows(my_rank){
            yield(row);
          }
        }
        //forall row in (l.id..(this.nr-1) by numLocales){
        //yield(row);
        //}
      }
  }
  
  /******************************************************************/
  iter these(param tag: iterKind) where tag == iterKind.standalone{
    //parallel class iterator
    forall row in this.rows(){
      for col in this.row_nz(row){
        yield(row, col);
      }
    }
  }
  /*coforall l in Locales do on l{
    forall row in (l.id..(this.nr-1) by numLocales){
    for col in this.row_nz(row){
    yield(row, col);
    }
    }
    }  
  */
  /******************************************************************/
  iter nzs(param tag: iterKind) ref where tag == iterKind.standalone{
    // parallel nonzero (column indices) iterator
    coforall l in Locales do on l{
        forall row in (l.id..(this.nr-1) by numLocales){
          for col in this.row_nz(row){
            yield(col);
          }
        }
      }
  }

  
  /******************************************************************/
  proc print(detail=1){
    if(detail == 1){
      writeln("MATRIX: nr = ", this.nr, " nc = ", this.nc, " nnz = ", this.nnz);
    }else if(detail > 1){
      writeln("****************************************************");
      writeln("MATRIX: nr = ", this.nr, " nc = ", this.nc, " nnz = ", this.nnz);
      writeln("first_row_per_task: ", this.first_row_per_task);
      writeln("num_rows_per_task: ", this.num_rows_per_task);
      writeln("row_stride_per_task: ", this.row_stride_per_task);
      writeln("OFFSET: ", this.offset);
      writeln("NONZERO: ", this.nonzero);
      
      if(detail == 3){
        for row in this.rows(){
          write("row ", row, ": ");
          for col in this.row_nz(row){
            write(col, ",");
          }
          writeln();
        }
      }
      writeln("****************************************************");
    }
  }
}


// We are generating an Erdos-Renyi random graph using
// a method inspired by the paper "Efficient Generation
// of Large Random Networks" by Batagelj and Brandes (2005).
// We observe that since the wait time between nonzeros
// is geometric and since the geometric distribution is memoryless,
// we can just start each row from scratch rather than doing it
// the way the paper does it (this may lead to more random
// draws, but it makes it easier to parallelize!).
// This routine generates the lower half of the adjacency matrix

proc generate_erdos_renyi_half(const n: int, const p:real, const seed=1208, const unit_diagonal=true){
  var th_per_locale = numTasksPerLocale;
  var rowcounts = newCyclicArr(0..#n, int);  
  rowcounts = 0;
  var lnnz = newCyclicArr(0..#numLocales, int);
  coforall l in Locales do on l{
      var tot: [0..#numTasksPerLocale] int;
      ref lrc = rowcounts.localSlice[rowcounts.localSubdomain()];      
      coforall t in 0..#numTasksPerLocale{
        var my_rank = l.id*numTasksPerLocale + t;
        var rndstream = new owned RandomStream(real, seed + my_rank, parSafe=false);
        var lq = log(1 - p);
        for row in task_local_inds(rowcounts, t, numTasksPerLocale){
          var cnt = if unit_diagonal then 1 else 0;
          var col = 0;
          var bound = row-1;
          while(col < bound){
            var r = rndstream.getNext();
            col += floor(log(1 - r)/lq): int;
            if(col < bound){
              cnt += 1;
              col += 1;
            }
          }
          lrc[row] = cnt;
          tot[t] += cnt;
        }        
      }
      lnnz[l.id] = (+ reduce tot);
    }
  
  var A = new owned Sparsemat(n, n, lnnz);
  A.get_offsets_from_rowcounts(rowcounts);

  coforall l in Locales do on l{
      const local_nonzero_inds = A.nonzero.localSubdomain();
      ref lnonzero = A.nonzero.localSlice[local_nonzero_inds];
      ref loffset = A.offset.localSlice[A.offset.localSubdomain()];
      coforall t in 0..#numTasksPerLocale{
        var my_rank = l.id*numTasksPerLocale + t;
        var rndstream = new owned RandomStream(real, seed + my_rank, parSafe=false);
        var lq = log(1 - p);
        for row in task_local_inds(rowcounts, t, numTasksPerLocale){
          var nzpos = loffset[row];
          var col = 0;
          var bound = row-1;
          if(unit_diagonal){
            lnonzero[nzpos] = row; //diagonal element
            nzpos += local_nonzero_inds.stride;
          }
          while(col < bound){
            var r = rndstream.getNext();
            col += floor(log(1 - r)/lq): int;
            if(col < bound){
              lnonzero[nzpos] = col;
              nzpos += local_nonzero_inds.stride;
              col += 1;
            }          
          }
        }
      }
    }

  return(A);
  
}

/***************************************************************/
/* this routine generates a random permutaion of [0,1,...,n-1]
   using the dart-throwing algorithm from bale_classic.
*/
proc rand_permp(N:int, seed:int) : []int{
  var D = newCyclicDom(0..#(2*N));
  var target: [D] atomic int;
  var perm = newCyclicArr(0..#N, int);
  forall p in target{
    p.write(-1);
  }
  if(0){
    coforall l in Locales do on l{
        var th_per_locale = numTasksPerLocale;
        coforall t in (0..#th_per_locale){
          var my_rank = l.id*th_per_locale + t;
          var rndstream = new owned RandomStream(real, seed + my_rank, parSafe=false);
          var i = my_rank;
          var npes = th_per_locale*numLocales;
          while( i < N ){
            var r = (rndstream.getNext() * 2*N): int;
            if(target[r].compareAndSwap(-1, i)){
              i += npes;
            }
          }
        }
      }

    // count how many darts landed on each Locale
    var cnt: [0..#(numLocales+1)] int;
    coforall l in Locales do on l{
        var tmpcnt = 0;
        forall p in target.localSlice[target.localSubdomain()] with (+ reduce tmpcnt){
          if(p.read() > -1){
            tmpcnt += 1;
          }
        }
        cnt[l.id + 1] = tmpcnt;
      }

    var offset = (+ scan cnt);

    // this loop is only parallel over locales
    coforall l in Locales do on l{
        var pos = 0;
        for p in target.localSlice[target.localSubdomain()]{
          if (p.read() > -1){
            perm[offset[l.id] + pos] = p.read();
            pos += 1;
          }
        }
      }
  }else{
    // new try
    var thpl = numTasksPerLocale;
    var b = new Barrier(numLocales*thpl);
    var offset = newBlockArr(0..numLocales*thpl, int);
    coforall l in Locales do on l{
        var th_per_locale = numTasksPerLocale;
        coforall t in (0..#th_per_locale){
          var npes = numLocales*th_per_locale;
          var my_rank = l.id*th_per_locale + t;
          var rndstream = new owned RandomStream(real, seed + my_rank, parSafe=false);
          // throw the darts
          for i in my_rank..(N-1) by npes{
            do{ var r = (rndstream.getNext() * 2*N): int;
            }while(!target[r].compareAndSwap(-1, i));
          }
        
          b.barrier();

          //count the number of darts that landed in my task domain
          for i in ((l.id*th_per_locale + t)..(target.size-1) by npes){
            if(target[i].read() > -1){ offset[my_rank]+=1;}
          }
          
        }
      }
    
    offset = (+ scan offset);

    // redistribute the darts in a CYCLIC fashion to create the final permutation
    coforall l in Locales do on l{
        var th_per_locale = numTasksPerLocale;
        coforall t in (0..#th_per_locale){
          var my_rank = l.id*th_per_locale + t;
          var npes = numLocales*th_per_locale;
          
          var pos = 0;
          if(my_rank) {pos = offset[my_rank - 1];}
        
          for i in ((l.id*th_per_locale + t)..(target.size-1) by npes){
            var p = target[i].read();
            if(p > -1){
              perm[pos] = p;
              pos += 1;
            }
          }
        }
      }
  }
  if (!is_perm(perm)){
    writeln("ERROR!!!!! not a perm!");
  }
  return(perm);
}

/***************************************************************/
proc invert_perm(perm:[] int){
  var inv_perm = newCyclicArr(0..#perm.size, int);
  [i in 0..#perm.size] inv_perm[perm[i]] = i;
  return(inv_perm);
}
/***************************************************************/
/* verify that an array is a permutation. */
proc is_perm(perm) : bool{

  var N = perm.size;
  var flag = newCyclicArr(0..#N, int(8));
  flag = 0;
  var error = 0;
  forall i in perm with (+ reduce error){
    if((i >= N) || (i < 0)){
      error += 1;
    }else{
      flag[i] = 1;
    }
  }
  forall i in flag with (+ reduce error){
    if(i != 1){
      error +=1;
    }
  }
  if(error){
    return(false);
  }else{
    return(true);
  }
}

/***************************************************************/
// why can't I put ": Sparsemat" as a return type?
proc permute_matrix(A: Sparsemat, rperm: []int, cperm: []int){
  var th_per_locale = numTasksPerLocale;
  
  // get the rowcounts array for the original matrix
  var rowcounts_orig = A.getrowcounts();
  
  // get the permuted rowcounts
  var rowcounts = newCyclicArr(0..#A.nr, int);
  forall i in 0..#rperm.size{
    rowcounts[i] = rowcounts_orig[rperm[i]];
  }

  var lnnz: [0..#numLocales] int;
  coforall l in Locales do on l{
      lnnz[l.id] = (+ reduce rowcounts.localSlice[rowcounts.localSubdomain()]);
    }

  var Ap = new owned Sparsemat(A.nr, A.nc, lnnz);
  Ap.get_offsets_from_rowcounts(rowcounts);
  
  // shuffle rows  
  // I think this zipped iteration makes it so that Aprow
  // is local (even though it looks like the follower...)
  forall (Arow, Aprow) in zip(rperm, 0..#Ap.nr){
    var nzpos = Ap.offset[Aprow];
    for col in A.row_nz(Arow){
      Ap.nonzero[nzpos] = col;
      nzpos += Ap.nonzero.domain.stride;
    }
  }
  
  // invert cperm
  var cperm_inv = invert_perm(cperm);
  
  // relabel column indicies
  forall c in Ap.nzs(){
    c = cperm_inv[c];
  }
  //forall i in 0..#Ap.nonzero.size{
  //Ap.nonzero[i] = cperm_inv[Ap.nonzero[i]];
  //}
  
  return(Ap);
}