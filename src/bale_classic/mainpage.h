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
/*!
 * \defgroup libgetputgrp Libgetput library (libgetput)
 * This is a simple parallel macro library that exposes very few
 * functions to the users. The most important are the 'get' and 'put' functions.
 *
 * \defgroup exstackgrp exstack library (exstack)
 * Our first attempts at buffered communications librarys. Contains both synchronous
 * (exstack) and asynchronous (exstack2) modes.
 *  
 * \defgroup spmatgrp Sparse Matrix library (spmat)
 * This is a simple parallel sparse matrix library built on top of libgetput. It
 * implements a distributed compressed sparse row data structure.

 */

/*! \mainpage 
 *
 ******************************************************************************************
 
 \section IntroToBale Why bale?  
 The bale effort is, first and
 foremost, <b> a vehicle for discussion</b>.  Central to bale is a
 directory of "apps" that exhibit interesting communication patterns
 and programming demands. They allow us to test and exchange ideas
 with others. We hope the bale discussion can lead to improved
 programmer productivity (including existing and/or new programming
 models) and performance.

 - Part of bale is to keep revisiting the question of ``How would
 we \e like to write this code?''.
 
 - We think all of us want to write programs that are simple, elegant,
 and look like the algorithms they are implementing. In bale, we call
 this (obviously subjective) ideal the AGI model (``As God
 Intended''). For each app, we have included our current best effort
 at an AGI implementation -- this usually makes heavy use of single
 word gets and puts.
 
 - We want something worthy of being called AGI, but with an
 acceptable fraction of the ``gold standard'' (highest possible
 performing) implementation. bale purposely does not define
 "acceptable" or "performance". We want to shrink the gap and
 explore the continuim between high productivity and high
 performance. <b>bale is not meant to be a benchmark</b>.
 
 bale is released under the 3-clause BSD license to try to
 encourage others to join this discussion. Please send correspondence to bale@super.org.
 
 
\section WhatsInBale What is in bale?
  The bale archive is a collection of applications and support libraries that have been
  written to study programming models for HPC system as applied to irregular problems.

  The bale archive currently contains the following libraries and applications packages.

  Apps
  - \subpage histogram_page simple random puts (remote writes)
  - \subpage indexgather_page simple random gets (remote reads)
  - \subpage toposort_page a bfs algorithm to re-order the rows and columns of  "morally" upper triangular matrix
  - \subpage transpose_matrix_page sparse matrix manipulation that is interesting in it own right
  - \subpage permute_matrix_page sparse matrix manipulation that is interesting in it own right
  - \subpage randperm_page a parallel algorithm to find a large permutation in parallel
  - \subpage triangles_page counting triangles in a graph
  - \subpage sssp_page Single Source Shortest Path in a graph.
  - \subpage sparse_matrix_io_page Write and Read a sparse matrix in a binary form.

  Libraries
  - \ref libgetputgrp a set of macros and functions that can be built on UPC or SHMEM
   to supports basic PGAS functions (like get and put) as well as some atomic operations.  
   It allows for basic PGAS programming in a UPC vs SHMEM agnostic way.
  - \ref exstackgrp Our first attempt at a buffered communications library was called exstack.  
  It uses barriers and forces one to write in a bulk-synchronous programming style. 
  Our second attempt was exstack2, which allowed for more asynchronous programming. 
    Both models are included in this library.
  - \b convey: The successor to the exstack libraries. Conveyors should offer better performance and more flexibility.
  It can be used in either a bulk-synchronous or asynchronous style.
  - \ref spmatgrp A parallel sparse matrix library that supports the creation and manipulation of {0,1}-matrices,
  often call "positional" matrices.
  
  Also
  - INSTALL: Instructions on how to get bale up and running.
  - make_bale: A python script to make compiling bale easier.
  - run_apps.py: A python script to run tests on applications.

  \subsection History History

  A number of irregular problems have a natural PGAS description.
  However, this model often requires the use of fine-grain,
  low-latency, and sometimes atomic, references to remote memory.
  Most current and anticipated HPC platform achieve their optimal
  global communication bandwidth only by sending large packages.  This
  work began with the development of a library (exstack) to help
  relieve the programming burden of managing remote memory operation
  using buffered communication.
  
  Exstack uses a bulk synchronous model: All threads fill their
  outgoing buffers, then barrier and exchange the buffers, then
  barrier and empty their incoming buffers.  Using exstack did greatly
  increase performance---sometimes achieving near optimal network
  injection bandwidth.  However, asynchronous PGAS code with inserted
  blocks of bulk synchronization code is ugly, hard to follow and hard
  to maintain.  We also believe that programs with an irregular
  program flow, like toposort, incur a performance penalty if they are
  not asynchronous.

  We hoped we could improve exstack by making the process of
  exchanging the buffers asynchronous. This resulted in two separate
  and simultaneous efforts: exstack2 and conveyors (conveyors unifies
  the APIs for, and adds many features to, exstack and exstack2).  The
  results were interesting because they led to more questions than
  answers. Fairly quickly it became clear that conveyors were superior
  to exstack2 from a features perspective and we stopped development
  of exstack2. We include it in bale for historical purposes -- its
  performance can be quite good in some applications but there are
  some applications (especially in SHMEM) where it performs terribly
  and we have not put the energy into finding out why yet. There are
  also applications (such as randperm and sparse_matrix_io) which
  are hard to implement in exstack2. Features like steady-mode and
  giant-messages in coveyors makes some applications easier to program
  in conveyors than exstack/2.

  As the example applications in bale demonstrate, the
  use of such libraries is an unacceptable solution to the problem of
  programming irregular problems on future HPC systems because of the
  burden it places on the programmer.

  ******************************************************************************************
 

 */

