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

/*! \file exstack.h
 * \brief The header file for exstack library.
 * \ingroup exstackgrp
 */

#ifndef exstack_INCLUDED
#define exstack_INCLUDED /*!< std trick */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <stdarg.h>
#include <fcntl.h>
#include <assert.h>
#include <stddef.h>
#include <getopt.h>
#include<libgetput.h>

#if __UPC__
#include <upc_nb.h>
#endif

/*! \brief the struct for each of the individual buffers (or stacks) */
typedef struct exstack_buffer_t{
  int64_t count;        /*!< number of things in the buffer */
  char data[];          /*!< bytes to hold the buffer */
}exstack_buffer_t;

/*******************  Classic Exstack  ***************************************************/
/*!
\brief A structure to hold the exstack buffers and the flags used in the endgame.
\ingroup exstackgrp
*/
typedef struct exstack_t {
  int64_t bytes_per_stack;       /*!< bytes version of buf_cnt */
  SHARED char * snd_buf;         /*!< SHARED send buffers */
  SHARED char * rcv_buf;         /*!< SHARED receive buffers */
  exstack_buffer_t ** l_snd_buf; /*!< local pointer to a threads part of the send buffer */
  exstack_buffer_t ** l_rcv_buf; /*!< local pointer to a threads part of the receive buffer */
  char **fifo_ptr;               /*!< (internal) pointer to oldest work item in each rcv buffer */
  char **push_ptr;               /*!< (internal) pointer to next avail work item position in each snd buffer */
  int *put_order;                /*!< local pointer to random upc_memput ordering */
  uint64_t buf_cnt;              /*!< The number of structs to allocate for each send and recieve buffer */
  uint64_t pkg_size;             /*!< The size of each struct (work item) in bytes */
  uint64_t first_ne_rcv;         /*!< (internal) 1st possible non-empty rcv buffer */

  int64_t notify_done;           /*!< flag to say this pe has sent its done status. */
  SHARED int64_t *wait_done;     /*!< The shared array which will communicate done condition. */
  int64_t *l_wait_done;          /*!< The local pointer to wait_done. */
} exstack_t;

/* public routines */
exstack_t *exstack_init(int64_t buf_len_alloc , size_t wk_item_siz);
int64_t exstack_proceed(exstack_t *q , int done_cond);
int64_t exstack_push(exstack_t *q, void *push_item, int64_t pe);
void exstack_exchange(exstack_t *q );
int64_t exstack_pop_thread(exstack_t *q, void *pop_item, int64_t pe);
void exstack_unpop_thread(exstack_t *q , int64_t pe);
int64_t exstack_pop(exstack_t *q, void *pop_item ,  int64_t *from_pe);
void exstack_unpop(exstack_t *q);
void *exstack_pull(exstack_t *q, int64_t *from_pe);
void exstack_unpull(exstack_t *q);
int64_t exstack_min_headroom(exstack_t *q);
int64_t exstack_headroom(exstack_t *q , int64_t pe);
void exstack_clear(exstack_t *q );
void exstack_reset(exstack_t *q);

/**************************  Exstack2      ***************************************************/

/*****  Message queue format  ********/
/*
 * These macros work with the messages sent internal to exstack2 
 * A message is sent when one thread sends (does a upc_memput) a stack to another thread
 * Included in the message is the number of pkgs sent, the sending thread and 
 * whether or not this is the last stack to sent between said threads
 *
 * the format is ( (number of pkgs sent ) | ( the sending thread )  | (islast flag bit) )
 *
 * the following macros help build the message
 */
#define MSG_cnt_shift 32               /*!< 31 bits for the requesting pe */
#define MSG_pe_mask 0x00000000FFFFFFFF /*!< mask for those bits */ 

#define msg_pack( cnt, islast ) (((cnt) << MSG_cnt_shift) | (MYTHREAD << 1) | ((islast) & 0x1)) /*!< makes the message */
#define msg_cnt(msg)     ((msg) >> MSG_cnt_shift )                                              /*!< returns number of packages */
#define msg_pe(msg)     ( ((msg) & MSG_pe_mask) >> 1 )                                          /*!< returns who sent it */
#define msg_islast(msg) ((msg) & 0x1)                                                           /*!< returns whether or not they are done */


/*!
\brief A structure to hold the buffers and the flags used to control them  along with flags for the endgame.
\ingroup exstackgrp
*/
typedef struct exstack2_t { 
  
  /*******  Buffers  ***********/
  // Note: these are shared only so that we can use memget to transfer the buffers
  // any node-to-node transfer method would be fine
  int64_t buf_cnt;              /*!< the number of pkgs in a buffer */
  size_t pkg_size;              /*!< pkg_size = the number of bytes in a pkg (one element in a buffer) */
  SHARED char * s_rcv_buffer;   /*!< receive buffer space */
  SHARED char * s_snd_buffer;   /*!< send buffer space */
  char ** l_rcv_buffer;  /*!< 2D array of local pointers to view the receive buffers */ 
  char ** l_snd_buffer;  /*!< 2D array of local pointers to view the send buffers */ 
                                          /*!< locally... as if each thread has THREADS different send and rcv buffers */
  int64_t * push_cnt;     /*!< THREADS long array, number of pkgs pushed to each send buffer */
  char ** push_ptr;       /*!< THREADS long array, pointers to current position to push to in send buffers */
  int64_t * pop_cnt;      /*!< THREADS long array, number of pkgs popped in each received buffer */
  char ** pop_ptr;        /*!< THREADS long array, pointers to current position to pop from in receive buffers */
  int64_t pop_pe;         /*!< the pe that is currently being popped or pulled (-1 if no buffer has been received) */

  int64_t * flush_order;  /*!< THREADS+1 long array records successful flushs as a linked list */
  int64_t * flush_perm;   /*!< THREADS long array stores a random permutation of threads */
 
  /*******  Synchronization  ***********/
  SHARED int64_t * s_can_send; /*!< length THREADS on each thread (one for each of the send buffers). */
                                         /*!< 1 means safe to send this buffer i.e. the receiver is ready */
                                         /*!< 0 means buffer is not safe to send yet i.e. the receiver has not popped your last sent buffer */
  int64_t * l_can_send;        /*!< local pointer to the shared array of buffer states */
  // the sending thread sets l_can_send[pe] to 0 after sending a buffer to mark that the buffer can be used
  // when the popping pe finishes popping said buffer it globally resets it to 1
  // This is shared volatile as part of the barrier free synchronization

  /*******  Message Queue  ***********/
  int64_t msg_Q_mask;      /*!< the queue is a circular queue of length a power of 2 that is greater */
                           /* than (2*THREADS). We can do "mod" the queue length with this mask */
  SHARED int64_t * s_msg_queue; /*!< keeps track of the order of pull requests on each thread */
  int64_t * l_msg_queue;      /*!< local pointer to s_msg_queue */
  // This is shared volatile as part of the barrier free synchronization
  // The sender put the message on the receivers queue, the receiver uses it locally to pop the stacks

  SHARED int64_t * s_num_msgs;  /*!< keeps track of total number of pull requests received on this thread  */
                                /*!< This is the head of the msg_queue. It is updated with a fetch_and_add */
  int64_t * l_num_msgs;         /*!< my local pointer to s_num_msgs */
  int64_t num_popped;           /*!< keeps track of the total number of stacks popped.  This is the tail of the queue */

  int64_t * active_buffer_queue; /*!< array of indices into the msg_queue for active (=unpopped and recvd) buffers. */
  int64_t num_active_buffers;  /*!< current number of active buffers */
  int64_t current_active_index;/*!< our current index in the active_buffer_queue */
  int64_t num_made_active;     /*!< total number of active buffers we have seen */

  int all_done;         /*!< This exstack is all done (exstack2_proceed should return 0), useful in nested exstacks */
  int num_done_sending; /*!< This is a counter of the number of threads that are done sending to me. */
} exstack2_t;

/* public routines */
exstack2_t *exstack2_init(int64_t buf_cnt, size_t pkg_size);
int64_t exstack2_proceed(exstack2_t *Xstk2, int done_pushing);
int64_t exstack2_push(exstack2_t *Xstk2, void *pkg, int64_t pe);
int64_t exstack2_pop(exstack2_t *Xstk2, void *pkg, int64_t *from_pe); 
int64_t exstack2_send(exstack2_t * Xstk2, int64_t pe, int islast);
   void exstack2_unpop(exstack2_t *Xstk2);
  void *exstack2_pull(exstack2_t *Xstk2, int64_t *from_pe); 
   void exstack2_unpull(exstack2_t *Xstk2);
   void exstack2_reset(exstack2_t *Xstk2);
   void exstack2_clear(exstack2_t *Xstk2);

#endif 

