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
/*! \file exstack.upc
 * \brief A library to do bulk synchronous buffered communications in parallel programs.
 */

#include "exstack.h"
#include "libgetput.h"
//#include <upc_castable.h>

/*! \brief Initialize an exstack_t.
 * \param buf_cnt the number of package that will fit into a send or receive stack
 * \param pkg_size The number of bytes per package
 * \return pointer to the initialized exstack_t
 * \ingroup exstackgrp
 */
exstack_t * exstack_init( int64_t buf_cnt, size_t pkg_size)
{
  uint64_t th;
  int i , j , k , tmp ;

  if(!buf_cnt) return(NULL);
  exstack_t *XStk = calloc(1, sizeof(exstack_t));

  srand(THREADS*MYTHREAD + 1);                    // create  a different random
  XStk->put_order = calloc(THREADS, sizeof(int)); // memput sequence for each thread
  if(XStk->put_order == NULL) 
    return(NULL);
  for(i=0; i<THREADS; i++) 
    XStk->put_order[i] = i;
  for(j=THREADS-1; j>1; j--) {
    k               = rand()%(j+1);
    tmp             = XStk->put_order[k];
    XStk->put_order[k] = XStk->put_order[j];
    XStk->put_order[j] = tmp;
  }

  // get the size of an exstack_buffer_t struct (with buffer)
  // and round up to a multiple of an int64_t.
  XStk->bytes_per_stack = pkg_size*buf_cnt + sizeof(exstack_buffer_t);
  XStk->bytes_per_stack =((XStk->bytes_per_stack + sizeof(int64_t) - 1)/sizeof(int64_t))*sizeof(int64_t);
  /* each thread has a send and receive buffer for every other thread   */
  XStk->snd_buf = lgp_all_alloc(XStk->bytes_per_stack*THREADS*THREADS,sizeof(char));
  if(XStk->snd_buf   == NULL) return(NULL);
  XStk->rcv_buf = lgp_all_alloc(XStk->bytes_per_stack*THREADS*THREADS,sizeof(char));
  if(XStk->rcv_buf   == NULL) return(NULL);
  
  XStk->l_snd_buf = calloc(THREADS,sizeof(exstack_buffer_t *));
  if(XStk->l_snd_buf == NULL) return(NULL);
  XStk->l_rcv_buf = calloc(THREADS,sizeof(exstack_buffer_t *));
  if(XStk->l_rcv_buf == NULL) return(NULL);
  XStk->fifo_ptr  = calloc(THREADS,sizeof(char *));
  if(XStk->fifo_ptr  == NULL) return(NULL);
  XStk->push_ptr  = calloc(THREADS,sizeof(char *));
  if(XStk->push_ptr  == NULL) return(NULL);

  XStk->l_snd_buf[0] = lgp_local_part(exstack_buffer_t,  XStk->snd_buf);
  XStk->l_rcv_buf[0] = lgp_local_part(exstack_buffer_t,  XStk->rcv_buf);
  for(th=0; th<THREADS; th++) {
    if(th > 0){
      XStk->l_snd_buf[th] = (exstack_buffer_t*)(((char*)XStk->l_snd_buf[th-1]) + XStk->bytes_per_stack);
      XStk->l_rcv_buf[th] = (exstack_buffer_t*)(((char*)XStk->l_rcv_buf[th-1]) + XStk->bytes_per_stack);
    }
    XStk->l_snd_buf[th]->count  = 0L;
    XStk->l_rcv_buf[th]->count  = 0L;
    XStk->push_ptr[th] = XStk->l_snd_buf[th]->data;
    XStk->fifo_ptr[th] = XStk->l_rcv_buf[th]->data;
  }

  XStk->buf_cnt               = buf_cnt;
  XStk->pkg_size              = pkg_size;
  XStk->first_ne_rcv          = 0;

  XStk->wait_done = lgp_all_alloc(THREADS*THREADS,sizeof(int64_t));
  if(XStk->wait_done == NULL) return(NULL);             
  
  XStk->l_wait_done = lgp_local_part(int64_t, XStk->wait_done);
  for(th=0; th<THREADS; th++) 
    XStk->l_wait_done[th] = 0;
  XStk->notify_done = 0;

  lgp_barrier();
  return XStk;
}


/*! \brief proceeds until all threads have
    reported that they have hit the done condition.
 * \param Xstk a pointer to the exstack struct
 * \param im_done a flag that indicates that this thread has finished
 * \return 0 if all threads have signaled done,  1 otherwise.
 * \ingroup exstackgrp
 */
int64_t exstack_proceed(exstack_t *Xstk , int im_done) {
  int i;
  int64_t breakout;
  if( im_done && !Xstk->notify_done ) {
    // If I am done pushing and I haven't notified everyone, do so now.
    for(i=0; i<THREADS; i++) 
      lgp_put_int64(Xstk->wait_done, MYTHREAD*THREADS + Xstk->put_order[i], 1L);
    Xstk->notify_done = 1;
  }
  lgp_barrier();
  
  breakout = 0;
  
  // See if everyone else is all done pushing.
  if( Xstk->notify_done ){ 
    breakout = 1L; 
    for(i=0; i<THREADS; i++)
      breakout &= Xstk->l_wait_done[i];
  }
  
  if(breakout){ 
    return(0); 
  }else{
    return(1);
  } 
}

/*! \brief push a pkg from a thread onto a send stack
 * \param Xstk pointer to an exstack_t
 * \param th_num the id of the thread to which we are sending these package
 * \param push_item a pointer to the package being pushed
 * \return headroom the amount of room (in packages) that was on the stack when called
 * Note: hence returns 0 if push fails and non-zero otherwise
 * \ingroup exstackgrp
 */
int64_t exstack_push(exstack_t *Xstk, void *push_item, int64_t th_num)
{
  uint64_t headroom = Xstk->buf_cnt - Xstk->l_snd_buf[th_num]->count;
  if( headroom ) {
    memcpy(Xstk->push_ptr[th_num], (char*)push_item, Xstk->pkg_size);
    Xstk->l_snd_buf[th_num]->count++;
    Xstk->push_ptr[th_num] += Xstk->pkg_size;
  }
  return(headroom);
}


/*! \brief Uses memputs to essentially do any an all_to_all of 
     send stacks to receive stacks
   Note: this is a variable length all_to_all, since we on send whats 
     needed from each stack
 * \param Xstk pointer to an exstack_t
 * \ingroup exstackgrp
 */
void exstack_exchange(exstack_t *Xstk )
{
  uint64_t ran, th, snd_buf_ct;

  // copy the contents of all non-empty send buffers from my thread
  // to the receive buffers of all other threads
  for(ran=0; ran<THREADS; ran++) {
    th = (uint64_t)Xstk->put_order[ran];
    snd_buf_ct = Xstk->l_snd_buf[th]->count;
    //printf("PE %d: sending %ld to %ld\n", MYTHREAD, snd_buf_ct, th);
    if(snd_buf_ct) {// only send if buffer has something in it
      lgp_memput(Xstk->rcv_buf,
                 Xstk->l_snd_buf[th],
                 snd_buf_ct*Xstk->pkg_size + sizeof(exstack_buffer_t),
                 Xstk->bytes_per_stack*MYTHREAD*THREADS + th);      
    }
  }

  lgp_barrier();
  //fprintf(stderr,"recvd %ld\n", Xstk->l_rcv_buf[th]->count);
  // clear send buffers for my thread and 
  // reset crnt_min_headroom and first_ne_rcv
  for(th=0; th<THREADS; th++) {
    Xstk->l_snd_buf[th]->count = 0L;
    Xstk->push_ptr[th] = Xstk->l_snd_buf[th]->data;
  }
  Xstk->first_ne_rcv = 0;

  // Note the barrier is after the memputs so that each thread can start as soon as its stacks are ready
  lgp_barrier();

  // set rcv buffer pointers after we know that all the memputs have completed
  for(th=0; th<THREADS; th++) {
    Xstk->fifo_ptr[th]  = Xstk->l_rcv_buf[th]->data;
  }
  
}

/*! \brief pops a package from specified thread
 * \param Xstk pointer to the exstack_t
 * \param pop_item pointer to object holding the package being pop'd
 * \param th_num thread id of the stack to be popped.
 * \return 1 on success, 0 if the stack was empty
 * Note: this uses the buffer as a fifo
 * \ingroup exstackgrp
 */
int64_t exstack_pop_thread(exstack_t *Xstk, void *pop_item, int64_t th_num)
{
  uint64_t i;

  if(Xstk->l_rcv_buf[th_num]->count == 0) 
    return(0);

  memcpy((char*)pop_item , Xstk->fifo_ptr[th_num] , Xstk->pkg_size);
  Xstk->fifo_ptr[th_num] += Xstk->pkg_size;
  Xstk->l_rcv_buf[th_num]->count--;
  return(1);
}

/*! \brief unpops an item after a exstack_pop_thread.
 * \param Xstk pointer to an exstack_t
 * \param th_num thread id of item that gets unpopped
 * Note: it is only safe to unpop from an exstack immediately after 
 * a successful pop.
 * \ingroup exstackgrp
 */
void exstack_unpop_thread(exstack_t *Xstk , int64_t th_num)
/**********************************************************************/
{
  Xstk->fifo_ptr[th_num] -= Xstk->pkg_size;
  Xstk->l_rcv_buf[th_num]->count++;
}


/*! \brief first-in-first-out pop from any thread, ie the next available package
 * \param Xstk pointer to an exstack_t
 * \param pop_item pointer to object holding the package being popped
 * \param from_th  pointer to thread id of item that gets popped or NULL if you don't care
 * \return 1 on success, 0 if all stacks are empty
 * \ingroup exstackgrp
 */
int64_t exstack_pop(exstack_t *Xstk, void *pop_item ,  int64_t *from_th)
{
  int64_t th;
  //fprintf(stderr,"I AM Calling exstack pop\n");
  for(th=Xstk->first_ne_rcv; th<THREADS; th++) {
    //fprintf(stderr,"calling exstack_pop count = %ld\n", Xstk->l_rcv_buf[th]->count);
    if(Xstk->l_rcv_buf[th]->count == 0L) continue;

    // non-empty rcv buffer found: pop work item:
    memcpy((char*)pop_item , Xstk->fifo_ptr[th] , Xstk->pkg_size);
    Xstk->fifo_ptr[th] += Xstk->pkg_size;
    Xstk->l_rcv_buf[th]->count--;
    Xstk->first_ne_rcv = th;
    if( from_th != NULL )
       *from_th = th;
    return(1);
  }
  //fprintf(stderr,"I AM ALL DONE\n");
  // all receive buffers empty:  get ready for next exstack_memcpy and  return failure
  Xstk->first_ne_rcv = 0;
  for(th=0; th<THREADS; th++) {
    Xstk->fifo_ptr[th] = Xstk->l_rcv_buf[th]->data;
  }
  return (0);
}

/*! \brief unpops an item after a exstack_pop
 * \param Xstk pointer to an exstack_t
 * Note, in order to unpop a pop_any 
 * use the from_thread filled in by the that pop command.
 * Note: it is only safe to unpop from an exstack immediately after 
 * a successful pop.
 * \ingroup exstackgrp
 */
void exstack_unpop(exstack_t *Xstk)
/**********************************************************************/
{
  // ONLY A GUESS AT WHAT MIGHT WORK 
  Xstk->fifo_ptr[Xstk->first_ne_rcv] -= Xstk->pkg_size;
  Xstk->l_rcv_buf[Xstk->first_ne_rcv]->count++;
}

/*! \brief first-in-first-out pop from any thread, ie the next available package
 * \param Xstk pointer to an exstack_t
 * \param from_th  pointer to thread id of item that gets popped or NULL if you don't care
 * \return a ptr to the popped item with the exstack, or on failure
 * \ingroup exstackgrp
 */
void *exstack_pull(exstack_t *Xstk, int64_t *from_th)
{
  char *ret;
  int64_t th;
  
  for(th=Xstk->first_ne_rcv; th<THREADS; th++) {
    //printf("TH %d: %ld items from %ld\n", MYTHREAD, Xstk->l_rcv_buf[th]->count, th);
    if(Xstk->l_rcv_buf[th]->count == 0L) continue;

    // non-empty rcv buffer found: pop work item:
    ret = Xstk->fifo_ptr[th];
    Xstk->fifo_ptr[th] += Xstk->pkg_size;
    Xstk->l_rcv_buf[th]->count--;
    Xstk->first_ne_rcv = th;
    if( from_th != NULL )
       *from_th = th;
    
    return((void*)ret);
  }

  // all receive buffers empty:  get ready for next exstack_memcpy and  return failure
  Xstk->first_ne_rcv = 0;
  for(th=0; th<THREADS; th++) {
    Xstk->fifo_ptr[th] = Xstk->l_rcv_buf[th]->data;
  }
  return (NULL);
}

/*! \brief unpulls an item after a exstack_pull
 * \param Xstk pointer to an exstack_t
 * Note: it is only safe to unpull from an exstack immediately after 
 * a successful pull.
 * \ingroup exstackgrp
 */
void exstack_unpull(exstack_t *Xstk)
/**********************************************************************/
{
  // ONLY A GUESS AT WHAT MIGHT WORK 
  Xstk->fifo_ptr[Xstk->first_ne_rcv] -= Xstk->pkg_size;
  Xstk->l_rcv_buf[Xstk->first_ne_rcv]->count++;
}



/*! \brief returns the minimum headroom over all of my send buffers
 * \param Xstk is the pointer to the exstack
 * \return the min headroom (in number of packages) across all send buffers
 * \ingroup exstackgrp
 */
int64_t exstack_min_headroom(exstack_t *Xstk)
{
  uint64_t t, ret, snd_buf_ct_ptr;
  // first 64 bits of the buffer hold the current buffer count, regardless of the package size.
  ret = 0;
  for(t=0; t<THREADS; t++) {
    snd_buf_ct_ptr = Xstk->l_snd_buf[t]->count;   
    ret += (Xstk->buf_cnt - snd_buf_ct_ptr); // the headroom in packages going to MYTHREAD

  }
  return(ret);      // min across all my send buffers
}

/*! \brief returns the headroom for a send buffer
 * \param Xstk is the pointer to the exstack
 * \param th_num is the desired thread number
 * \return the headroom (in number of packages) for the send buffer for th_num
 * \ingroup exstackgrp
 */
int64_t exstack_headroom(exstack_t *Xstk , int64_t th_num)
{
  int64_t snd_buf_ct_ptr = Xstk->l_snd_buf[th_num]->count;
  return( Xstk->buf_cnt - snd_buf_ct_ptr);
}


/*! \brief frees all memory allocated by exstack_init
 * \param Xstk is the pointer to the exstack stack
 * \ingroup exstackgrp
 */
void exstack_clear(exstack_t *Xstk )
{
  int64_t i;
  lgp_all_free(Xstk->wait_done);
  lgp_all_free(Xstk->snd_buf);
  lgp_all_free(Xstk->rcv_buf);
  free(Xstk->put_order);
  free(Xstk->l_snd_buf);
  free(Xstk->l_rcv_buf);
  free(Xstk->fifo_ptr );
  free(Xstk->push_ptr );
}

/*! \brief resets the exstack buffers and resets the end game 
 * \param Xstk is the pointer to the exstack stack
 * \ingroup exstackgrp
 */
void exstack_reset(exstack_t * Xstk){
  
  for(int th=0; th<THREADS; th++) {
    Xstk->l_snd_buf[th]->count  = 0L;
    Xstk->l_rcv_buf[th]->count  = 0L;
    // set the push and pop ptrs to the item in the buffers
    Xstk->push_ptr[th] = Xstk->l_snd_buf[th]->data;
    Xstk->fifo_ptr[th] = Xstk->l_rcv_buf[th]->data;
  }
  Xstk->first_ne_rcv = 0;
  for(int th=0; th<THREADS; th++) 
    Xstk->l_wait_done[th] = 0;
  Xstk->notify_done = 0;
  lgp_barrier();
}

