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
#include <shmem.h>
#include <stdio.h>
#if HAVE_CONFIG_H
#include "config.h"
#endif
#include "mpp2shmem.h"

int
main(void)
{
  shmem_init();
  int* data = mpp_alloc(sizeof(int));
  shmem_barrier_all();

  if (shmem_my_pe() == 0) {
    void* friend = shmem_ptr(data, 1);
    *data = (friend == NULL);
    if (friend == NULL)
      fputs("ERROR: shmem_ptr() is not working.  Either fix your environment to\n"
            "make shmem_ptr() work, or configure --without-shmem-ptr.\n", stderr);
  }

  shmem_barrier_all();
  int status = shmem_int_g(data, 0);
  shmem_barrier_all();
  shmem_finalize();
  return status;
}
