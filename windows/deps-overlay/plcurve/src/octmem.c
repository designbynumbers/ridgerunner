/*
 *  Octree memory management.
 *
 *    reserve_octmem() either registers the space handed to it in mem or,
 *    if mem is NULL, reserves the desired amount of space (perhaps in addition
 *    to space already reserved).
 *
 *    free_octmem() frees any reserved octmem and cleans up pointers and such.
 *
 *    oct_alloc() doles out the space.
 *
 *  $Id: octmem.c,v 1.8 2006-04-18 19:14:37 ashted Exp $
 *
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of liboctrope.
   
liboctrope is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

liboctrope is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with liboctrope; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_MALLOC_H
  #include <malloc.h>
#endif

#include "octrope.h"
#include "octrope_tls.h"

/* Per-thread arenas so concurrent octrope callers do not share bump state. */
static OCTROPE_TLS void *octmem = NULL;
static OCTROPE_TLS void *top = NULL;
static OCTROPE_TLS int   provided = 0; /* Memory was or was not provided by the *
                                        * calling program                       */
static OCTROPE_TLS int   octmem_size = 0;
static OCTROPE_TLS void *aux_mem = NULL;
static OCTROPE_TLS void *aux_mem_top = NULL;
static OCTROPE_TLS int   aux_mem_size = 0;
static OCTROPE_TLS int   malloced_aux = 0; 

void reserve_octmem(const int bytes_desired,void *mem,const int memsize) {

  if (mem != NULL && memsize >= bytes_desired) {
    top = octmem = mem;
    octmem_size = memsize;
    provided = 1;
  } else {
    top = octmem = (void *)(malloc(bytes_desired));
    if (octmem == NULL) {
      fprintf(stderr,
        "Unable to reserve %d bytes in reserve_octmem.\n",bytes_desired);
      exit(-1);
    }
    octmem_size = bytes_desired;
    provided = 0;
  }
#ifdef DEBUG
  if (octrope_debug_level() > 0) {
    printf("Total memory is: %d\n",octmem_size);
  }
#endif
}

void reserve_aux_mem(const int bytes_desired) {
  
  if (bytes_desired + (top - octmem) <= octmem_size) {
    /* reserve the space at the top of allocated memory */
    aux_mem_top = aux_mem = octmem + octmem_size - bytes_desired;
    aux_mem_size = bytes_desired;
    octmem_size -= bytes_desired;
  } else {
    aux_mem_top = aux_mem = (void *)(malloc(bytes_desired));
    if (aux_mem == NULL) {
      fprintf(stderr,
        "Unable to reserve %d bytes in reserve_aux_mem.\n",bytes_desired);
      exit(-1);
    } 
    malloced_aux = 1;
    aux_mem_size = bytes_desired;
  }
#ifdef DEBUG
  if (octrope_debug_level() > 0) {
    printf("Total memory is now: %d\n",octmem_size + aux_mem_size);
  }
#endif
}
    
void free_octmem() {
#ifdef DEBUG
  int size1,size2;

  if (octrope_debug_level() > 0) {
    size1 = octmem_size + (octmem - top);
    size2 = aux_mem_size + (aux_mem - aux_mem_top);
    printf("Bytes unused: %d + %d = %d\n",size1,size2,size1+size2);
  }
#endif
  if (octmem != NULL && provided == 0) {
    free(octmem);
  }
  if (malloced_aux) {
    free(aux_mem);
  }
  octmem = top = NULL;
  octmem_size = 0;
  provided = 0;
  aux_mem = aux_mem_top = NULL;
  aux_mem_size = 0;
  malloced_aux = 0;
}

void *oct_alloc(int nelem,int size) {
  int   bytes_requested;
  void *retptr;

  bytes_requested = nelem*size;
  if (bytes_requested + (top - octmem) > octmem_size) {
    return NULL;
  }
  retptr = top;
  top += bytes_requested;
  return retptr;
}

void *aux_alloc(int nelem,int size) {
  int   bytes_requested;
  void *retptr;

  bytes_requested = nelem*size;
  if (bytes_requested + (aux_mem_top - aux_mem) > aux_mem_size) {
    return NULL;
  }
  retptr = aux_mem_top;
  aux_mem_top += bytes_requested;
  return retptr;
}
