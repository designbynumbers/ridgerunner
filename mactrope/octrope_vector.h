/*
 * 
 * Prototypes for routines in vector.c.
 *
 * $Id: octrope_vector.h,v 1.1 2005-03-10 02:46:02 michaelp Exp $
 *
 */

/* Copyright 2004 The University of Georgia */

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

#ifndef __OCTROPE_VECTOR_H
#define __OCTROPE_VECTOR_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct octrope_vector_type {    /* A point in 3-space. */
  double c[3];
} octrope_vector;

/*
 * Prototypes for vector routines. 
 *
 */

octrope_vector octrope_vplus(octrope_vector A,octrope_vector B);
octrope_vector octrope_vminus(octrope_vector A,octrope_vector B);
octrope_vector octrope_cross(octrope_vector A,octrope_vector B);
octrope_vector octrope_scalarmult(double x,octrope_vector A);

/*
 * Macro replacements (requires some reprogramming)
 *
 */

#define octrope_dot(A,B)    ((A).c[0]*(B).c[0] + (A).c[1]*(B).c[1] + (A).c[2]*(B).c[2])
#define octrope_norm(A)     sqrt(octrope_dot((A),(A)))
#define octrope_vadd(A,B)  \
  (A).c[0] += (B).c[0]; (A).c[1] += (B).c[1]; (A).c[2] += (B).c[2];
#define octrope_vsub(A,B)  \
  (A).c[0] -= (B).c[0]; (A).c[1] -= (B).c[1]; (A).c[2] -= (B).c[2];
#define octrope_vsmult(s,V) (V).c[0] *= s; (V).c[1] *= s; (V).c[2] *= s;
  /* Add a multiple of B to A */
#define octrope_vmadd(A,s,B)  (A).c[0] += (s)*(B).c[0]; \
                              (A).c[1] += (s)*(B).c[1]; \
                              (A).c[2] += (s)*(B).c[2];
  /* A = B + s(C-B)                               *
   * equivalent to                                *
   *   A = C; vsub(A,B); vsmult(s,A); vsadd(A,B); */
#define octrope_vweighted(A,s,B,C)  \
    (A).c[0] = (B).c[0] + s*((C).c[0] - (B).c[0]); \
    (A).c[1] = (B).c[1] + s*((C).c[1] - (B).c[1]); \
    (A).c[2] = (B).c[2] + s*((C).c[2] - (B).c[2]); 

#endif
