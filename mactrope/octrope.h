/*
 *
 * Data structures and prototypes necessary for using liboctrope.a
 *
 *  $Id: octrope.h,v 1.2 2005-03-10 02:46:02 michaelp Exp $
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

#ifndef __OCTROPE_H
#define __OCTROPE_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif 
  
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include "octrope_link.h"
#include "octrope_truefalse.h"
#ifdef HAVE_FLOAT_H
#include <float.h>
#endif
#ifndef DBL_MAX
#  ifdef FLT_MAX
#    define DBL_MAX FLT_MAX
#  else
#    define DBL_MAX 3.40282347e+38F
#  endif
#endif

  
typedef struct octrope_strut_type {
  int component[2];    /* On which component the strut starts/ends */
  int lead_vert[2];    /* Lead vertex of the edge on which it starts/ends */
  double position[2];  /* How far along the starting/ending edge (0 to 1) */
  double length;       /* How long it is */
} octrope_strut;
  
typedef struct octrope_mrloc_type {
  int component;       /* On which component the vertex is found */
  int vert;            /* The location of the vertex. */
  double mr;
} octrope_mrloc;
    
static int octrope_error_num;
static char octrope_error_str[80];

/************************ Now for the main routines **************************/

/*
 * Build a list of all the POCAs which are within epsilon of the shortest
 * POCA length (so long as they are at most max_strut_len long) and return the
 * number that are in that list.
 *
 * strutlist is the place to store the actual struts.  If strutlist is NULL,
 * only the number of struts is returned.  sl_size is the number of struts
 * which can be stored in strutlist--if strutlist is too small, 
 * octrope_struts() returns as many as will fit in the list.
 * shortest holds the length of the shortest POCA.  
 *
 * If mem is non-NULL and memsize > 0, ropelength will use that space as its
 * initial workspace (need not be zeroed).
 */
int octrope_struts(const octrope_link *L, 
                   const double strut_cutoff, const double epsilon, 
                   octrope_strut *strutlist, int sl_size,
                   double *shortest, void *mem, const int memsize);

/*
 * Calculate the "minrad" (as defined by Eric Rawdon) of a link.
 *
   If <min_rad_locs> is not NULL, and <mr_size> is not zero, we 
   assume that min_rad_locs is a list of octrope_mrloc structures.
   We then search the link for vertices with minrad < max_minrad.
   The ones we find (if any) are entered into the <min_rad_locs> array.

 */
double octrope_minrad(const octrope_link *L, 
                      const double mr_cutoff, const double epsilon,
                      octrope_mrloc *min_rad_locs,
                      int mr_size, int *num_min_rad_locs);

/*
 * Calculate the actual curve length of the knot or link.
 *
   Simply adds up the lengths of the various edges and returns a number.
 */
double octrope_curvelength(const octrope_link *L);

/*
 * Return the length of the shortest strut.
 */
double octrope_poca(const octrope_link *L, void *mem,const int memsize);

/*
 * Return the value of MinRad.
 */
double octrope_minradval(const octrope_link *L);

/*
 * Return the thickness of the knot, which is the minimum of
 *   1) All strut lengths
 *   2) 2*MinRad(L)/lambda
 *   2) 2*MinRad(L)/lambda
 *
 * If mem is non-NULL and memsize > 0, ropelength will use that space as its
 * initial workspace (need not be zeroed).
 */

double octrope_thickness(const octrope_link *L,const double lambda,
                          void *mem,const int memsize);

/*
 * Return the curve length divided by the thickness (see above).
 * 
 * If mem is non-NULL and memsize > 0, ropelength will use that space as its
 * initial workspace (need not be zeroed).
 */

double octrope_ropelength(const octrope_link *L,const double lambda,
                          void *mem,const int memsize);

/* 

  The main "combined" call allows the user to obtain all of the information
  determined by the octrope library in a single call. 

  This call computes both the ropelength and the strut list, as well as
  finding the locations at which minrad is in control of ropelength. 

  A brief explanation of the variables follows:

     <strutlist>    will contain a list of all struts with length less than
                    <strut_cutoff> or, if strut_cutoff=0, within <epsilon> of
                    the shortest strut. The size of the list is passed in
                    <sl_size>, and the number of entries returned in
                    <num_struts>.

     <lambda>       is the stiffness of the rope. One can view lambda as an a
                    priori lower bound on the diameter of curvature of the unit
                    thickness rope. 

     <min_rad_locs> will contain a list of all the vertices at which min_rad is
                    at the controlling value. As with the strut list, <mr_size>
                    gives the size of the list, and the number of such entries
                    found is passed back in <num_min_rad_locs>. 

     <mem> and 
     <memsize>      can be used to pass octrope a memory buffer to use for the
                    computation.  The memory does not have to be zeroed ahead
                    of time. If <memsize> is 0, or <mem> is NULL, octrope will
                    allocate its own memory. 

  Octrope returns -1 on failure, or 1 for success. 

*/

void octrope(const octrope_link *L,      /* The knot or link */

             /* Here we have the lambda for the thickness calculation and
                places to put the ropelength and thickness values.  */
             const double lambda,         
             double *ropelength,
             double *thickness_ret,

             /* The actual length of the curve(s) in the knot/link */
             double *curve_len_ret,

             /* Places to store MinRad value and length of shortest strut */
             double *min_rad_ret, 
             double *min_strut_ret, 

             /* mr_epsilon is the tolerance which decides whether a given
              * MinRad location is saved or not.  If MinRad locations are being
              * saved, all four of the following need to have reasonable
              * values.  Otherwise, minrad_locs and num_min_rad_locs should be
              * NULL and mr_size should be 0. 
              */
             const double mr_cutoff,
             const double mr_epsilon, 
             octrope_mrloc *min_rad_locs, 
             const int mr_size,
             int *num_min_rad_locs, 

             /* strut_epsilon serves as mr_epsilon and above for deciding which
              * which struts are stored in the strutlist.  As above, all 5
              * parameters need to be reasonable if struts are to be saved.
              * Otherwise, strutlist and num_strut_ret should be NULL and
              * sl_size should be 0. 
              */
             const double strut_cutoff,
             const double strut_epsilon, 
             octrope_strut *strutlist, 
             const int sl_size,
             int *num_strut_ret,

             void *mem,                          
             const int memsize);                
/* 
 * Estimate the amount of memory (in bytes) needed to call struts() or
 * ropelength() on a link with num_edges edges.
 *
 */
int octrope_est_mem(const int num_edges);

/*
 * Set the octree's maximum number of levels.  Calling with levels equal to 0
 * returns to the library default calculation.
 *
 */
void octrope_set_levels(int levels);

/*
 * Set debugging level 0 through 9, higher is more output
 *
 */
void octrope_set_debug(int level); 
int octrope_debug_level();

/*
 * Translate the strut-format information into two vectors which are the
 * ends of the strut in question.  S points to the strut, se should 
 * be the address of a vector[2].
 *
 */
void octrope_strut_ends(const octrope_link *L, const octrope_strut *S, 
                        octrope_vector se[2]);

#if (__cplusplus || c_plusplus)
};
#endif
#endif

