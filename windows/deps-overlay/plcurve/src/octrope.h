/*
 *
 * Data structures and prototypes necessary for using liboctrope.a
 *
 *  $Id: octrope.h,v 1.25 2007-06-14 21:04:06 cantarel Exp $
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
  
#include <plCurve.h>
  
typedef struct octrope_strut_type {
  int component[2];    /* On which component the strut starts/ends */
  int lead_vert[2];    /* Lead vertex of the edge on which it starts/ends */
  double position[2];  /* How far along the starting/ending edge (0 = at lead vert, 1 = at trailing vert) */
  double length;       /* How long it is */
  double compression;  /* A Lagrange multiplier for this active constraint */
                       /* This is only set by external programs. */
} octrope_strut;
  
typedef struct octrope_mrloc_type {
  int component;       /* On which component the vertex is found */
  int vert;            /* The location of the "primary" vertex, at which minrad is measured. */
  int svert;           /* The "secondary" vertex. The edgelength in minrad is the length of edge (svert,vert). */
                       /* The secondary vertex should be either vert+1 or vert-1.*/

  double mr;
  double compression;  /* Again, a Lagrange multiplier */
                       /* Only set by external programs. */
} octrope_mrloc;
    
extern int octrope_error_num;
extern char octrope_error_str[1024];

/************************ Now for the main routines **************************/

/*
 * Build a list of POCAs. If strut_cutoff > 0, we return a list of up to 
 * sl_size POCAs with length < strut_cutoff. In this case, epsilon is ignored.
 * If there are more than sl_size such POCAs, the selection of POCAs returned
 * is essentially random. In particular, there is no guarantee that the shortest
 * POCAs will be returned.

 * If strut_cutoff = 0, then we return a list of all POCAs which are within 
 * epsilon of the shortest POCA length and return the number that are in that list.
 * Again, if too many POCAs are in this group, the selection of POCAs returned is
 * random.
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
int octrope_struts(plCurve *L, 
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
double octrope_minrad(plCurve *L, 
                      const double mr_cutoff, const double epsilon,
                      octrope_mrloc *min_rad_locs,
                      int mr_size, int *num_min_rad_locs);

/*
 * Calculate the actual curve length of the knot or link.
 *
   Simply adds up the lengths of the various edges and returns a number.
 */
double octrope_curvelength(const plCurve *L);

/*
 * Return the length of the shortest strut.
 */
double octrope_poca(plCurve *L, void *mem,const int memsize);

/*
 * Return the value of MinRad.
 */
double octrope_minradval(plCurve *L);

/*
 * Return the thickness of the knot, which is the minimum of
 *   1) All strut lengths/2
 *   2) MinRad(L)/lambda
 *
 * If mem is non-NULL and memsize > 0, ropelength will use that space as its
 * initial workspace (need not be zeroed).
 */

double octrope_thickness(plCurve *L,void *mem,const int memsize,
                         const double lambda);

/*
 * Return the curve length divided by the thickness (see above).
 * 
 * If mem is non-NULL and memsize > 0, ropelength will use that space as its
 * initial workspace (need not be zeroed).
 */

double octrope_ropelength(plCurve *L,void *mem,const int memsize,
                          const double lambda);

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

void octrope(plCurve *L,      /* The knot or link */

             /* Here are places to put the ropelength and thickness values.*/
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
             const int memsize,

             /* The lambda for thickness calculations was moved down here to
              * get folks to notice that thickness has changed definitions. 
              */
             const double lambda);                
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
 * OpenMP thread count for octrope_struts edge search (no-op without OpenMP).
 * Pass nthreads > 0 to set; octrope_get_threads returns omp_get_max_threads().
 */
void octrope_set_threads(int nthreads);
int octrope_get_threads(void);

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
void octrope_strut_ends(const plCurve *L, const octrope_strut *S, 
                        plc_vector se[2]);

/* 
 * Read a list of octrope_strut and octrope_mrloc records stored on disk. 
 */

void octrope_strutfile_read(plCurve *L,
                            int *n_struts,octrope_strut **strutlist,
                            int *n_mrloc,octrope_mrloc **mrlist,
                            FILE *file);

/*
 * Write a list of octrope_strut and octrope_mrloc records to disk.
 */

void octrope_strutfile_write(int n_struts,octrope_strut *strutlist,
                             int n_mrloc,octrope_mrloc *mrlist,
                             FILE *file);
/*
 * Print out or write to buffer the version of octrope used. 
 */

void octrope_version(char *version, size_t strlen);

#if (__cplusplus || c_plusplus)
};
#endif
#endif

