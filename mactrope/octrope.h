/*
 *
 * Data structures and prototypes necessary for using liboctrope.a
 *
 *  $Id: octrope.h,v 1.1 2004-08-16 22:35:30 michaelp Exp $
 *
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
#include "link.h"
#include "truefalse.h"
  
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
    

/************************ Now for the main routines **************************/

/*
 * Return the length of the link divided by the minimum of
 *   1) All strut lengths/2
 *   2) MinRad(L)/lambda
 * 
 * If mem is non-NULL and memsize > 0, ropelength will use that space as its
 * initial workspace (need not be zeroed).
 */

  double octrope_ropelength(const octrope_link *L,const double lambda,
                          double *min_rad_ret, double *min_strut_ret,
                          double *curve_len_ret, void *mem,const int memsize);

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
  int octrope_struts(const octrope_link *L, const double epsilon, 
                   const double max_strut_len, 
                   octrope_strut *strutlist, int sl_size,
                   double *shortest, void *mem, const int memsize);

/* 

  The main "combined" call allows the user to obtain all of the information
  determined by the octrope library in a single call. 

  This call computes both the ropelength and the strut list, as well as
  finding the locations at which minrad is in control of ropelength. 

  A brief explanation of the variables follows:

     <strutlist> will contain a list of all struts with length 
                 within <epsilon> of <shortest> and less than 
                 <max_strut_len>. The size of the list is passed
		 in <sl_size>, and the number of entries returned
		 in <num_struts>.

     <lambda>    is the stiffness of the rope. One can view lambda
                 as an a priori lower bound on the diameter of curvature
                 of the unit thickness rope. 

     <min_rad_locs> will contain a list of all the vertices at which 
                    min_rad is at the controlling value. As with the
                 strut list, <mr_size> gives the size of the list, and 
                 the number of such entries found is passed back in 
                 <num_min_rad_locs>. 

     <mem> and <memsize> can be used to pass octrope a memory buffer
                         to use for the computation. It does not have
                 to be zeroed. If <memsize> is 0, or <mem> is NULL, 
                 octrope will allocate its own memory. 

  Octrope returns -1 on failure, or 1 for success. 

*/


  int octrope (const octrope_link *L, const double epsilon, 
		const double max_strut_len, 
		const double lambda, 
		double *ropelength,
		
		octrope_mrloc *min_rad_locs, int mr_size,
		double *min_rad_ret, int *num_min_rad_locs, 
		
		octrope_strut *strutlist, int sl_size,
		double *shortest, int *num_struts,

      	        double *curve_len_ret,
		
		void *mem, 
		const int memsize);


/*
 * Calculate the "minrad" (as defined by Eric Rawdon) of a link.
 *
   If <min_rad_locs> is not NULL, and <mr_size> is not zero, we 
   assume that min_rad_locs is a list of octrope_mrloc structures.
   We then search the link for vertices with minrad < the value <target>.
   The ones we find (if any) are entered into the <min_rad_locs> array.

 */
  double octrope_minrad(const octrope_link *L, octrope_mrloc *min_rad_locs,
			int mr_size, int *num_min_rad_locs, double target);

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

