/*
 *
 * Data structures and prototypes for octrope_links
 *
 *  $Id: octrope_link.h,v 1.1 2005-03-10 02:46:02 michaelp Exp $
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

#ifndef __OCTROPE_LINK_H
#define __OCTROPE_LINK_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include <octrope_truefalse.h>
#include <octrope_vector.h>


typedef struct octrope_color_type {
  double r;
  double g;
  double b;
  double alpha;
} octrope_color;

typedef struct octrope_pline_type {
  int             acyclic;   /* This is an "open" pline (with distinct ends) */
  int             nv;        /* Number of vertices */
  int             cc;        /* Color count (number of colors) */
  octrope_vector *vt;        /* Actual vertices */
  octrope_color  *clr;       /* Colors */
  /***** Need a way to specify constraints on endpoints here *****/
} octrope_pline;

typedef struct octrope_link_type {	
  int nc;			/* Number of components */
  octrope_pline *cp;            /* Components */
} octrope_link;

/* 
 * Prototypes for routines to deal with links.  More in-depth documentation is
 * available in link.c.
 *
 */

/* Build a new link (with associated plines) */
octrope_link *octrope_link_new(int components, 
                               const int *nv, 
                               const int *acyclic,
                               const int *cc);

/* Free the link (and plines) */
void          octrope_link_free(octrope_link *L);

/* Read link data from a file */
octrope_link *octrope_link_read(FILE *infile);

/* Write link data to a file */
int           octrope_link_write(FILE *outfile, const octrope_link *L);

/* Fix the "hidden vertices" for easy handling of closed components */
void          octrope_link_fix_wrap(const octrope_link *L);

/* Count the edges in a link (correctly handling open/closed) */
int           octrope_link_edges(const octrope_link *L);

/* Copy a link */
octrope_link *octrope_link_copy(const octrope_link *L);

#if (__cplusplus || c_plusplus)
};
#endif
#endif
