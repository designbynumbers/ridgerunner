/*
 *
 * Data structures and prototypes for octrope_links
 *
 *  $Id: link.h,v 1.1 2004-08-16 22:35:30 michaelp Exp $
 *
 */

#ifndef __LINK_H
#define __LINK_H

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdio.h>
#include "truefalse.h"

typedef struct octrope_vector_type {    /* A point in 3-space. */
  double c[3];
} octrope_vector;

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
