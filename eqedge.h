
/*

   eqedge.h : Header file for demo equilateralization code using LAPACK. 

*/

#ifndef _H_eqedge
#define _H_eqedge

/************************* Include Files **************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>


#ifdef __APPLE__
#include <vecLib/vBLAS.h>
#else
#include <gsl/gsl_cblas.h>
#endif // __APPLE__


#include"link.h"
#include"vector.h"

/************************* Global Variables ***********************/

int VERBOSITY;
int GRAPHICS;
FILE *gclpipe;

/************************* Function Prototypes ********************/

/* Utility Routines */

void our_matrix_write(double val, double *A, int LDA, int i, int j);
double maxError(octrope_link *L, int comp, double target_length);

int    octrope_pline_edges(octrope_pline *P);
double octrope_pline_length(octrope_pline *P); 
double octrope_link_short_edge(octrope_link *L);
double octrope_link_long_edge(octrope_link *L);

octrope_link* octrope_fixlength( octrope_link* inLink );

octrope_link* octrope_fixlength( octrope_link* inLink );

octrope_vector octrope_link_tangent_vector(octrope_link *L,int comp,int vert);
octrope_vector octrope_link_edge_dir(octrope_link *L,int comp,int edge);

void   octrope_link_draw(FILE *outfile, octrope_link *L);

/* Display Routines */

void	init_display();
void	shutdown_display();
void   refresh_display(octrope_link *L);
void	export_struts(octrope_link* inLink, octrope_strut* inStruts, int inSize, double* compressions);
void	exportVect( const octrope_vector* dl, octrope_link* link, const char* fname );

/* Computational Routines. */

double *tangential_rigidity_matrix(octrope_link *L, int comp, int *LDA, int *N);
double *edgelength_error(octrope_link *L,int comp, double target_length);

/* External interface. */

void lapack_eqedge(octrope_link *L, int nsteps);

/* Demos. */

void segment_demo();
void ellipse_demo();


#endif // _H_eqedge
