
/*

   eqedge.h : Header file for demo equilateralization code using LAPACK. 

*/

/************************* Include Files **************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>

/*#include"cblas.h"*/
#ifdef __APPLE__
#include <vecLib/vBLAS.h>
#include <vecLib/clapack.h>
#else
#include "gsl_cblas.h"
#endif // __APPLE__

#include "octrope_link.h"
#include "octrope_vector.h"
#include "stepper.h"

/************************* Global Variables ***********************/

int VERBOSITY;
int GRAPHICS;
FILE *gclpipe;

/************************* Function Prototypes ********************/

/* Utility Routines */

void our_matrix_write(double val, double *A, int LDA, int i, int j);
double maxError(octrope_link *L, int comp, search_state* inState);

void			link_scale( octrope_link* inLink, double factor );
octrope_link* octrope_fixlength( octrope_link* inLink );
octrope_link*	octrope_equilize_density( octrope_link* inLink, search_state* inState );
octrope_link* octrope_double_edges( octrope_link* inLink );
octrope_link*	octrope_double_component( octrope_link* inLink, int comp );

int    octrope_pline_edges(octrope_pline *P);
double octrope_pline_length(octrope_pline *P); 
double octrope_link_short_edge(octrope_link *L);
double octrope_link_long_edge(octrope_link *L);

octrope_vector octrope_link_tangent_vector(octrope_link *L,int comp,int vert);
octrope_vector octrope_link_edge_dir(octrope_link *L,int comp,int edge);

void   octrope_link_draw(FILE *outfile, octrope_link *L);

/* Display Routines */

void	init_display();
void	shutdown_display();
void   refresh_display(octrope_link *L);
void	export_struts(octrope_link* inLink, octrope_strut* inStruts, int inSize, double* compressions, double time);
void	exportVect( const octrope_vector* dl, octrope_link* link, const char* fname );

/* Computational Routines. */

double *tangential_rigidity_matrix(octrope_link *L, int comp, int *LDA, int *N);
double *edgelength_error(octrope_link *L,int comp, search_state* inState);

/* External interface. */

void lapack_eqedge(octrope_link *L, int nsteps, search_state* inState);

/* Demos. */

void segment_demo();
void ellipse_demo();




