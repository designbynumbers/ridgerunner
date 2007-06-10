/*

   ridgerunner.h :

   Part of main ridgerunner project. Defines function prototypes and
   global variables for entire ridgerunner project. 

*/

#ifndef _H_ridgerunner
#define _H_ridgerunner

#include "portability.h"
#include "errors.h"
#include "stepper.h"
#include "settings.h"

#include "plCurve.h"
#include "octrope.h"
#include "libtsnnls/tsnnls.h"

/************************* Global Variables ***********************/

int VERBOSITY;
int GRAPHICS;
FILE *gclpipe;

/************************* Function Prototypes ********************/

/* Utility Routines */

void our_matrix_write(double val, double *A, int LDA, int i, int j);
double maxError(plCurve *L, int comp, search_state* inState);

plCurve*  octrope_fixlength( plCurve* inLink );
plCurve*  octrope_double_edges( plCurve* inLink );

double	  plCurve_torsion( plCurve* inLink, FILE* outPlot );

double plCurve_short_edge(plCurve *L);
double plCurve_long_edge(plCurve *L);

plc_vector plCurve_tangent_vector(plCurve *L,int comp,int vert);
plc_vector plCurve_edge_dir(plCurve *L,int comp,int edge);

void   plCurve_draw(FILE *outfile, plCurve *L);

/* Display Routines */

void	init_display();
void	shutdown_display();
void    refresh_display(plCurve *L);
void	export_ted(plCurve* inLink, octrope_strut* strutSet, 
		   int inSize, octrope_mrloc* minradSet, int minradLocs, 
		   double* compressions, search_state* inState);

void	export_struts(plCurve* inLink, octrope_strut* inStruts, 
		      int inSize, double* compressions, search_state* inState);
void	exportVect( const plc_vector* dl, plCurve* link, const char* fname );


#endif
