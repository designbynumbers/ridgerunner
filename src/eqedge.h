
/*

   eqedge.h

   Header file for demo equilateralization code using
   LAPACK.  

*/

/************************* Include Files **************************/

#include "portability.h"
#include "plCurve.h"
#include "plc_vector.h"
#include "stepper.h"

/************************* Global Variables ***********************/

int VERBOSITY;
int GRAPHICS;
FILE *gclpipe;

/************************* Function Prototypes ********************/

/* Utility Routines */

void our_matrix_write(double val, double *A, int LDA, int i, int j);
double maxError(plCurve *L, int comp, search_state* inState);

void	       link_scale( plCurve* inLink, double factor );
plCurve*  octrope_fixlength( plCurve* inLink );
plCurve*  octrope_equilize_density( plCurve* inLink, search_state* inState );
plCurve*  octrope_double_edges( plCurve* inLink );
plCurve*  octrope_double_component( plCurve* inLink, int comp );

double	       plCurve_torsion( plCurve* inLink, FILE* outPlot );

int    plc_strand_edges(plc_strand *P);
double plc_strand_length(plc_strand *P); 
double plCurve_short_edge(plCurve *L);
double plCurve_long_edge(plCurve *L);

plc_vector plCurve_tangent_vector(plCurve *L,int comp,int vert);
plc_vector plCurve_edge_dir(plCurve *L,int comp,int edge);

void   plCurve_draw(FILE *outfile, plCurve *L);

/* Display Routines */

void	init_display();
void	shutdown_display();
void   refresh_display(plCurve *L);
void	export_ted(plCurve* inLink, octrope_strut* strutSet, 
			int inSize, octrope_mrloc* minradSet, int minradLocs, 
			double* compressions, search_state* inState);
			
void	export_struts(plCurve* inLink, octrope_strut* inStruts, int inSize, double* compressions, search_state* inState);
void	exportVect( const plc_vector* dl, plCurve* link, const char* fname );

/* Computational Routines. */

double *tangential_rigidity_matrix(plCurve *L, int comp, int *LDA, int *N);
double *edgelength_error(plCurve *L,int comp, search_state* inState);

/* External interface. */

void lapack_eqedge(plCurve *L, int nsteps, search_state* inState);

/* Demos. */

void segment_demo();
void ellipse_demo();




