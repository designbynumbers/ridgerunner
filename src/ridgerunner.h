/*

   ridgerunner.h :

   Part of main ridgerunner project. Defines function prototypes and
   global variables for entire ridgerunner project. 

*/

#ifndef _H_ridgerunner
#define _H_ridgerunner

#include "portability.h"
#include "errors.h"

#include "plCurve.h"
#include "octrope.h"
#include "libtsnnls/tsnnls.h"
#include "argtable2.h"

#define kStepScale 0.01
//#define kMinStepSize 1e-5
//#define kMaxStepSize 1e-3

#define kMinStepSize 1e-8
//#define kMaxStepSize 1e-2


/************************* Global Variables ***********************/

int VERBOSITY;
int GRAPHICS;
FILE *gclpipe,*gLogfile;

int gVerboseFiling = 0;
int gSuppressOutput = 0;
int gQuiet = 0;
int gSurfaceBuilding = 0;
int gAvoidTmpConflicts = 0;
int gPaperInfoInTmp = 0;
int gFastCorrectionSteps = 0;
double gLambda = 1.0;   /* lambda-stiffness of rope */
int gMaxCorrectionAttempts = 25;

/************* Defined Data Types ***********************/

#define LOG_FLUSH_INTERVAL 10

enum GraphTypes
{
  kLength = 0,
  kRopelength,
  kStrutCount,
  kStepSize,
  kThickness,
  kMinrad, 
  kResidual,
  kMaxOverMin,
  kRcond,	 // reciprocal condition number as estimated by LAPACK
		 // -- note this slows things down a bit
  kWallTime,		// the wall time since the beginning of run,
			// wall time vs evolution time is useful to
			// see where we're spending computation time
  kMaxVertexForce,
  kCorrectionStepsNeeded, // the number of correction steps required to converge,
  kEQVariance,	// the variance of the set of (edge lengths - the average)
  
  kTotalLogTypes
};

/* "State" a long data structure for holding everything related to the 
   current state of the computation. */

typedef struct
{
  /* Stopping criteria. */
  
  long	 maxItrs;     // maximum number of steps to perform before quit
  double stop20;     // must decrease by this amount per 20 steps to continue
  double residualThreshold;	// threshold for residual stopping, will use if nonzero

  short	saveConvergence;    // saving convergence info?
  short	movie;		    // are we generating movie frames?
  short	fancyVisualization; // fancy visualization stuff w/ geomview, gnuplot, and tube
  FILE*	fancyPipe;	    // geomview pipe if we're fancy
  
  char	fname[1024];
  char  finalfilename[1024];
  char  finalstrutname[1024];
  char  logfilename[1024];
  char  vectprefix[1024];
  char  fprefix[1024];
  
  /* Correction step information. */

  double correctionStepDefault;	// default correction step size
 
  double  overstepTol;		        // the amount we are willing to
                                        // overstep before correcting strut length
  double  minradOverstepTol;
  
  double  tube_radius;	 // the user specified tube radius of the link
                         // we enforce the constraint 
   
                         //     thickness >= tube_radius

                         // throughout the run.

  int     last_cstep_attempts; // # of attempts required to 
                               // converge in last round of 
                               // correction stepping.

  int     last_step_attempts; // # of attempts required to find correct bsearch
                              // stepsize in last curvature step. 

  /* Geometric data about the current state of the link. */

  double  minrad;	 // Rawdon's minimum radius of curvature
  double  ropelength;	 // the ropelength of the link
  double  length;	 // the polygonal length of the link
  double  shortest;	 // the shortest strut on last firstVariation call
  double  thickness;     // the overall thickness of the link

  double  factor;	 // curvature scaling factor

  /* Data about the current run. */

  unsigned int   steps;	 // the number of steps we've taken

  double  stepSize;	 // the current step size of our stated run
  double  maxStepSize;	 // as high as it's allowed to get
  
  double  time;
  double  cstep_time;
  double  lastMaxMin;	  // max/min of last step
  
  int*	compOffsets;	  
  
  /* At many points in the program we'll convert back and forth from a
     "flat" representation of a link, where a single array of
     plc_vectors will represent all vertices of all components of the
     link in order, or a single array of doubles will represent all of
     the components of all these vectors. The "compOffsets" array
     contains pointers into such an array _of plc_vectors_, where 

     compOffset[i] is the 0th vertex of component i. */

  int*	conserveLength;	  // array of size link->nc that specifies 
                          // whether or not to treat
                          // the given component's edges as bars 
                          // in the rigidity matrix
  
  octrope_strut* lastStepStruts; // the struts from the last successful step, 
                                 // or null at start
  int		lastStepStrutCount;

  octrope_mrloc* lastStepMRlist; // minrad locs from last completed step
  int		lastStepMinradStrutCount;
  
  int		totalVerts;
  int		totalSides;
  
  int		tsnnls_evaluations;
  
  double*       sideLengths;
  double	avgSideLength;
  
  double*	perVertexResidual;   // if we're recording paper info... 
                                     // we do this and a whole lot of other 
                                     // time consuming stuff
  
  double  rcond; // rcond from the last step, nonzero only 
                 // if graphing[kRcond] == 1 && struts > 0
  
  double  ofvNorm;	// ofvB L2
  double  ofvResidual;	// residual of the ofv lsqr call
  
  int	  graphing[kTotalLogTypes];
  
  double  maxPush;	// the maximum of all the strut force compression sums
  double  avgDvdtMag;   // avg of all dvdt norms
  double  residual;	// residual of last snnls call
  
  double  oldLength;
  double  oldLengthTime;
  /*  double  checkDelta; */   // the amount of time to wait to 
  // check for checkThreshold progress for stopping
  /*  double  checkThreshold; */ // the amount to check for */
  
  /* This have been replaced by stop20 and maxItrs */

  double  minminrad;
  
  double  eqMultiplier;	// scale of eq force, increased as things get less and less eq
  double  eqVariance;	// if graphing eq variance, this will be set to 
                        // the variance of the set of edge length difference 
                        // from the average

  double  eqAvgDiff;

  FILE    *logfiles[128]; /* The logfiles hold the various data that can be recorded.*/
  char    *logfilenames[128]; /* These buffers hold the names of the log files. */ 

#ifdef HAVE_TIME
  time_t  start_time;
#endif

} search_state;

typedef struct
{
  double  rescaleThreshold;
} RSettings;

extern RSettings gRidgeSettings;

/****************************************************************/
/* Function Prototypes                                          */
/****************************************************************/

/* stepper.c */

void bsearch_stepper( plCurve** inLink, search_state* inState );
/*
 * The main loop of the program. Calls bsearch_step until one of
 
 1) inState->maxItrs steps are performed
 2) the change in ropelength over the last 20 steps is less than inState->stop20
 3) inState->residual < inState->residualThreshold
 
 *
 */

void correct_thickness(plCurve *inLink,search_state *inState);
/* Newton's method thickness correction algorithm */

void update_runtime_logs(search_state *state);
/* Writes data on current run to various log files. */

void update_vect_directory(plCurve * const link, const search_state *state);

void free_search_state(search_state *inState);

double maxovermin( plCurve* inLink, search_state* inState );
void updateSideLengths( plCurve* inLink, search_state* inState );
void reloadDump( double* A, int rows, int cols, double* x, double* b );

void our_matrix_write(double val, double *A, int LDA, int i, int j);
double maxError(plCurve *L, int comp, search_state* inState);

plCurve*  octrope_fixlength( plCurve* inLink );
plCurve*  plCurve_fixresolution( plCurve* inLink, double newres );

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

void    strut_vectfile_write(plCurve *inLink, octrope_strut *strutlist, 
			     int strutCount, FILE *fp);

void	export_ted(plCurve* inLink, octrope_strut* strutSet, 
		   int inSize, octrope_mrloc* minradSet, int minradLocs, 
		   double* compressions, search_state* inState);

void	export_struts(plCurve* inLink, octrope_strut* inStruts, 
		      int inSize, double* compressions, search_state* inState);
void	exportVect( const plc_vector* dl, plCurve* link, const char* fname );

void    preptmpname( char* outName, const char* inName, search_state* inState );

/* Error Handling Routines. */

void FatalError(char *debugmsg,const char *file,int line);
void dumpAxb_full( search_state *inState, 
		   double* A, int rows, int cols, 
		   double* x, double* b );
void dumpAxb_sparse( search_state *inState, taucs_ccs_matrix* A, double* x, double* b );
void dumpVertsStruts(plCurve* link, octrope_strut* strutSet, int strutCount);
void checkDuplicates( octrope_strut* struts, int num );
void collapseStruts( octrope_strut** struts, int* count );

FILE *fopen_or_die(const char *filename,const char *mode,
		   const char *file,const int line); 

#endif
