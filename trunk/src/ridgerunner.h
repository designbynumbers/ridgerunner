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

#define kStepScale 0.01
//#define kMinStepSize 1e-5
//#define kMaxStepSize 1e-3

#define kMinStepSize 1e-8
//#define kMaxStepSize 1e-2


/************************* Global Variables ***********************/

int VERBOSITY;
int GRAPHICS;
FILE *gclpipe;

int gVerboseFiling = 0;
int gSuppressOutput = 0;
int gQuiet = 0;
int gSurfaceBuilding = 0;
int gAvoidTmpConflicts = 0;
int gPaperInfoInTmp = 0;
int gFastCorrectionSteps = 0;
double gLambda = 1.0;   /* lambda-stiffness of rope */

/************* Defined Data Types ***********************/

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
  kConvergence,   // the ropelength vs the number of steps,
  kEQVariance,	// the variance of the set of (edge lengths - the average)
  
  kTotalGraphTypes
};

/* "State" a long data structure for holding everything related to the 
   current state of the computation. */

typedef struct
{
  long	maxItrs;     // maximum number of steps to perform before quit
  
  double stop20;     // must decrease by this amount per 20 steps to continue
  
  short	saveConvergence;    // saving convergence info?
  short	movie;		    // are we generating movie frames?
  short	fancyVisualization; // fancy visualization stuff w/ geomview, gnuplot, and tube
  FILE*	fancyPipe;	    // geomview pipe if we're fancy
  
  char	fname[1024];
  char  finalfilename[1024];
  char  finalstrutname[1024];
  char  logfilename[1024];
  char  vectprefix[1024];
  
  double  correctionStepDefault;	// default correction step size
  double  overstepTol;		        // the amount we are willing to
                                        // overstep before correcting strut length
  double  minradOverstepTol;
  
  double  injrad;	 // the user specified injectivity radius of the link
  double  minrad;	 // Rawdon's minimum radius of curvature
  double  ropelength;	 // the ropelength of the link
  double  length;	 // the polygonal length of the link
  double  shortest;	 // the shortest strut on last firstVariation call
  
  double  factor;	 // curvature scaling factor
  
  double  stepSize;	 // the current step size of our stated run
  unsigned int   steps;	 // the number of steps we've taken
  double  maxStepSize;	 // as high as it's allowed to get
  
  double  time;
  double  cstep_time;
  
  double  eqThreshold;	  // the max/min edge value at which to eq. 0 means never eq.
  double  lastMaxMin;	  // max/min of last step
  
  int*	compOffsets;	  // we need to know the component offsets 
  // to index struts when there are multiple components
  
  int*	conserveLength;	  // array of size link->nc that specifies 
                          // whether or not to treat
                          // the given component's edges as bars 
                          // in the rigidity matrix
  
  octrope_strut* lastStepStruts; // the struts from the last successful step, 
                                 // or null at start
  int		lastStepStrutCount;
  int		lastStepMinradStrutCount;
  
  int		totalVerts;
  int		totalSides;
  
  int		tsnnls_evaluations;
  int		curvature_step;
  int		eq_step;
  
  double*       sideLengths;
  double	avgSideLength;
  
  double*	perVertexResidual;   // if we're recording paper info... 
                                     // we do this and a whole lot of other 
                                     // time consuming stuff
  
  double  rcond; // rcond from the last step, nonzero only 
                 // if graphing[kRcond] == 1 && struts > 0
  
  double	ofvNorm;	// ofvB L2
  double	ofvResidual;	// residual of the ofv lsqr call
  
  int		graphing[kTotalGraphTypes];
  
  double  maxPush;	// the maximum of all the strut force compression sums
  double  avgDvdtMag;   // avg of all dvdt norms
  double  residual;	// residual of last snnls call
  
  double  oldLength;
  double  oldLengthTime;
  double  checkDelta;	// the amount of time to wait to 
                        // check for checkThreshold progress for stopping
  double  checkThreshold; // the amount to check for 
  
  double  minminrad;
  
  double  eqMultiplier;	// scale of eq force, increased as things get less and less eq
  double  eqVariance;	// if graphing eq variance, this will be set to 
                        // the variance of the set of edge length difference 
                        // from the average

  double  eqAvgDiff;
  double  residualThreshold;	// threshold for residual stopping, will use if nonzero
  int	  refineUntil;	// keep running until discretization refineuntil x ropelength
  int	  ignore_minrad; // control minrad

  FILE    (*logfiles)[128]; /* The logfiles hold the various data that can be recorded.*/

  char    *logfilenames[32]; 
  int     nlogs;

  struct rusage	frameStart;
  
} search_state;

typedef struct
{
  double  rescaleThreshold;
} RSettings;

extern RSettings gRidgeSettings;

/*
 * Performs inMaxSteps SUCCESSFUL gradient steps on *inLink and stores the result in 
 * **inLink. Keep in mind that the initial *inLink will be freed in this 
 * process. Successful gradient steps are those that either move the curve by 
 * dLen with no changes to the strut set or that create only a single strut.
 *
 * An equilateralization flow is also performed after each successful step to keep
 * the link approximately equilateral.
 *
 */

void bsearch_stepper( plCurve** inLink, search_state* inState );

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
void	export_ted(plCurve* inLink, octrope_strut* strutSet, 
		   int inSize, octrope_mrloc* minradSet, int minradLocs, 
		   double* compressions, search_state* inState);

void	export_struts(plCurve* inLink, octrope_strut* inStruts, 
		      int inSize, double* compressions, search_state* inState);
void	exportVect( const plc_vector* dl, plCurve* link, const char* fname );

void preptmpname( char* outName, const char* inName, search_state* inState );

#endif
