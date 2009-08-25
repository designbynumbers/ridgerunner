/*

   ridgerunner.h :

   Part of main ridgerunner project. Defines function prototypes and
   global variables for entire ridgerunner project. 

*/

#ifndef _H_ridgerunner
#define _H_ridgerunner

#include "ncurses.h"         

/* Ncurses rudely stomps the "bool" type without checking to see if it's already */
/* defined (say, by stdbool.h). As a result, it must appear before we include    */
/* portability (and hence stdbool). */

#include "portability.h"
#include "errors.h"

#include "plCurve.h"
#include "octrope.h"
#include "libtsnnls/tsnnls.h"
#include "libtsnnls/lsqr.h"
#include "argtable2.h"


#define DEBUG 1 
/* Turn on all the asserts in the code. */

// #define CURSES_DISPLAY 1

#define kStepScale 0.01
//#define kMinStepSize 1e-5
//#define kMaxStepSize 1e-3

#define kMinStepSize 1e-8
//#define kMaxStepSize 1e-2


/************************* Global Variables ***********************/

extern int VERBOSITY;
extern FILE *gLogfile;
extern FILE *gclpipe;

extern int gSuppressOutput;
extern int gQuiet;
extern double gLambda;                  /* lambda-stiffness of rope */
extern int gMaxCorrectionAttempts;
extern int gNoRcond;  /* Used to suppress rcond calls on systems with buggy ATLAS. */
extern int gLsqrLogging; /* Used to enable lsqr output logging. */
extern int gNoTimeWarp; /* If true, turns off "timewarp" acceleration of free vertices */
extern int gEqIt; /* If true, runs the equilaterization code every time max/min > 3. */
extern int gAnimationStepper;

#ifdef CURSES_DISPLAY

WINDOW *gLogwin;

#endif

extern int gNumTubeColors;
extern plc_color gTubeColors[5];
extern plc_color gStraightSegColor;
extern plc_color gKinkColor;
extern plc_color gHelixColor;



/************* Defined Data Types ***********************/

#define LOG_FLUSH_INTERVAL 1000

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
  klsqrlog,
  kMemused,  // total amount of memory allocated
  kEffectiveness, // a measure of how successful we are in avoiding error
  
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
  int    stopTime;    // Maximum time for the run, in seconds.

  double moviefactor;           // multiplier determines how often to save movie frames. 
  int    maxmovieframes;        // maximum number of movie frames to save.

  short	saveConvergence;    // saving convergence info?
  short	movie;		    // are we generating movie frames?
  short	fancyVisualization; // fancy visualization stuff w/ geomview, gnuplot, and tube
  FILE*	fancyPipe;	    // geomview pipe if we're fancy
  
  char	fname[1024];           // full filename of input file:      ../data/3.1.vect
  char  basename[1024];        // name of input file w/o extension: 3.1  
  char  workingfilename[1024]; // name for intermediate output files
  char  workingstrutname[1024];// name for intermediate strut files
  char  logfilename[1024];
  char  vectprefix[1024];
  char  fprefix[1024];
  char  snapprefix[1024];      // prefix for "snapshots" of given steps
  
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

  int     cstep_count;

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
  
  octrope_strut* lastStepStruts; // the struts from the last successful step or cstep, 
                                 // or null at start
  int		lastStepStrutCount;

  octrope_mrloc* lastStepMRlist; // minrad locs from last completed step or cstep
  int		lastStepMinradStrutCount;

  double        lastStepPocaEffectiveness;  // The error introduced in pocas by the last 
  double        lastStepMREffectiveness;    // step, as a fraction of the error that would
                                            // have been introduced by taking a step of the 
                                            // same (L^2) size in the dLen direction.
  int		totalVerts;
  int		totalSides;
  
  int		tsnnls_evaluations;
  int           octrope_calls;
  
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

  int     snapinterval; // save a complete "snapshot" of the computation every 
                        // snapinterval steps.

  int     maxlogsize;   // Maximum logfile size (entries).
  int     loginterval;  // Log data every nth step.

  FILE    *logfiles[128];     /* The logfiles hold the various data that can be recorded.*/
  char    *logfilenames[128]; /* These buffers hold the names of the log files. */ 

  int     ndisplay;           /* Number of data items to display and data to display */
  enum GraphTypes runtime_display[128]; 

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

int correct_thickness(plCurve *inLink,search_state *inState);
/* Newton's method thickness correction algorithm */

void correct_constraints(plCurve *inLink,search_state *inState);
/* Correct the position of the link so that constraints are obeyed. */

taucs_ccs_matrix *buildRigidityMatrix(plCurve *inLink,search_state *inState);
/* Creates rigidity matrix corresponding to current set of struts, kinks, constraints. */

int dlenPos( plCurve *inLink,int cmp,int vt);
/* Converts a position on a link to a "flat" address in the dlen buffer. */

void dlenForce( plc_vector* ioDL, plCurve* inLink, search_state* inState );
/* Adds gradient of length to the buffer ioDL. */

void eqForce( plc_vector* dlen, plCurve* inLink, search_state* inState );
/* Adds an "equilateralization force" to the buffer dlen. */

void spinForce( plc_vector* dlen, plCurve* inLink, search_state* inState );
/* Adds a "spin force" to the buffer dLen. */

void specialForce( plc_vector* dlen, plCurve* inLink, search_state* inState );
/* A stub, used in future versions to add other forces to dlen. */

void constraintForce( plc_vector* dlen, plCurve* inLink, search_state* inState );
/* Alters the force to make sure that it does not attempt to violate constraints. */

plc_vector *resolveForce( plc_vector* dl, plCurve* inLink, search_state* inState);
/* Uses rigidity matrix to resolve the force dl over struts, kinks, and constraints. */

void open_runtime_logs(search_state *state, char opentype);
/* Opens the runtime logs */

void close_runtime_logs(search_state *state);
/* Closes the runtime logs. */

void compress_runtime_logs(search_state *state);
/* Compresses the various runtime logs to save disk space during a long run, 
   if we need to save space. */

void update_runtime_logs(search_state *state);
/* Writes data on current run to various log files. */

void update_vect_directory(plCurve * const link, search_state *state);
void compress_vectdir(const search_state *inState);

void free_search_state(search_state *inState);

void reloadDump( double* A, int rows, int cols, double* x, double* b );

void our_matrix_write(double val, double *A, int LDA, int i, int j);
double maxError(plCurve *L, int comp, search_state* inState);

/* linklib_additions.c */

plCurve*  octrope_fixlength( plCurve* inLink );
plCurve*  plCurve_fixresolution( plCurve* inLink, double newres );

double	  plCurve_torsion( plCurve* inLink, FILE* outPlot );

double plCurve_short_edge(plCurve *L);
double plCurve_long_edge(plCurve *L);

plc_vector plCurve_tangent_vector(plCurve *L,int comp,int vert);
plc_vector plCurve_edge_dir(plCurve *L,int comp,int edge);

int plCurve_score_constraints(plCurve *inLink);

void   plCurve_draw(FILE *outfile, plCurve *L);
void   plc_color_curve(plCurve *L,plc_color col);

/* Display Routines */

void init_runtime_display(search_state *inState);
void update_runtime_display(plCurve *inLink,search_state *inState);
void close_runtime_display();

void init_gv_display();
void shutdown_gv_display();
void refresh_display(plCurve *L);

void strut_vectfile_write(plCurve *inLink, octrope_strut *strutlist, 
			  int strutCount, FILE *fp);

void export_struts(plCurve* inLink, octrope_strut* inStruts, 
		   int inSize, double* compressions, search_state* inState);

plCurve *vectorfield_to_plCurve(plc_vector *vf, plCurve *inLink);

/* Error Handling Routines. */

void FatalError(char *debugmsg,const char *file,int line);
void NonFatalError(char *debugmsg,const char *file,int line);

void dumpAxb_full( search_state *inState, 
		   double* A, int rows, int cols, 
		   double* x, double* b );
void dumpAxb_sparse( search_state *inState, taucs_ccs_matrix* A, 
		     double* x, double* b );
void dumpVertsStruts(plCurve* link, octrope_strut* strutSet, int strutCount);

void dumpLink( plCurve *inLink, search_state *inState, char *dumpname);
void dumpNamedLink( plCurve *inLink, search_state *inState, char *tag, char *dumpname);
void dumpStruts( plCurve *inLink, search_state *inState, char *dumpname);
void dumpDvdt( plc_vector* dvdt, plCurve *inLink, search_state *inState );
void dumpdLen( plc_vector* dLen, plCurve *inLink, search_state *inState );

void snapshot( plCurve *inLink,
	       plc_vector *dVdt,plc_vector *dlen,
	       search_state *inState );

double rigidityEntry( taucs_ccs_matrix* A, int strutCount, 
		      int minradLocs, int row, int col );
void checkDuplicates( octrope_strut* struts, int num );
void collapseStruts( octrope_strut** struts, int* count );

FILE *fopen_or_die(const char *filename,const char *mode,
		   const char *file,const int line); 
int   system_or_die(char *cmdline,const char *file,int line);
void *malloc_or_die(size_t size, const char *file, const int line);
void  remove_or_die(char *filename,const char *file, const int line);
void  rename_or_die(char *oldname,char *newname,const char *file, const int line);
int   mkstemp_or_die(char *template,const char *file, const int line);
FILE *fdopen_or_die(int fd, const char *opentype, const char *file, const int line);
DIR  *opendir_or_die(char *dirname,const char *file, const int line);

void  logprintf(char *format, ... );
/* Prints to stdout and to the system log. */

/* Code for handling "free" regions of the curve */

int strut_free_vertices( plCurve* inLink, double tube_radius, bool *freeFlag);
/* Fills a "flat" buffer freeFlag, expected to be plc_num_verts(inLink) in size, 
   with "true" if the vertex is no closer than 4 vertices to a strut and "false" otherwise.
   The buffer refers to vertices on the plCurve inlink in dictionary order on (cp,vt). */

void highlight_curve( plCurve *L, search_state *state );
/* Highlight straight segments, kinks, and other understood regions of the curve. */

void accelerate_free_vertices( plc_vector *dLen, plCurve *L, search_state *inState);
/* Scales up the dLen force on vertices that are not close to any vertex with a strut. */

#endif
