/*
 *  stepper.h
 *  ridgerunner
 *
 *  Created by Michael Piatek on Fri May 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "octrope.h"

#ifndef _H_stepper
#define _H_stepper

#define kStepScale 0.01
#define kMinStepSize 1e-4

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
	kRcond,			// reciprocal condition number as estimated by LAPACK -- note this slows things down a bit
	kWallTime,		// the wall time since the beginnging of run, wall time vs evolution time is useful to see where we're spending computation time
	kMaxVertexForce,
	kConvergence,   // the ropelength vs the number of steps
	
	kTotalGraphTypes
};

typedef struct
{
	double  injrad;				// the user specified injectivity radius of the link
	double  minrad;				// Rawdon's minimum radius of curvature
	double  ropelength;			// the ropelength of the link
	double  length;				// the polygonal length of the link
	double  shortest;			// the shortest strut on last firstVariation call
	
	double  factor;				// curvature scaling factor
	
	double  stepSize;			// the current step size of our stated run
	unsigned int   steps;		// the number of steps we've taken
	double  maxStepSize;		// as high as it's allowed to get
	
	double  time;
	
	double  eqThreshold;		// the max/min edge value at which to eq. 0 means never eq.
	
	int*	compOffsets;		// we need to know the component offsets 
								// to index struts when there are multiple components
	
	int*	conserveLength;		// array of size link->nc that specifies whether or not to treat
								// the given component's edges as bars in the rigidity matrix
	
	octrope_strut* lastStepStruts; // the struts from the last successful step, or null at start
	int		lastStepStrutCount;
	int		lastStepMinradStrutCount;
	
	int		totalVerts;
	int		totalSides;
	
	int		tsnnls_evaluations;
	int		curvature_step;
	
	double* sideLengths;
	
	double  rcond;				// rcond from the last step, nonzero only if graphing[kRcond] == 1 && struts > 0
	
	int		graphing[kTotalGraphTypes];
	
	double  maxPush;		// the maximum of all the strut force compression sums
	double  avgDvdtMag;   // avg of all dvdt norms
	double  residual;		// residual of last snnls call
	
} search_state;

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
void bsearch_stepper( octrope_link** inLink, unsigned int inMaxSteps, search_state* inState );

double				maxovermin( octrope_link* inLink );

void	updateSideLengths( octrope_link* inLink, search_state* inState );

void reloadDump( double* A, int rows, int cols, double* x, double* b );

#endif // _H_stepper

