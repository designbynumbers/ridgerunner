/*
 *  stepper.c
 *  ridgerunner
 *
 *  Created by Michael Piatek on Fri May 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "ridgerunner.h"

/* Function prototypes */

plCurve*  bsearch_step( plCurve* inLink, search_state* inState );
void	  step( plCurve* inLink, double stepSize, plc_vector* dVdt, 
		search_state* inState );
void	  firstVariation( plc_vector* inOutDvdt, plCurve* inLink, search_state* inState,
			  octrope_strut** outStruts, int* outStrutsCount, int dlenStep);
void	  computeCompressPush( plCurve* inLink, octrope_strut* strutSet,
			       octrope_mrloc* minradSet, int strutCount, 
			       int minradLocs );

extern void sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod );	  	
extern void fast_lsqr_mult( long mode, dvec* x, dvec* y, void* prod );		

// lesser utility functions
int    equalStruts( const octrope_strut* s1, const octrope_strut* s2 );
void   normalizeStruts( plc_vector* strutDirections, octrope_strut* strutSet, 
			plCurve* inLink, int strutCount );

void   placeVertexBars( double* A, plCurve* inLink, int contactStruts, 
			int totalBarVerts, int totalBars, search_state* inState );

/* Global variables */

int displayEveryFrame = 0;

extern int gQuiet;
extern int gSurfaceBuilding;
extern int gPaperInfoInTmp;
extern double gLambda;

static void
swap_rigidity_verts( taucs_ccs_matrix* A, int v1, int v2, int col )
{
  /* here, verts are given in offset from 0 of the 4 or 3 per strut */
  double vals[3];
  double inds[3];
  int j;
	
  for( j=0; j<3; j++ ) {

    vals[j] = A->values.d[A->colptr[col] + 3*v1 + j];
    inds[j] = A->rowind[A->colptr[col] + 3*v1 + j];
 
  }
  
  for( j=0; j<3; j++ ) {

    A->values.d[A->colptr[col] + 3*v1 + j] = A->values.d[A->colptr[col] + 3*v2 + j];
    A->rowind[A->colptr[col] + 3*v1 + j] = A->rowind[A->colptr[col] + 3*v2 + j];
  
  }
	
  for( j=0; j<3; j++ ) {

    A->values.d[A->colptr[col] + 3*v2 + j] = vals[j];
    A->rowind[A->colptr[col] + 3*v2 + j] = inds[j];
  
  }
}

static void
taucs_enforce_ccs_sort( taucs_ccs_matrix* A, int strutCount )
{
  int cItr, rightBound, i;
  for( cItr=0; cItr<A->n; cItr++ ) {

    int col = cItr;
    // we need to move things in blocks of three.
    if( cItr < strutCount ) {

      for( rightBound = 3; rightBound > 0; rightBound-- ) {

	for( i=0; i<rightBound; i++ ) {

	  if(A->rowind[A->colptr[col]+3*i] > A->rowind[A->colptr[col]+3*(i+1)])
	    swap_rigidity_verts(A, i, i+1, col);
	}
      }
    }// < strutCount? (thickness strut)
    else {
      // bubble sort with three vertices.
      if( A->rowind[A->colptr[col]] > A->rowind[A->colptr[col]+3] )
	swap_rigidity_verts(A, 0, 1, col);
      if( A->rowind[A->colptr[col]+3] > A->rowind[A->colptr[col]+6] )
	swap_rigidity_verts(A, 1, 2, col);
      if( A->rowind[A->colptr[col]] > A->rowind[A->colptr[col]+3] )
	swap_rigidity_verts(A, 0, 1, col);
    }
    
  } // for over columns
}

// note to self -- get rid of THIS later too
static double
rigidityEntry( taucs_ccs_matrix* A, int strutCount, int minradLocs, int row, int col )
{
  int vItr;
  
  // the most common case will be a 0, so get that out of the way quick --
  // movd this check to firstVariation for speed purposes
  
  // also note that this exploits the t_snnls enforced ordering of row
  // entries and doesn't work for generic ccs matrices
  /*	if( row < A->rowind[A->colptr[col]] || row > A->rowind[A->colptr[col+1]-1] )
	return 0;
  */
  // now we look.
  if( col < strutCount ) {

    // we can manually unroll this loop -- data dependencies in here
    // result in large slowdown according to shark.
    /*		for( vItr=A->colptr[col]; vItr<A->colptr[col+1]; vItr++ )
		{
		if( A->rowind[vItr] == row )
		return A->values.d[vItr];
		}
    */
    /* Remember that if this column represents a strut, there at most
       12 rows in the column (less only if the strut is vertex-vertex,
       or vertex-edge).  So we need only search through 12 possible
       rowind values looking for the desired row. 

       If there are too few rows in this column, the search could in 
       principle wrap around and return a row from the next column.*/

    vItr=A->colptr[col];

    if( A->rowind[vItr+0] == row ) return A->values.d[vItr+0];
    if( A->rowind[vItr+1] == row ) return A->values.d[vItr+1];
    if( A->rowind[vItr+2] == row ) return A->values.d[vItr+2];
    if( A->rowind[vItr+3] == row ) return A->values.d[vItr+3];
    
    if( A->rowind[vItr+4] == row ) return A->values.d[vItr+4];
    if( A->rowind[vItr+5] == row ) return A->values.d[vItr+5];
    if( A->rowind[vItr+6] == row ) return A->values.d[vItr+6];
    if( A->rowind[vItr+7] == row ) return A->values.d[vItr+7];
    
    if( A->rowind[vItr+8] == row ) return A->values.d[vItr+8];
    if( A->rowind[vItr+9] == row ) return A->values.d[vItr+9];
    if( A->rowind[vItr+10] == row ) return A->values.d[vItr+10];
    if( A->rowind[vItr+11] == row ) return A->values.d[vItr+11];
 
  } else {

    /* This is not a strut constraint, so we're not sure how many 
       columns to search. Check them all just to be sure. */
    
    for( vItr=A->colptr[col]; vItr<A->colptr[col+1]; vItr++ ) {

      if( A->rowind[vItr] == row )
	return A->values.d[vItr];
    }
  }

  /* We should never reach this point in the code. */

  char errmsg[1024];

  sprintf(errmsg,
	  "ridgerunner: illegal call to rigidityEntry\n"
	  "             tried for entry (%d,%d) of %d x %d matrix A\n"
	  "             strutcount %d\n"
	  "             minradLocs %d\n",
	  row,col,A->n,A->m,strutCount,minradLocs);

  FatalError(errmsg, __FILE__ , __LINE__ );

  return 0;
}

int gOutputFlag = 1;
int gConditionCheck = 0;
extern int gFastCorrectionSteps;

#define kOutputItrs 1
#define SECS(tv)        (tv.tv_sec + tv.tv_usec / 1000000.0)

int	gEQevents = 0;


extern int gSuppressOutput;
extern int gVerboseFiling;

int	gCorrectionAttempts = 0;

void 
bsearch_stepper( plCurve** inLink, search_state* inState )
{
  unsigned int i, offset = 0;
  int cItr, vItr;
  int lastEq = 0;
  unsigned int cSteps = 0;
  
  double maxmaxmin = 0;
  double minthickness = 500;
  double nextMovieOutput = 0.0;
  
  gConditionCheck = 0;
  int firstRun = 1, secondRun = 0;
  double rop_20_itrs_ago = {DBL_MAX};
  
#ifdef HAVE_CLOCK  
  
  clock_t startTime;
  startTime = clock();

#endif
  
  inState->steps = 0;
  int stepItr;
  
  for( stepItr=0; /* Main loop, incorporates stopping criteria. */
       (stepItr < inState->maxItrs) && 
	 (rop_20_itrs_ago - inState->ropelength > inState->stop20) &&
	 (inState->residual > inState->residualThreshold);
       stepItr++ ) {
    
    int lastSet;      
    lastSet = inState->lastStepStrutCount;
    
    /* Decide whether to initiate thickness correction. */
    
    if( (inState->shortest < (2*inState->tube_radius*(1-inState->overstepTol))) ||
	(inState->minrad < gLambda*inState->minradOverstepTol) ) {

      correct_thickness(link,inState);

    }        
  
    /************************************************************************/
    
    *inLink = curvature_step(*inLink, inState);
    
    /************************************************************************/
    
    inState->steps++;
    
    /* Manage data output */
      
    if( (stepItr%kOutputItrs) == 0 && gQuiet == 0 ) { /* Note "k" means constant */
      
      update_runtime_display(inState);  
      
    }

    update_runtime_logs(inState);

    if (inState->time >= nextMovieOutput) {

      update_vect_directory(inState);
      nextMovieOutput += 0.041666666667;

    }

    /* Last, we update "inState" to keep track of running variables. */

    if( inState->shortest < minthickness )
      minthickness = inState->shortest;

    inState->oldLength = inState->length;
    inState->oldLengthTime = inState->cstep_time;

    double maxmin;
    maxmin = maxovermin(*inLink, inState);
    inState->lastMaxMin = maxmin;
    if( maxmin > maxmaxmin )
      maxmaxmin = maxmin;
    
  } 

  /* We have now terminated. The final output files will be written 
     in ridgerunner_main.c */

}


void update_vect_directory(const plCurve *link, const search_state *state)

     /* If it is time for the next vect output, go ahead and write 
	another file to the appropriate directory. */
{
  char tmpfilename[1024];
  char tmpfullname[1024];
  FILE *outfile;
  
  sprintf(tmpfilename,"%s%s.%7d.vect",state->vectprefix,state->fname,state->steps);
  outfile = fopen_or_die(tmpfilename,"w", __FILE__ , __LINE__ );
  plc_write(outfile,link);
  fclose(outfile);
  
}

void update_runtime_logs(search_state *state)

     /* We now update the various data logs for the run. */
     /* These logs are flushed every LOG_FLUSH_INTERVAL steps, */
     /* where this is defined in ridgerunner.h. */
{
  
  fprintf(state->logfiles[kLength],"%g \n",state->length);
  fprintf(state->logfiles[kRopelength],"%g \n",state->ropelength);
  fprintf(state->logfiles[kStrutCount],"%d %d\n",
	  state->lastStepStruts,state->lastStepMinradStrutCount);
  fprintf(state->logfiles[kStepSize],"%g \n",state->stepSize);
  fprintf(state->logfiles[kThickness],"%g \n",state->thickness);
  fprintf(state->logfiles[kMinrad],"%g \n",state->minrad);
  fprintf(state->logfiles[kResidual],"%g \n",state->residual);
  fprintf(state->logfiles[kMaxOverMin],"%g \n",state->lastMaxMin);
  fprintf(state->logfiles[kRcond],"%g \n",state->rcond);

#ifdef HAVE_CLOCK
  fprintf(state->logfiles[kWallTime],"%g \n",clock()/CLK_TCK);
#endif

  fprintf(state->logfiles[kMaxVertexForce],"%g \n",maxPush);
  fprintf(state->logfiles[kCorrectionStepsNeeded],"%d \n",state->last_cstep_attempts);
  fprintf(state->logfiles[kEQVariance],"%d \n",state->eqVariance);

  if (state->steps%LOG_FLUSH_INTERVAL == 0) {

    for(i=0;i<kTotalGraphTypes;i++) fflush(state->logfiles[i]);

  }

}

/********************************************************************/

/*    correction_step and bsearch_stepper components                */

/********************************************************************/

void
step( plCurve* inLink, double stepSize, plc_vector* dVdt, search_state* inState )
{
  int cItr, vItr, dVdtItr;
  for( cItr=0, dVdtItr=0; cItr<inLink->nc; cItr++ ) {

    for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++, dVdtItr++ ) {

      /* Since our eq strategy isn't perfect, and we have non-eq theory
       * we use it here. Scale force at point by the deviation from the 
       * average edge length of the averge of adjacent edges
       */
      double avg, scale=1.0;
      
      if( inState->curvature_step != 0 ) {

	avg = (inState->sideLengths[inState->compOffsets[cItr] + vItr] + 
	       inState->sideLengths[inState->compOffsets[cItr] + vItr+1])/2.0;
	scale = avg / inState->avgSideLength;
			
      }

      plc_M_vmadd(inLink->cp[cItr].vt[vItr],scale*stepSize,dVdt[dVdtItr]);

    }
  }
}

      
void correct_thickness(plCurve *link,search_state *inState) 

     /* Newton's method correction algorithm for thickness. */
     
     /* Attempts to fix things so that the shortest strut and 
	worst minrad vertex are within "greenZone" of obeying 
	the constraints. 

	May take most of runtime if the number of verts is large,
	since we use the (slow) Stanford LSQR code to compute
	direction for Newton steps. */
     

     /* INCOMPLETE FOR NOW! */

{
  double greenZone = 
    2*inState->tube_radius*(1 - 0.5*inState->overstepTol);
  double mrgreenZone = 
    gLambda*inState->tube_radius*(1-0.5*inState->minradOverstepTol);
  
  if( (inState->curvature_step == 0 && inState->shortest < greenZone) ||
      (inState->curvature_step == 0 && inState->minrad < mrgreenZone) ) {
    
    // we haven't finished yet, record
    gCorrectionAttempts++;	    
    inState->curvature_step = 0;
    
      } else {
	
	if (inState->curvature_step == 0) { 
	  /* We are finishing a round of correction steps. 
	     Record the number of steps in the round in a log file. */
	  
	  fprintf(inState->logfiles[kCorrectionStepsNeeded],
		  "%d\n",gCorrectionAttempts);
	}
	
	inState->curvature_step = 1;
	
      }


  else {
      step(workerLink, correctionStepSize, dVdt, inState);
  }
  
}


int
equalStruts( const octrope_strut* s1, const octrope_strut* s2 )
{
  if( s1->component[0] == s2->component[0] &&
      s1->component[1] == s2->component[1] &&
      s1->lead_vert[0] == s2->lead_vert[0] &&
      s1->lead_vert[1] == s2->lead_vert[1] &&
      fabs(s1->position[0]-s2->position[0]) < 0.05 &&
      fabs(s1->position[1]-s2->position[1]) < 0.05 )
    {
      return 1;
    }
  return 0;
}


static void
normalizeVects( plc_vector* dvdt, int size )

     /* Scale the 3*size vector dvdt to norm 10 in R^{3*size}. */

{
  double sum = 0;
  int i;

  for( i=0; i<size; i++ ) { sum += plc_M_dot(dvdt[i],dvdt[i]); } 
  sum = sqrt(sum);
  for( i=0; i<size; i++ ) { plc_M_scale_vect(10/sum,dvdt[i]) }

}

plCurve*
bsearch_step( plCurve* inLink, search_state* inState )
{	
  int stepAttempts = 0;
  int dvdtItr;
  double  lastDCSD, lastMR, eps = inState->stepSize*inState->stepSize;
  plCurve* workerLink;
	
  double ERROR_BOUND = 1e-5;
  double MR_ERROR_BOUND = 1e-5;
  
  // create initial vector field for which we want to move along in this step
  plc_vector* dVdt; 

  dVdt = calloc(inState->totalVerts, sizeof(plc_vector));
  dlenForce(dVdt,inLink,inState);
  eqForce(dVdt,inLink,inState);

  resolveForce(dVdt,inLink,inState); /* Built from the bones of firstVariation. */

  /* Now we loop to find largest stepsize which doesn't violate constraints. */

  if( inState->shortest > 2*inState->tube_radius )
    inState->shortest = 2*inState->tube_radius;
			
  lastDCSD = inState->shortest;
  lastMR = inState->minrad;
  
  stepAttempts = 0;
  workerLink = NULL;
  double curr_error = 0, mr_error=0, newpoca = 0, newmr=0, 
    correctionStepSize=inState->correctionStepDefault;
  short	improvedNorm = 0;
	
  normalizeVects(dVdt, inState->totalVerts);
    
  double oldLength = inState->length;

  /* We're going to have to call octrope every time we go through this 
     loop in order to compute the level of error we have so far. In order
     to speed these calls up (slightly), we allocate memory once, in advance. */

  void *octmem;
  int  octmem_size;

  octmem_size = octrope_est_mem(plc_num_edges(inLink));
  octmem = malloc(sizeof(char)*octmem_size);
  fatalifnull_(octmem);
    
  double newthi,newrop,newlen;

  do {			
    
    stepAttempts++;
		
    // move along dVdt
	
    if( workerLink != NULL ) plc_free(workerLink);
    workerLink = plc_copy(inLink); 

    step(workerLink, inState->stepSize, dVdt, inState);

    /* We now compute the error in minRad and poca.  Since we need to
       update inState later if we're accepting this configuration, we
       just go ahead and make a full-on octrope call.
    
       This is pretty much a wash, speed-wise. I don't know whether
       octrope computes length in order to do an octrope_poca call, so
       there might be some speed savings in calling octrope_poca and
       octrope_minrad separately.

       On the other hand, having done so you need to either assemble
       these into thickness yourself after the loop (likely to be
       buggy) or call octrope again in order to do the assembly and
       length calculation (slow). 

       I'm probably overthinking this. */
  
    octrope(inLink,&newrop,&newthi,&newlen,&newmr,&newpoca,
	    0,0,NULL,0,NULL,0,0,NULL,0,NULL,
	    octmem,octmem_size);
    
    curr_error = (newpoca < 2*inState->tube_radius) ? max(lastDCSD-newpoca,0) : 0;
    mr_error = (newmr < gLambda*tube->radius) ? max(lastMR-newmr,0) : 0;

    /**********************************************************************************/
	
    /* This is where step size is being adjusted! */
    
    /* This algorithm may or may not be a good idea. It's an
       "optimistic" stepper, which concludes after each successful
       step that it should attempt to double step size.  If step size
       is basically constant, this results in one failure per step.
	   
       Further, it doesn't really try to converge on a nice stepsize,
       but rather just maintains a strict "powers of 2" hierarchy. It
       might be interesting to look this up and see about other
       possible designs here. */
	
    if( curr_error < ERROR_BOUND && mr_error < MR_ERROR_BOUND )
      inState->stepSize *= 2;
    else
      inState->stepSize /= 2; 

    /***************************************************/	
  
  } while( ((curr_error > ERROR_BOUND ) || (mr_error > MR_ERROR_BOUND)) &&
	   inState->stepSize < inState->maxStepSize 
	   && inState->stepSize > kMinStepSize );

  /* Notice that this loop condition restricts the maximum stepsize,
     even if the error is small, and causes us to accept a step which
     may cause a lot of error if the stepsize gets too small (we're
     hoping to fix things on the next correction step). */

  plc_free(inLink);		/* Commit to the step by changing inLink */
  inLink = workerLink;

  /* Now we update "inState" to reflect the changes in inLink. */
  
  inState->last_step_attempts = stepAttempts++;  
  
  inState->minrad = newmr;
  inState->shortest = newpoca;
  inState->length = newlen;
  inState->ropelength = newrop;
  inState->thickness = newthi;
  
  if( inState->stepSize < kMinStepSize )
    inState->stepSize = kMinStepSize;
  
  inState->cstep_time += inState->stepSize;
  inState->time += inState->stepSize;
   
  if( inState->stepSize > inState->maxStepSize ) inState->stepSize = inState->maxStepSize;
  
  // grab average dvdt now that we have finished mungering it
  inState->avgDvdtMag=0;
  for( dvdtItr=0; dvdtItr<inState->totalVerts; dvdtItr++ ) {

    inState->avgDvdtMag += plc_M_norm(dVdt[dvdtItr]);
  
  }
  inState->avgDvdtMag /= inState->totalVerts;
  
  // we should make sure stepsize isn't > avgDvdt^2 as that's the minrad control bound
  if( inState->stepSize > inState->avgDvdtMag*inState->avgDvdtMag && 
      inState->curvature_step != 0 &&
      inState->avgDvdtMag*inState->avgDvdtMag > kMinStepSize) 
    /* last this keeps us from zeroing in dVdt is really small */
    {
      
    inState->stepSize = inState->avgDvdtMag*inState->avgDvdtMag;
    
    }

  // also shouldn't be > 10% of edgelength
  if( inState->stepSize > inState->length/inState->totalVerts*.1 )
    inState->stepSize = inState->length/inState->totalVerts*.1;
  
  free(dVdt);
  free(octmem);
  
  return inLink;
}

/******************************************************************/
/*                 Building the rigidity matrix                   */
/******************************************************************/

static void
placeMinradStruts2Sparse( taucs_ccs_matrix* rigidityA, plCurve* inLink, 
			  octrope_mrloc* minradStruts, 
			  int minradLocs, search_state* inState, 
			  int contactStruts )

     /* Using the information provided by octrope in the octrope_mrloc array,
	and assuming that rigidityA already contains "contactstruts" 12 entry
	columns, loop through the mrloc locations, computing the gradient of 
	minrad at each vertex and adding the appropriate 9-entry column to the 
        matrix rigidityA at each step. 

	We assert that rigidityA contains at least contactstruts + minradlocs
	columns to begin with, and that rigidityA->colptr is set correctly for
        contactStruts 12-row columns and minradLocs 9-row columns. This _does_
        happen in firstVariation. */

{
  int mItr;
  int totalStruts = contactStruts + minradLocs;
  char errmsg[1024];
  bool ok,sortFlag = {FALSE};
  
  for( mItr=0; mItr<minradLocs; mItr++ ) {

    /* For each vertex at minimum minrad radius... compute the gradient of mr */

    plc_vector B, A, cross, As, Bs, Cs;
    double	bmag, amag;
    double value, angle;
    double kappa, dot, prevLen, thisLen;
    plc_vector  prevSide, thisSide, N, fancyL, fancyM, fancyN;
    
    int vItr = minradStruts[mItr].vert;
    int cItr = minradStruts[mItr].component;
		
    prevSide = plc_vect_diff(inLink->cp[cItr].vt[vItr],inLink->cp[cItr].vt[vItr-1]);
    thisSide = plc_vect_diff(inLink->cp[cItr].vt[vItr+1],inLink->cp[cItr].vt[vItr]);
		
    /* dot = plc_M_dot(prevSide, thisSide); */

    prevLen = plc_M_norm(prevSide);
    thisLen = plc_M_norm(thisSide);

    // B = b-v = thisSide. 

    B = thisSide;
    A = plc_scale_vect(-1,prevSide);

    bmag = plc_M_norm(B);
    amag = plc_M_norm(A);
		
    /* value = dot/(prevLen*thisLen); */
    angle = plc_angle(prevSide,thisSide);

    /* We now check that angle is high enough for the following 
       stuff to work. */

    if (angle < 1e-12 || angle > PI - 1e-12) {
      
      sprintf(errmsg,
	      "ridgerunner: Can't compute minrad gradient when edges are\n"
	      "             almost colinear. Angle between edges is %g.\n",
	      angle);

      FatalError(errmsg, __FILE__ , __LINE__ );

    }
    		
    if( thisLen < prevLen )
      {
	// says... maple?
	kappa = -bmag/(2-2*cos(angle));
	
	// BxA

	plc_M_cross(N,B,A);
	N = plc_normalize_vect(N,&ok);
	assert(ok);
	
	double Lconst, Mconst, Nconst;
	
	Lconst = (1/(2*tan(angle/2) * bmag));	
        fancyL = plc_scale_vect(Lconst,B);

	Mconst = kappa*(1/(amag*amag));
	// A x N

	plc_M_cross(cross,A,N);
	fancyM = plc_scale_vect(Mconst,cross);
	
	Nconst = kappa*(1/(bmag*bmag));
	// N x B

	plc_M_cross(cross,N,B);
	fancyN = plc_scale_vect(Nconst,cross);
	
	As = fancyM;
		
	Bs.c[0] = -fancyM.c[0] - fancyN.c[0] - fancyL.c[0];
	Bs.c[1] = -fancyM.c[1] - fancyN.c[1] - fancyL.c[1];
	Bs.c[2] = -fancyM.c[2] - fancyN.c[2] - fancyL.c[2];

	/* This is really the fastest way to accomplish this operation. The plCurve */
	/* option would be the code: */

	/* Bs = plc_scale_vect(-1,plc_vect_sum(plc_vect_sum(fancyM,fancyN),fancyL)); */

	/* which is really pretty awful. */
	
	Cs = plc_vect_sum(fancyN,fancyL);
	
      }
    else
      {
	// says... maple?
	kappa = -amag/(2-2*cos(angle));
	
	// BxA

	plc_M_cross(N,B,A);
	N = plc_normalize_vect(N,&ok);
	assert(ok);
	
	double Lconst, Mconst, Nconst;
	
	Lconst = (1/(2*tan(angle/2) * amag));

	fancyL = plc_scale_vect(Lconst,A);
	Mconst = kappa*(1/(amag*amag));

	// A x N

	plc_M_cross(cross,A,N);
	fancyM = plc_scale_vect(Mconst,cross);
	
	Nconst = kappa*(1/(bmag*bmag));
	// N x B
	plc_M_cross(cross,N,B);
	fancyN = plc_scale_vect(Nconst,cross);

	As = plc_vect_sum(fancyM,fancyL);
	
	Bs.c[0] = -fancyM.c[0] - fancyN.c[0] - fancyL.c[0];
	Bs.c[1] = -fancyM.c[1] - fancyN.c[1] - fancyL.c[1];
	Bs.c[2] = -fancyM.c[2] - fancyN.c[2] - fancyL.c[2];

	/* See above. This is better than using plCurve. */
	
	Cs = fancyN;

      }
    
    norm = sqrt(plc_M_dot(As,As) + plc_M_dot(Bs,Bs) + plc_M_dot(Cs,Cs));
    
    // if we're doing fast tsnnls based correction steps, the normalization 
    // doesn't really matter
    if( inState->curvature_step != 0 || gFastCorrectionSteps != 0 ) {

      plc_M_scale_vect(1/norm,As);
      plc_M_scale_vect(1/norm,Bs);
      plc_M_scale_vect(1/norm,Cs);
    
    }

    /* We have now computed the mr gradient, and it's time to fill in the 
       appropriate locations in rigidityA. */
    
    // temporarily increment the strut's verts based on their component interactions
    // we undo this change at the end of the for loop in case the user
    // wants to keep the strut set
    
    int aVert, bVert, cVert, entry;
    aVert = minradStruts[mItr].vert-1;
    bVert = minradStruts[mItr].vert;
    cVert = minradStruts[mItr].vert+1;
    
    aVert += inState->compOffsets[minradStruts[mItr].component];
    bVert += inState->compOffsets[minradStruts[mItr].component];
    cVert += inState->compOffsets[minradStruts[mItr].component];

    /* We start by placing vector As. "entry" will record the */
    /* row location of As.c[0] in rigidity matrix. Since there */
    /* is no wraparound addressing in the rigidity matrix, we */
    /* must deal with the situation if aVert = -1. */

    if ( aVert == -1 ) {

      if (!inLink->cp[minradStruts[mItr].component].open) {

	aVert = inLink->cp[minradStruts[mItr].component].nv-1;
	sortFlag = TRUE; 

	/* This screws up assumption that rowind values are increasing, and 
	   we'll have to fix that later on in the construction process. */

      } else { 

	/* We should not have an mrloc at vert 0 or nv-1 of an open component. */

	sprintf(errmsg,
		"ridgerunner: Should not have mrloc at vert 0 of an open\n"
		"             polyline.\n");
	FatalError(errmsg, __FILE__ , __LINE__ );

      }

    }

    /* This is the code that this fragment replaced: */

    /*if( minradStruts[mItr].vert == 0 &&
	(inLink->cp[minradStruts[mItr].component].acyclic == 0) ) {
      
	entry = 3*inLink->cp[minradStruts[mItr].component].nv-1 + 
	inState->compOffsets[minradStruts[mItr].component];
	
	  (I think this is a bug-- entry is used below as a row number,
	   so we shouldn't involve compOffsets.)
    
	   } else {
	   
	   entry = 3*aVert;
	   
	   }
    */

    assert(0 <= entry < rigidityA->m-2);
    assert(contactStruts+mItr < rigidityA->n);

    int col;

    colptr = rigidityA->colptr[contactsStruts+mItr];  
    /* Find rowind and values offset for column contactStruts + mItr]; */
		
    rigidityA->values.d[colptr+0] = As.c[0];
    rigidityA->values.d[colptr+1] = As.c[1];
    rigidityA->values.d[colptr+2] = As.c[2];
    
    rigidityA->rowind[colptr+0] = entry + 0;
    rigidityA->rowind[colptr+1] = entry + 1;
    rigidityA->rowind[colptr+2] = entry + 2;

    /* We now fill in the appropriate entries for b. We know that
       bVert should be a "normal" vertex, so there's no need to worry
       about open/closed plCurves. */
    
    entry = (3*bVert);
    assert(0 <= entry < rigidityA->m-2);

    rigidityA->values.d[colptr+3] = Bs.c[0];
    rigidityA->values.d[colptr+4] = Bs.c[1];
    rigidityA->values.d[colptr+5] = Bs.c[2];	
    
    rigidityA->rowind[colptr+3] = entry + 0;
    rigidityA->rowind[colptr+4] = entry + 1;
    rigidityA->rowind[colptr+5] = entry + 2;

    /* We now fill in the cVert values. Again, cVert might be too large
       if we have a minradloc on the n-1 vertex of a closed plCurve. So
       we fix cVert if needed. */

    if ( cVert == inLink->cp[minradStruts[mItr].component].nv ) {

      if (!inLink->cp[minradStruts[mItr].component].open) {

	cVert = 0;
	sortFlag = TRUE; 

	/* Again, this screws up the assumption that the rowind values 
	   will be increasing. We'll sort later to fix that. */

      } else { 

	/* We should not have an mrloc at vert 0 or nv-1 of an open component. */

	sprintf(errmsg,
		"ridgerunner: Should not have mrloc at vert nv-1 of an open\n"
		"             polyline.\n");
	FatalError(errmsg, __FILE__ , __LINE__ );

      }

    }

    entry = 3*cVert;
    assert(0 <= entry < rigidityA->m-2);
        
    rigidityA->values.d[colptr+6] = Cs.c[0];
    rigidityA->values.d[colptr+7] = Cs.c[1];
    rigidityA->values.d[colptr+8] = Cs.c[2];		
    
    rigidityA->rowind[colptr+6] = entry + 0;
    rigidityA->rowind[colptr+7] = entry + 1;
    rigidityA->rowind[colptr+8] = entry + 2;

  } // for over minrad struts

  /* We may have a few bad columns where we filled in rowind in the wrong order. */
  /* If so, sort to fix it. */

  if (sortFlag) { 

    taucs_enforce_ccs_sort(rigidityA, contactStruts);

  }
  
}

static void
placeContactStrutsSparse( taucs_ccs_matrix* A, plCurve* inLink, 
			  octrope_strut* strutSet, int strutCount, 
			  search_state* inState, int minradStruts )

     /* Part of the buildrigiditymatrix call. We assume when we're
        building the matrix that the strut columns are placed FIRST,
        and the minrad columns are placed after this. */
     
{
  int sItr, totalStruts;
  
  if( strutCount == 0 )
    return;
  
  // in constructing the ridigity matrix, we will need the struts as viewed 
  // as force vectors on the edges, so we create normalized vectors for each strut
  // here
  plc_vector* strutDirections = (plc_vector*)calloc(strutCount, sizeof(plc_vector));
  fatalifnull_(strutDirections);

  normalizeStruts( strutDirections, strutSet, inLink, strutCount );
  
  totalStruts = minradStruts + strutCount;
  
  for( sItr=0; sItr<strutCount; sItr++ ) {

    int		entry;
		
    // temporarily increment the strut's verts based on their component interactions
    // we undo this change at the end of the for loop in case the user
    // wants to keep the strut set
    strutSet[sItr].lead_vert[0] += inState->compOffsets[strutSet[sItr].component[0]];
    strutSet[sItr].lead_vert[1] += inState->compOffsets[strutSet[sItr].component[1]];
		
    // entry is the offset in A which begin this strut's influce
    // it corresponds to the x influence on the lead_vert[0]th vertex
    // after this line, entry+1 is y, +2: z.
    entry = 3*strutSet[sItr].lead_vert[0];
	
    // the strut information includes the position from the strut.lead_vert
    // so we assign "1-position[0]" of the force to the lead vert and "position[0]"
    // of the force to lead_vert+1
  
    // our column is 12*sItr from the start, and we need to set our rowInds
    A->values.d[12*sItr+0] = (1-strutSet[sItr].position[0]) * strutDirections[sItr].c[0];
    A->values.d[12*sItr+1] = (1-strutSet[sItr].position[0]) * strutDirections[sItr].c[1];
    A->values.d[12*sItr+2] = (1-strutSet[sItr].position[0]) * strutDirections[sItr].c[2];
			
    A->rowind[12*sItr+0] = entry + 0;
    A->rowind[12*sItr+1] = A->rowind[12*sItr+0] + 1;
    A->rowind[12*sItr+2] = A->rowind[12*sItr+0] + 2;
			
    // now for the next vertex, receiving "position[0]" of the force, this is 
    // potential wrapping case
    if( (strutSet[sItr].lead_vert[0]-inState->compOffsets[strutSet[sItr].component[0]]) 
	== (inLink->cp[strutSet[sItr].component[0]].nv-1) &&
	(inLink->cp[strutSet[sItr].component[0]].acyclic == 0) )
      {
	entry = 3*inState->compOffsets[strutSet[sItr].component[0]];
      }
    else
      {
	entry = 3*(strutSet[sItr].lead_vert[0]+1);
      }
  
    A->values.d[12*sItr+3] = (strutSet[sItr].position[0]) * strutDirections[sItr].c[0];
    A->values.d[12*sItr+4] = (strutSet[sItr].position[0]) * strutDirections[sItr].c[1];
    A->values.d[12*sItr+5] = (strutSet[sItr].position[0]) * strutDirections[sItr].c[2];
    
    A->rowind[12*sItr+3] = entry + 0;
    A->rowind[12*sItr+4] = A->rowind[12*sItr+3] + 1;
    A->rowind[12*sItr+5] = A->rowind[12*sItr+3] + 2;
		
    // we do the same thing at the opposite end of the strut, except now the 
    // force is negated
    entry = 3*strutSet[sItr].lead_vert[1];
					
    A->values.d[12*sItr+6] = (1-strutSet[sItr].position[1]) * -strutDirections[sItr].c[0];
    A->values.d[12*sItr+7] = (1-strutSet[sItr].position[1]) * -strutDirections[sItr].c[1];
    A->values.d[12*sItr+8] = (1-strutSet[sItr].position[1]) * -strutDirections[sItr].c[2];
    
    A->rowind[12*sItr+6] = entry + 0;
    A->rowind[12*sItr+7] = A->rowind[12*sItr+6] + 1;
    A->rowind[12*sItr+8] = A->rowind[12*sItr+6] + 2;
    
    if( (strutSet[sItr].lead_vert[1]-inState->compOffsets[strutSet[sItr].component[1]]) 
	== (inLink->cp[strutSet[sItr].component[1]].nv-1) &&
	(inLink->cp[strutSet[sItr].component[1]].acyclic == 0) ) {

      entry = 3*inState->compOffsets[strutSet[sItr].component[1]];
    
    } else {

      entry = 3*(strutSet[sItr].lead_vert[1]+1);
      
    }
    /*
      A[entry] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[0];
      A[entry+totalStruts] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[1];
      A[entry+(2*totalStruts)] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[2];
    */
    A->values.d[12*sItr+9] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[0];
    A->values.d[12*sItr+10] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[1];
    A->values.d[12*sItr+11] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[2];
    A->rowind[12*sItr+9] = entry + 0;
    A->rowind[12*sItr+10] = A->rowind[12*sItr+9] + 1;
    A->rowind[12*sItr+11] = A->rowind[12*sItr+9] + 2;
    
    strutSet[sItr].lead_vert[0] -= inState->compOffsets[strutSet[sItr].component[0]];
    strutSet[sItr].lead_vert[1] -= inState->compOffsets[strutSet[sItr].component[1]];
  }
	
  free(strutDirections);
}


void taucs_ccs_matrix *buildRigidityMatrix(plCurve *inLink,search_state *inState)

     /* Procedure calls octrope and uses the results to allocate and
	build a rigidity matrix for the curve, calling
	placeMinradStruts, placeContactStruts, and
	placeConstraintStruts. */

{
  int strutStorageSize = 0;
  int minradStorageSize = 0;

  int strutCount = 0;
  int minradLocs = 0;
  int constraintCount = 0;

  octrope_strut *strutSet = NULL;
  octrope_mrloc *minradSet = NULL;

  char dumpname[1024],errmsg[1024];

  /* We need to start by allocating some memory for the strutSet. It
     is clear that we can kill the problem of finding too many minrad
     struts by setting minradStorageSize to the number of verts. We
     believe (with no proof) that choosing strutStorageSize = 6*the
     number of edges in the link we will similarly overestimate the #
     of struts. We will double-check the results anyway. */

  /* We note that this process could be marginally more efficient if
     we did this malloc "once and for all" at the start of
     computation. */

  strutStorageSize = 6*plc_verts(inLink);
  minradStorageSize = plc_verts(inLink);

  strutSet = malloc(sizeof(octrope_strut)*strutStorageSize);
  minradSet = malloc(sizeof(octrope_mrloc)*minradStorageSize);

  fatalifnull_(strutSet);
  fatalifnull_(minradSet); 

  octrope(inLink, 
		
	  &inState->ropelength,
	  &inState->thickness,		
		
	  &inState->length,		
	  &inState->minrad,
	  &inState->shortest,
		
	  // minrad struts
	  gLambda*inState->tube_radius,  /* Cutoff reports all mrlocs less than this */
	  0,                             /* This value will be ignored. */
	  minradSet, 
	  minradStorageSize,
	  &minradLocs,
		
	  // strut info
	  2*inState->tube_radius, /* Cutoff reports struts shorter than this. */
	  0,                      /* This epsilon value also ignored. */
	  strutSet,
	  strutStorageSize,
	  &strutCount,
		
	  NULL, 0,

	  gLambda);               /* The global "stiffness" parameter. */  

  /* We now need to make sure that we didn't exceed the size of the strut
     and/or minrad buffers. */

  if (strutCount >= strutStorageSize || minradLocs >= minradStorageSize ||
      strutCount < 0 || minradLocs < 0) {

    dumpLink(inLink,inState,dumpname);
    sprintf(errmsg,
	    "ridgerunner: octrope found %d struts and %d mrlocs on\n"
	    "             the %d vertex link (dumped to %s), too close to\n"
	    "             minradStorageSize of %d or strutStorageSize %d.\n",
	    strutCount,minradLocs,plc_verts(inLink),dumpname,
	    minradStorageSize,strutStorageSize);
    FatalError(errmsg, __FILE__ , __LINE__ );

  }

  inState->lastStepStrutCount = strutCount;
  inState->lastStepMinradStrutCount = minradLocs;

  constraintCount = plCurve_score_constraints(inLink);

  /* We now build the rigidity matrix. */
  
  /* A maps from the strut space to the space of variations of
   * vertices, and gives the force at each vertex resulting from a
   * compressive force pushing _out_ from each strut.
   *
   * If the strut strikes in the middle of an edge, we apply its force
   * to both endpoints of the edge, divided according to the position
   * of the end of the strut along the edge.
   *
   * Each row of A corresponds to a single component of a single
   * _vertex_ of the overall picture.  The entries in the row
   * corresponding to each strut that pushes on the vertex are that
   * component of the unit vector pointing _out_ from that strut's
   * endpoint at the edge incident to the given vertex.
   */
				 	
  int nnz = 12*strutCount + 9*minradLocs + 3*constraintCount; // we KNOW this
  taucs_ccs_matrix *cleanA;
  
  cleanA = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  fatalifnull_(cleanA);

  cleanA->n = strutCount+minradLocs+constraintCount;
  cleanA->m = 3*inState->totalVerts;
  cleanA->flags = TAUCS_DOUBLE;
  
  cleanA->colptr = (int*)malloc(sizeof(int)*(cleanA->n+1));
  cleanA->rowind = (int*)malloc(sizeof(int)*nnz);
  cleanA->values.d = (double*)malloc(sizeof(taucs_double)*nnz);

  fatalifnull_(cleanA->colptr);
  fatalifnull_(cleanA->rowind);
  fatalifnull_(cleanA->values.d);
  
  // just as we know nnz, we know the column positions. 
  // struts involve 4 vertices, 
  // minrad struts involve 3 vertices,
  // constraint struts involve 1 vertex.

  cleanA->colptr[0] = 0;

  for( sItr=1; sItr<=strutCount; sItr++ ) /* The end condition isn't clear here-- */
    cleanA->colptr[sItr] = sItr*12;       /* do we ignore column 0? does tsnnls? */
  
  for(sItr=strutCount+1; sItr<=strutCount+minradLocs; sItr++ )
    cleanA->colptr[sItr] = cleanA->colptr[sItr-1]+9;

  for(sItr=strutCount+minradLocs+1; 
      sItr<=strutCount+minradLocs+constraintCount;sItr++) {

    cleanA->colptr[sItr] = cleanA->colptr[sItr-1]+3;

  }
  

  
  
}

/************************************************************************/





void
updateSideLengths( plCurve* inLink, search_state* inState )
{
  int cItr, vItr;
  
  inState->totalSides = 0;
  for( cItr=0; cItr<inLink->nc; cItr++ ) {

    inState->totalSides += inLink->cp[cItr].nv;
  }
  
  if( inState->sideLengths != NULL )
    free(inState->sideLengths);
  
  inState->sideLengths = (double*)malloc(inState->totalSides*sizeof(double));
  fatalifnull_(inSide->sideLengths);
  
  int tot = 0;
  for( cItr=0; cItr<inLink->nc; cItr++ ) {

    for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ ) {
         
      inState->sideLengths[inState->compOffsets[cItr] + vItr] = 
	plc_M_distance(inLink->cp[cItr].vt[vItr],
		       inLink->cp[cItr].vt[vItr+1]) /*;*/
	
      inState->avgSideLength += plc_M_norm(side);
      tot++;
    }
  }
  inState->avgSideLength /= (double)tot;
}

static void
spinForce( plc_vector* dlen, plCurve* inLink, search_state* inState )

     /* This adds to dlen a tangential force causing the tube to spin
	as it is tightening. This doesn't seem to have any useful
	numerical effect. */

{
  int cItr, vItr;
  int *edges = malloc(sizeof(int)*inLink->nc);

  // note to self - this can be made much faster by not being dumb
  
  plc_fix_wrap(inLink);
  plc_edges(inLink,edges);
  
  for( cItr=0; cItr<inLink->nc; cItr++ ) {

    if (!inLink->cp[cItr].acyclic) {  
      
      /* We can only spin _closed_ components. */

      // the first thing to do is grab edge lengths
      plc_vector* sides;
      plc_vector* adjustments;
      
      sides = (plc_vector*)malloc(sizeof(plc_vector)*edges[cItr]);
      adjustments = (plc_vector*)calloc(edges[cItr], sizeof(plc_vector));
      
      for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ ) {
	
	plc_vector s1, s2;
	bool ok;
	char errmsg[1024],dumpname[1024];

	sides[vItr] = plc_normalize_vect(
					 plc_vect_diff(inLink->cp[cItr].vt[vItr+1],
						       inLink->cp[cItr].vt[vItr]),
					 &ok);

	if (!ok) { 

	  dumpLink(inLink,inState,dumpname);
	  sprintf("ridgerunner: spinforce can't normalize side %d of component %d\n"
		  "             of link %s.\n",vItr,cItr,dumpname);
	  FatalError(errmsg, __FILE__ , __LINE__ );

	}
	
      }

      double spinFactor = 0.5;
      // fix zero and compute the tangential change necessary for the rest of the vertices
      for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ ) {
	
	adjustments[vItr] = plc_scale_vect(spinFactor,sides[vItr]);

	dlen[((vItr)) + inState->compOffsets[cItr]].c[0] += adjustments[(vItr)].c[0];
	dlen[((vItr)) + inState->compOffsets[cItr]].c[1] += adjustments[(vItr)].c[1];
	dlen[((vItr)) + inState->compOffsets[cItr]].c[2] += adjustments[(vItr)].c[2];
      }
      
      free(sides);
      free(adjustments);
    
    }

  }
  
  free(edges);

  char fname[1024];
  sprintf("adjustedDL.vect");
  exportVect( dlen, inLink, fname );

}

void
specialForce( plc_vector* dlen, plCurve* inLink, search_state* inState )

     /* This stub is left here in case we want to experiment with knot
	tightening with respect to some external body force (like
	gravity or something). */

{ 

}

static void
eqForce( plc_vector* dlen, plCurve* inLink, search_state* inState )
{
  int cItr, vItr;
  int *edges = malloc(sizeof(int)*inLink->nc);
  fatalifnull_(edges);
	
  // note to self - this can be made much faster by not being dumb
	
  plc_fix_wrap(inLink);
		
  double* diffFromAvg;
  double* averages;

  diffFromAvg = (double*)malloc(sizeof(double)*inState->totalSides);
  fatalifnull_(diffFromAvg);
  averages = (double*)malloc(sizeof(double)*inLink->nc);
  fatalifnull_(averages);

  plc_edges(inLink,edges);
  
  for( cItr=0; cItr<inLink->nc; cItr++ ) {

    // the first thing to do is grab edge lengths
 
    double* lengths;
    plc_vector* adjustments;
    
    double  averageLength;
    double  usedLength;
    double  scaleFactor;
    
    lengths = (double*)malloc(sizeof(double)*edges[cItr]);
    fatalifnull_(lengths);
    adjustments = (plc_vector*)calloc(edges[cItr], sizeof(plc_vector));
    fatalifnull_(adjustments);

    averageLength = 0;
    
    for( vItr=0; vItr<edges[cItr]; vItr++ ) {

      lengths[vItr] = plc_M_distance(inLink->cp[cItr].vt[vItr+1],
				     inLink->cp[cItr].vt[vItr]) /*;*/
            	
      averageLength += lengths[vItr];
    
    }
    
    averageLength /= edges[cItr];
    averages[cItr] = averageLength;
     
    for( vItr=0; vItr<edges[cItr]; vItr++ ) {
     
      diffFromAvg[vItr+inState->compOffsets[cItr]] = 
	fabs(lengths[vItr] - averageLength)/averageLength;	
      
    }
    
    usedLength = 0;
    
    // fix zero and compute the tangential change necessary for the rest of the vertices
    for( vItr=1; vItr<edges[cItr]-1; vItr++ ) {

      usedLength += lengths[vItr];
      scaleFactor = (vItr)*averageLength - usedLength;
      // scale the side by the amount we need to move
      scaleFactor *= inState->eqMultiplier;
      
      /* "we change the eq code so that it moves in the direction
	 perpendicular to the gradient of length (that is, the
	 average of the two tangent vectors) rather than in the
	 direction of the tangent vector to the "left" edge at each
	 vertex.  This might make a small difference in eq
	 performance (it certainly shouldn't hurt) and it makes the
	 theory nicer."
      */

      bool ok;
      
      adjustments[vItr] = plc_scale_vect(scaleFactor,
					 plc_mean_tangent(inLink,cItr,vItr,&ok));

      if (!ok) {

	char dumpname[1024], errmsg[1024];
	
	dumpLink(inLink,inState,dumpname);
	sprintf(errmsg,
		"ridgerunner: eqForce can't compute tangent at vertex %d of comp %d \n"
		"             of link %s.\n",vItr,cItr,dumpname);
	FatalError(errmsg, __FILE__ , __LINE__ );
	
      }

      plc_M_add_vect(dlen[vItr + inState->compOffsets[cItr]],
		     adjustments[vItr]) /*;*/

    }
    
    free(lengths);
    free(adjustments);

  }

  free(edges);

  /* Now we compute some statistics about how well we're EQing. */
  
  double sum = 0, isum, sampleAvg=0;
  int N = 0;
  
  // remember you made these all percent diffs 
  
  for( cItr=0; cItr<inLink->nc; cItr++ ) {
      for( vItr=0; vItr<inState->totalSides; vItr++ ) {

	sampleAvg += diffFromAvg[vItr+inState->compOffsets[cItr]];
	N++;
      }
  }
  sampleAvg /= (double)N;
  
  for( cItr=0; cItr<inLink->nc; cItr++ ) {
    for( vItr=0; vItr<inState->totalSides; vItr++ ) {
      
      isum = diffFromAvg[vItr+inState->compOffsets[cItr]] - sampleAvg;
      sum += isum*isum;

    }
  }
  
  inState->eqAvgDiff = sampleAvg;
  inState->eqVariance = (1.0/(double)N)*sum;
  
  free(diffFromAvg);
  free(averages);
}

static double*
stanford_lsqr( taucs_ccs_matrix* sparseA, double* minusDL, double* residual )

     /* This presents a problem. By default, we don't link a copy of
	lsqr directly into the RR build, preferring to get our linear
	algebra through the tsnnls package. On the other hand, this
	functionality is not directly exposed in the tsnnls build either. */

     /* The difference is just that we return the norm of the residual
	of the lsqr solution, wheras the t_lsqr function doesn't. It's
	probably going to be best to replace this with a t_lsqr call
	anyway and work around any problems introduced that way. */

     /* We will make a decision on this later, once we understand 
	the code below. */

{ /* if there are no constrained struts, t_snnls won't actually work,
     so use SOL LSQR */

  lsqr_input   *lsqr_in;
  lsqr_output  *lsqr_out;
  lsqr_work    *lsqr_work;
  lsqr_func    *lsqr_func;
  int bItr;
  double*      result;
  
  alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, sparseA->m, sparseA->n );
  
  /* we let lsqr() itself handle the 0 values in this structure */
  lsqr_in->num_rows = sparseA->m;
  lsqr_in->num_cols = sparseA->n;
  lsqr_in->damp_val = 0;
  lsqr_in->rel_mat_err = kZeroThreshold;
  lsqr_in->rel_rhs_err = 0;
  lsqr_in->cond_lim = 0;
  lsqr_in->max_iter = lsqr_in->num_rows + lsqr_in->num_cols + 5000;
  lsqr_in->lsqr_fp_out = NULL;	
 
  for( bItr=0; bItr<sparseA->m; bItr++ )
    {
      lsqr_in->rhs_vec->elements[bItr] = minusDL[bItr];
    }
  /* Here we set the initial solution vector guess, which is 
   * a simple 1-vector. You might want to adjust this value for fine-tuning
   * t_snnls() for your application
   */
  for( bItr=0; bItr<sparseA->n; bItr++ )
    {
      lsqr_in->sol_vec->elements[bItr] = 0; 
    }
  
  /* This is a function pointer to the matrix-vector multiplier */
  lsqr_func->mat_vec_prod = fast_lsqr_mult;
  
  lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, sparseA );
  
  result = (double*)malloc(sizeof(double)*sparseA->n);
  for( bItr=0; bItr<sparseA->n; bItr++ ) // not really bItr here, but hey
    result[bItr] = lsqr_out->sol_vec->elements[bItr];
  
  if( residual != NULL )
    *residual = lsqr_out->resid_norm;
  
  free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );
  
  return result;
}

int gFeasibleThreshold = 0;
int gDeferredStrutExport = 0;
int gFoo = 0;

void
resolveForce( plc_vector* dl, plCurve* inLink, search_state* inState)

     /* This function resolves the given force dl over the struts of
	inLink, building a rigidity matrix along the way and calling
	octrope and tsnnls to do their jobs. 
     */

{
  /*
   * this function creates dVdt, the vector field under which we should actually be 
   * moving at this point. We resolve dlen curvature force over the strut field as 
   * returned by liboctrope.
   */
  
  plc_vector*	  dVdt = NULL;
  int		  strutStorageSize = 0;
  int					strutCount = 0;
  octrope_strut*		strutSet = NULL;
  octrope_mrloc*		minradSet = NULL;
  double				thickness = 2*inState->tube_radius;
  int					cItr, sItr; // loop iterators over components, verticies, and struts
  int					minradLocs;
  int					dlItr;
  double				dummyThick;
  
  double*				A = NULL; // the rigidity matrix
	taucs_ccs_matrix*	cleanA = NULL;
	
	strutStorageSize = (inState->lastStepStrutCount != 0) ? (2*inState->lastStepStrutCount) : 10;
	
	if( outStrutsCount != NULL )
		*outStrutsCount = 0;
	if( outStruts != NULL )
		*outStruts = NULL;
	
	// grab dlen
	if( dl == NULL )
		dl = calloc(inState->totalVerts, sizeof(double));
	
	if( dlenStep != 0 )
		dlenForce(dl, inLink, inState);
	
//	dl = (plc_vector*)calloc(inState->totalVerts, sizeof(struct plc_vector_type));
		
	// acquire the strut set, we probably haven't increased the strut count much 
	// beyond inState.lastStrutCount, so it's probably safe to assume that we can fit
	// the strut set in 2*inState.lastStrutCount, but maybe not, so we loop over 
	// strut finder until we have stored less than we possibly can.
	strutSet = NULL;
	do
	{
	
		if( strutSet != NULL )
			free(strutSet);
		if( minradSet != NULL )
			free(minradSet);
			
		strutSet = (octrope_strut*)calloc(strutStorageSize, sizeof(octrope_strut));
		minradSet = (octrope_mrloc*)malloc(strutStorageSize * sizeof(octrope_mrloc));
			
		if( gSurfaceBuilding == 0 )
		{
			octrope(	inLink, 

						1,	// factor 1
						&inState->ropelength,
						&dummyThick,		
						
						&inState->length,
						
						&inState->minrad,
						&inState->shortest,

						// minrad struts
						0.5,
						100,
						minradSet, 
						strutStorageSize,
						&minradLocs,
						
						// strut info
						thickness,
						100,
						strutSet,
						strutStorageSize,
						&strutCount,
						
						NULL, 0 );
		}
		else
		{
			// use epsilon to get a more complete strut set
			octrope(	inLink, 

						1,	// factor 1
						&inState->ropelength,
						&dummyThick,		
						
						&inState->length,
						
						&inState->minrad,
						&inState->shortest,

						// minrad struts
						0.5,
						100,
						minradSet, 
						strutStorageSize,
						&minradLocs,
						
						// strut info
						0,
				  //	  0.001,
						2*inState->overstepTol,
						strutSet,
						strutStorageSize,
						&strutCount,						
						NULL, 0 );
		}
					
		if( inState->ignore_minrad )
			minradLocs = 0;
					
	//	minradLocs = 0; // these actually aren't helpful
		gFeasibleThreshold = minradLocs + strutCount;
		
		strutStorageSize += 10; // increase if we have to repeat this loop
	} while( (strutCount == (strutStorageSize-10)) || (minradLocs==(strutStorageSize-10)) );
	
	strutStorageSize -= 10; // return it to its actual value
			
	//gFeasibleThreshold = strutCount;
	
	//collapseStruts(&strutSet, &strutStorageSize);
	
	inState->lastStepStrutCount = strutCount;
	inState->lastStepMinradStrutCount = minradLocs;
	
	int barVerts = 0;
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		if( inState->conserveLength[cItr] != 0 )
		{
			barVerts += inLink->cp[cItr].nv;
			minradLocs += plc_strand_edges(&inLink->cp[cItr]);
		}
	}

	if( dlenStep != 0 )
		eqForce(dl, inLink, inState);

//	if( dlenStep != 0 || inState->eq_step != 0 )
//	if( inState->eq_step != 0 )
//	{
//		eqForce(dl, inLink, inState);
	//	spinForce(dl, inLink, inState);
//	}
//	specialForce(dl, inLink, inState);
	
	// this actually doesn't DO anything except print minrad state right now (if gOutputFlag)
	//minradForce(dl, inLink, inState);

	//specialForce(dl, inLink, inState);
	
	// if there are no struts, dVdt is completely controlled by dlen as there 
	// are no self contact constraints to worry about
	if( strutCount+minradLocs == 0 )
	{
		// we don't need this space if there are no struts
		free(strutSet);
		free(minradSet);
	
		// if there are no struts...
		//inState->shortest = 2*inState->tube_radius;

		//dVdt = dl;
		
		inState->avgDvdtMag = 0;
		for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
		{
			inState->avgDvdtMag += plc_M_norm(dl[dlItr]);
		}
		inState->avgDvdtMag /= inState->totalVerts;
	}
	else
	{
		double* compressions = NULL;
		double* minusDL;
		double* totalPushes = NULL;
		int		Acols = 0;
		
		if( outStrutsCount != NULL )
			*outStrutsCount = strutCount;
		if( outStruts != NULL )
			*outStruts = strutSet;
			
		/* if there are struts, we need to construct A, the rigidity matrix.
		 * A maps from the strut space to the space of variations
		 * of vertices, and gives the force at each vertex resulting from a compressive force pushing _out_
		 * from each strut. 
		 *
		 * If the strut strikes in the middle of an edge, we apply its force to both endpoints of the 
		 * edge, divided according to the position of the end of the strut along the edge.
		 *
		 * Each row of A corresponds to a single component of a single _vertex_ of the overall picture.
		 * The entries in the row corresponding to each strut that pushes on the vertex are that component 
	 	 * of the unit vector pointing _out_ from that strut's endpoint at the edge incident to the
	 	 * given vertex.
		 */
				 
	/*	for( sItr=0; sItr<strutCount; sItr++ )
		{
			printf( "strut %d norm: %lf %lf %lf\n", sItr, strutDirections[sItr].c[0], 
				strutDirections[sItr].c[1], strutDirections[sItr].c[2] );
		}
	*/	
		// A is size:   rows - 3*totalVerts
		//				cols - strutCount		
		
//		taucs_ccs_matrix* sparseA = NULL;
		int nnz = 12*strutCount + 9*minradLocs; // we KNOW this
		
		  cleanA = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
		  cleanA->n = strutCount+minradLocs;
		  cleanA->m = 3*inState->totalVerts;
		  cleanA->flags = TAUCS_DOUBLE;
		  
		  cleanA->colptr = (int*)malloc(sizeof(int)*(cleanA->n+1));
		  cleanA->rowind = (int*)malloc(sizeof(int)*nnz);
		  cleanA->values.d = (double*)malloc(sizeof(taucs_double)*nnz);
		  
		  // just as we know nnz, we know the column positions. struts only
		  // involve 4 vertices and minrad struts 3.
		cleanA->colptr[0] = 0;
		  for( sItr=1; sItr<=strutCount; sItr++ )
			cleanA->colptr[sItr] = sItr*12;
		
		for(sItr=strutCount+1; sItr<=strutCount+minradLocs; sItr++ )
			cleanA->colptr[sItr] = cleanA->colptr[sItr-1]+9;
						
		Acols = (3*inState->totalVerts)*(strutCount+minradLocs);
//		A = (double*)calloc(Acols, sizeof(double)); // calloc zeros A
		// first, place contact struts
//		placeContactStruts(A, inLink, strutSet, strutCount, inState, minradLocs);
		placeContactStrutsSparse(cleanA, inLink, strutSet, strutCount, inState, minradLocs);
		
//		placeVertexBars(A, inLink, strutCount, barVerts, minradLocs, inState);
		
		// and then minrad struts
		if( minradLocs > 0 )
		{
			placeMinradStruts2Sparse(cleanA, inLink, minradSet, minradLocs, inState, strutCount);
	//		placeMinradStruts2(A, inLink, minradSet, minradLocs, inState, strutCount);
		}
		
		/* We now try to cancel as much of the motion as possible with strut forces. 
		 * This involves finding the best nonnegative partial solution to the system 
		 *
		 *                               AX = -dl. 
		 *
		 * Of course, we can't solve this equation in general (unless the knot is critical!)
		 * so we settle for the closest partial solution in a least-squares sense. 
		 * A further description of this function is contained in taucs_snnls.c
		 */
	/*	if( minradLocs > 0 )
		{
			taucs_ccs_matrix* sparseA = taucs_construct_sorted_ccs_matrix(A, strutCount+minradLocs, 3*inState->totalVerts);
			taucs_print_ccs_matrix(sparseA);
			taucs_ccs_free(sparseA);
		}
	*/	
			
		
		// fill this in now if we're not correcting, otherwise we might 
		// have to flip strut gradients for those already in the green zone
	/*	if( dlenStep != 0 || inState->eq_step != 0 )
			sparseA = taucs_construct_sorted_ccs_matrix(A, strutCount+minradLocs, 3*inState->totalVerts);
	*/		
		int*	greenZoneStruts;
		int		greenZoneCount = 0;
			
		// construct minusDL -- this is a column vector of size totalVerts
		// we must operate using strictly doubles to interoperate with taucs_snnls
		minusDL = (double*)calloc((3*inState->totalVerts), sizeof(double));
		if( dlenStep != 0 || inState->eq_step != 0 )
		{
			for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
			{
				minusDL[dlItr*3+0] = -dl[dlItr].c[0];
				minusDL[dlItr*3+1] = -dl[dlItr].c[1];
				minusDL[dlItr*3+2] = -dl[dlItr].c[2];
			}
		}
		else
		{
			if( gFastCorrectionSteps == 0 )
			{
				taucs_ccs_matrix* sparseAT = NULL;
				double*	ofvB = (double*)calloc(strutCount+minradLocs, sizeof(double));
			
			/*	sparseA = taucs_construct_sorted_ccs_matrix(A, strutCount+minradLocs, 3*inState->totalVerts);
				sparseAT = taucs_ccs_transpose(sparseA);
				taucs_ccs_free(sparseA);
			*/
				greenZoneStruts = (int*)malloc(sizeof(int)*strutCount);
			
			//	if(  inState->shortest > thickness - (thickness*inState->overstepTol)*0.5 )
				{
					for( sItr=strutCount; sItr<strutCount+minradLocs; sItr++ )
					{
						double minradSecondGreenZone = ((thickness/2.0)-inState->minminrad) / 2.0;
						minradSecondGreenZone = (thickness/2.0)-minradSecondGreenZone;
						// we're in green green, and we should shorten rather than lengthen, but not so much that
						// we dip below the min again
						if( minradSet[sItr-strutCount].mr > minradSecondGreenZone )
						{
					//		ofvB[sItr] = minradSet[sItr-strutCount].mr - minradSecondGreenZone;
					//		greenZoneMR[greenZoneMRCount++] = sItr; // keep in mind this is offset by strutCount
						}
						else // we want to lengthen to the second green zone
							ofvB[sItr] = minradSecondGreenZone-minradSet[sItr-strutCount].mr;
					}
				}
							
			//	if( inState->shortest < thickness - (thickness*inState->overstepTol) )
				{
					for( sItr=0; sItr<strutCount; sItr++ )
					{
						double secondGreenZone = thickness - (thickness*inState->overstepTol)*.25;
						double greenZone = thickness - (thickness*inState->overstepTol)*0.5;
						
						// this only works if there aren't minrad struts for 
						// deep numerical reasons. talk to jason or me.
						if( strutSet[sItr].length > secondGreenZone && minradLocs == 0 )
						{
							// shorten and allow shrinkage
							ofvB[sItr] = strutSet[sItr].length-(secondGreenZone);
							// if we don't in some way protect existing struts, they might 
							// be destroyed during correction, so we'll record everyone
							// who is currently in the green zone and flip their 
							// gradients in the rigidity matrix for this step
							greenZoneStruts[greenZoneCount++] = sItr;
						}
						else if( strutSet[sItr].length < greenZone )
						{
							double target = thickness;
							ofvB[sItr] = secondGreenZone-strutSet[sItr].length;
						}
					}
				}
				
				int sIndex, mItr, spotItr;
				for( sIndex=0; sIndex<greenZoneCount; sIndex++ )
				{
					// flip the entries in this guy's column
					int totalStruts = strutCount+minradLocs;
					sItr = greenZoneStruts[sIndex];
					int entry = (totalStruts*3*strutSet[sItr].lead_vert[0])+sItr;
					
					// this is actually MUCH easier since we don't have to compute the rows
					
					for(spotItr=0; spotItr<12; spotItr++ )
						cleanA->values.d[12*sItr+spotItr] *= -1;
				
			/*		A[entry] *= -1;
					A[entry+totalStruts] *= -1;
					A[entry+(2*totalStruts)] *= -1;
					
					if( (strutSet[sItr].lead_vert[0]-inState->compOffsets[strutSet[sItr].component[0]]) == (inLink->cp[strutSet[sItr].component[0]].nv-1) &&
						(inLink->cp[strutSet[sItr].component[0]].acyclic == 0) )
					{
						entry = (totalStruts*3*inState->compOffsets[strutSet[sItr].component[0]]);
					}
					else
					{
						entry = (totalStruts*3*(strutSet[sItr].lead_vert[0]+1))+sItr;
					}
					A[entry] *= -1;
					A[entry+totalStruts] *= -1;
					A[entry+(2*totalStruts)] *= -1;
					
					// other end
					entry = (totalStruts*3*strutSet[sItr].lead_vert[1])+sItr;
					A[entry] *= -1;
					A[entry+totalStruts] *= -1;
					A[entry+(2*totalStruts)] *= -1;
					if( (strutSet[sItr].lead_vert[1]-inState->compOffsets[strutSet[sItr].component[1]]) == (inLink->cp[strutSet[sItr].component[1]].nv-1) &&
						(inLink->cp[strutSet[sItr].component[1]].acyclic == 0) )
					{
						entry = (totalStruts*3*inState->compOffsets[strutSet[sItr].component[1]]);
					}
					else
					{
						entry = (totalStruts*3*(strutSet[sItr].lead_vert[1]+1))+sItr;
					}
					A[entry] *= -1;
					A[entry+totalStruts] *= -1;
					A[entry+(2*totalStruts)] *= -1;*/
				}
							
				// done with these
				free(greenZoneStruts);
				
				// now we can build since we've flipped gradients in full A.
			//	sparseA = taucs_construct_sorted_ccs_matrix(A, strutCount+minradLocs, 3*inState->totalVerts);
			//	sparseAT = taucs_ccs_transpose(sparseA);
		
				sparseAT = taucs_ccs_transpose(cleanA);
		
				inState->ofvNorm = 0;
				for( sItr=0; sItr<strutCount+minradLocs; sItr++ )
				{
					inState->ofvNorm += ofvB[sItr]*ofvB[sItr];
				}
				inState->ofvNorm = sqrt(inState->ofvNorm);
							
				/*
				 * Fact: The rigidity matrix A (compressions) = (resulting motions of verts).
				 * So it's also true that
				 * 
				 *		  A^T (a motion of verts) = (resulting change in edge length)
				 *
				 *	Now we _have_ a desired change in edge lengths, namely (1 - l_i), (1-l_j)
				 *	and all that. Call it b. So we really ought to compute velocity for a
				 *	correction step by
				 *
				 *		  A^T v = b.
				 *
				 *	Is this equation always solvable? Probably not. So we'll take the results
				 *	of
				 *
				 *		 lsqr(A^T,b) = overstep fixing velocity (OFV)
				 *
				 */
				
				double* ofv;
				
				/* Jason note: we might have to change the lsqr initial guess later
				 * me: remember the funky motion thing?
				 *
				 * If we have to change it: 
				 *	want: linear combo of gradient of active constraints 
				 *			coeffs are k's computed just like k in test
				 *			(the change you want in that constraint divided
				 *			by square of norm o' gradient)
				 */
				ofv = stanford_lsqr(sparseAT, ofvB, &inState->ofvResidual);

			/* test
				double dx, k, norm;
				double diff[5000];
				dx = 0.5 - octrope_minradval(inLink);
				norm = 0;
				int nItr;
				for( nItr=0; nItr<3*inState->totalVerts; nItr++ )
					norm += A[nItr]*A[nItr];
				norm = sqrt(norm);
				k = dx / (norm*norm);
				for( nItr=0; nItr<3*inState->totalVerts; nItr++ )
				{
					diff[nItr] = ofv[nItr] - k*A[nItr];
			//		printf( "%e\n", diff[nItr] );
				}
				
				norm = 0;
				for( nItr=0; nItr<3*inState->totalVerts; nItr++ )
					norm += ofv[nItr]*ofv[nItr];
				norm = sqrt(norm);
				printf( "|ofv|: %e\n", norm );
		
				norm = 0;
				for( nItr=0; nItr<3*inState->totalVerts; nItr++ )
					norm += A[nItr]*A[nItr];
				norm = sqrt(norm);
				printf( "k*|A|: %e\n", k*norm );
				
				norm = 0;
				for( nItr=0; nItr<3*inState->totalVerts; nItr++ )
					norm += diff[nItr]*diff[nItr];
				norm = sqrt(norm);
				printf( "|diff|: %e\n", norm );
				end o' test		*/ 

				
				/*plc_vector* ofvVec;
				ofvVec = (plc_vector*)malloc(sizeof(plc_vector)*inState->totalVerts);
				int i;
				for( i=0; i<inState->totalVerts; i++ )
				{
					ofvVec[i].c[0] = ofv[3*i+0];
					ofvVec[i].c[1] = ofv[3*i+1];
					ofvVec[i].c[2] = ofv[3*i+2];
				}
				exportVect( ofvVec, inLink, "/tmp/ofv.vect" );
				free(ofvVec);
				*/
			/*	for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
				{
					plc_vector  vert;
					double  norm;
					vert.c[0] = ofv[3*dlItr+0];
					vert.c[1] = ofv[3*dlItr+1];
					vert.c[2] = ofv[3*dlItr+2];
					norm = plc_M_norm(vert);
					if( norm != 0 )
					{
						ofv[3*dlItr+0] /= norm;
						ofv[3*dlItr+1] /= norm;
						ofv[3*dlItr+2] /= norm;
					
					/*	ofv[3*dlItr+0] *= 0.1;
						ofv[3*dlItr+1] *= 0.1;
						ofv[3*dlItr+2] *= 0.1;
					*
					}
					else
					{
						ofv[3*dlItr+0] = 0;
						ofv[3*dlItr+1] = 0;
						ofv[3*dlItr+2] = 0;
					}
				}*/
				
				memcpy(minusDL, ofv, sizeof(double)*3*inState->totalVerts);
				
		/*		{
					int cItr, vItr, dVdtItr;
					for( cItr=0, dVdtItr=0; cItr<inLink->nc; cItr++ )
					{
						for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++, dVdtItr++ )
						{
							inLink->cp[cItr].vt[vItr].c[0] += inState->correctionStepDefault*ofv[3*dVdtItr+0];
							inLink->cp[cItr].vt[vItr].c[1] += inState->correctionStepDefault*ofv[3*dVdtItr+1];
							inLink->cp[cItr].vt[vItr].c[2] += inState->correctionStepDefault*ofv[3*dVdtItr+2];
						}
					}
				}
		*/		
							
				free(ofvB);
				free(ofv);
				taucs_ccs_free(sparseAT);
				
		/*		// extra stuff due to changed method
				taucs_ccs_free(sparseA);
				free(A);
				free(minusDL);
				
				free(minradSet);
				
				if( outStruts == NULL )
					free(strutSet);
					
				if( inState->graphing[kRcond] != 0 )
				{
					inState->rcond = taucs_rcond(sparseA);
					
					if( gPaperInfoInTmp )
					{
						char fname[512];
						preptmpname(fname, "rcond", inState);
						FILE* rcondF = fopen(fname,"a");
						fprintf( rcondF, "%d %e\n", strutCount+minradLocs, inState->rcond );
						fclose(rcondF);
					}
				}

				return;
		*/		
				// now we need to normalize for tsnnls again -- oh wait, we're using 
				// just lsqr! perhaps won't need this if lsqr alone works...
			
				// this doesn't seem to work
			
			/*	for( mItr=0; mItr<minradLocs; mItr++ )
				{
					double norm;
					norm = 0;
					int rItr;
					// column is mItr, so just get the norm of all entries and normalize
					// since only implicated vertices are in column
					for( rItr=sparseA->colptr[mItr+strutCount]; rItr<sparseA->colptr[mItr+strutCount+1]; rItr++ )
					{
						norm += sparseA->values.d[rItr]*sparseA->values.d[rItr];
					}
					norm = sqrt(norm);
					
					for( rItr=sparseA->colptr[mItr+strutCount]; rItr<sparseA->colptr[mItr+strutCount+1]; rItr++ )
					{
						 sparseA->values.d[rItr] /= norm;
					}
				}
			*/	
			/*	for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
				{
					plc_vector  vert;
					double  norm;
					vert.c[0] = minusDL[3*dlItr+0];
					vert.c[1] = minusDL[3*dlItr+1];
					vert.c[2] = minusDL[3*dlItr+2];
					norm = plc_M_norm(vert);
					if( norm != 0 )
					{
						minusDL[3*dlItr+0] /= norm;
						minusDL[3*dlItr+1] /= norm;
						minusDL[3*dlItr+2] /= norm;
					
						minusDL[3*dlItr+0] *= 0.1;
						minusDL[3*dlItr+1] *= 0.1;
						minusDL[3*dlItr+2] *= 0.1;
					}
					else
					{
						minusDL[3*dlItr+0] = 0;
						minusDL[3*dlItr+1] = 0;
						minusDL[3*dlItr+2] = 0;
					}
				}*/
			
		/*		if( gOutputFlag == 1 )
				{
					plc_vector  debug[500];
					for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
					{
						debug[dlItr].c[0] = 0.5*minusDL[3*dlItr+0];
						debug[dlItr].c[1] = 0.5*minusDL[3*dlItr+1];
						debug[dlItr].c[2] = 0.5*minusDL[3*dlItr+2];
					}
					exportVect(debug, inLink, "/tmp/minusDL.vect");
				}
		*/	
			} // gFastCorrectionSteps
			else
			{
				int rItr;
				// do this. follow grads.
				double	scaleFactor;
				double secondGreenZone = thickness - (thickness*inState->overstepTol)*.25;
				double greenZone = thickness - (thickness*inState->overstepTol)*0.5;
				
				//sparseA = taucs_construct_sorted_ccs_matrix(A, strutCount+minradLocs, 3*inState->totalVerts);	
				
				// struts are rows
				for( sItr=0; sItr<strutCount; sItr++ )
				{
					if( strutSet[sItr].length < greenZone )
					{
						scaleFactor = thickness-strutSet[sItr].length;
						// all we need to do is copy over the column in an appropriately 
						// scaled way.
						for( rItr=cleanA->colptr[sItr]; rItr<cleanA->colptr[sItr+1]; rItr++ )
						{
							minusDL[cleanA->rowind[rItr]] = scaleFactor*cleanA->values.d[rItr];
						}
					} // length < greenZone
				}// for all struts
				
				// exact same thing for MR struts
				for( sItr=strutCount; sItr<strutCount+minradLocs; sItr++ )
				{
					double minradSecondGreenZone = ((thickness/2.0)-inState->minminrad) / 2.0;
					minradSecondGreenZone = (thickness/2.0)-minradSecondGreenZone;
					// we're in green green, and we should shorten rather than lengthen, but not so much that
					// we dip below the min again
					if( minradSet[sItr-strutCount].mr < minradSecondGreenZone )
					{
						scaleFactor = ((thickness/2.0)-minradSet[sItr-strutCount].mr);
						
						for( rItr=cleanA->colptr[sItr]; rItr<cleanA->colptr[sItr+1]; rItr++ )
						{
							minusDL[cleanA->rowind[rItr]] += scaleFactor*cleanA->values.d[rItr];
						}
					}
				}

			} // doing fast corrections?
		}
			
			// if we're only doing 1 itr and we're in paper recording mode, 
			// save the rigidity matrix is /tmp/ for the MATLAB
			// condition number calculation.
			if( inState->maxItrs == 1 && gPaperInfoInTmp != 0 )
			{
				FILE* matrixF = fopen("/tmp/rigidity","w");
				taucs_print_ccs_matrix(cleanA, matrixF);
				fclose(matrixF);
			}
										
		// solve AX = -dl, x is strut compressions
		
			// tsnnls routines are sentsitive to sorting of rows in ccs matrices, 
			// so we enforce the rules here at low cost.
			taucs_enforce_ccs_sort(cleanA, strutCount);
		
		//	cleanA = taucs_construct_sorted_ccs_matrix(taucs_convert_ccs_to_doubles(cleanA),
		//				strutCount+minradLocs, 3*inState->totalVerts);
			
			// here, we compare
			int rItr;
		/*	for( cItr=0; cItr<strutCount+minradLocs; cItr++ )
			{
				for( rItr=0; rItr<sparseA->colptr[cItr+1]-sparseA->colptr[cItr]; rItr++ )
				{
					if( sparseA->rowind[sparseA->colptr[cItr]+rItr] != cleanA->rowind[cleanA->colptr[cItr]+rItr] )
					{
						printf("row indices disagree\n");
						exit(-1);
					}
					if( sparseA->values.d[sparseA->colptr[cItr]+rItr] != cleanA->values.d[cleanA->colptr[cItr]+rItr] )
					{
						printf("values disagree\n");
						exit(-1);
					}
				}
			}
		*/	

		/*	FILE* foo = fopen("foo","w");
			taucs_print_ccs_matrix(cleanA, foo);
			fclose(foo);
		*/
			if( gPaperInfoInTmp == 0 )
			{
				compressions = t_snnls(cleanA, minusDL, &inState->residual, NULL, 2, 0);
			}
			else
			{
				compressions = t_snnls(cleanA, minusDL, &inState->residual, inState->perVertexResidual, 2, 0);
			}
			
			
			if( compressions == NULL )
			{
				printf( "****** NULL compressions!\n" );
				fprintf( stderr, "****** NULL compressions!\n" );
				exit(-1);
				return;
			}
			
	//		for( foo=0; foo<sparseA->n; foo++ )
	//			printf( "(%d) %lf ", minradSet[foo].vert, compressions[foo] );
	//		printf( "\n" );
	//	else
	//	{
	/*		int *F;
			int i;
			F = malloc(sizeof(int)*sparseA->n);
			for( i=0; i<sparseA->n; i++ )
				F[i] = i;
			compressions = t_snnlslsqr(sparseA, minusDL, taucs_ccs_aprime_times_a(sparseA), F, 1, 1);
			free(F);
	*/	
			/* if there are no constrained struts, t_snnls won't actually work, so use SOL LSQR *
			lsqr_input   *lsqr_in;
			lsqr_output  *lsqr_out;
			lsqr_work    *lsqr_work;
			lsqr_func    *lsqr_func;
			int bItr;
								
			alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, sparseA->m, sparseA->n );
			
			/* we let lsqr() itself handle the 0 values in this structure *
			lsqr_in->num_rows = sparseA->m;
			lsqr_in->num_cols = sparseA->n;
			lsqr_in->damp_val = 0;
			lsqr_in->rel_mat_err = kZeroThreshold;
			lsqr_in->rel_rhs_err = 0;
			lsqr_in->cond_lim = 0;
			lsqr_in->max_iter = lsqr_in->num_rows + lsqr_in->num_cols + 5000;
			lsqr_in->lsqr_fp_out = NULL;	
			for( bItr=0; bItr<sparseA->m; bItr++ )
			{
				lsqr_in->rhs_vec->elements[bItr] = minusDL[bItr];
			}
			/* Here we set the initial solution vector guess, which is 
			 * a simple 1-vector. You might want to adjust this value for fine-tuning
			 * t_snnls() for your application
			 *
			for( bItr=0; bItr<sparseA->n; bItr++ )
			{
				lsqr_in->sol_vec->elements[bItr] = 1; 
			}
			
			/* This is a function pointer to the matrix-vector multiplier *
			lsqr_func->mat_vec_prod = sparse_lsqr_mult;
			
			lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, sparseA );
			
			compressions = (double*)malloc(sizeof(double)*sparseA->n);
			for( bItr=0; bItr<sparseA->n; bItr++ ) // not really bItr here, but hey
				compressions[bItr] = lsqr_out->sol_vec->elements[bItr];
			
			free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );*/
	//	}
		
		if( dlenStep == 0 )
			gDeferredStrutExport = 1;
		
		if( (gOutputFlag == 1 && dlenStep != 0) || (gDeferredStrutExport != 0 && dlenStep != 0) || gVerboseFiling )
		{
			export_struts(inLink, strutSet, strutCount, compressions, inState);
			export_ted(inLink, strutSet, strutCount, minradSet, minradLocs, compressions, inState);
			gDeferredStrutExport = 0;
		}

		// if we are graphing rcond, we should record it here
		if( inState->graphing[kRcond] != 0 )
		{
	//		taucs_ccs_matrix* apda = taucs_ccs_aprime_times_a(sparseA);
			inState->rcond = taucs_rcond(cleanA);
			
			if( gPaperInfoInTmp )
			{
				char fname[512];
				preptmpname(fname, "rcond", inState);
				FILE* rcondF = fopen(fname,"a");
				fprintf( rcondF, "%d %e\n", strutCount+minradLocs, inState->rcond );
				fclose(rcondF);
			}
			
	//		taucs_ccs_free(apda);
		}
		
		inState->tsnnls_evaluations++;
		// dump results in reload.m format for benchmarking
		if( compressions == NULL )
		{
			// probably ill-conditioned, try predictor-corrector algorithm
		//	compressions = predictor_corrector(A, 3*inState->totalVerts, strutCount, minusDL);
			exit(-1);
		}
				
		// we scale compressions by their thickness overrun to create an 
		// incentive to return to thickness 1 that CANNOT be resolved by the 
		// least squares solver
	/*	for( sItr=0; sItr<strutCount; sItr++ )
		{
			if( strutSet[sItr].length < (thickness-(0.0001*thickness)) )
			{
				compressions[sItr] *= 1.01*(thickness/strutSet[sItr].length*thickness/strutSet[sItr].length);
			}
		}
	*/		
		// debug code
	/*	for( dlItr=0; dlItr<strutCount; dlItr++ )
			printf( "strut %d compression: %lf\n", dlItr, compressions[dlItr] );
	*/	
		
		dVdt = (plc_vector*)calloc(inState->totalVerts, sizeof(plc_vector));
		// dVdt = dl + A*compressions is this for loop
		for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
		{
			double partialMult = 0;
			int totalStruts = strutCount + minradLocs;
			
			for( sItr=0; sItr<totalStruts; sItr++ )
			{
				if( !(3*dlItr < cleanA->rowind[cleanA->colptr[sItr]] || 3*dlItr > cleanA->rowind[cleanA->colptr[sItr+1]-1] ) )
				{
					if( sItr < strutCount )
					{
					//	partialMult += compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+0];
						partialMult += compressions[sItr]*rigidityEntry(cleanA, strutCount, minradLocs, 3*dlItr, sItr);
					}
					else // we're in minrad land (-- actually this means nothing now)
					{
					//	partialMult += 1*compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+0];
						partialMult += 1*compressions[sItr]*rigidityEntry(cleanA, strutCount, minradLocs, 3*dlItr, sItr);
					}
				}
			}
			dVdt[dlItr].c[0] = dl[dlItr].c[0] + partialMult;
		
			partialMult=0;
			for( sItr=0; sItr<totalStruts; sItr++ )
			{
				if( !(3*dlItr < cleanA->rowind[cleanA->colptr[sItr]] || 3*dlItr > cleanA->rowind[cleanA->colptr[sItr+1]-1] ) )
				{
					if( sItr < strutCount )
					{
					//	partialMult += compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+totalStruts];
						partialMult += compressions[sItr]*rigidityEntry(cleanA, strutCount, minradLocs, 3*dlItr + 1, sItr);
					}
					else
					{
						partialMult += 1*compressions[sItr]*rigidityEntry(cleanA, strutCount, minradLocs, 3*dlItr + 1, sItr);
				//		partialMult += 1*compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+totalStruts];
					}
				}
			}
			dVdt[dlItr].c[1] = dl[dlItr].c[1] + partialMult;
			
			partialMult = 0;
			for( sItr=0; sItr<totalStruts; sItr++ )
			{
				if( !(3*dlItr < cleanA->rowind[cleanA->colptr[sItr]] || 3*dlItr > cleanA->rowind[cleanA->colptr[sItr+1]-1] ) )
				{
					if( sItr < strutCount )
					{
				//		partialMult += compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+(2*totalStruts)];
						partialMult += compressions[sItr]*rigidityEntry(cleanA, strutCount, minradLocs, 3*dlItr + 2, sItr);
					}
					else
					{
				//		partialMult += 1*compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+(2*totalStruts)];
						partialMult += compressions[sItr]*rigidityEntry(cleanA, strutCount, minradLocs, 3*dlItr + 2, sItr);
					}
				}
			}
			dVdt[dlItr].c[2] = dl[dlItr].c[2] + partialMult;
		}
		
		if( gOutputFlag )
		{
			int next0, next1, pItr;
			
			computeCompressPush( inLink, strutSet, minradSet, strutCount, minradLocs );
		
			totalPushes = (double*)calloc(inState->totalVerts, sizeof(double));
			inState->maxPush = 0;
			// note to self - we're not currently including minrad struts here
			for( sItr=0; sItr<(strutCount); sItr++ )
			{
				if( inLink->cp[strutSet[sItr].component[0]].acyclic != 0 )
					next0 = strutSet[sItr].lead_vert[0] + 1;
				else
					next0 = ((strutSet[sItr].lead_vert[0] + 1) % inLink->cp[strutSet[sItr].component[0]].nv);
			
				if( inLink->cp[strutSet[sItr].component[1]].acyclic != 0 )
					next1 = strutSet[sItr].lead_vert[1] + 1;
				else
					next1 = ((strutSet[sItr].lead_vert[1] + 1) % inLink->cp[strutSet[sItr].component[1]].nv);
			
				strutSet[sItr].lead_vert[0] += inState->compOffsets[strutSet[sItr].component[0]];
				strutSet[sItr].lead_vert[1] += inState->compOffsets[strutSet[sItr].component[1]];
				next0 += inState->compOffsets[strutSet[sItr].component[0]];
				next1 += inState->compOffsets[strutSet[sItr].component[1]];
				
				totalPushes[strutSet[sItr].lead_vert[0]] += compressions[sItr]*(1-strutSet[sItr].position[0]);
				totalPushes[next0] += compressions[sItr]*(strutSet[sItr].position[0]);
				totalPushes[strutSet[sItr].lead_vert[1]] += compressions[sItr]*(1-strutSet[sItr].position[1]);
				totalPushes[next1] += compressions[sItr]*(strutSet[sItr].position[1]);
			
				strutSet[sItr].lead_vert[0] -= inState->compOffsets[strutSet[sItr].component[0]];
				strutSet[sItr].lead_vert[1] -= inState->compOffsets[strutSet[sItr].component[1]];
			}
			
			for( pItr=0; pItr<inState->totalVerts; pItr++ )
			{
				if( totalPushes[pItr] > inState->maxPush )
					inState->maxPush = totalPushes[pItr];
			}
			
			if( gOutputFlag /*&& dlenStep != 0*/ )
			{
				char	fname[1024];
				preptmpname(fname,"compress_push.vect",inState);
				export_pushed_edges(inLink, inState, totalPushes, fname, 0);
				if( inState->fancyVisualization != 0 )
				{
					fprintf(inState->fancyPipe, "(delete g0)\n");
					fprintf(inState->fancyPipe, "(load %s)\n", fname);
					fprintf(inState->fancyPipe, "(look g0)\n");
					fflush(inState->fancyPipe);
				}
			}
			free(totalPushes);
		}
		
		// clean up
		taucs_ccs_free(cleanA);
		//free(A);
		free(compressions);
		// we only free this if the user wasn't interested in keeping it
		if( outStruts == NULL )
			free(strutSet);
		free(minradSet);
		// we free this here and not outside the else{} since dVdt _IS_ dl if
		// strutCount == 0, which is the other case.
	/*	if( gOutputFlag == 1 )
		{
			exportVect(dl, inLink, "/tmp/dl.vect");
		}*/
		//free(dl);
		// throw dVdt back into dl, which on output, now includes dVdt
		for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
		{
			dl[dlItr].c[0] = dVdt[dlItr].c[0];
			dl[dlItr].c[1] = dVdt[dlItr].c[1];
			dl[dlItr].c[2] = dVdt[dlItr].c[2];
		}
		
		// let's try just moving along the ofv
/*		if( dlenStep == 0 )
		{
			// minusDL has ofv at this point
			for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
			{
				dl[dlItr].c[0] = minusDL[3*dlItr+0];
				dl[dlItr].c[1] = minusDL[3*dlItr+1];
				dl[dlItr].c[2] = minusDL[3*dlItr+2];
			}
		}
*/		
		if( gOutputFlag == 1 )
		{
			char	fname[1024];
			preptmpname(fname,"dVdt.vect",inState);
			exportVect(dVdt, inLink, fname);
		}
		
		free(minusDL);
		free(dVdt);
	}
		
//	return dVdt;
}

#pragma mark -

void
normalizeStruts( plc_vector* strutDirections, octrope_strut* strutSet, plCurve* inLink, int strutCount )
{
	int sItr;
	for( sItr=0; sItr<strutCount; sItr++ )
	{
		plc_vector  points[2];
		double norm = 0;
		octrope_strut_ends( inLink, &strutSet[sItr], points );
		
		// the normalized difference of pointOne, pointTwo is the strut force vector
		strutDirections[sItr].c[0] = points[0].c[0] - points[1].c[0];
		strutDirections[sItr].c[1] = points[0].c[1] - points[1].c[1];
		strutDirections[sItr].c[2] = points[0].c[2] - points[1].c[2];
		
		norm = plc_M_norm(strutDirections[sItr]);
		strutDirections[sItr].c[0] /= norm;
		strutDirections[sItr].c[1] /= norm;
		strutDirections[sItr].c[2] /= norm;
	}

}

double
maxovermin( plCurve* inLink, search_state* inState )
{
	int cItr, vItr;
	double max = 0, min = DBL_MAX;
	plc_vector s1, s2;
	double max_maxovermin = 0, len=0;
	int edges, totalEdges=0;
	int maxVert, maxComp, minVert, minComp;
	
	inState->avgSideLength = 0;
	int tot = 0;
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		max = 0;
		min = DBL_MAX;
		
		edges = plc_strand_edges(&inLink->cp[cItr]);
		totalEdges += edges;
		
		for( vItr=0; vItr<edges; vItr++ )
		{
			int cvertex;
			double norm;
			s1.c[0] = inLink->cp[cItr].vt[vItr].c[0];
			s1.c[1] = inLink->cp[cItr].vt[vItr].c[1];
			s1.c[2] = inLink->cp[cItr].vt[vItr].c[2];
			
			cvertex = vItr+1;
			
			// wrap around case
			if( vItr == (inLink->cp[cItr].nv-1) && inLink->cp[cItr].acyclic == 0 )
				cvertex = 0;
			
			s2.c[0] = inLink->cp[cItr].vt[cvertex].c[0];
			s2.c[1] = inLink->cp[cItr].vt[cvertex].c[1];
			s2.c[2] = inLink->cp[cItr].vt[cvertex].c[2];
			
			plc_M_sub_vect(s1, s2);
			norm = plc_M_norm(s1);
			
			inState->sideLengths[inState->compOffsets[cItr] + vItr] = norm;
			inState->avgSideLength += norm;
			tot++;
			
			len += norm;
		
			if( norm < min )
			{
				minVert = vItr;
				minComp = cItr;
				min = norm;
			}
			if( norm > max )
			{
				maxVert = vItr;
				maxComp = cItr;
				max = norm;
			}
		}
		if( (max/min) > max_maxovermin )
			max_maxovermin = (max/min);
	}
	
	if( gQuiet == 0 )
	{
		printf( "max edge: %f (%f/%f) at %d on %d / min edge: %f (%f/%f) at %d on %d\n", 
			max, fabs(max-(len/totalEdges)), (fabs(max-(len/totalEdges))/(len/totalEdges))*100,
			maxVert, maxComp,
			min, fabs(min-(len/totalEdges)), (fabs(min-(len/totalEdges))/(len/totalEdges))*100,
			minVert, minComp );
	}
	
	//printf( "max/min: %lf\n", max_maxovermin );
	
	inState->avgSideLength /= (double)tot;
	
	return max_maxovermin;
}

void
computeCompressPush( plCurve* inLink, octrope_strut* strutSet,
				octrope_mrloc* minradSet, int strutCount, int minradLocs )
{
	
}
