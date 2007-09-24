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

struct ccs_sort {
  
  int col;
  int row;
  double val;

};

int compare_ccs_sort (const void *a, const void *b) 

     /* Sort in dictionary order on (col,row). */

{

  struct ccs_sort *A = (struct ccs_sort *)a;
  struct ccs_sort *B = (struct ccs_sort *)b;

  if (A->col != B->col) {

    return A->col - B->col; 

  } else {

    return A->row - B->row;

  }

}
  

static void
taucs_enforce_ccs_sort( taucs_ccs_matrix* A )

     /* Note to self: need to r/w this to be more general. 
	We shouldn't need strutCount in order to do this job. */
{

  int nnz = A->colptr[A->n];
  struct ccs_sort *sortbuf = malloc(nnz*sizeof(struct ccs_sort));
  int colItr,eltItr;

  /* Load the buffer of ccs_sort structures. */

  for( colItr = 0; colItr < A->n; colItr++) {

    for ( eltItr = A->colptr[colItr]; eltItr < A->colptr[colItr+1]; eltItr++) {

      sortbuf[eltItr].col = colItr;
      sortbuf[eltItr].row = A->rowind[eltItr];
      sortbuf[eltItr].val = A->values.d[eltItr];

    }

  }

  /* Now sort. */

  qsort(sortbuf,nnz,sizeof(struct ccs_sort),compare_ccs_sort);
  
  /* Now rewrite the rowind and value arrays in this order. Since columns trump
     rows in the sort, the elements in each column should be in the same slots 
     they were in to start. Hence there is no need to update A->colptr. */
  
  for (eltItr = 0; eltItr < nnz; eltItr++ ) {

    A->rowind[eltItr] = sortbuf[eltItr].row;
    A->values.d[eltItr] = sortbuf[eltItr].val;

  }

  free(sortbuf);

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

void   *gOctmem;
int     gOctmem_size;

void 
bsearch_stepper( plCurve** inLink, search_state* inState )
{
  double maxmaxmin = 0;
  double minthickness = 500;
  double nextMovieOutput = 0.0;
  
  gConditionCheck = 0;
  double rop_20_itrs_ago = {DBL_MAX};
  double oldrops[20];
  int rItr;
  
#ifdef HAVE_CLOCK  
  
  clock_t startTime;
  startTime = clock();

#endif
  
  inState->steps = 0;
  int stepItr;

  /* We allocate a global buffer for all the octrope calls inside
     the stepper loop. Note that this means we are assuming that 
     the number of verts stays more or less constant inside here. */

  gOctmem_size = octrope_est_mem(plc_num_edges(*inLink));
  gOctmem = malloc_or_die(sizeof(char)*gOctmem_size, __FILE__ , __LINE__ );
  
  for( stepItr=0; /* Main loop, incorporates stopping criteria. */
       (stepItr < inState->maxItrs) && 
	 (rop_20_itrs_ago - inState->ropelength > inState->stop20) &&
	 (inState->residual > inState->residualThreshold);
       stepItr++ ) {
    
    int lastSet;      
    lastSet = inState->lastStepStrutCount;
    
    /* Decide whether to initiate thickness correction. */
    
    if( (inState->shortest < (2*inState->tube_radius*(1-inState->overstepTol))) ||
	(inState->minrad < gLambda*inState->tube_radius*(1-inState->minradOverstepTol))) {

      correct_thickness(*inLink,inState);

    }        
  
    /************************************************************************/
    
    *inLink = bsearch_step(*inLink, inState);
    
    /************************************************************************/
    
    inState->steps++;
    
    /* Manage data output */
      
    if( (stepItr%kOutputItrs) == 0 && gQuiet == 0 ) { 
      /* Note "k" means constant */
      
      update_runtime_display(*inLink,inState);  
      
    }

    if ( inState->steps%inState->loginterval == 0 ) {

      update_runtime_logs(inState);

    }

    if (inState->time >= nextMovieOutput) {

      update_vect_directory(*inLink,inState);
      nextMovieOutput += 0.041666666667*inState->moviefactor;
      
    } 

    /* Last, we update "inState" to keep track of running variables. */

    if( inState->shortest < minthickness )
      minthickness = inState->shortest;

    inState->oldLength = inState->length;
    inState->oldLengthTime = inState->cstep_time;

    double maxmin;
    maxmin = plCurve_long_edge(*inLink)/plCurve_short_edge(*inLink);
    inState->lastMaxMin = maxmin;
    if( maxmin > maxmaxmin )
      maxmaxmin = maxmin;
   
    /* And we update our list of old ropelength values. */

    for(rItr=19;rItr>0;rItr--) { oldrops[rItr] = oldrops[rItr-1]; }
    oldrops[0] = inState->ropelength;

    if (stepItr > 20) { rop_20_itrs_ago = oldrops[19]; }
    
  } 

  free(gOctmem);

  /* We have now terminated. The final output files will be written 
     in ridgerunner_main.c. However, we log the reason for termination. */

  logprintf("ridgerunner: run complete. Terminated because \n");

  if (!(stepItr < inState->maxItrs)) {
    
    logprintf("ridgerunner: reached maximum number of steps (%d).\n",
	      stepItr);

  }

  if (!(rop_20_itrs_ago - inState->ropelength > inState->stop20)) {

    logprintf("ridgerunner: change in rop over last 20 iterations %g > stop20 = %g.\n",
	      rop_20_itrs_ago - inState->ropelength, inState->stop20);

  }

  if (!(inState->residual > inState->residualThreshold)) {

    logprintf("ridgerunner: residual %g < residualThreshold %g\n",
	      inState->residual, inState->residualThreshold);

  }

}


void update_vect_directory(plCurve * const inLink, search_state *inState)

     /* If it is time for the next vect output, go ahead and write 
	another file to the appropriate directory. */
{
  char tmpfilename[1024];
  FILE *outfile;
  static int frames_written_since_compression = 0;
  
  sprintf(tmpfilename,"%s.%07d.vect",
	  inState->vectprefix,inState->steps);
  outfile = fopen_or_die(tmpfilename,"w", __FILE__ , __LINE__ );
  plc_write(outfile,inLink);
  fclose(outfile);
  
  frames_written_since_compression++;

  if (frames_written_since_compression == inState->maxmovieframes/2) { 

    compress_vectdir(inState);
    frames_written_since_compression = 0;
    inState->moviefactor *= 2;

  }
      
}

void compress_vectdir(const search_state *inState)

/* Deletes half of the files in the vect directory. */

{
  int  j = 0;
  char dirname[1024];
  DIR  *dirstream;
  struct dirent *ent;
  char fullname[2048];
  int deleted = 0;

  sprintf(dirname,"./%s.rr/vectfiles",inState->basename);  
  dirstream = opendir_or_die(dirname, __FILE__ , __LINE__ );

  while ( (ent = readdir(dirstream)) ) {  /* Note! The single = is not a bug. */
  
    if (strstr(ent->d_name,"vect") != NULL) { /* Make sure we have a vect file */

      j++;

      if (j % 2 == 0) { /* Delete the odd vect files. */

	sprintf(fullname,"%s/%s",dirname,ent->d_name);
	remove_or_die(fullname, __FILE__ , __LINE__ );
	deleted++;

      }

    }

  }

  closedir(dirstream);

  logprintf("Compressed vect file directory by deleting %d files.\n",
	    deleted);

}   
  

void open_runtime_logs(search_state *state,char code)

{
  char tmpfilename[1024];
  char opentype[10];
  int  i;

  sprintf(opentype,"%c",code);

  for(i=0;i<kTotalLogTypes;i++) {

    sprintf(tmpfilename,"./%s.rr/logfiles/%s.dat",
	    state->basename,
	    state->logfilenames[i]);

    state->logfiles[i] = fopen_or_die(tmpfilename,opentype, __FILE__ , __LINE__ );

  }

}

void close_runtime_logs(search_state *state)

{
  int i;

  for(i=0;i<kTotalLogTypes;i++) {

    fclose(state->logfiles[i]);

  }

}

void compress_runtime_logs(search_state *state)

/* We check the runtime logs against our maximum logfilesize. If the
   logs appear to be getting close to the maximum, we thin them by 
   deleting every other entry. */

{
  int i,j;

  char *linebuf;
  size_t linebufsize;

  char tmpfilename[4096];
  char logfilename[4096];
  
  int  tmpdes;
  FILE *tmpfile;

  for(i=0;i<kTotalLogTypes;i++) { 
    
    /* We begin by rebuilding the filename of this logfile. */
    
    sprintf(logfilename,"./%s.rr/logfiles/%s.dat",
	    state->basename,
	    state->logfilenames[i]);
    
    /* There is a special case to handle here. "lsqroutput" is not 
       a standard line-based logfile. The best we can do is delete 
       what we've got and start over. */
    
    if (strcmp(state->logfilenames[i],"lsqroutput") == 0) {
      
      fclose(state->logfiles[i]);
      state->logfiles[i] = fopen_or_die(logfilename,"w", __FILE__ , __LINE__ );
      
    } else { 
      
      /* We are actually going to try to preserve some data from the
	 start of the run. We now open a new file to read the data
	 we're going to save into. */
      
      sprintf(tmpfilename,"%sXXXXXX",logfilename);
      
      tmpdes = mkstemp_or_die(tmpfilename, __FILE__ , __LINE__ );
      tmpfile = fdopen_or_die(tmpdes,"w", __FILE__ , __LINE__ );
      
      /* We now close and reopen the original log file to rewind it. */
      
      fclose(state->logfiles[i]);
      state->logfiles[i] = fopen_or_die(logfilename,"r", __FILE__ , __LINE__ );
      
      /* We now copy the data from file to file. */
      
      linebuf = malloc_or_die(4096*sizeof(char), __FILE__ , __LINE__ );
      linebufsize = 4096*sizeof(char);
      
      for(j=0;fgets(linebuf,linebufsize,state->logfiles[i]) != NULL;j++) {
	
	if (j % 2 == 0) { fprintf(tmpfile,"%s",linebuf); }
	
      }
	
      free(linebuf);
	
      /* Finally, we delete the old file, replace it with the new, smaller file, */
      /* and reopen the new logfile for appending. */
      
      fclose(tmpfile);
      fclose(state->logfiles[i]);
      remove_or_die(logfilename, __FILE__ , __LINE__);
      rename_or_die(tmpfilename,logfilename, __FILE__ , __LINE__ );
      state->logfiles[i] = fopen_or_die(logfilename,"a", __FILE__ , __LINE__ );
   	
    }

  }

  /* Now we report our progress in the log. */
  
  if (VERBOSITY > 0) {
   
    logprintf("Compressed logfiles at step %d = (%d * %d).\n",
	      state->steps,state->loginterval,state->maxlogsize);
  
  }
}

void update_runtime_logs(search_state *state)

     /* We now update the various data logs for the run. */
     /* These logs are flushed every LOG_FLUSH_INTERVAL steps, */
     /* where this is defined in ridgerunner.h. */
{
  int i;
  static int logged_cstep_count = 0;
  static int log_entries_since_compression = 0;

  fprintf(state->logfiles[kLength],"%d %2.8g \n",state->steps,state->length);
  fprintf(state->logfiles[kRopelength],"%d %2.8g \n",state->steps,state->ropelength);
  fprintf(state->logfiles[kStrutCount],"%d %d %d\n",state->steps,
	  state->lastStepStrutCount,state->lastStepMinradStrutCount);
  fprintf(state->logfiles[kStepSize],"%d %g \n",state->steps,state->stepSize);
  fprintf(state->logfiles[kThickness],"%d %2.7g \n",state->steps,state->thickness);
  fprintf(state->logfiles[kMinrad],"%d %2.7g \n",state->steps,state->minrad);
  fprintf(state->logfiles[kResidual],"%d %g \n",state->steps,state->residual);
  fprintf(state->logfiles[kMaxOverMin],"%d %g \n",state->steps,state->lastMaxMin);
  fprintf(state->logfiles[kRcond],"%d %g \n",state->steps,state->rcond);

  log_entries_since_compression++;

#ifdef HAVE_TIME
  time_t now;
  int hrs,min,sec;
  double elapsed;

  now = time(NULL);
  elapsed = difftime(now,state->start_time);

  hrs = floor(elapsed/3600.0);
  elapsed -= 3600*hrs;
  
  min = floor(elapsed/60.0);
  elapsed -= 60*min;

  sec = floor(elapsed);

  fprintf(state->logfiles[kWallTime],"%d %d:%d:%d\n",state->steps,hrs,min,sec);
#else
  fprintf(state->logfiles[kWallTime],"Could not link with 'time()' at compile.\n");
#endif

  fprintf(state->logfiles[kMaxVertexForce],"%d %g \n",state->steps,state->maxPush);

  /* We only update the cstep log when we have recently corrected. */
  /* This log file marks everything by step number as well. */

  if (state->cstep_count != logged_cstep_count) {

    fprintf(state->logfiles[kCorrectionStepsNeeded],"%d %d\n",
	    state->last_cstep_attempts,state->steps);

    logged_cstep_count = state->cstep_count;

  }

  fprintf(state->logfiles[kEQVariance],"%d %g \n",state->steps,state->eqVariance);

  if (log_entries_since_compression%LOG_FLUSH_INTERVAL == 0) {

    for(i=0;i<kTotalLogTypes;i++) fflush(state->logfiles[i]);

  }

  if (state->steps == (state->loginterval*state->maxlogsize)) {

    compress_runtime_logs(state); /* Drop the size of each log by half.   */
    state->loginterval *= 2;      /* Log half as often from now on. */
    log_entries_since_compression = 0;

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

      plc_M_vmadd(inLink->cp[cItr].vt[vItr],stepSize,dVdt[dVdtItr])

    }
  }

  plc_fix_wrap(inLink);		

  /* Very important to call here! We have just modified the entries of a plCurve */
  /* directly, and such a call must _always_ be followed by fix_wrap. */

}

  
void convert_force_to_velocity(plc_vector *F, plCurve *inLink,search_state *inState)

     /* There is a difference between _force_ fields and _velocity_
	fields for non-equilateral polygons. Roughly, the curvature
	_force_ should be smaller when we deal with smaller edges,
	since there will be more struts, and the added-up force will
	be the same.  On the other hand, the _velocity_ of the
	curvature flow should be the same, regardless of
	discretization.

	This procedure goes ahead and rescales a force field to a 
	velocity field for a non-eq curve. If we decided to do this,
        it should be added to bsearch_step. */

{
  int vItr,cItr,Fitr=0;
  double avg, scale=1.0;

  inState->avgSideLength = plc_arclength(inLink,NULL)/(double)(plc_edges(inLink,NULL));
      
  for(cItr=0;cItr < inLink->nc;cItr++) {

    for(vItr=0;vItr < inLink->cp[cItr].nv;vItr++,Fitr++) {

      avg = 

	plc_distance(inLink->cp[cItr].vt[vItr-1],
		     inLink->cp[cItr].vt[vItr]) +
	plc_distance(inLink->cp[cItr].vt[vItr],
		     inLink->cp[cItr].vt[vItr+1]);

      avg /= ((vItr == 0 || vItr == inLink->cp[cItr].nv-1) && inLink->cp[cItr].open) ?
	1 : 2;
      
      scale = inState->avgSideLength/avg;

      plc_M_scale_vect(scale,F[Fitr])

    }
    
  }

}

/***********************************************************************/
/*                     Error correction code                           */
/***********************************************************************/ 

extern void sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod );
/* This is linked in from tsnnls. */

double*
stanford_lsqr( search_state *inState, 
	       taucs_ccs_matrix* sparseA, double* minusDL)
{
  /* if there are no constrained struts, t_snnls won't actually work, so use SOL LSQR */
  lsqr_input   *lsqr_in;
  lsqr_output  *lsqr_out;
  lsqr_work    *lsqr_work;
  lsqr_func    *lsqr_func;
  int bItr;
  double*		result;
  char errmsg[1024];
  
  alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, sparseA->m, sparseA->n );
  
  /* we let lsqr() itself handle the 0 values in this structure */
  lsqr_in->num_rows = sparseA->m;
  lsqr_in->num_cols = sparseA->n;
  lsqr_in->damp_val = 0;
  lsqr_in->rel_mat_err = kZeroThreshold;
  lsqr_in->rel_rhs_err = 0;
  lsqr_in->cond_lim = 1/(10*sqrt(DBL_EPSILON));
  lsqr_in->max_iter = 4*lsqr_in->num_cols;  /* Suggested by lsqr docs */
  lsqr_in->lsqr_fp_out = inState->logfiles[klsqrlog];	

  for( bItr=0; bItr<sparseA->m; bItr++ ) {

    lsqr_in->rhs_vec->elements[bItr] = minusDL[bItr];
  
  }
  
  /* Here we set the initial solution vector guess, which is 
   * a vector of zeros. We might want to adjust this value later.
   */

  for( bItr=0; bItr<sparseA->n; bItr++ ) {

    lsqr_in->sol_vec->elements[bItr] = 0; 
  
  }
	
  /* This is a function pointer to the matrix-vector multiplier */
  lsqr_func->mat_vec_prod = sparse_lsqr_mult;
  
  
  /**************************************************************/

  lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, sparseA );

  /**************************************************************/
  

  /* We now check the results, and return. */

  if (lsqr_out->term_flag == 3 || lsqr_out->term_flag == 6) {

    dumpAxb_sparse(inState,sparseA,NULL,minusDL);
    sprintf(errmsg,
	    "ridgerunner: lsqr failed in correction stepper because the\n"
	    "             condition number of %d x %d matrix A^T > %g.\n"
	    "             The matrix A^T has been dumped to A.dat and \n"
	    "             the vector of desired corrections dumped to b.dat.\n"
	    "\n"
	    "             Check the log file 'lsqrlog' for more output.\n",
	    sparseA->m,sparseA->n,lsqr_in->cond_lim);
    FatalError(errmsg, __FILE__ , __LINE__ );

  }

  if (lsqr_out->term_flag == 7) {

    dumpAxb_sparse(inState,sparseA,NULL,minusDL);
    sprintf(errmsg,
	    "ridgerunner: lsqr failed to converge in correction stepper\n"
	    "             on %d x %d matrix A^T after %ld iterations.\n"
	    "             The matrix A^T has been dumped to A.dat and \n"
	    "             the vector of desired corrections dumped to b.dat.\n"
	    "\n"
	    "             Check the log file 'lsqrlog' for more output.\n" ,
	    sparseA->m,sparseA->n,lsqr_in->max_iter);
    FatalError(errmsg, __FILE__ , __LINE__ );

  }
  
  result = (double*)malloc(sizeof(double)*sparseA->n);
  fatalifnull_(result);
  
  for( bItr=0; bItr<sparseA->n; bItr++ ) // not really bItr here, but hey
    result[bItr] = lsqr_out->sol_vec->elements[bItr];
  
  inState->ofvResidual = lsqr_out->resid_norm;
  inState->ofvNorm = lsqr_out->sol_norm;
  
  free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );
  
  return result;
}



double *thicknessError(plCurve *inLink,
		       int strutCount, octrope_strut *strutList,
		       int minradCount, octrope_mrloc *mrList,
		       search_state *inState) 

     /* We use the strutList and mrList passed in the header, ignoring
	the versions stored in inState, which may not apply to this link. */

{
  int constraintCount = plCurve_score_constraints(inLink);
  double* C = (double*)(malloc((strutCount+minradCount+constraintCount)*sizeof(double)));
  
  /* We now create the C (corrections) vector giving the desired adjustment 
     to each strut and minrad length. */
  
  int sItr;
  
  for( sItr=0; sItr < inState->lastStepStrutCount ; sItr++) {
    
    double secondGreenZone = 
      2*inState->tube_radius*(1 - inState->overstepTol*.25);
    double greenZone = 
      2*inState->tube_radius*(1 - inState->overstepTol*0.5);
    
    // this only works if there aren't minrad struts for 
    // deep numerical reasons. talk to jason or me.
    if( strutList[sItr].length > secondGreenZone && minradCount == 0 ) {
      
      // Actually make these struts _shorter_ (!) to hold them in place.
      C[sItr] = -(strutList[sItr].length-(secondGreenZone));
      
    } else if( strutList[sItr].length < greenZone ) {
      
      C[sItr] = secondGreenZone-strutList[sItr].length;
      
    } else { /* length is between greenZone and secondGreenZone */
      
      C[sItr] = 0; /* We're happy, request no adjustment. */
      
    }
    
  }
  
  /* Now we deal with the minRad struts */
  
  for( sItr=strutCount; 
       sItr<strutCount+minradCount;
       sItr++ ) {
    
    double minradSecondGreenZone  /* Note used to be set by minminrad? */
      = gLambda*(inState->tube_radius)*(1 - inState->minradOverstepTol*0.25);

    double minradGreenZone 
      = gLambda*(inState->tube_radius)*(1 - inState->minradOverstepTol*0.5);
    
    /* At one point, we had a similar scheme as with struts: attempt
       to get everything between minradGreenZone and
       minradSecondGreenZone, making things worse if need be. This
       doesn't seem to have worked for numerical reasons (maybe b/c
       the minrad gradients are so much larger than the strut force
       gradients?) so we don't actually _do_ anything different here. */
    
    if( mrList[sItr-strutCount].mr > minradSecondGreenZone ) {
      
      /*ofvB[sItr] = minradSet[sItr-strutCount].mr - minradSecondGreenZone;
	greenZoneMR[greenZoneMRCount++] = sItr; // keep in mind this is offset */
      
      C[sItr] = 0;
      
    } else if (mrList[sItr-strutCount].mr < minradGreenZone) {
      
      C[sItr] = minradSecondGreenZone - mrList[sItr-strutCount].mr;
      
    } else { /* Between minradGreenZone and minradSecondGreenZone */
      
      C[sItr] = 0;
      
    }
    
  }
  
  /* Now we deal with the constraints. This gets a little bit
     interesting, since the values to which the constraints are
     supposed to be set have to be deduced on the fly from the
     constraint data structure in the plCurve inState. */
  
  plc_constraint *thisCst;
  int vItr;
  plc_vector errVect;
  char errmsg[1024];
  
  for( sItr = strutCount + minradCount,
	 thisCst = inLink->cst;
       thisCst != NULL;
       thisCst = thisCst->next ) {
    
    for (vItr = thisCst->vert; vItr < thisCst->vert+thisCst->num_verts; vItr++) {
      
      if (thisCst->kind == unconstrained) {
	
      } else if (thisCst->kind == fixed) {
	
	errVect = plc_vect_diff(inLink->cp[thisCst->cmp].vt[vItr],thisCst->vect[0]);
	
	C[sItr++] = -errVect.c[0];
	C[sItr++] = -errVect.c[1];
	C[sItr++] = -errVect.c[2];
	
      } else if (thisCst->kind == plane) {
	
	C[sItr++] = -(plc_dot_prod(inLink->cp[thisCst->cmp].vt[vItr],thisCst->vect[0]) - 
		      plc_dot_prod(thisCst->vect[1],thisCst->vect[0]));
	
      } else if (thisCst->kind == line) {

	sprintf(errmsg,
		"ridgerunner: line constraints are not implemented,"
		"              and should have been detected before now.");
	FatalError(errmsg, __FILE__ , __LINE__ );
	
      } else {
	
	sprintf(errmsg,
		"ridgerunner: Unknown constraint type %d in inLink.\n",
		thisCst->kind);
	FatalError(errmsg, __FILE__ , __LINE__ );
	
      }
      
    }
    
  }

  return C;

}

double l2norm(double *V, int N)

     /* Computes the l2 norm of a vector of size N */

{

  double val = 0;
  int i;

  for(i=0;i<N;i++) { val += V[i]*V[i]; }
  return sqrt(val);

}

void correct_thickness(plCurve *inLink,search_state *inState) 

     /* Newton's method correction algorithm for thickness. */
     
     /* Attempts to fix things so that the shortest strut and 
	worst minrad vertex are within "greenZone" of obeying 
	the constraints. 

	May take most of runtime if the number of verts is large,
	since we use the (slow) Stanford LSQR code to compute
	direction for Newton steps. */
     
{
  double stepSize = 1.0;
  double currentError = 0, newError;
  double greenZone = 
    2*inState->tube_radius*(1 - 0.5*inState->overstepTol);
  double mrgreenZone = 
    gLambda*inState->tube_radius*(1-0.5*inState->minradOverstepTol);
  int Csize = inState->lastStepStrutCount+inState->lastStepMinradStrutCount+
    plCurve_score_constraints(inLink);
 

  /* The "green zone" is half of the error allowed by the overstepTol
     variables. We initiate correction if we are less than the
     overstepTol amount, but terminate correction if we are less than
     half this amount. This keeps the strut set relatively stable. */

  /* To compute the error after a trial step is taken, we'll have to
     call octrope. This will require some setup. */
  
  plCurve *workerLink = NULL;
  double  *C = NULL,*workerC = NULL;
  
  int      strutCount;
  int      minradCount;
  
  octrope_strut *strutSet = NULL;
  octrope_mrloc *minradSet = NULL;
  
  int strutStorageSize = 6*inState->totalVerts;
  int minradStorageSize = inState->totalVerts;

  if (VERBOSITY >= 5) { logprintf("Correction step...\n"); }
  
  strutSet = (octrope_strut *)(malloc(sizeof(octrope_strut)*strutStorageSize));
  minradSet = (octrope_mrloc *)(malloc(sizeof(octrope_mrloc)*minradStorageSize));
  
  fatalifnull_(strutSet); 
  fatalifnull_(minradSet);
  
  double rop, thi, len, minrad, shortest;				  

  inState->cstep_count++;

  /* This is the main loop for Newton's method. */

  for(gCorrectionAttempts=0;  
      ((inState->shortest < greenZone) ||
	(inState->minrad < mrgreenZone)) && gCorrectionAttempts < gMaxCorrectionAttempts;
      gCorrectionAttempts++) {

    taucs_ccs_matrix *sparseA = NULL;
    taucs_ccs_matrix *sparseAT = NULL;
    
    /* We start by building a rigidity matrix, and taking its' transpose */

    sparseA = buildRigidityMatrix(inLink,inState);  /* Note: strut set in inState updates */
    sparseAT = taucs_ccs_transpose(sparseA);

    C = thicknessError(inLink,
		       inState->lastStepStrutCount,inState->lastStepStruts,
		       inState->lastStepMinradStrutCount,inState->lastStepMRlist,
		       inState);

    currentError = l2norm(C,Csize);

    /* We have now constructed the correction vector C and recorded its' norm for 
       stepping purposes. The next step is to call lsqr. */
	  
    /*
     * Fact: The rigidity matrix A (compressions) = (resulting motions of verts).
     * So it's also true that
     * 
     *		  A^T (a motion of verts) = (resulting change in edge length)
     *
     *	Now we _have_ a desired change in edge lengths, namely C.  So
     *	we really ought to compute velocity for a correction step by
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
    plc_vector *ofv_vect;
	  
    /* Jason note: we might have to change the lsqr initial guess later
     * me: remember the funky motion thing?
     *
     * If we have to change it: 
     *	want: linear combo of gradient of active constraints 
     *	      coeffs are k's computed just like k in test
     *	      (the change you want in that constraint divided
     *	       by square of norm o' gradient)
     */
    
    if (VERBOSITY >= 10) { logprintf("\tCalling stanford_lsqr..."); }

    ofv = stanford_lsqr(inState,sparseAT,C); 

    if (VERBOSITY >= 10) { logprintf("ok\n"); }

    /* This will die if it doesn't work, so we don't need error checking. */

    /* Now convert to plc_vector format for step's benefit. */
    /* This could be done faster, but more dangerously, by recasting */
    /* the ofv pointer. */

    ofv_vect = (plc_vector *)(malloc(inState->totalVerts*sizeof(plc_vector)));
    fatalifnull_(ofv_vect);

    int i;

    for(i=0;i<inState->totalVerts;i++) {

      ofv_vect[i].c[0] = ofv[3*i];
      ofv_vect[i].c[1] = ofv[3*i+1];
      ofv_vect[i].c[2] = ofv[3*i+2];

    }
    
    /* Now we have computed the Newton direction and are ready to take a step. */

    double stepSize = 2;
    double alpha = 1e-4; /* 10^-4 is a standard alpha for this kind of thing */

    do {

      if (VERBOSITY >= 10) { logprintf("Attempting cstep at size %2.8g.\n",stepSize/2); }

      stepSize /= 2;

      if (workerLink != NULL) plc_free(workerLink);
      workerLink = plc_copy(inLink);
      
      step(workerLink,stepSize,ofv_vect,inState);

      /* We now have to compute the error in the proposed configuration. */

      inState->octrope_calls++;

      octrope(workerLink, 
	      
	      &rop,
	      &thi,		
		
	      &len,		
	      &minrad,
	      &shortest,
		
	      // minrad struts
	      gLambda*inState->tube_radius,  /* Cutoff reports all mrlocs less than this */
	      0,                             /* This value will be ignored. */
	      minradSet, 
	      minradStorageSize,
	      &minradCount,
		
	      // strut info
	      2*inState->tube_radius, /* Cutoff reports struts shorter than this. */
	      0,                      /* This epsilon value also ignored. */
	      strutSet,
	      strutStorageSize,
	      &strutCount,
		
	      NULL, 0,                /* Let octrope allocate its' own memory. */

	      gLambda);               /* The global "stiffness" parameter. */  
      
      /* We now need to make sure that we didn't exceed the size of the strut
	 and/or minrad buffers. */
      
      if (strutCount >= strutStorageSize || minradCount >= minradStorageSize ||
	  strutCount < 0 || minradCount < 0) {

	char dumpname[1024];
	char errmsg[1024];
	
	dumpLink(inLink,inState,dumpname);
	sprintf(errmsg,
		"ridgerunner: octrope found %d struts and %d mrlocs on\n"
		"             the %d vertex link (dumped to %s), too close to\n"
		"             minradStorageSize of %d or strutStorageSize %d.\n",
		strutCount,minradCount,inState->totalVerts,dumpname,
		minradStorageSize,strutStorageSize);
	FatalError(errmsg, __FILE__ , __LINE__ );
	
      }

      /* Now we can compute the new error. */

      if (workerC != NULL) { free(workerC); }

      workerC = thicknessError(workerLink,
			       strutCount,strutSet,
			       minradCount,minradSet,
			       inState);

      newError = l2norm(workerC,Csize);

    } while (newError > (1 - alpha*stepSize)*currentError);

    /* This "sufficient decrease" condition comes from Kelley,
       "Solving nonlinear equations with Newton's method", p. 13. 

       We have now chosen a stepSize, stepped, and accepted the
       results.  Before we return to the top of the main for-loop and
       compute a new Newton direction, we must move the temp data in
       workerLink (and friends) to current data in inLink and
       inState. */

    step(inLink,stepSize,ofv_vect,inState); /* Trust me-- this was better than copying. */

    /* Now we update inState. */

    assert(inState->lastStepStruts != NULL && inState->lastStepMRlist != NULL);
   
    free(inState->lastStepStruts);
    free(inState->lastStepMRlist);
    
    inState->lastStepStruts = (octrope_strut *)(malloc(strutCount*sizeof(octrope_strut)));
    inState->lastStepMRlist = (octrope_mrloc *)(malloc(minradCount*sizeof(octrope_mrloc)));

    fatalifnull_(inState->lastStepStruts);
    fatalifnull_(inState->lastStepMRlist);

    int sItr;

    for(sItr=0;sItr<strutCount;sItr++) { 
      
      inState->lastStepStruts[sItr] = strutSet[sItr];

    }

    for(sItr=0;sItr<minradCount;sItr++) {

      inState->lastStepMRlist[sItr] = minradSet[sItr];

    }
    
    inState->length = len;
    inState->thickness = thi;
    inState->ropelength = rop;
    inState->shortest = shortest;
    inState->minrad = minrad;

    /* Finally, we free memory used in this step. */

    free(C);
    free(ofv);
    free(ofv_vect);

    taucs_ccs_free(sparseA);
    taucs_ccs_free(sparseAT);
    
  }

  /* We have now either converged or reached gMaxCorrectionAttempts. */

  if (gCorrectionAttempts >= gMaxCorrectionAttempts) {

    char dumpname[1024], errmsg[1024];

    dumpLink(inLink,inState,dumpname);
    sprintf(errmsg,
	    "ridgerunner: Newton's method error correction solver has failed\n"
	    "             to converge after %d iterations (maximum of %d). \n"
	    "             Current iterate dumped to %s. At this point\n"
	    "\n"
	    "             shortest strut has length %g \n"
	    "             minrad values is %g \n"
	    "\n"
	    "             last ofvnorm was %g \n"
	    "             stepsize is %g \n"
	    "\n"
	    "             Can set --MaxCorrectionAttempts=XXX from command\n"
	    "             line and try again or debug solver problem.\n",
	    gCorrectionAttempts,gMaxCorrectionAttempts,dumpname,
	    inState->shortest, inState->minrad, inState->ofvNorm,
	    stepSize);

    FatalError(errmsg, __FILE__ , __LINE__ );

  }

  inState->last_cstep_attempts = gCorrectionAttempts; 

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


void
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
  double  lastDCSD, lastMR;
  plCurve* workerLink;
	
  double ERROR_BOUND = 1e-5;
  double MR_ERROR_BOUND = 1e-5;
  
  // create initial vector field for which we want to move along in this step
  plc_vector  *dVdt;
  plc_vector  *dLen; 

  if (VERBOSITY >= 5) {  /* Verbose or higher */ 

    logprintf("Bsearch step %d.\n",inState->steps);

  }

  dLen = (plc_vector *)(calloc(inState->totalVerts, sizeof(plc_vector)));
  /* Note: this must be calloc since the xxxForce procedures add to given buffer. */
  fatalifnull_(dLen);
    
  dlenForce(dLen,inLink,inState);
  eqForce(dLen,inLink,inState);
  specialForce(dLen,inLink,inState);

  dVdt = resolveForce(dLen,inLink,inState); 
  /* Built from the bones of firstVariation. */

  free(dLen);

  /* Now we loop to find largest stepsize which doesn't violate constraints. */

  if( inState->shortest > 2*inState->tube_radius )
    inState->shortest = 2*inState->tube_radius;
			
  lastDCSD = inState->shortest;
  lastMR = inState->minrad;
  
  stepAttempts = 0;
  workerLink = NULL;
  double curr_error = 0, mr_error=0, newpoca = 0, newmr=0;
	
  /* normalizeVects(dVdt, inState->totalVerts); */

  /* We found that this was too aggressive as a step policy and tended
     to overrun the correction stepper. It does give very fast performance,
     though. */
    
  /* We're going to have to call octrope every time we go through this 
     loop in order to compute the level of error we have so far. 

     We will use the globals "gOctmem" and "gOctmem_size" for memory. */
    
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

       I'm probably overthinking this. --Jason. */
  
    octrope(workerLink,&newrop,&newthi,&newlen,&newmr,&newpoca,
	    0,0,NULL,0,NULL,0,0,NULL,0,NULL,
	    gOctmem,gOctmem_size,gLambda);

    inState->octrope_calls++;
    
    curr_error = (newpoca < 2*inState->tube_radius) ? max(lastDCSD-newpoca,0) : 0;
    mr_error = (newmr < gLambda*inState->tube_radius) ? max(lastMR-newmr,0) : 0;

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
      inState->avgDvdtMag*inState->avgDvdtMag > kMinStepSize) 
    /* last this keeps us from zeroing in dVdt is really small */
    {
      
    inState->stepSize = inState->avgDvdtMag*inState->avgDvdtMag;
    
    }

  // also shouldn't be > 10% of edgelength
  if( inState->stepSize > inState->length/inState->totalVerts*.1 )
    inState->stepSize = inState->length/inState->totalVerts*.1;
  
  free(dVdt);
  
  return inLink;
}

/******************************************************************/
/*                 Building the rigidity matrix                   */
/******************************************************************/

int rownum(search_state *inState, plCurve *inLink, int cmp, int vert, int coord)

     /* Computes the row number in the rigidity matrix corresponding
	to the given data. Handles wraparound correctly, and does
	error checking. We should use this to get stuff into rows of
	vectors and the matrix A exclusively, rather than try to
	recode the computation at various places in the program. */

{
  char errmsg[1024];
  int rnum;

  /* Basic sanity check before we begin. */

  if (vert < -1 || vert > inLink->cp[cmp].nv) {

    sprintf(errmsg,
	    "ridgerunner: rownum reference to vert %d of %d vertex cmp %d"
	    "             is illegal. \n",vert,inLink->cp[cmp].nv,cmp);
    FatalError(errmsg, __FILE__ , __LINE__ );
    
  }

  /* Handle wraparound manually */

  if (vert == -1) { 

    if (inLink->cp[cmp].open) {

      sprintf(errmsg,"ridgerunner: rownum reference to vertex -1 of open component.\n");
      FatalError(errmsg, __FILE__ , __LINE__ );

    } else {

      vert = inLink->cp[cmp].nv-1;

    }

  }

  if (vert == inLink->cp[cmp].nv) {

    if (inLink->cp[cmp].open) {

      sprintf(errmsg,"ridgerunner: rownum reference to vertex nv of open component.\n");
      FatalError(errmsg, __FILE__ , __LINE__ );

    } else {

      vert = 0;

    }

  }

  /* Now compute the row number. Recall that compOffset counts the
     number of _vectors_ which the 0-th vertex of this component is
     offset from the start, so we must multiply by three to get the
     number of _rows_. */

  rnum = 3*inState->compOffsets[cmp] + 3*vert + coord;
  return rnum;

}

static void
placeContactStruts ( taucs_ccs_matrix* A, plCurve* inLink, 
		     octrope_strut* strutSet, int strutCount,  
		     search_state* inState )

     /* Part of the buildrigiditymatrix call. We assume when we're
        building the matrix that the strut columns are placed FIRST,
        and the minrad columns are placed after this. */
     
{
  int sItr;
  
  if( strutCount == 0 )
    return;
  
  // in constructing the ridigity matrix, we will need the struts as viewed 
  // as force vectors on the edges, so we create normalized vectors for each strut
  // here
  plc_vector* strutDirections = (plc_vector*)calloc(strutCount, sizeof(plc_vector));
  fatalifnull_(strutDirections);

  normalizeStruts( strutDirections, strutSet, inLink, strutCount );
  
  for( sItr=0; sItr<strutCount; sItr++ ) {

    int		entry;
				
    // entry is the offset in A which begin this strut's influce
    // it corresponds to the x influence on the lead_vert[0]th vertex
    // after this line, entry+1 is y, +2: z.
    entry = rownum(inState,inLink,
		   strutSet[sItr].component[0],strutSet[sItr].lead_vert[0],0);

    // the strut information includes the position from the strut.lead_vert
    // so we assign "1-position[0]" of the force to the lead vert and "position[0]"
    // of the force to lead_vert+1
  
    // our column is 12*sItr from the start, and we need to set our rowInds
    A->values.d[12*sItr+0] = (1-strutSet[sItr].position[0]) * strutDirections[sItr].c[0];
    A->values.d[12*sItr+1] = (1-strutSet[sItr].position[0]) * strutDirections[sItr].c[1];
    A->values.d[12*sItr+2] = (1-strutSet[sItr].position[0]) * strutDirections[sItr].c[2];
			
    A->rowind[12*sItr+0] = entry + 0;
    A->rowind[12*sItr+1] = entry + 1;
    A->rowind[12*sItr+2] = entry + 2;
			
    // now for the next vertex, receiving "position[0]" of the force, 

    entry = rownum(inState,inLink,
		   strutSet[sItr].component[0],strutSet[sItr].lead_vert[0]+1,0);
      
    A->values.d[12*sItr+3] = (strutSet[sItr].position[0]) * strutDirections[sItr].c[0];
    A->values.d[12*sItr+4] = (strutSet[sItr].position[0]) * strutDirections[sItr].c[1];
    A->values.d[12*sItr+5] = (strutSet[sItr].position[0]) * strutDirections[sItr].c[2];
    
    A->rowind[12*sItr+3] = entry + 0;
    A->rowind[12*sItr+4] = entry + 1;
    A->rowind[12*sItr+5] = entry + 2;
		
    // we do the same thing at the opposite end of the strut, except now the 
    // force is negated
    entry = rownum(inState,inLink,
		   strutSet[sItr].component[1],strutSet[sItr].lead_vert[1],0);
					
    A->values.d[12*sItr+6] = (1-strutSet[sItr].position[1]) * -strutDirections[sItr].c[0];
    A->values.d[12*sItr+7] = (1-strutSet[sItr].position[1]) * -strutDirections[sItr].c[1];
    A->values.d[12*sItr+8] = (1-strutSet[sItr].position[1]) * -strutDirections[sItr].c[2];
    
    A->rowind[12*sItr+6] = entry + 0;
    A->rowind[12*sItr+7] = entry + 1;
    A->rowind[12*sItr+8] = entry + 2;

    entry = rownum(inState,inLink,
		   strutSet[sItr].component[1],strutSet[sItr].lead_vert[1]+1,0);
    
    A->values.d[12*sItr+9] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[0];
    A->values.d[12*sItr+10] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[1];
    A->values.d[12*sItr+11] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[2];
    
    A->rowind[12*sItr+9]  = entry + 0;
    A->rowind[12*sItr+10] = entry + 1;
    A->rowind[12*sItr+11] = entry + 2;
    
  }
	
  free(strutDirections);
}

static void
placeMinradStruts( taucs_ccs_matrix* rigidityA, plCurve* inLink, 
		   octrope_mrloc* minradStruts, 
		   int minradLocs, int startColumn, search_state* inState)

     /* Using the information provided by octrope in the octrope_mrloc array,
	and assuming that rigidityA already contains "contactstruts" 12 entry
	columns, loop through the mrloc locations, computing the gradient of 
	minrad at each vertex and adding the appropriate 9-entry column to the 
        matrix rigidityA at each step. 

	We assert that rigidityA contains at least contactstruts + minradlocs
	columns to begin with, and that rigidityA->colptr is set correctly for
        contactStruts 12-row columns and minradLocs 9-row columns. This _does_
        happen in buildRigidityMatrix. */

{
  int mItr;
  char errmsg[1024];
  
  for( mItr=0; mItr<minradLocs; mItr++ ) {

    /* For each vertex at minimum minrad radius... compute the gradient of mr */

    plc_vector B, A, cross, As, Bs, Cs;
    double bmag, amag, norm;
    double angle;
    double kappa, prevLen, thisLen;
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
    bool ok = true;
    angle = plc_angle(prevSide,thisSide,&ok);

    if (!ok) {

      sprintf(errmsg,
	      "ridgerunner: Couldn't compute turning angle at vertex %d of cmp %d\n"
	      "             of polyline. \n",
	      vItr,cItr);
      FatalError(errmsg, __FILE__ , __LINE__ );

    }

    /* We now check that angle is high enough for the following 
       stuff to work. */

    double PI = 3.14159265358979323846;

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
    
    plc_M_scale_vect(1/norm,As);
    plc_M_scale_vect(1/norm,Bs);
    plc_M_scale_vect(1/norm,Cs);
    
    /* We have now computed the mr gradient, and it's time to fill in the 
       appropriate locations in rigidityA. */
    
    // temporarily increment the strut's verts based on their component interactions
    // we undo this change at the end of the for loop in case the user
    // wants to keep the strut set

    int colptr;

    colptr = rigidityA->colptr[startColumn+mItr];  
    
    int aVert, bVert, cVert, entry;
    aVert = minradStruts[mItr].vert-1;
    bVert = minradStruts[mItr].vert;
    cVert = minradStruts[mItr].vert+1;

    entry = rownum(inState,inLink,minradStruts[mItr].component,aVert,0);
		
    rigidityA->values.d[colptr+0] = As.c[0];
    rigidityA->values.d[colptr+1] = As.c[1];
    rigidityA->values.d[colptr+2] = As.c[2];
    
    rigidityA->rowind[colptr+0] = entry + 0;
    rigidityA->rowind[colptr+1] = entry + 1;
    rigidityA->rowind[colptr+2] = entry + 2;

    /* We now fill in the appropriate entries for b. */

    entry = rownum(inState,inLink,minradStruts[mItr].component,bVert,0);
    assert(0 <= entry && entry < rigidityA->m-2);

    rigidityA->values.d[colptr+3] = Bs.c[0];
    rigidityA->values.d[colptr+4] = Bs.c[1];
    rigidityA->values.d[colptr+5] = Bs.c[2];	
    
    rigidityA->rowind[colptr+3] = entry + 0;
    rigidityA->rowind[colptr+4] = entry + 1;
    rigidityA->rowind[colptr+5] = entry + 2;

    /* We now fill in the cVert values. */

    entry = rownum(inState,inLink,minradStruts[mItr].component,cVert,0);
    assert(0 <= entry && entry < rigidityA->m-2);
        
    rigidityA->values.d[colptr+6] = Cs.c[0];
    rigidityA->values.d[colptr+7] = Cs.c[1];
    rigidityA->values.d[colptr+8] = Cs.c[2];		
    
    rigidityA->rowind[colptr+6] = entry + 0;
    rigidityA->rowind[colptr+7] = entry + 1;
    rigidityA->rowind[colptr+8] = entry + 2;

  } // for over minrad struts

}


void placeConstraintStruts(taucs_ccs_matrix *rigidityA, plCurve *inLink, 
			   int startColumn, search_state *inState )

     /* Places "constraint" struts on the vertices in inLink which are
	subject to fixed, line, or plane constraints. Note that we use
	plc_build_frame to generate line constraints-- there are some
	stability issues here if the lines are in particularly
	unfortunate directions. */

     /* All of the placeXXXX functions assume that
	rigidityA->colptr[i] is allocated correctly and that rigidityA
	has enough space in the values array to hold all of the
	constraints that we're placing. */

{
  plc_constraint *thisConst;
  int vItr, cstItr,entry;
  char errmsg[1024];
  int i;

  cstItr = startColumn; /* The first column of constraints. */

  for( thisConst = inLink->cst; thisConst != NULL; thisConst = thisConst->next ) {

    for ( vItr = thisConst->vert; vItr < thisConst->vert+thisConst->num_verts; vItr++) {

      if ( thisConst->kind == unconstrained ) {

      } else if ( thisConst->kind == plane ) {

	/* The constraint gradient is the normal vector to the plane. */

	rigidityA->values.d[rigidityA->colptr[cstItr]+0] = thisConst->vect[0].c[0];
	rigidityA->values.d[rigidityA->colptr[cstItr]+1] = thisConst->vect[0].c[1];
	rigidityA->values.d[rigidityA->colptr[cstItr]+2] = thisConst->vect[0].c[2];
	
	entry = rownum(inState,inLink,thisConst->cmp,vItr,0);

	rigidityA->rowind[rigidityA->colptr[cstItr]+0] = entry + 0;
	rigidityA->rowind[rigidityA->colptr[cstItr]+1] = entry + 1;
	rigidityA->rowind[rigidityA->colptr[cstItr]+2] = entry + 2;

	cstItr++;

      } else if ( thisConst->kind == line ) {

	/* We should not have any constraints of this type. */

	sprintf(errmsg,
		"ridgerunner: Illegal constraint of type 'line' found.\n");
	FatalError(errmsg, __FILE__ , __LINE__ );

      } else if (thisConst->kind == fixed) {

	/* Here there are three gradients, one in each coordinate direction. */

	for (i=0;i<3;i++) {

	  /* Set the ith value to 1. We zero everything first so that
	     we can do this inside a loop (it's easier to debug than
	     repetitive code). */
	  
	  rigidityA->values.d[rigidityA->colptr[cstItr]+0] = 0;
	  rigidityA->values.d[rigidityA->colptr[cstItr]+1] = 0;
	  rigidityA->values.d[rigidityA->colptr[cstItr]+2] = 0;

	  rigidityA->values.d[rigidityA->colptr[cstItr]+i] = 1;
	  
	  entry = rownum(inState,inLink,thisConst->cmp,vItr,0);

	  rigidityA->rowind[rigidityA->colptr[cstItr]+0] = entry + 0;
	  rigidityA->rowind[rigidityA->colptr[cstItr]+1] = entry + 1;
	  rigidityA->rowind[rigidityA->colptr[cstItr]+2] = entry + 2;

	  cstItr++;

	}

      } else {  /* We should never get to this point. */

	sprintf(errmsg,
		"ridgerunner: Unknown constraint type %d in placeConstraints.\n",
		thisConst->kind);

	FatalError(errmsg, __FILE__ , __LINE__ );

      }

    }

  }

}

	
	

taucs_ccs_matrix *buildRigidityMatrix(plCurve *inLink,search_state *inState)

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

  int sItr;

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

  if (VERBOSITY >= 10) { logprintf("\tbuildRigidityMatrix..."); }

  strutStorageSize = 6*inState->totalVerts;
  minradStorageSize = inState->totalVerts;

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

  inState->octrope_calls++;

  /* We now need to make sure that we didn't exceed the size of the strut
     and/or minrad buffers. */

  if (strutCount >= strutStorageSize || minradLocs >= minradStorageSize ||
      strutCount < 0 || minradLocs < 0) {

    dumpLink(inLink,inState,dumpname);
    sprintf(errmsg,
	    "ridgerunner: octrope found %d struts and %d mrlocs on\n"
	    "             the %d vertex link (dumped to %s), too close to\n"
	    "             minradStorageSize of %d or strutStorageSize %d.\n",
	    strutCount,minradLocs,inState->totalVerts,dumpname,
	    minradStorageSize,strutStorageSize);
    FatalError(errmsg, __FILE__ , __LINE__ );

  }

  /* We now swap in the strutSet and minRad set into inState. */
  
  inState->lastStepStrutCount = strutCount;
  inState->lastStepMinradStrutCount = minradLocs;
  
  if (inState->lastStepStruts != NULL) { free(inState->lastStepStruts); }
  if (inState->lastStepMRlist != NULL) { free(inState->lastStepMRlist); }

  inState->lastStepStruts = strutSet;
  inState->lastStepMRlist = minradSet;

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
				 	
  constraintCount = plCurve_score_constraints(inLink);
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

  // Keep in mind that a taucs_ccs_matrix has one _more_ colptr entries
  // than columns-- the last holds the number of entries in the values array.

  cleanA->colptr[0] = 0;

  for( sItr=1; sItr<=strutCount; sItr++ ) {
    cleanA->colptr[sItr] = sItr*12; }      
  
  for(sItr=strutCount+1; sItr<=strutCount+minradLocs; sItr++ ) {
    cleanA->colptr[sItr] = cleanA->colptr[sItr-1]+9; }

  for(sItr=strutCount+minradLocs+1; 
      sItr<=strutCount+minradLocs+constraintCount;sItr++) {

    cleanA->colptr[sItr] = cleanA->colptr[sItr-1]+3;

  }
  
  // We now double-check that the pointers here are all legal

  assert(cleanA->colptr[sItr-1] == nnz);

  // We now fill the matrix A. 

  placeContactStruts( cleanA, inLink, strutSet, strutCount, inState );
  placeMinradStruts( cleanA, inLink, minradSet, minradLocs, strutCount, inState );
  placeConstraintStruts( cleanA, inLink, strutCount+minradLocs, inState );

  taucs_enforce_ccs_sort(cleanA);

  // We may want to put some kind of strut-exporting code here if we want 
  // intermediate strut sets. 

  if (VERBOSITY >= 10) { logprintf(" ok\n"); }
  
  return cleanA;

}

void freeRigidityMatrix(taucs_ccs_matrix **A)

     /* Frees memory associated with taucs_ccs_matrix A */
     /* and sets the pointer to NULL to prevent further use. */
     /* Can be called multiple times without danger. */
{

  if (A != NULL) {

    free((*A)->colptr);
    free((*A)->rowind);
    free((*A)->values.d);

    free((*A));
    *A = NULL;

  }
}
    
/************************************************************************/

void
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

    if (!inLink->cp[cItr].open) {  
      
      /* We can only spin _closed_ components. */

      // the first thing to do is grab edge lengths
      plc_vector* sides;
      plc_vector* adjustments;
      
      sides = (plc_vector*)malloc(sizeof(plc_vector)*edges[cItr]);
      adjustments = (plc_vector*)calloc(edges[cItr], sizeof(plc_vector));
      
      for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ ) {
	
	bool ok = true;
	char errmsg[1024],dumpname[1024];

	sides[vItr] = plc_normalize_vect(
					 plc_vect_diff(inLink->cp[cItr].vt[vItr+1],
						       inLink->cp[cItr].vt[vItr]),
					 &ok);

	if (!ok) { 

	  dumpLink(inLink,inState,dumpname);
	  sprintf(errmsg,
		  "ridgerunner: spinforce can't normalize side %d of component %d\n"
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

}

void
specialForce( plc_vector* dlen, plCurve* inLink, search_state* inState )

     /* This stub is left here in case we want to experiment with knot
	tightening with respect to some external body force (like
	gravity or something). */

{ 

}

void
eqForce( plc_vector* dlen, plCurve* inLink, search_state* inState )

{
  int cItr, vItr;
  int *edges = malloc(sizeof(int)*inLink->nc);
  fatalifnull_(edges);
  double *lengths = malloc(sizeof(double)*inLink->nc);
  fatalifnull_(lengths);

  double lenUsed,goalUsed;
  int dlItr;
  double varianceSum = 0;
  double scaleFactor;
  bool ok = true;

  plc_vector eqF;

  plc_edges(inLink,edges);
  plc_arclength(inLink,lengths);

  for(cItr=0,dlItr=0;cItr<inLink->nc;cItr++) {

    double avg_edge_length, this_edge_length;
    
    avg_edge_length = lengths[cItr]/edges[cItr];

    for(lenUsed=0,vItr=0;vItr<inLink->cp[cItr].nv;vItr++,dlItr++) {

      goalUsed = vItr*avg_edge_length;
      scaleFactor = inState->eqMultiplier*(goalUsed - lenUsed);
      
      eqF = plc_scale_vect(scaleFactor,plc_mean_tangent(inLink,cItr,vItr,&ok));

      if (!ok) {

	char dumpname[1024], errmsg[1024];
	
	dumpLink(inLink,inState,dumpname);
	sprintf(errmsg,
		"ridgerunner: eqForce can't compute tangent at vertex %d of comp %d \n"
		"             of link %s.\n",vItr,cItr,dumpname);
	FatalError(errmsg, __FILE__ , __LINE__ );
	
      }

      plc_M_add_vect(dlen[dlItr],eqF);

      /* Now update lenUsed to prepare for the next iteration. */

      this_edge_length = plc_distance(inLink->cp[cItr].vt[vItr],
				      inLink->cp[cItr].vt[vItr+1]);

      lenUsed += this_edge_length;
      varianceSum += (this_edge_length - avg_edge_length)*
	(this_edge_length - avg_edge_length);

    }

  }

  inState->eqVariance = varianceSum/inState->totalVerts;

  free(edges);
  free(lengths);

}

int gFeasibleThreshold = 0;
int gDeferredStrutExport = 0;
int gFoo = 0;

plc_vector 
*resolveForce( plc_vector* dl, plCurve* inLink, search_state* inState)

     /* This function resolves the given force dl over the struts of
	inLink, building a rigidity matrix along the way and calling
	octrope and tsnnls to do their jobs. The resulting force is
	returned in a new buffer.  
     */

{
  
  plc_vector*	     dVdt = malloc(sizeof(plc_vector)*inState->totalVerts);
  taucs_ccs_matrix*  A = NULL;
  int		     vItr, sItr; 
  // loop iterators over components, vertices, and struts
  int		     dlItr;
  char               errmsg[1024],vectdumpname[1024],strutdumpname[1024];

  double* compressions = NULL;
  double* minusDL = NULL;
  double *dVdtflat;

  /* We start with a little error checking. */
	
  assert(dl != NULL && inLink != NULL && inState != NULL);
  fatalifnull_(dVdt);

  if (VERBOSITY >= 10) { logprintf("\tresolveForce...\n"); }
  
  /* Now we get to work. */

  A = buildRigidityMatrix(inLink,inState);

  if (A->n == 0) {  /* There are no constraints yet. 
		       Copy dl to dvdt and quit. */		
    
    for(vItr = 0;vItr < inState->totalVerts;vItr++) { dVdt[vItr] = dl[vItr]; }
    
  } else { 			/* There is some linear algebra to do. */

    /* We now try to cancel as much of the motion as possible with strut
       forces.  This involves finding the best nonnegative partial
       solution to the system
       
       AX = -dl.
       
       Of course, we can't solve this equation in general (unless the
       knot is critical!)  so we settle for the closest partial solution
       in a least-squares sense.  A further description of this function
       is contained in taucs_snnls.c
       
    */
    
    // construct minusDL -- this is a column vector of size totalVerts
    // we must operate using strictly doubles to interoperate with taucs_snnls
    
    minusDL = (double*)malloc((3*inState->totalVerts)*sizeof(double));
    fatalifnull_(minusDL);
    
    for( dlItr=0; dlItr<inState->totalVerts; dlItr++ ) {
      
      minusDL[dlItr*3+0] = -dl[dlItr].c[0];
      minusDL[dlItr*3+1] = -dl[dlItr].c[1];
      minusDL[dlItr*3+2] = -dl[dlItr].c[2];
      
    }
    
    // solve AX = -dl, x is strut compressions. We set inRelErrTolerance to 2
    // to make sure that we don't try to do a final lsqr step inside t_snnls
    // (this is way too slow to do right now). We also set inPrintErrorWarnings
    // to TRUE (1), since we want to know if the linear algebra goes bad. 

    // Note to self: We really should implement some kind of fail-safe on t_snnls,
    // which causes the code to time out after some number of minutes. 

    if (inState->steps == 55) {

      printf("ridgerunner: Step 55 dumper enabled. Kill this debugging"
	     "code at %s : %d.\n", __FILE__ , __LINE__ );
      dumpAxb_sparse( inState, A, NULL, minusDL);

    }

    if (VERBOSITY >= 10) { logprintf("\tCalling t_snnls..."); }
    
    compressions = t_snnls(A, minusDL, &inState->residual, 2, 1);
    
    if (VERBOSITY >= 10) { logprintf("ok\n"); }

    if (compressions == NULL) {  

      /* Here's where we expect to fail. Put recovery code in here when designed. */

      dumpAxb_sparse(inState,A,NULL,minusDL);
      dumpStruts(inLink,inState,strutdumpname);
      dumpLink(inLink,inState,vectdumpname);

      sprintf(errmsg,
	      "ridgerunner: Linear algebra failure.\n"
	      "             Dumped link to %s.\n"
	      "             Dumped matrix to A.dat. \n"
	      "             Dumped minusDL to b.dat.\n"
	      "             Dumped struts to %s.\n",
	      vectdumpname,strutdumpname);

      FatalError(errmsg, __FILE__ , __LINE__ );

    }

    /* We have survived the linear algebra step! We record our victory in inState. */

    for(sItr=0;sItr<inState->lastStepStrutCount;sItr++) {

      inState->lastStepStruts[sItr].compression = compressions[sItr];

    }

    for(sItr=0;sItr<inState->lastStepMinradStrutCount;sItr++) {

      inState->lastStepMRlist[sItr].compression = 
	compressions[sItr + inState->lastStepStrutCount];

    }

    /* Note: what to do with constraint strut "compressions"? They aren't recorded! */

    inState->rcond = taucs_rcond(A);
    inState->tsnnls_evaluations++;

    /* We now build dVdt = dl + cleanA*compressions; */

    dVdtflat = (double *)malloc(inState->totalVerts * 3 * sizeof(double));
    fatalifnull_(dVdtflat);

    ourtaucs_ccs_times_vec(A,compressions,dVdtflat);

    for(dlItr=0; dlItr<inState->totalVerts; dlItr++) {

      dVdt[dlItr].c[0] = dVdtflat[3*dlItr+0];
      dVdt[dlItr].c[1] = dVdtflat[3*dlItr+1];
      dVdt[dlItr].c[2] = dVdtflat[3*dlItr+2];
      
    }

    free(dVdtflat);

    /* We have now computed the compression force A*(compressions)
       pushing on each vertex and stored it in dVdt. (There should
       probably be some way to write this to the vect file directory
       if desired.)

       We now add dl to dVdt to finish our computation

       dVdt = A*compressions + dl.

    */
    
    for(dlItr=0;dlItr<inState->totalVerts;dlItr++) {

      plc_M_add_vect(dVdt[dlItr],dl[dlItr]);

    }

    free(minusDL);
    free(compressions);

  }

  /* At this point, we take a snapshot of the computation every state.snapinterval */

  if (inState->steps % inState->snapinterval == 0) {
    
    snapshot(inLink,dVdt,dl,inState);

  }

  /* We have now generated dVdt. Go ahead and free memory and then return it. */

  taucs_ccs_free(A);

  return dVdt;

}

void
normalizeStruts( plc_vector* strutDirections, octrope_strut* strutSet, 
		 plCurve* inLink, int strutCount )
{
  int sItr;
  bool ok = true;
  char errmsg[1024];

  for( sItr=0; sItr<strutCount; sItr++ ) {

    plc_vector  points[2];
    octrope_strut_ends( inLink, &strutSet[sItr], points );
    
    // the normalized difference of pointOne, pointTwo is the strut force vector
    strutDirections[sItr] = plc_normalize_vect(plc_vect_diff(points[0],points[1]),&ok);
    
    if (!ok) {

      sprintf(errmsg,
	      "ridgerunner: Could not normalize strut %d of %d struts on %d vertex"
	      "             plCurve inLink.\n",
	      sItr,strutCount,plc_num_verts(inLink));
      FatalError(errmsg, __FILE__ , __LINE__ );

    }
    
  }
  
}

void
computeCompressPush( plCurve* inLink, octrope_strut* strutSet,
				octrope_mrloc* minradSet, int strutCount, int minradLocs )
{
	
}
