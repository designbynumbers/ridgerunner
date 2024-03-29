/*
 *  stepper.c
 *  ridgerunner
 *
 *  Created by Michael Piatek on Fri May 28 2004.

Copyright Jason Cantarella.

This file is part of ridgerunner. ridgerunner is free software: you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

ridgerunner is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.  You should have received a copy of the GNU General
Public License along with ridgerunner. If not, see
<https://www.gnu.org/licenses/>.

*/


#include "ridgerunner.h"

/* Function prototypes */

plCurve*  bsearch_step( plCurve* inLink, search_state* inState );
void	  step( plCurve* inLink, double stepSize, plc_vector* dVdt );
void	  firstVariation( plc_vector* inOutDvdt, plCurve* inLink, search_state* inState,
			  octrope_strut** outStruts, int* outStrutsCount, int dlenStep);
void	  computeCompressPush( plCurve* inLink, octrope_strut* strutSet,
			       octrope_mrloc* minradSet, int strutCount, 
			       int minradLocs );
plCurve*  steepest_descent_step( plCurve *inLink, search_state *inState);

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
  double minthickness = 500;
  double nextMovieOutput = 0.0;
  
  gConditionCheck = 0;
  double rop_20_itrs_ago = {DBL_MAX};
  double oldrops[20];
  int rItr;
  plCurve *tempLink;
  
#ifdef HAVE_TIME  
  
  time_t startTime,now;
  startTime = time(NULL);
  now = startTime;

#else

  int startTime,now;
  startTime = 0;
  now = 0;

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
	 (inState->residual > inState->residualThreshold) &&
	 (now - startTime < inState->stopTime);
       stepItr++ ) {
    
    int lastSet;      
    lastSet = inState->lastStepStrutCount;
    assert(lastSet >= 0);

    double maxmin;
    maxmin = plCurve_long_edge(*inLink)/plCurve_short_edge(*inLink);
    inState->lastMaxMin = maxmin;
   
    if( gEqIt && (maxmin > 3.0 || plc_num_verts(*inLink)/inState->ropelength < 1.5)) {
      
      int i;
    
      if (plc_num_verts(*inLink)/inState->ropelength < 1.5) { 

	logprintf("Resolution = %g < 1.5. ",plc_num_verts(*inLink)/inState->ropelength);

	/* Ridgerunner won't work reliably below resolution 2. We need to spline up to more verts. */
      
	double     *length = malloc(sizeof(double)*(*inLink)->nc);
	int        *nv = malloc(sizeof(int)*(*inLink)->nc);
	double      totlen,thi;
	int         j;

	assert(nv != NULL);
	assert(length != NULL);
	
	totlen = plc_arclength((*inLink),length);
	assert(totlen > 0);

	thi = inState->thickness;
	
	for(j=0;j<(*inLink)->nc;j++) { length[j] /= thi; }
	totlen /= thi;
	
	if (thi < 1e-10) { 
	  
	  FatalError("ridgerunner: Thickness of curve is %g, which is too small to respline.\n",__FILE__,__LINE__);
	  
	}
	
	/* Computes an effective resolution */
	
	double res;
	res = 2.0;
	
	int totvt = 0;
	assert(totvt==0); /* This is a debugging variable. We just put in the check to prevent compiler warnings here. */
	
	for(j=0;j<(*inLink)->nc;j++) {nv[j] = ceil(res*length[j]); totvt += nv[j];}  
	
	plc_spline *spline;
	bool ok;
   
	spline = plc_convert_to_spline((*inLink),&ok);
	if (!ok) { FatalError("ridgerunner: Couldn't spline to resolution 2.0.\n",__FILE__,__LINE__); }
	plc_free(*inLink);
	(*inLink) = plc_convert_from_spline(spline,nv);

	free(nv);
	free(length);

	logprintf("After splining, resolution = %g.\n",plc_num_verts(*inLink)/octrope_ropelength(*inLink,NULL,0,gLambda));

	/* Now we've changed the size of the link, so we need to update inState. */

	inState->totalVerts = plc_num_verts(*inLink);
	 
      } else {  /* We can use the ordinary fixlength method to distribute the vertices that we've got */
      
	for(i=0;i<3;i++) {
	
	  tempLink = *inLink;
	  *inLink = octrope_fixlength(tempLink);
	  plc_free(tempLink);
	  
	}

      }
	
      free(inState->newDir); inState->newDir = NULL; /* Change the link, frag the direction. */

      
      if (maxmin > 3.0) {

	logprintf("Max edgelength/min edgelength = %g > 3. After equilateralization, max/min = %g.\n",maxmin,
		  plCurve_long_edge(*inLink)/plCurve_short_edge(*inLink));

      } 

      plc_scale(*inLink,(inState->tube_radius + 0.01)/octrope_thickness(*inLink,NULL,0,gLambda));
      logprintf("Rescaled to thickness %g.\n",octrope_thickness(*inLink,NULL,0,gLambda));
      
    }
    
    /* Decide whether to initiate thickness correction. */

    if (!gSONO) { /* SONO mode is self-correcting, so we don't run the correction stepper at any point. */
    
      if( (inState->shortest < (2*inState->tube_radius*(1-inState->overstepTol))) ||
	  (inState->minrad < gLambda*inState->tube_radius*(1-inState->minradOverstepTol))) {
	
	if (gAnimationStepper) { 
	  
	  /* The correct_thickness Newton's method produces a very small correction which looks good
	     in animation. It is _really_ slow for highres knots, and can fail completely when coupled
	     with the steepest_descent stepper, which is too aggressive for correct_thickness to work. */
	  
	  if (gMangleMode) {

	     plc_scale(*inLink,(inState->tube_radius)/octrope_thickness(*inLink,NULL,0,gLambda));

	  } else if (!correct_thickness(*inLink,inState)) {
	    
	    if (inState->oktoscale) {
	      
	      plc_scale(*inLink,(inState->tube_radius)/octrope_thickness(*inLink,NULL,0,gLambda));
	      logprintf("Correction stepper failed. Rescaled manually.\n");
	      
	    } else {
	      
	      char errmsg[1024],dumpname[1024];
	      
	      dumpLink(*inLink,inState,dumpname);
	      sprintf(errmsg,
		      "Correction stepper failed, and we cannot rescale due to constraints of type \"fixed\".\n"
		      "Dumped link to %s.\n"
		      "Terminating run.\n",
		      dumpname);
	      FatalError(errmsg,__FILE__,__LINE__);
	      
	    }

	  }
	  
	} else {        	  
	  if (inState->oktoscale) { 
	    if (gTryNewton) {
	      if (!correct_thickness(*inLink,inState)) {
		plc_scale(*inLink,(inState->tube_radius)/octrope_thickness(*inLink,NULL,0,gLambda));
	      }
	    } else {
	      plc_scale(*inLink,(inState->tube_radius)/octrope_thickness(*inLink,NULL,0,gLambda));
	    }	    
	  } else { 
	    if (!correct_thickness(*inLink,inState)) { 
	      NonFatalError("Newton error correction failed. Will try to continue run anyway.\n",__FILE__,__LINE__);
	    }
	  }
	  
	  free(inState->newDir); inState->newDir = NULL; 
	  inState->score = stepScore(*inLink,inState,NULL,0); /* We have changed the curve, so reset the score */

	  /* We have changed the curve, so attempt to resymmetrize. Harmless if we have no symmetry group. */

	  plc_symmetrize(*inLink);
	  
	}
	
      }   

    }
  
    /************************************************************************/
    
    if (gAnimationStepper) { *inLink = bsearch_step(*inLink, inState); } 
    else if (gSONO) { *inLink = sono_step(*inLink,inState); }
    else { *inLink = steepest_descent_step(*inLink,inState); }

    correct_constraints(*inLink,inState); // This is lightweight, so just keep us honest.
    plc_symmetrize(*inLink); // Likewise: if we are symmetrizing, this keeps us honest. 

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
   
    /* And we update our list of old ropelength values. */

    for(rItr=19;rItr>0;rItr--) { oldrops[rItr] = oldrops[rItr-1]; }
    oldrops[0] = inState->ropelength;

    if (stepItr > 20) { rop_20_itrs_ago = oldrops[19]; }

    #ifdef HAVE_TIME

    now = time(NULL);

    #endif
    
  } 

  if (inState->newDir != NULL) { free(inState->newDir); inState->newDir = NULL; }

  free(gOctmem);

  /* We have now terminated. The final output files will be written 
     in ridgerunner_main.c. However, we log the reason for termination. */

  logprintf("ridgerunner: run complete. Terminated because \n");

  if (!(stepItr < inState->maxItrs)) {
    
    logprintf("ridgerunner: reached maximum number of steps (%d).\n",
	      stepItr);

  }

  if (!(rop_20_itrs_ago - inState->ropelength > inState->stop20)) {

    logprintf("ridgerunner: change in rop over last 20 iterations %g < stop20 = %g.\n",
	      rop_20_itrs_ago - inState->ropelength, inState->stop20);

  }

  if (!(inState->residual > inState->residualThreshold)) {

    logprintf("ridgerunner: residual %g < residualThreshold %g\n",
	      inState->residual, inState->residualThreshold);

  }

  if (!(now - startTime < inState->stopTime)) {

    logprintf("ridgerunner: elapsed time %d (sec) > stopTime %d (sec).\n",
	      now - startTime, inState->stopTime);

  }

  /* If needed, we now make a final log entry, reflecting the condition when we quit. */

  if ( inState->steps%inState->loginterval != 0 ) {

    update_runtime_logs(inState);
    
  }

}


void update_vect_directory(plCurve * const inLink, search_state *inState)

     /* If it is time for the next vect output, go ahead and write 
	another file to the appropriate directory. */
{
  char tmpfilename[1024];
  FILE *outfile;
  static int frames_written_since_compression = 0;

  if (!gSuppressOutput) {
  
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
  fprintf(state->logfiles[kEffectiveness],"%d %g %g\n",state->steps,
	  state->lastStepPocaEffectiveness,state->lastStepMREffectiveness);
  fprintf(state->logfiles[kScore],"%d %g\n",state->steps,state->score);

#ifdef HAVE_MALLINFO
  struct mallinfo m;
  m = mallinfo();
  fprintf(state->logfiles[kMemused],"%d %g M\n",state->steps,((double)(m.uordblks + m.hblkhd))/(1024.0*1024.0) ); // Total mem alloc'd by malloc
#else
  fprintf(state->logfiles[kMemused],"%d (this system does not have mallinfo)\n",state->steps);
#endif

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
step( plCurve* inLink, double stepSize, plc_vector* dVdt )
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

/* Try calling lsqr on this data. If it doesn't work, throw a warning, log things,
   and return NULL. */

{
  /* if there are no constrained struts, t_snnls won't actually work, so use SOL LSQR */
  lsqr_input   *lsqr_in;
  lsqr_output  *lsqr_out;
  lsqr_work    *lsqr_work;
  lsqr_func    *lsqr_func;
  int bItr;
  double*		result;
  char errmsg[1024],logfilename[4096];
  
  alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, sparseA->m, sparseA->n );
  
  /* we let lsqr() itself handle the 0 values in this structure */
  lsqr_in->num_rows = sparseA->m;
  lsqr_in->num_cols = sparseA->n;
  lsqr_in->damp_val = 0;
  lsqr_in->rel_mat_err = 0; // was kZeroThreshold
  lsqr_in->rel_rhs_err = 0;
  lsqr_in->cond_lim = 1/(10*sqrt(DBL_EPSILON));
  lsqr_in->max_iter = 4*lsqr_in->num_cols;  /* Suggested by lsqr docs */

  /* We only want to store the last output in the log. So we need to close and 
     reopen the lsqr logfile here. */

  if (gLsqrLogging) {

    sprintf(logfilename,"./%s.rr/logfiles/%s.dat",
	    inState->basename,
	    inState->logfilenames[klsqrlog]);

    fclose(inState->logfiles[klsqrlog]);
    inState->logfiles[klsqrlog] = fopen_or_die(logfilename,"w",__FILE__,__LINE__);
    
    lsqr_in->lsqr_fp_out = inState->logfiles[klsqrlog];	

  } else {

    lsqr_in->lsqr_fp_out = NULL;

  }

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
	    "             The matrix A^T has been dumped to A.sparse and \n"
	    "             the vector of desired corrections dumped to b.mat.\n"
	    "\n"
	    "             Check the log file 'lsqroutput' for more output.\n",
	    sparseA->m,sparseA->n,lsqr_in->cond_lim);
    NonFatalError(errmsg, __FILE__ , __LINE__ );
    return NULL;

  }

  if (lsqr_out->term_flag == 7) {

    dumpAxb_sparse(inState,sparseA,NULL,minusDL);
    sprintf(errmsg,
	    "ridgerunner: lsqr failed to converge in correction stepper\n"
	    "             on %d x %d matrix A^T after %ld iterations.\n"
	    "             The matrix A^T has been dumped to A.sparse and \n"
	    "             the vector of desired corrections dumped to b.mat.\n"
	    "\n"
	    "             Check the log file 'lsqroutput' for more output.\n" ,
	    sparseA->m,sparseA->n,lsqr_in->max_iter);
    NonFatalError(errmsg, __FILE__ , __LINE__ );
    return NULL;
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
		       search_state *inState,int *errorSize) 

     /* We use the strutList and mrList passed in the header, ignoring
	the versions stored in inState, which may not apply to this link. */

     // We don't score constraint error here... it's handled elsewhere.

{
  int constraintCount = 0; //plCurve_score_constraints(inLink);
  double* C = (double*)(calloc((strutCount+minradCount+constraintCount),sizeof(double)));
  *errorSize = strutCount+minradCount+constraintCount;

  /* We now create the C (corrections) vector giving the desired adjustment 
     to each strut and minrad length. */
  
  int sItr;
  
  for( sItr=0; sItr < strutCount ; sItr++) {  /* Changed from instate->lastStepStrutCount */
    
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
  
  /* plc_constraint *thisCst; */
/*   int vItr; */
/*   plc_vector errVect; */
/*   char errmsg[1024]; */
  
/*   for( sItr = strutCount + minradCount, */
/* 	 thisCst = inLink->cst; */
/*        thisCst != NULL; */
/*        thisCst = thisCst->next ) { */
    
/*     for (vItr = thisCst->vert; vItr < thisCst->vert+thisCst->num_verts; vItr++) { */
      
/*       if (thisCst->kind == unconstrained) { */
	
/*       } else if (thisCst->kind == fixed) { */
	
/* 	errVect = plc_vect_diff(inLink->cp[thisCst->cmp].vt[vItr],thisCst->vect[0]); */
	
/* 	C[sItr++] = -errVect.c[0]; */
/* 	C[sItr++] = -errVect.c[1]; */
/* 	C[sItr++] = -errVect.c[2]; */
	
/*       } else if (thisCst->kind == plane) { */
	
/* 	C[sItr++] = -(plc_dot_prod(inLink->cp[thisCst->cmp].vt[vItr],thisCst->vect[0]) -  */
/* 		      plc_dot_prod(thisCst->vect[1],thisCst->vect[0])); */
	
/*       } else if (thisCst->kind == line) { */

/* 	sprintf(errmsg, */
/* 		"ridgerunner: line constraints are not implemented," */
/* 		"              and should have been detected before now."); */
/* 	FatalError(errmsg, __FILE__ , __LINE__ ); */
	
/*       } else { */
	
/* 	sprintf(errmsg, */
/* 		"ridgerunner: Unknown constraint type %d in inLink.\n", */
/* 		thisCst->kind); */
/* 	FatalError(errmsg, __FILE__ , __LINE__ ); */
	
/*       } */
      
/*     } */
    
/*   } */

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

double plc_l2norm(plc_vector *V,int n)
/* Computes the l^2 norm of an array of vectors of size n */
{
  double sum=0;
  int i;

  for(i=0;i<n;i++) { sum += plc_M_dot(V[i],V[i]); }
  return sqrt(sum);
}


void correct_constraints(plCurve *inLink,search_state *inState) 

     /* Make sure that the curve inLink obeys constraints, if any. */

{

  if (inLink->cst == NULL) { return; } // We can save some time if there are no constraints.

  plc_fix_cst(inLink);

  /* We have now modified the link. We have several responsibilities to discharge. */
  /* Basically, we have to update inState to match the new length/strutSet/etc.... */

  int      strutCount;
  int      minradCount;
  
  octrope_strut *strutSet = NULL;
  octrope_mrloc *minradSet = NULL;
  
  int strutStorageSize = 6*inState->totalVerts;
  int minradStorageSize = 2*inState->totalVerts;

  strutSet = (octrope_strut *)(malloc(sizeof(octrope_strut)*strutStorageSize));
  minradSet = (octrope_mrloc *)(malloc(sizeof(octrope_mrloc)*minradStorageSize));
  
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
	    "correct_constraints: octrope found %d struts and %d mrlocs on\n"
	    "                     the %d vertex link (dumped to %s), too close to\n"
	    "                     minradStorageSize of %d or strutStorageSize %d.\n",
	    strutCount,minradCount,inState->totalVerts,dumpname,
	    minradStorageSize,strutStorageSize);
    FatalError(errmsg, __FILE__ , __LINE__ );
    
  }
  
  /* Now we come up with new buffers for the strutSet and mrLocSet and copy over */
  
  if (inState->lastStepStruts != NULL) { free(inState->lastStepStruts); }
  if (inState->lastStepMRlist != NULL) { free(inState->lastStepMRlist); }
  
  if (strutCount > 0) {
    inState->lastStepStruts = (octrope_strut *)(malloc((strutCount)*sizeof(octrope_strut)));
    fatalifnull_(inState->lastStepStruts);
  } else {
    inState->lastStepStruts = NULL;
  }

  if (minradCount > 0) {
    inState->lastStepMRlist = (octrope_mrloc *)(malloc((minradCount)*sizeof(octrope_mrloc)));
    fatalifnull_(inState->lastStepMRlist);
  } else {
    inState->lastStepMRlist = NULL;
  }

  int sItr;
  
  for(sItr=0;sItr<strutCount;sItr++) { 
    
    inState->lastStepStruts[sItr] = strutSet[sItr];
    
  }
  
  for(sItr=0;sItr<minradCount;sItr++) {
    
    inState->lastStepMRlist[sItr] = minradSet[sItr];
    
  }
  
  free(strutSet);
  free(minradSet);
  
  inState->score = stepScore(inLink,inState,NULL,0); 
  /* Update the score after constraint enforcement */

}

double predict_deltarop(plCurve *inLink,plc_vector *stepDir,double stepSize,double tube_radius,double lambda)

/* Uses the rigidity matrix to predict the change in ropelength if we
   take the step stepSize * stepDir. This is primarily a debugging check. */

{
  taucs_ccs_matrix *sparseA = NULL;
 
  sparseA = buildRigidityMatrix(inLink,tube_radius,lambda,NULL);  /* Note: called "stateless", so doesn't update state. */

  double *step;

  step = (double*)(calloc(3*plc_num_verts(inLink),sizeof(double)));
  fatalifnull_(step);
  
  int i;
  for(i=0;i<plc_num_verts(inLink);i++) { 
   
    step[3*i]     = stepSize*stepDir[i].c[0]; 
    step[(3*i)+1] = stepSize*stepDir[i].c[1];
    step[(3*i)+2] = stepSize*stepDir[i].c[2];
  
  }

  /* Now we take the matrix product. The output buffer has to be allocated in advance. 
     It should be equal in size to the number of columns in sparseA. */

  double *B;
  B = calloc(sparseA->n,sizeof(double));
  fatalifnull_(B);
  taucs_transpose_vec_times_matrix_nosub(step,sparseA,B);

  /* We now analyze the resulting B output. This is a little
     complicated.  This should predict the change in individual strut
     lengths, but remember that the final ropelength change is
     controlled by the length of the shortest strut (and not all of
     these struts are the same length). To determine the actual 
     results, we need to rebuild the strut and minrad sets. */

  double newrop,newthi,newmr,newpoca,newlen;
  octrope_strut struts[7000];
  octrope_mrloc mrlocs[7000];
  int nmr,nstrut;
  
  octrope(inLink,&newrop,&newthi,&newlen,&newmr,&newpoca,
	  lambda*tube_radius,0,mrlocs,7000,&nmr,
	  2*tube_radius,0,struts,7000,&nstrut,
	  NULL,0,lambda);
  
  double predictedThickness=100000,poca = 100000;
  int ptloc = -1,pocaloc = -1;
  int sitr,mritr;

  assert(ptloc == -1);
  assert(pocaloc == -1);

  for(i=0,sitr = 0;sitr < nstrut;i++,sitr++) {

    if ((struts[sitr].length + B[i])/2.0 < predictedThickness) {

      predictedThickness = (struts[sitr].length + B[i])/2.0;
      ptloc = i;

    }

    if (struts[sitr].length < poca) {

      poca = struts[sitr].length;
      pocaloc = sitr;

    }

  }

  for(mritr=0;mritr < nmr;mritr++,i++) {

    if ((mrlocs[mritr].mr + B[i])/lambda < predictedThickness) {

      predictedThickness = (mrlocs[mritr].mr + B[i])/lambda;
      ptloc = i;

    }
  
  }

  double dThidt;
  dThidt = predictedThickness - newthi;

  /* The prediction is that 

   (d/dt) Len/Thi = ((d/dt) Len * Thi - Len * (d/dt) Thi)/(Thi * Thi).

   we have computed (d/dt) Thi, which is expected to be leastDelta. But 
   we need the other numbers to make a prediction. To figure out the change in length, 
   we need the vector field dLen. Here we go:

  */

  plc_vector *dLen;

  dLen = (plc_vector *)(calloc(plc_num_verts(inLink), sizeof(plc_vector)));
  /* Note: this must be calloc since the xxxForce procedures add to given buffer. */
  fatalifnull_(dLen);
  dlenForce(dLen,inLink,NULL);

  double dlendt = 0;

  for(i=0;i<plc_num_verts(inLink);i++) {

    dlendt += -stepSize * plc_dot_prod(stepDir[i],dLen[i]);

  }

  free(dLen);

  /*

   As this is debugging code, we just do it the easy way. 

  */

  double Len, Thi;

  Len = octrope_curvelength(inLink);
  Thi = octrope_thickness(inLink,NULL,0,gLambda);

  double prediction;

  /*   (d/dt) Len/Thi = ((d/dt) Len * Thi - Len * (d/dt) Thi)/(Thi * Thi). */

  prediction = ((dlendt * Thi) - (Len * dThidt))/(Thi * Thi);

  /* Now we clean up after ourselves */

  free(B);
  free(step);
  taucs_ccs_free(sparseA);

  return prediction;

}

int correct_thickness(plCurve *inLink,search_state *inState) 

     /* Newton's method correction algorithm for thickness. */
     
     /* Attempts to fix things so that the shortest strut and 
	worst minrad vertex are within "greenZone" of obeying 
	the constraints. 

	May take most of runtime if the number of verts is large,
	since we use the (slow) Stanford LSQR code to compute
	direction for Newton steps. */

     /* If the Newton code fails, will return FALSE. */

     // We take out the constraint handling code here. 
     
{
  double stepSize = 1.0;
  double currentError = 0, newError;
  double greenZone = 
    2*inState->tube_radius*(1 - 0.5*inState->overstepTol);
  double mrgreenZone = 
    gLambda*inState->tube_radius*(1-0.5*inState->minradOverstepTol);
  int Csize = inState->lastStepStrutCount+inState->lastStepMinradStrutCount
    /*+plCurve_score_constraints(inLink)*/;
  int workerCsize;
 
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
  int minradStorageSize = 2*inState->totalVerts;

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
	(inState->minrad < mrgreenZone)) && 
	gCorrectionAttempts < gMaxCorrectionAttempts;
      gCorrectionAttempts++) {

    taucs_ccs_matrix *sparseA = NULL;
    taucs_ccs_matrix *sparseAT = NULL;
    
    /* We start by building a rigidity matrix, and taking its' transpose */

    sparseA = buildRigidityMatrix(inLink,inState->tube_radius,gLambda,inState);  /* Note: strut set in inState updates */
    sparseAT = taucs_ccs_transpose(sparseA);

    C = thicknessError(inLink,
		       inState->lastStepStrutCount,inState->lastStepStruts,
		       inState->lastStepMinradStrutCount,inState->lastStepMRlist,
		       inState,&Csize);

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

    if (ofv == NULL) { 
      
      char dumpname[1024];
      dumpLink(inLink,inState,dumpname);
      logprintf("stanford_lsqr failed.\n File dumped to: %s.\n",dumpname); 
      return FALSE;}

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
    int maxsplits = 18;
    int splits = 0;

    do {
      
      splits++;

      if (VERBOSITY >= 10) { logprintf("Attempting cstep at size %2.8g.\n",stepSize/2); }

      stepSize /= 2;

      if (workerLink != NULL) plc_free(workerLink);
      workerLink = plc_copy(inLink);
      
      step(workerLink,stepSize,ofv_vect);

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
	
	dumpLink(workerLink,inState,dumpname);
	sprintf(errmsg,
		"ridgerunner: octrope found %d struts and %d mrlocs on\n"
		"             the %d vertex link (dumped to %s), too close to\n"
		"             minradStorageSize of %d or strutStorageSize %d.\n",
		strutCount,minradCount,inState->totalVerts,dumpname,
		minradStorageSize,strutStorageSize);
	FatalError(errmsg, __FILE__ , __LINE__ );
	
      }

      /* Now we can compute the new error. */

      workerC = thicknessError(workerLink,          /* This is new memory, allocated by thickErr */
			       strutCount,strutSet,
			       minradCount,minradSet,
			       inState,&workerCsize);

      newError = l2norm(workerC,workerCsize);

      free(workerC); /* We free it here to keep things tidy. */
      workerC = NULL;
      workerCsize = 0;

    } while (newError > (1 - alpha*stepSize)*currentError && splits < maxsplits);

    if (splits == maxsplits) {

      char dumpname[1024];
      char errmsg[1024];
          
      dumpLink(inLink,inState,dumpname);   
      sprintf(errmsg,
	      "correct_thickness: Couldn't get error to meet the sufficient\n"
	      "decrease condition in %d splits.\n"
	      "\n"
	      "newError = %g, alpha = %g, stepSize = %g, currentError = %g.\n"
	      "so\n"
	      "newError vs (1 - alpha*stepSize)*currentError is\n"
	      "%g vs %g.\n"
	      "\nDumping this link to %s\n",
	      maxsplits,newError,alpha,stepSize,currentError,
	      newError,(1 - alpha*stepSize)*currentError,
	      dumpname);
      
      NonFatalError(errmsg, __FILE__ , __LINE__ );

      /* The original link has not yet been modified, so it's ok to just return. */
      
      free(C);
      free(ofv);
      free(ofv_vect);
      free(strutSet);
      free(minradSet);
      
      taucs_ccs_free(sparseA);
      taucs_ccs_free(sparseAT);

      plc_free(workerLink);

      return 0;
 
    }

    /* This "sufficient decrease" condition comes from Kelley,
       "Solving nonlinear equations with Newton's method", p. 13. 

       We have now chosen a stepSize, stepped, and accepted the
       results.  Before we return to the top of the main for-loop and
       compute a new Newton direction, we must move the temp data in
       workerLink (and friends) to current data in inLink and
       inState. */

    step(inLink,stepSize,ofv_vect); 
    /* Trust me-- this was better than copying. */

    /* We now update the lastStepStrut and MRstrut information to  */
    /* reflect the fact that we've changed the link. Rather than calling */
    /* octrope again, we depend on the fact that the 

       strutCount, strutSet

       and

       minradCount, minradSet

       buffers were computed from the same step on the same initial 
       configuration, so we can copy them into the lastStep buffers in 
       inState. 

       Of course, the first step is to free the buffers if (and only if)
       they were actually allocated in the first place. */

    if (inState->lastStepStruts != NULL) {

      free(inState->lastStepStruts);
      inState->lastStepStruts = NULL;
      inState->lastStepStrutCount = 0;

    }

    if (inState->lastStepMRlist != NULL) {

      free(inState->lastStepMRlist);
      inState->lastStepMRlist = NULL;
      inState->lastStepMinradStrutCount = 0;

    }
    
    /* We don't know whether either buffer actually contains entries. So we 
       won't try to allocate and work with the pointers unless we are allocating
       a nonzero buffer. (In principle, if malloc is POSIX-compliant, we should
       be ok even if called with a size of 0, but for portability reasons, we 
       don't want to depend on this.) */

    if (strutCount > 0) { 

      inState->lastStepStruts = 
	(octrope_strut *)(malloc_or_die(strutCount*sizeof(octrope_strut),
					__FILE__, __LINE__ ));

      int sItr;

      for(sItr=0;sItr<strutCount;sItr++) { 
	
	inState->lastStepStruts[sItr] = strutSet[sItr];

      }

    } 

    inState->lastStepStrutCount = strutCount;

    /* Now handle MRlist similarly */

    if (minradCount > 0) {

      inState->lastStepMRlist 
	= (octrope_mrloc *)(malloc_or_die(minradCount*sizeof(octrope_mrloc),
					  __FILE__ , __LINE__ ));
      int mrItr;
      
      for(mrItr=0;mrItr<minradCount;mrItr++) {
	
	inState->lastStepMRlist[mrItr] = minradSet[mrItr];
	
      }

    }

    inState->lastStepMinradStrutCount = minradCount;

    /* Now update the rest of the geometric information in inState. */
    
    inState->length = len;
    inState->thickness = thi;
    inState->ropelength = rop;
    inState->shortest = shortest;
    inState->minrad = minrad;

    /* Finally, we free memory used in this step. */

    free(C);
    free(ofv);
    free(ofv_vect);
    if (workerC != NULL) { free(workerC); workerC = NULL; }

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
	    "             minrad value for current link is %g \n"
	    "\n"
	    "             last ofvnorm was %g \n"
	    "             stepsize is %g \n"
	    "\n"
	    "             Can set --MaxCorrectionAttempts=XXX from command\n"
	    "             line and try again or debug solver problem.\n",
	    gCorrectionAttempts,gMaxCorrectionAttempts,dumpname,
	    inState->shortest, inState->minrad, inState->ofvNorm,
	    stepSize);

    NonFatalError(errmsg, __FILE__ , __LINE__ );
    plc_free(workerLink);
    free(strutSet);
    free(minradSet);

    return 0;

  }

  inState->last_cstep_attempts = gCorrectionAttempts; 
  plc_free(workerLink);

  free(strutSet);
  free(minradSet);

  return 1; /* We have succeeded. */

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

double trialStep( plCurve *inLink, search_state *inState, 
		  plc_vector *stepDir, double stepSize, 
		  double *mr_error,double *curr_error) 
{
  /* Procedure computes the ropelength involved in stepping in 
     the direction stepDir for size stepSize. The procedure
     is nondestructive to inLink, and frees all the memory
     that it allocates in the computation. */

  /* This is essentially debugging code which is used in 
     computing a measure of step effectiveness. */

  /* This procedure should not change inState, aside from 
     incrementing octrope_calls. */

  plCurve *workerLink;
  double newrop, newthi, newlen, newmr, newpoca;

  workerLink = plc_copy(inLink); 

  step(workerLink, stepSize, stepDir);

  octrope(workerLink,&newrop,&newthi,&newlen,&newmr,&newpoca,
	  0,0,NULL,0,NULL,0,0,NULL,0,NULL,
	  gOctmem,gOctmem_size,gLambda);

  inState->octrope_calls++;
    
  *curr_error = (newpoca < 2*inState->tube_radius) ? max(inState->shortest-newpoca,0) : 0;
  *mr_error = (newmr < gLambda*inState->tube_radius) ? max(inState->minrad-newmr,0) : 0;

  plc_free(workerLink);

  return newrop;
  
}

plc_vector *inputForce( plCurve *inLink, double tube_radius, double eqMultiplier, double lambda,search_state *inState )

/* Builds the entire "input" force, including the dLen term, any timewarp force, eqForce, constraintForce or specialForce. */

{
  plc_vector *dLen;

  dLen = (plc_vector *)(calloc(plc_num_verts(inLink), sizeof(plc_vector)));
  /* Note: this must be calloc since the xxxForce procedures add to given buffer. */
  fatalifnull_(dLen);
    
  dlenForce(dLen,inLink,inState);
  if (!gNoTimeWarp) { accelerate_free_vertices(dLen,inLink,tube_radius); }  // This feature tries to straighten free sections faster
  eqForce(dLen,inLink,eqMultiplier,inState);
  specialForce(dLen,inLink,inState);
  if (gSpinForce) { spinForce(dLen,inLink,inState); }
  constraintForce(dLen,inLink); // Make sure that dLen doesn't try to violate constraints.

  return dLen;

}

plc_vector *stepDirection( plCurve *inLink, double tube_radius, double eqMultiplier, double lambda,search_state *inState) 

/* Compute the step direction by building the rigidity matrix and resolving against constraints */
/* Can be called "stateless" by passing NULL for inState. */

{
  // create initial vector field for which we want to move along in this step
  plc_vector  *dVdt;
  plc_vector  *dLen; 
  
  dLen = inputForce(inLink,tube_radius,eqMultiplier,lambda,inState);
  dVdt = resolveForce(dLen,inLink,tube_radius,lambda,inState); 
  free(dLen);
    
  return dVdt;
}

double stepScore(plCurve *inLink, search_state *inState, plc_vector *stepDir, double stepSize)
 
/* Scoring a step is actually a little difficult, since the score is not strictly in terms of ropelength.
   When there are no struts, the score should simply be the length. If there _are_ struts, then the score
   is ropelength. This isn't hard to do, but we encapsulate it in order to make the stepper code readable.
   
   Note that you can score the current configuration by passing NULL as the stepDir and 0 as the 
   stepsize. 
*/
{ 
  plCurve *workerLink;
  double newrop, newthi, newlen, newmr, newpoca;
  double score;
 
  workerLink = plc_copy(inLink); 
  
  if (stepSize != 0 && stepDir != NULL) {

    step(workerLink, stepSize, stepDir);

  }
    
  octrope(workerLink,&newrop,&newthi,&newlen,&newmr,&newpoca,
	  0,0,NULL,0,NULL,0,0,NULL,0,NULL,
	  gOctmem,gOctmem_size,gLambda);

  inState->octrope_calls++;
  
  if (newmr < gLambda * inState->tube_radius || newpoca < 2*inState->tube_radius) {

    score=newrop;

  } else {

    score=2*newlen;

  }

  //if (!gNoTimeWarp) {

  //  score += 2*(9*strut_free_length(workerLink,inState));
    /* The portion of the curve which is strut free is weighted 10x more heavily. 
       We've already counted it once. But we need to add the extra weight. */

  //}

  plc_free(workerLink);

  return score;

}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10

double brent_step(double ax,double bx,double cx, 
		  plCurve *inLink, search_state *inState, plc_vector *dVdt, 
		  double tol,double *xmin)

/* Given a bracketing triple of step sizes ax, bx, and cx, finds the stepsize of minimum score (between ax and cx) 
   to within tol and returns it in xmin, returning the lowest score itself as the return value of the function. */

{
  int iter;
  double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;
  double d=0;
  
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;

  fw=fv=fx=stepScore(inLink,inState,dVdt,x);

  for (iter=1;iter<=ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
	d=p/q;
	u=x+d;
	if (u-a < tol2 || b-u < tol2)
	  d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));

    fu=stepScore(inLink,inState,dVdt,u);

    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
	SHFT(fv,fw,fx,fu)
	} else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      } else if (fu <= fv || v == x || v == w) {
	v=u;
	fv=fu;
      }
    }
  }
  char ErrMsg[1024];
  sprintf(ErrMsg,"brent_step failed to converge on a minimum after 100 iterations.\n");
  FatalError(ErrMsg,__FILE__,__LINE__);

  *xmin=x;
  return fx;
}

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	    plCurve *inLink,search_state *inState,plc_vector *dVdt)
{
  double ulim,u,r,q,fu,dum;
  
  *fa=stepScore(inLink,inState,dVdt,*ax);
  *fb=stepScore(inLink,inState,dVdt,*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
      }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=stepScore(inLink,inState,dVdt,*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=stepScore(inLink,inState,dVdt,u);
      if (fu < *fc) {
	*ax=(*bx);
	*bx=u;
	*fa=(*fb);
	*fb=fu;
	return;
      } else if (fu > *fb) {
	*cx=u;
	*fc=fu;
	return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=stepScore(inLink,inState,dVdt,u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=stepScore(inLink,inState,dVdt,u);
      if (fu < *fc) {
	SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
	  SHFT(*fb,*fc,fu,stepScore(inLink,inState,dVdt,u))
	  }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=stepScore(inLink,inState,dVdt,u);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=stepScore(inLink,inState,dVdt,u);
    }
    SHFT(*ax,*bx,*cx,u)
      SHFT(*fa,*fb,*fc,fu)
      }
}

plCurve *doStep( plCurve *inLink, plc_vector *dVdt, double InitialStepSize, search_state *inState)

{
  /* We attempt to actually take the step, starting at InitialStepSize and reducing the size until 
     we reach a point where the linear algebra works on the next step. We set the inState variable
     stepDir if we succeed. We do not return if we fail. */

  plCurve *workerLink;
  plc_vector *newGrad = NULL;
  double StepSize;

  for(StepSize = 2*InitialStepSize;newGrad==NULL;) {

    StepSize /= 2.0;  /* We change StepSize here so that it's correct when we exit the loop */

    workerLink = plc_copy(inLink);
    step(workerLink, StepSize, dVdt);
    newGrad = stepDirection( workerLink, inState->tube_radius, inState->eqMultiplier, gLambda, inState);
    
    if (StepSize < 1e-11*inState->maxStepSize) { 

      /* We couldn't get the linear algebra to work at any step size (note that we scale this relative
         to the maximum allowed step size). */

      char errmsg[1024],dumpname[1024];

      dumpLink(inLink,inState,dumpname);
      sprintf(errmsg,
	      "Linear algebra fails even at StepSize = %g.\n"
	      "Dumping link to %s.\n"
	      ,StepSize,dumpname);
      FatalError(errmsg,__FILE__,__LINE__);
      
    }

  }

  /* We've managed to actually take the step! */

  inState->stepSize = StepSize;

  double newrop,newthi,newlen,newmr,newpoca;

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
  
  /* Now we update "inState" to reflect the changes in inLink. */

   int dvdtItr;
  // grab average dvdt now that we have finished mungering it
  inState->avgDvdtMag=0;
  for( dvdtItr=0; dvdtItr<inState->totalVerts; dvdtItr++ ) {
    inState->avgDvdtMag += plc_M_norm(dVdt[dvdtItr]);
  }
  inState->avgDvdtMag /= inState->totalVerts;

  if (gConjugateGradient) {

    /* We now update the conjugate gradient direction. 
       
       We have the last (h) step direction (inState->newDir), the last gradient (inState->lastGrad),
       and the current negative gradient (newGrad). We try to assemble these into the new search direction
       using the Polack-Ribiere formula:
       
       h = newGrad - (lastGrad - newgrad).(newgrad)/(lastGrad.lastGrad) newDir;
       
    */
    
    plc_vector *newDir;
    newDir = calloc(inState->totalVerts,sizeof(plc_vector));
    
    double num,denom;
    int dotItr;
    
    for(num=0,denom=0,dotItr = 0;dotItr < inState->totalVerts;dotItr++) {
      
      num += (plc_dot_prod(plc_vect_diff(inState->lastGrad[dotItr],newGrad[dotItr]),newGrad[dotItr]));
      denom += plc_dot_prod(inState->lastGrad[dotItr],inState->lastGrad[dotItr]);
      
    }
    
    for(dotItr=0;dotItr<inState->totalVerts;dotItr++) {
      
      newDir[dotItr] = plc_vlincomb(1,newGrad[dotItr],-num/denom,inState->newDir[dotItr]);
      
    }

    /* Now here's a serious question: should we resolve newDir against the struts? In practice, it's a disaster. */
    /* In principle, if the strut set is stable, we don't need to. However, if the strut set is unstable, 
       this might make a real difference. The question is nontrivial because it essentially doubles the 
       step time (we've already made a linear algebra call when we got newGrad). */
    
    free(inState->newDir);
    inState->newDir = newDir;

    free(inState->lastGrad);
    inState->lastGrad = newGrad;
    
  } else { /* Use ordinary gradient descent */

    free(inState->newDir);
    inState->newDir = newGrad;

  }

  inState->last_step_attempts = 0; /* This is no longer recorded, as it doesn't really make sense. */
  inState->minrad = newmr;
  inState->shortest = newpoca;
  inState->length = newlen;
  inState->ropelength = newrop; 
  inState->thickness = newthi;
  
  inState->cstep_time += inState->stepSize;
  inState->time += inState->stepSize;

  return workerLink;

}

double brentSearch(plCurve *inLink,search_state *inState,plc_vector *dVdt,double *best_step)
 /* Find a bracketing triple and do our best to find a stepsize. */
 /* Function evaluations are octrope calls, so we try to be efficient. Essentially, this
     is a combination of the Numerical Recipes routines mnbrak and brent. */

 /* The terrifying thing about the steepest descent stepper is that there is NO UPPER BOUND
    on stepsize. This means that this probably won't be good for moviemaking. On the other hand,
    it may do very well in the endgame. */

{
  double ax = 0, bx = 1.0*1e-4 /* fabs(inState->stepSize) */, cx = 2.0*1e-4 /*fabs(inState->stepSize) */;
  double fa, fb, fc;

  mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,inLink,inState,dVdt);
  
  /* We have now found a bracketing triple axv, bxv, cxv with corresponding step scores fav, fbv, fcv. */
  /* Our goal will be to pass them to brent and take this opportunity to zero in on the best possible step. */

  if (fb > fa || fb > fc) {

    char ErrMsg[1024];
    sprintf(ErrMsg,"sd_step: mnbrak returned triple of step sizes (%g,%g,%g) with values (%g,%g,%g) which do not appear to bracket a min.\n",
	    ax,bx,cx,fa,fb,fc);
    NonFatalError(ErrMsg,__FILE__,__LINE__);

  }

  /* We now invoke the brent code. Again, we lift this from a web source. We're trying to converge on the best 
     step size, so relatively high (1e-2) error is acceptable, especially since function evaluations (octrope calls) 
     are quite expensive. */

  double best_score;
  best_score = brent_step(ax,bx,cx,inLink,inState,dVdt,1e-2,best_step);

  return best_score;
}

plCurve* 
sono_step( plCurve *inLink, search_state *inState)

/* We try a version of the stepper which takes a large number of SONO steps in a streamlined
   way, then attempts to compute the usual stuff as an afterthought. */

{
  plc_vector *dVdt;
  int itr;
  plCurve *workerA;
  double rop,thi,len,mr,poca;
  double newrop, newthi, newlen, newmr, newpoca;

  workerA = plc_copy(inLink);

  for (itr=0;itr<1000;itr++) {

    octrope(workerA,&rop,&thi,&len,&mr,&poca,
	    0,0,NULL,0,NULL,0,0,NULL,0,NULL,
	    gOctmem,gOctmem_size,gLambda);

    dVdt = inputForce(workerA,inState->tube_radius,inState->eqMultiplier,gLambda,inState); 
    step(workerA, 0.001, dVdt);    

    octrope(workerA,&newrop,&newthi,&newlen,&newmr,&newpoca,
	    0,0,NULL,0,NULL,0,0,NULL,0,NULL,
	    gOctmem,gOctmem_size,gLambda);
    
    plc_scale(workerA,(inState->tube_radius)/newthi);
    
    octrope(workerA,&rop,&thi,&len,&mr,&poca,
	    0,0,NULL,0,NULL,0,0,NULL,0,NULL,
	    gOctmem,gOctmem_size,gLambda);
    	    
  }
  
  plc_free(inLink);  
  return workerA;

}


plCurve* 
steepest_descent_step( plCurve *inLink, search_state *inState)

/* We will have to make this a little more modular in order to compute step quality */
/* and consider the effect of moving simply in the dLen/scale direction instead of  */
/* resolving. */

{
  double best_step,best_score;
  plc_vector *dVdt;
  int itr;

  if (inState->newDir != NULL) {  /* We may have cached a step direction. If we did, use it! */
    dVdt = inState->newDir;
  } else {
    dVdt = stepDirection(inLink,inState->tube_radius,inState->eqMultiplier,gLambda, inState);

    if (dVdt == NULL) {

      char errMsg[1024],dumpname[1024];

      dumpLink(inLink,inState,dumpname);
      sprintf(errMsg,
	      "Could not compute a step direction at the start of a steepest descent step.\n"
	      "Either this link cannot be run with -c or correction stepping has resulted\n"
	      "in a degenerate configuration.\n"
	      "Dumped link to %s.\n"
	      "Terminating run.\n",dumpname);
      FatalError(errMsg, __FILE__, __LINE__ );

    }

    /* We are starting (or restarting) the cached-direction process. We reset the CG apparatus as well,
       copying both initial directions to the negative gradient. */

    inState->newDir = calloc(inState->totalVerts,sizeof(plc_vector));
    if (inState->lastGrad != NULL) {free(inState->lastGrad);}
    inState->lastGrad = calloc(inState->totalVerts,sizeof(plc_vector));

    for(itr=0;itr < inState->totalVerts;itr++) { inState->lastGrad[itr] = dVdt[itr]; inState->newDir[itr] = dVdt[itr];}

  }

  double best_dVdt_score, best_dVdt_step;
  best_dVdt_score = brentSearch(inLink,inState,dVdt,&best_dVdt_step);

  /* Now we are going to loop to figure out the best step in the current direction. */
 
  if (best_dVdt_step < inState->minStep) { 

    if (inState->residual > 1) { /* We only try alternate step directions if residual is low */

      best_step = inState->minStep;
      best_score = stepScore(inLink, inState, dVdt, inState->minStep);

    } else {

      /* We suspect "unhealthy" step behavior. We are going to try a dVdt step instead. */

      plc_vector *dLen;
      dLen = inputForce(inLink,inState->tube_radius,inState->eqMultiplier,gLambda,inState);
      
      double best_dLen_score,best_dLen_step;
      best_dLen_score = brentSearch(inLink,inState,dLen,&best_dLen_step);
      
      if (best_dLen_step < inState->minStep) {
	
	best_dLen_step = inState->minStep;
	best_dLen_score = stepScore(inLink,inState,dLen,inState->minStep);
	
      }
      
      if (best_dLen_score < best_dVdt_score) {
	
	int i,nv;
	for(i=0,nv=plc_num_verts(inLink);i<nv;i++) { dVdt[i] = dLen[i]; }
	best_step = best_dLen_step;
	best_score = best_dLen_score;

	logprintf("(dLen step)");
	
      } else {
	
	best_step = inState->minStep; 
	best_score = stepScore(inLink, inState, dVdt, inState->minStep);

      }

      free(dLen);

    }

  } else { /* Normal step was of an ok size. We check to make sure that it's not too big. */

    if (best_dVdt_step > inState->maxStepSize) { 

      best_step = inState->maxStepSize;
      best_score = stepScore(inLink,inState,dVdt,inState->maxStepSize);

    } else {

      best_step = best_dVdt_step;
      best_score = best_dVdt_score;

    }

  }
  
  plCurve *workerLink;
  workerLink = doStep(inLink,dVdt,best_step,inState);
  inState->score = best_score;
  
  plc_free(inLink);  
  return workerLink;

}


plCurve*
bsearch_step( plCurve* inLink, search_state* inState )

/* This is the default stepper design. It attempts to increase stepsize as long the induced
   error per step is less than the fixed constants ERROR_BOUND and MR_ERROR_BOUND. This has 
   some interesting consequences; among them, ropelength can actually go UP in an accepted
   step. The output is a new link. We are expected to destroy inLink during the function. 
*/

{	
  int stepAttempts = 0;
  int dvdtItr;
  double  lastDCSD, lastMR;
  plCurve* workerLink;
	
  double ERROR_BOUND = 1e-5;
  double MR_ERROR_BOUND = 1e-5;

  /* Be conservative if we can't afford to crash the error-correction stepper. */
  if (!inState->oktoscale) { ERROR_BOUND *= 0.1; MR_ERROR_BOUND *= 0.1; }  
  
  // create initial vector field for which we want to move along in this step
  plc_vector  *dVdt;
  plc_vector  *dLen; 

  double maxDvDtnorm = 0;
  int i;

  if (VERBOSITY >= 5) {  /* Verbose or higher */ 

    logprintf("Bsearch step %d.\n",inState->steps);

  }

  dLen = inputForce( inLink, inState->tube_radius, inState->eqMultiplier,gLambda,inState );
  dVdt = resolveForce(dLen,inLink,inState->tube_radius,gLambda,inState); 
  /* Built from the bones of firstVariation. */ 

  /* We compute the maximum size of a vector in dVdt to help with debugging code. */

  for(i=0;i<inState->totalVerts;i++) 
    { if (plc_norm(dVdt[i]) > maxDvDtnorm) { maxDvDtnorm = plc_norm(dVdt[i]); } }

  double dVdtNorm; 
  dVdtNorm = plc_l2norm(dVdt,inState->totalVerts);

  double dLenNorm;
  dLenNorm = plc_l2norm(dLen,inState->totalVerts);
  
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
  double stepTaken;

  do {			
    
    stepAttempts++;
		
    // move along dVdt
	
    if( workerLink != NULL ) { plc_free(workerLink); }
    workerLink = plc_copy(inLink); 

    step(workerLink, inState->stepSize, dVdt);
    stepTaken = inState->stepSize; // Record the step actually taken this time.

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

    /* Now in principle, we just moved the vertices of the curve by no
       more than maxDvDtnorm*inState->stepSize.  So the error in
       thickness should be no more than twice that. If that's wrong,
       then one of workerLink and inLink has a poca that's being
       computed incorrectly. */

    if (curr_error > 5*maxDvDtnorm*inState->stepSize) {

      char dname[1024];

      logprintf("bsearch_step: Detected a possible error (in octrope?) at step %d.\n",inState->steps);
      logprintf("              5 * maxDvDtnorm (%g) * inState->stepSize (%g) = %g < curr_error (%g).\n",
		maxDvDtnorm, inState->stepSize, 5*maxDvDtnorm*inState->stepSize, curr_error);

      logprintf("              Dumping pair of links involved to filenames ",inState->steps);
      dumpNamedLink(inLink,inState,"a",dname);
      logprintf("%s and ",dname);
      dumpNamedLink(workerLink,inState,"b",dname);
      logprintf("%s.\n",dname);   

    }

    if (gLambda > 0 && mr_error > 0) {
    
      /* We can also compute the expected change in minrad from a
	 motion of this size. This is a little more complicated. We
	 summarize the computation as follows. If E is the minimum
	 edge length at vertex i,
	 
	 minrad = E/2 tan(theta/2), or 2 minrad = E/tan(theta/2), or tan(theta/2) = E/2minrad, or 
	 
      */
      
      double mintheta,edgelen;
      
      edgelen = newlen/plc_num_edges(workerLink);
      mintheta = 2*atan(edgelen/(2*gLambda*inState->tube_radius));
      
      /*
	
	where theta is the turning angle at i. Now the derivative of
	this function is
	
	d/d(theta) minrad = E/(2 cos theta - 2).
	
	This is an increasing and negative function, at least on 0
	.. pi. So how much can minrad drop on a given angle step?
	Well, at the minimum turning angle of minturn, the derivative
	of minrad should be pretty close to
	
      */
      
      double dminrad;
      dminrad = edgelen/(2*cos(mintheta) - 2);
      
      /* Now what's the maximum change in angle that we expect from
	 this much motion of the edges?  A computation reveals that
	 each edge should change direction by at most */
      
      double maxDeltaAngle;
      maxDeltaAngle = asin(maxDvDtnorm*inState->stepSize/(edgelen/2.0));
      
      /* So the total change in minrad should be no more than twice
	 this times the dminrad bound...  we allow up to 300% of this
	 error just to be on the safe side. We get */
      
      if (mr_error > fabs(3*dminrad*2*maxDeltaAngle)) {
	
	char dname[1024];
	
	logprintf("bsearch_step: Detected a possible error in minrad-controlled stepping at step %d.\n",inState->steps);
	logprintf("              fabs(3*dminrad (%g) * 2 * maxDeltaAngle (%g)) = %g < mr_error (%g).\n",
		  dminrad, maxDeltaAngle, fabs(3*dminrad*2*maxDeltaAngle), mr_error);
	logprintf("              Dumping pair of links involved to filenames ",inState->steps);
	dumpNamedLink(inLink,inState,"a",dname);
	logprintf("%s and ",dname);
	dumpNamedLink(workerLink,inState,"b",dname);
	logprintf("%s.\n",dname);   

      }
    
    }

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

  /* At this point, we have settled on an acceptable step size. */
  /* We now check the effectiveness of our step direction scheme by
       comparing the error resulting from this step to the error that
       would have resulting if we had just stepped in the dLen
       direction. We have to choose the comparison step size
       carefully, since dLen and dVdt have different norms. */
  
  double comp_curr_error, comp_mr_error;

  trialStep(inLink,inState,dLen,(dVdtNorm/dLenNorm)*stepTaken,
	    &comp_mr_error,&comp_curr_error);
  
  /* Now we compute effectiveness in terms of mr_error and curr_err. */

  if (comp_curr_error < 1e-15) { 
    inState->lastStepPocaEffectiveness = (curr_error > 0) ? -1 : 1;
  } else {
   inState->lastStepPocaEffectiveness = 1 - curr_error/comp_curr_error;
  }

  if (comp_mr_error < 1e-15) { 
    inState->lastStepMREffectiveness = (mr_error > 0) ? -1 : 1;
  } else {
    inState->lastStepMREffectiveness = 1 - mr_error/comp_mr_error;
  }

  free(dLen); /* We are now done with dLen, so let's get rid of it before we forget to do so. */
  
  if (curr_error > 100*ERROR_BOUND) {

    char dumpName[1024];

    dumpLink(workerLink,inState,dumpName);

    logprintf("bsearch_step: Couldn't reduce error (%g) to < 100 * ERROR_BOUND (%g), even at \n"
	      "              stepsize %g. Dumping workerLink to %s.\n",
	      curr_error,ERROR_BOUND,inState->stepSize,dumpName);

    dumpDvdt(dVdt,inLink,inState);
    exit(1);

  }

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
 
  /* Experimental Version of RR, July 31, 2009. We're trying to run stuff
     without this dvdt bound. This may or may not work. */
 
  /* IF we have minrad struts, then we need to make sure that we
     should make sure stepsize isn't > avgDvdt^2 as that's the minrad
     control bound. However, we don't want to do this if we don't
     currently have minrad struts because it makes everything very
     slow. We don't just check minrad, because stiffness will impact
     things as well. */

 /*  if( inState->stepSize > inState->avgDvdtMag*inState->avgDvdtMag &&  */
/*       inState->avgDvdtMag*inState->avgDvdtMag > kMinStepSize &&  */
/*       inState->lastStepMinradStrutCount > 0)  */
/*     /\* last this keeps us from zeroing if dVdt is really small *\/ */
/*     { */
      
/*       inState->stepSize = inState->avgDvdtMag*inState->avgDvdtMag; */
    
/*     } */

  // the final motion of any vertex shouldn't be > 10% of edgelength

  if( inState->stepSize*inState->avgDvdtMag > (inState->length/inState->totalVerts)*(0.1) ) {
    inState->stepSize = ((inState->length/inState->totalVerts)*.01)/(inState->avgDvdtMag);
  }

  free(dVdt);
  
  return inLink;
}

/******************************************************************/
/*                 Building the rigidity matrix                   */
/******************************************************************/

void cvc(plCurve *inLink,int rownum,int *cmp,int *vert,int *coord)

/* Computes the component, vector and coordinate corresponding to a given row
   in the rigidity matrix, using the componentOffsets in inState. cvc should be 
   the inverse of rownum. */

{
  char errmsg[1024];
  
  /* Basic sanity check before we begin. */

  if (rownum < 0 || rownum >= plc_num_verts(inLink)*3) {

    sprintf(errmsg,
	    "ridgerunner: reference to row %d of rigidity matrix for polyline"
	    "             with %d vertices (and %d rows in rmatrix) is illegal.\n",
	    rownum,plc_num_verts(inLink),plc_num_verts(inLink)*3);
    FatalError(errmsg, __FILE__ , __LINE__ );
    
  }

  /* We now convert. The idea is that we start subtracting entire components until 
     we go negative, then use that point to figure out which component we're in.
     We know that we _will_ go negative because rownum < plc_num_verts(inLink)*3. */

  int cp=0;
  for(cp=0;rownum >= 0;rownum -= 3*inLink->cp[cp].nv,cp++);
  *cmp = cp-1; // On the last trip through the loop, we incremented cp.
  rownum += 3*inLink->cp[*cmp].nv; /* Now we have just the remainder. */
  
  *vert = (int)(floor(rownum/3.0));
  *coord = rownum % 3;

  /* /\* The following is legacy code which used the compOffsets array. *\/ */

/*   /\* We now convert. Remember that the compOffsets are in terms of a list of VECTORS, */
/*      so they should be multiplied by 3 when we are offsetting in the list of doubles */
/*      referred to by rownum. *\/ */

/*   int cp = 0; */
/*   for(cp=0;rownum >= 3*inState->compOffsets[cp] && cp < inLink->nc;cp++); // Search for the right component. */
/*   *cmp = cp-1; */

/*   *vert = (int)(floor((rownum - 3*inState->compOffsets[*cmp])/(3.0)));             // rownum - 3*compOffsets is the coordinate number in this cmp */
/*   *coord = rownum - (3*inState->compOffsets[*cmp] + 3*(*vert));                    // could also be rownum - 3*compOffsets % 3. */

}
  

int rownum(plCurve *inLink, int cmp, int vert, int coord)

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

  /* Now compute the row number. */

  rnum = 0;
  int i;
  for(i=0;i<cmp;i++) { rnum += 3*inLink->cp[i].nv; }
  rnum += 3*vert + coord;
  
  /* The old version of this used compOffset. Here's the legacy code:

     Recall that compOffset counts the
     number of _vectors_ which the 0-th vertex of this component is
     offset from the start, so we must multiply by three to get the
     number of _rows_. */

  /* rnum = 3*inState->compOffsets[cmp] + 3*vert + coord; */

  return rnum;

}

void placeContactStruts ( taucs_ccs_matrix* A, plCurve* inLink, 
			  octrope_strut* strutSet, int strutCount)

     /* Part of the buildrigiditymatrix call. We assume when we're
        building the matrix that the strut columns are placed FIRST,
        and the minrad columns are placed after this. */
     
{
  int sItr;
  
  if( strutCount == 0 )
    return;
  
  // in constructing the rigidity matrix, we will need the struts as viewed 
  // as force vectors on the edges, so we create normalized vectors for each strut
  // here

  /* We follow Proposition 3.3 in the ridgerunner paper, remembering that  */
  /* ridgerunner, for historical reasons, uses the diameter form of ropelength */

  for( sItr=0; sItr<strutCount; sItr++ ) {

    int		entry;
    plc_vector  points[2],pmq;
    double      alpha,beta;
    double      strutlen;

    octrope_strut_ends( inLink, &strutSet[sItr], points );
    pmq = plc_scale_vect(0.5,plc_vect_diff(points[0],points[1]));
    alpha = 1 - strutSet[sItr].position[0];
    beta  = 1 - strutSet[sItr].position[1];
    strutlen = plc_distance(points[0],points[1]);

    // Here we compute (1/strutlen) (p - q). 
				
    // entry is the offset in A which begin this strut's influce
    // it corresponds to the x influence on the lead_vert[0]th vertex
    // after this line, entry+1 is y, +2: z.
    entry = rownum(inLink,
		   strutSet[sItr].component[0],strutSet[sItr].lead_vert[0],0);

    // the strut information includes the position from the strut.lead_vert
    // so we assign "1-position[0]" of the force to the lead vert and "position[0]"
    // of the force to lead_vert+1
  
    // our column is 12*sItr from the start, and we need to set our rowInds
    A->values.d[12*sItr+0] = alpha * pmq.c[0] / strutlen;
    A->values.d[12*sItr+1] = alpha * pmq.c[1] / strutlen;
    A->values.d[12*sItr+2] = alpha * pmq.c[2] / strutlen;
			
    A->rowind[12*sItr+0] = entry + 0;
    A->rowind[12*sItr+1] = entry + 1;
    A->rowind[12*sItr+2] = entry + 2;
			
    // now for the next vertex, receiving "position[0]" of the force, 

    entry = rownum(inLink,
		   strutSet[sItr].component[0],strutSet[sItr].lead_vert[0]+1,0);
      
    A->values.d[12*sItr+3] = (1 - alpha) * pmq.c[0] / strutlen;
    A->values.d[12*sItr+4] = (1 - alpha) * pmq.c[1] / strutlen;
    A->values.d[12*sItr+5] = (1 - alpha) * pmq.c[2] / strutlen;
    
    A->rowind[12*sItr+3] = entry + 0;
    A->rowind[12*sItr+4] = entry + 1;
    A->rowind[12*sItr+5] = entry + 2;
		
    // we do the same thing at the opposite end of the strut, except now the 
    // force is negated
    entry = rownum(inLink,
		   strutSet[sItr].component[1],strutSet[sItr].lead_vert[1],0);
					
    A->values.d[12*sItr+6] = beta * (-pmq.c[0]) / strutlen;
    A->values.d[12*sItr+7] = beta * (-pmq.c[1]) / strutlen; 
    A->values.d[12*sItr+8] = beta * (-pmq.c[2]) / strutlen;
    
    A->rowind[12*sItr+6] = entry + 0;
    A->rowind[12*sItr+7] = entry + 1;
    A->rowind[12*sItr+8] = entry + 2;

    entry = rownum(inLink,
		   strutSet[sItr].component[1],strutSet[sItr].lead_vert[1]+1,0);
    
    A->values.d[12*sItr+9]  = (1 - beta) * (-pmq.c[0]) / strutlen;
    A->values.d[12*sItr+10] = (1 - beta) * (-pmq.c[1]) / strutlen;
    A->values.d[12*sItr+11] = (1 - beta) * (-pmq.c[2]) / strutlen;
    
    A->rowind[12*sItr+9]  = entry + 0;
    A->rowind[12*sItr+10] = entry + 1;
    A->rowind[12*sItr+11] = entry + 2;
    
  }
	
}

void compute_minrad_gradient(plCurve *inLink,octrope_mrloc mrloc,plc_vector *AAs,plc_vector *BBs,plc_vector *CCs)

/* Returns the gradient of minrad corresponding to mrloc. This gradient has three vectors which 
   describe the motion of the vertex before the target (As), the target (Bs), and the vertex after the
   target (Cs). */

{
  plc_vector B, A, cross, As, Bs, Cs;
  double bmag, amag;
  double angle;
  double kappa, prevLen, thisLen;
  plc_vector  prevSide, thisSide, N, fancyL, fancyM, fancyN;
 		
  prevSide = plc_vect_diff(inLink->cp[mrloc.component].vt[mrloc.vert],inLink->cp[mrloc.component].vt[mrloc.vert-1]);
  thisSide = plc_vect_diff(inLink->cp[mrloc.component].vt[mrloc.vert+1],inLink->cp[mrloc.component].vt[mrloc.vert]);
  
  /* dot = plc_M_dot(prevSide, thisSide); */
  
  prevLen = plc_M_norm(prevSide);
  thisLen = plc_M_norm(thisSide);

  assert(prevLen >= 0.0);  /* These are debugging variables */
  assert(thisLen >= 0.0);
  
  // B = b-v = thisSide. 
  
  B = thisSide;
  A = plc_scale_vect(-1,prevSide);
  
  bmag = plc_M_norm(B);
  amag = plc_M_norm(A);
  
  /* value = dot/(prevLen*thisLen); */
  bool ok = true;
  angle = plc_angle(prevSide,thisSide,&ok);
  
  if (!ok) {
    
    char errmsg[1024];
    sprintf(errmsg,
	    "ridgerunner: Couldn't compute turning angle at vertex %d of cmp %d\n"
	    "             of polyline. \n",
	    mrloc.vert,mrloc.component);
    FatalError(errmsg, __FILE__ , __LINE__ );
    
  }
  
  /* We now check that angle is high enough for the following 
     stuff to work. */
  
  double PI = 3.14159265358979323846;
  
  if (angle < 1e-12 || angle > PI - 1e-12) {
    
    char errmsg[1024];
    sprintf(errmsg,
	    "ridgerunner: Can't compute minrad gradient when edges are\n"
	    "             almost colinear. Angle between edges is %g.\n",
	    angle);
    
    FatalError(errmsg, __FILE__ , __LINE__ );
    
  }
  
  if( mrloc.svert > mrloc.vert ) { // We are computing the gradient of "plus minrad". 
    
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
    
  } else { // We are computing the gradient of "minus minrad". 
    
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

  *AAs = As; *BBs = Bs; *CCs = Cs;
  
}
  
static void
placeMinradStruts( taucs_ccs_matrix* rigidityA, plCurve* inLink, 
		   octrope_mrloc* minradStruts, 
		   int minradLocs, int startColumn )

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
    
  for( mItr=0; mItr<minradLocs; mItr++ ) {

    plc_vector As, Bs, Cs;
    double norm;

    compute_minrad_gradient(inLink,minradStruts[mItr],&As,&Bs,&Cs); 
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

    entry = rownum(inLink,minradStruts[mItr].component,aVert,0);
		
    rigidityA->values.d[colptr+0] = As.c[0];
    rigidityA->values.d[colptr+1] = As.c[1];
    rigidityA->values.d[colptr+2] = As.c[2];
    
    rigidityA->rowind[colptr+0] = entry + 0;
    rigidityA->rowind[colptr+1] = entry + 1;
    rigidityA->rowind[colptr+2] = entry + 2;

    /* We now fill in the appropriate entries for b. */

    entry = rownum(inLink,minradStruts[mItr].component,bVert,0);
    assert(0 <= entry && entry < rigidityA->m-2);

    rigidityA->values.d[colptr+3] = Bs.c[0];
    rigidityA->values.d[colptr+4] = Bs.c[1];
    rigidityA->values.d[colptr+5] = Bs.c[2];	
    
    rigidityA->rowind[colptr+3] = entry + 0;
    rigidityA->rowind[colptr+4] = entry + 1;
    rigidityA->rowind[colptr+5] = entry + 2;

    /* We now fill in the cVert values. */

    entry = rownum(inLink,minradStruts[mItr].component,cVert,0);
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
	
	entry = rownum(inLink,thisConst->cmp,vItr,0);

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
	  
	  entry = rownum(inLink,thisConst->cmp,vItr,0);

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
	

taucs_ccs_matrix *buildRigidityMatrix(plCurve *inLink,double tube_radius,double lambda,search_state *inState)

     /* Procedure calls octrope and uses the results to allocate and
	build a rigidity matrix for the curve, calling
	placeMinradStruts, placeContactStruts, and
	placeConstraintStruts. The procedure can be called "stateless" by 
        passing NULL for search_state. */

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
     struts by setting minradStorageSize to twice the number of verts. We
     believe (with no proof) that choosing strutStorageSize = 6*the
     number of edges in the link we will similarly overestimate the #
     of struts. We will double-check the results anyway. */

  /* We note that this process could be marginally more efficient if
     we did this malloc "once and for all" at the start of
     computation. */

  if (VERBOSITY >= 10) { logprintf("\tbuildRigidityMatrix..."); }

  strutStorageSize = 6*plc_num_verts(inLink);
  minradStorageSize = 2*plc_num_verts(inLink);

  strutSet = malloc(sizeof(octrope_strut)*strutStorageSize);
  minradSet = malloc(sizeof(octrope_mrloc)*minradStorageSize);

  fatalifnull_(strutSet);
  fatalifnull_(minradSet); 

  double rop,thi,len,mr,shortest;

  octrope(inLink, 
		
	  &rop,
	  &thi,		
		
	  &len,		
	  &mr,
	  &shortest,
		
	  // minrad struts
	  lambda*tube_radius,  /* Cutoff reports all mrlocs less than this */
	  0,                             /* This value will be ignored. */
	  minradSet, 
	  minradStorageSize,
	  &minradLocs,
		
	  // strut info
	  2*tube_radius, /* Cutoff reports struts shorter than this. */
	  0,                      /* This epsilon value also ignored. */
	  strutSet,
	  strutStorageSize,
	  &strutCount,
		
	  NULL, 0,

	  lambda);               /* The global "stiffness" parameter. */  

  /* We now need to correct if stiffness == 0, because in this case, we want to throw 
     out the minradlocs buffer whatever its contents are.  */

  if (lambda == 0) {

    free(minradSet);
    minradSet = NULL;
    minradLocs = 0;

  }

  if (inState != NULL) {  /* Update the trailing cache in inState. */

    inState->ropelength = rop;
    inState->thickness = thi;
    inState->length = len;
    inState->minrad = mr;
    inState->shortest = shortest;
    inState->octrope_calls++;

  /* We now need to make sure that we didn't exceed the size of the strut
     and/or minrad buffers. */

  }

  if (strutCount >= strutStorageSize || minradLocs >= minradStorageSize ||
      strutCount < 0 || minradLocs < 0) {

    if (inState != NULL) { dumpLink(inLink,inState,dumpname); }
    sprintf(errmsg,
	    "ridgerunner: octrope found %d struts and %d mrlocs on\n"
	    "             the %d vertex link (dumped to %s), too close to\n"
	    "             minradStorageSize of %d or strutStorageSize %d.\n",
	    strutCount,minradLocs,inState->totalVerts,dumpname,
	    minradStorageSize,strutStorageSize);
    FatalError(errmsg, __FILE__ , __LINE__ );

  }

  if (inState != NULL) {

    /* We now swap in the strutSet and minRad set into inState. */
    
    inState->lastStepStrutCount = strutCount;
    inState->lastStepMinradStrutCount = minradLocs;
    
    if (inState->lastStepStruts != NULL) { free(inState->lastStepStruts); }
    if (inState->lastStepMRlist != NULL) { free(inState->lastStepMRlist); }
    
    if (strutCount > 0) {
      inState->lastStepStruts = strutSet;
    } else {
      free(strutSet);
      inState->lastStepStruts = NULL;
    }

    if (minradLocs > 0 || minradSet == NULL) {
      inState->lastStepMRlist = minradSet;
    } else {
      free(minradSet);
      inState->lastStepMRlist = NULL;
    }
      
  }
    
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

  // At the moment, we handle constraints outside the linear algebraic 
  // formalism to simplify things. Thus, we don't add constraints ever.
				 	
  constraintCount = 0; // plCurve_score_constraints(inLink);
  int nnz = 12*strutCount + 9*minradLocs + 3*constraintCount; // we KNOW this
  taucs_ccs_matrix *cleanA;
  
  cleanA = taucs_ccs_new(3*plc_num_verts(inLink),strutCount+minradLocs+constraintCount,nnz); /* (rows,columns,nnz) */
  
/* cleanA = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix)); */
/*   fatalifnull_(cleanA); */

/*   cleanA->n = strutCount+minradLocs+constraintCount; */
/*   cleanA->m = 3*inState->totalVerts; */
/*   cleanA->flags = TAUCS_DOUBLE; */
  
/*   cleanA->colptr = (int*)malloc(sizeof(int)*(cleanA->n+1)); */
/*   cleanA->rowind = (int*)malloc(sizeof(int)*nnz); */
/*   cleanA->values.d = (double*)malloc(sizeof(taucs_double)*nnz); */

/*   fatalifnull_(cleanA->colptr); */
/*   fatalifnull_(cleanA->rowind); */
/*   fatalifnull_(cleanA->values.d); */
  
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

  placeContactStruts( cleanA, inLink, strutSet, strutCount);
  placeMinradStruts( cleanA, inLink, minradSet, minradLocs, strutCount );
  //placeConstraintStruts( cleanA, inLink, strutCount+minradLocs, inState );

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

/* void */
/* spinForce( plc_vector* dlen, plCurve* inLink, search_state* inState ) */

/*      /\* This adds to dlen a tangential force causing the tube to spin */
/* 	as it is tightening. This doesn't seem to have any useful */
/* 	numerical effect. *\/ */

/* { */
/*   int cItr, vItr; */
/*   int *edges = malloc(sizeof(int)*inLink->nc); */

/*   // note to self - this can be made much faster by not being dumb */
  
/*   plc_fix_wrap(inLink); */
/*   plc_edges(inLink,edges); */
  
/*   for( cItr=0; cItr<inLink->nc; cItr++ ) { */

/*     if (!inLink->cp[cItr].open) {   */
      
/*       /\* We can only spin _closed_ components. *\/ */

/*       // the first thing to do is grab edge lengths */
/*       plc_vector* sides; */
/*       plc_vector* adjustments; */
      
/*       sides = (plc_vector*)malloc(sizeof(plc_vector)*edges[cItr]); */
/*       adjustments = (plc_vector*)calloc(edges[cItr], sizeof(plc_vector)); */
      
/*       for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ ) { */
	
/* 	bool ok = true; */
/* 	char errmsg[1024],dumpname[1024]; */

/* 	sides[vItr] = plc_normalize_vect( */
/* 					 plc_vect_diff(inLink->cp[cItr].vt[vItr+1], */
/* 						       inLink->cp[cItr].vt[vItr]), */
/* 					 &ok); */

/* 	if (!ok) {  */

/* 	  dumpLink(inLink,inState,dumpname); */
/* 	  sprintf(errmsg, */
/* 		  "ridgerunner: spinforce can't normalize side %d of component %d\n" */
/* 		  "             of link %s.\n",vItr,cItr,dumpname); */
/* 	  FatalError(errmsg, __FILE__ , __LINE__ ); */

/* 	} */
	
/*       } */

/*       double spinFactor = 0.5; */
/*       // fix zero and compute the tangential change necessary for the rest of the vertices */
/*       for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ ) { */
	
/* 	adjustments[vItr] = plc_scale_vect(spinFactor,sides[vItr]); */

/* 	dlen[((vItr)) + inState->compOffsets[cItr]].c[0] += adjustments[(vItr)].c[0]; */
/* 	dlen[((vItr)) + inState->compOffsets[cItr]].c[1] += adjustments[(vItr)].c[1]; */
/* 	dlen[((vItr)) + inState->compOffsets[cItr]].c[2] += adjustments[(vItr)].c[2]; */
/*       } */
      
/*       free(sides); */
/*       free(adjustments); */
    
/*     } */

/*   } */
  
/*   free(edges); */

/* } */

void
specialForce( plc_vector* dlen, plCurve* inLink, search_state* inState )

     /* This stub is left here in case we want to experiment with knot
	tightening with respect to some external body force (like
	gravity or something). */

{ 
  int dlp;
  int cp,vt;

  if (gMangleMode) {

    if (inState->steps % inState->steps_in_mode == 0) {  

      /* Time to pick a new axis and radius for the torus rotation randomly. */

      plc_vector com;
      double     diameter;

      com = plc_center_of_mass(inLink);
      diameter = plc_pointset_diameter(inLink);

      inState->torusrotate_axis = plc_random_vect();
      inState->torusrotate_center = plc_vect_sum(plc_scale_vect(0.3*diameter*rand()/RAND_MAX,plc_random_vect()),com);
      inState->torusrotate_radius = fabs(0.2*diameter + 0.3*diameter*rand()/RAND_MAX);

      logprintf("----------------------------------------------\n");
      logprintf("\nChanged Mangle settings to mode torusrotate.\n");
      logprintf("\tRotating around axis (%g,%g,%g).\n",plc_M_clist(inState->torusrotate_axis));
      logprintf("\tRotating with center (%g,%g,%g).\n",plc_M_clist(inState->torusrotate_center));
      logprintf("\tRotating with distance to focus %g.\n",inState->torusrotate_radius);

      /* what we need to do now is translate the entire link so that the center is at the origin */

      for(cp=0;cp<inLink->nc;cp++) {

	for(vt=0;vt<inLink->cp[cp].nv;vt++) {

	  plc_M_sub_vect(inLink->cp[cp].vt[vt],inState->torusrotate_center);

	}

      }

      plc_fix_wrap(inLink);

      /* now we rotate to make the axis the z - axis. */

      plc_random_rotate(inLink,inState->torusrotate_axis);

    }

    for(cp=0;cp<inLink->nc;cp++) {

      for(vt=0;vt<inLink->cp[cp].nv;vt++) {

	dlp = dlenPos(inLink,cp,vt);

	/* The vector field is computed to be tangent to the tau = constant tori
	   in toroidal coordinates. Our strategy is to compute the cylindrical 
	   radius of the center ring of the torus containing (x,y,z) and use it
	   to come up with the direction of the field. 

	   See the Mathematica notebook (MangleModePics.nb) in rrpaper/talks
	   for a more complete explanation of the code below.

	*/
	
	plc_vector v;
	double tau,d12,d22,rho,a,torusC;
	
	a = inState->torusrotate_radius;
	v = inLink->cp[cp].vt[vt];

	rho = sqrt(v.c[0]*v.c[0] + v.c[1]*v.c[1]);
	d12 = (rho + a)*(rho + a) + v.c[2]*v.c[2];
	d22 = (rho - a)*(rho - a) + v.c[2]*v.c[2];

	if (d22 > 0.01) { 

	  tau = log(d12/d22); 
	  torusC = a / tanh(tau);
	  
	  
	  dlen[dlp] = plc_vect_sum(dlen[dlp],
				   plc_build_vect(-v.c[2]*v.c[0]/rho,-v.c[2]*v.c[1],(rho - torusC)));
	  
	} /* else we're on the central circle and don't move */
				   
      }
	 
    }

  }

}

void
spinForce ( plc_vector *dlen, plCurve *inLink, search_state *inState) 

/* Adds a small component to the force in the tangent direction to the curve */

{

  int cmp;
  int vt;

  for(cmp=0;cmp < inLink->nc; cmp++) {

    if (!inLink->cp[cmp].open) {

      for(vt=0;vt < inLink->cp[cmp].nv;vt++) {
	
	int dlp;
	bool ok;
	
	dlp = dlenPos(inLink,cmp,vt);
	dlen[dlp] = plc_vect_sum(dlen[dlp],plc_scale_vect(0.01,plc_mean_tangent(inLink,cmp,vt,&ok)));
	
	if (!ok) {
	  
	  char errMsg[1024];
	  char dumpName[1024];
	  dumpLink(inLink,inState,dumpName);
	  
	  sprintf(errMsg,"ridgerunner: Couldn't compute tangent vector to vert %d of component %d of inlink.\n"
		  "             Dumped offending link to %s.\n",cmp,vt,dumpName);
	  
	  FatalError(errMsg,__FILE__,__LINE__);
	  
	}
	
      }
      
    }

  }
    
}			       

int dlenPos(plCurve *inLink,int cmp,int vt) 

     /* Converts a Link, cmp, vert triple to a position on the dlen vector. */

{
  int dlp=0,cmpItr;

  for(cmpItr=0;cmpItr<cmp;cmpItr++) {

    dlp += inLink->cp[cmpItr].nv;

  }

  dlp += vt;

  return dlp;

}

void
constraintForce( plc_vector* dlen, plCurve* inLink )

     /* Make sure that the force obeys any active constraints. */

{ 
  plc_constraint *thisCst;
  int vt,i;

  for(thisCst = inLink->cst;thisCst != NULL;thisCst=thisCst->next) {

    for(vt=thisCst->vert,i=0;i<thisCst->num_verts;vt++,i++) {

      if (thisCst->kind == fixed) {

	dlen[dlenPos(inLink,thisCst->cmp,vt)] = plc_build_vect(0,0,0);

      } else if (thisCst->kind == line) {

	plc_vector T,*dlv;
	bool ok;

	T = plc_normalize_vect(thisCst->vect[0],&ok); assert(ok);
	dlv = &(dlen[dlenPos(inLink,thisCst->cmp,vt)]);

	*dlv = plc_scale_vect(plc_dot_prod(*dlv,T),T);

      } else if (thisCst->kind == plane) {

	plc_vector N,*dlv;
	bool ok;

	N = plc_normalize_vect(thisCst->vect[0],&ok); assert(ok);
	dlv = &(dlen[dlenPos(inLink,thisCst->cmp,vt)]);

	*dlv = plc_vect_diff(*dlv,plc_scale_vect(plc_dot_prod(*dlv,N),N));

      } else {

	char errmsg[1024];

	sprintf(errmsg,"constraintForce: Fatal error. Unknown constraint type detected in inLink.\n");
	FatalError(errmsg,__FILE__,__LINE__);

      }

    }

  }

}

void
eqForce( plc_vector* dlen, plCurve* inLink, double eqMultiplier, search_state* inState )

/* Can be called "stateless" by passing NULL for inState. */

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
      scaleFactor = eqMultiplier*(goalUsed - lenUsed);
      
      eqF = plc_scale_vect(scaleFactor,plc_mean_tangent(inLink,cItr,vItr,&ok));

      if (!ok) {

	char dumpname[1024], errmsg[1024];
	
	if (inState != NULL) { dumpLink(inLink,inState,dumpname); }
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

  if (inState != NULL) {inState->eqVariance = varianceSum/inState->totalVerts;}

  free(edges);
  free(lengths);

}


taucs_ccs_matrix *taucs_ccs_delete_column(taucs_ccs_matrix *A,int col)
/* Constructs the matrix with col deleted. We expect col to be between 0 and A->n-1. */
{
  taucs_ccs_matrix *Ap;

  Ap = calloc(1,sizeof(taucs_ccs_matrix));

  if (col > A->n-1 || col < 0) { 

    char errcode[1024];

    sprintf(errcode,"taucs_ccs_delete_column: can't delete column %d of %d column matrix A.\n",col,A->n);
    FatalError(errcode,__FILE__,__LINE__);

  }

  Ap->n = A->n-1;
  Ap->m = A->m;
  Ap->flags = A->flags;

  Ap->colptr = calloc(Ap->n+1,sizeof(int));
  Ap->rowind = calloc(A->colptr[A->n],sizeof(int));      /* These are actually too big. But that's ok. */
  Ap->values.d = calloc(A->colptr[A->n],sizeof(double)); 

  /* First, we copy the column pointers into the new colptr array. */

  int i,shft,colstart,colstop;
  for(i=0;i<col;i++) { Ap->colptr[i] = A->colptr[i]; }
  colstart = A->colptr[col]; colstop = A->colptr[col+1];
  shft = colstop - colstart;  /* How many entries in rowind are we deleting? */
  for(;i<=Ap->n;i++) { Ap->colptr[i] = A->colptr[i+1] - shft; }

  /* Now we copy the row data and values. When we reach the break point, we skip ahead by shft entries. */

  for(i=0;i<colstart;i++) { Ap->rowind[i] = A->rowind[i]; Ap->values.d[i] = A->values.d[i]; }
  for(;i<Ap->colptr[Ap->n];i++) { Ap->rowind[i] = A->rowind[i+shft]; Ap->values.d[i] = A->values.d[i+shft]; }
  
  return Ap;

}

void vgscan_taucs_ccs_matrix(taucs_ccs_matrix *A)

/* Look at the data in A and make a (trivial) decision based on every item. 
   Intended to catch valgrind "conditional move or jump depends on uninitialized value" errors. */

{
  int i = 1;

  if (A->n > 0) { i += 2; }
  if (A->m > 0) { i += 2; }
  if (A->flags > 0) { i += 2; }
  
  int j;

  for(j=0;j<A->n+1;j++) {

    if (A->colptr[j] > 0) { i += 2; }
    
  }

  for (j=0;j<A->colptr[A->n];j++) {

    if (A->rowind[j] > 0) { i += 2; }
    if (A->values.d[j] > 0) { i += 2; }

  }

  assert(i>=0);  

}


int gFeasibleThreshold = 0;
int gDeferredStrutExport = 0;
int gFoo = 0;

plc_vector 
*resolveForce( plc_vector* dl, plCurve* inLink, double tube_radius, double lambda, search_state* inState)

     /* This function resolves the given force dl over the struts of
	inLink, building a rigidity matrix along the way and calling
	octrope and tsnnls to do their jobs. The resulting force is
	returned in a new buffer.  

	This function can be called "stateless" by passing inState = NULL. 

	If the linear algebra fails, we return NULL.
     */

{
  
  plc_vector*	     dVdt = calloc(plc_num_verts(inLink),sizeof(plc_vector));
  taucs_ccs_matrix*  A = NULL;
  int		     vItr, sItr; 
  // loop iterators over components, vertices, and struts
  int		     dlItr;
  //char               errmsg[1024],vectdumpname[1024],strutdumpname[1024];

  double* compressions = NULL;
  double* minusDL = NULL;
  double *dVdtflat;

  /* We start with a little error checking. */
	
  assert(dl != NULL && inLink != NULL && tube_radius > 0 && lambda >= 0);
  fatalifnull_(dVdt);

  if (VERBOSITY >= 10) { logprintf("\tresolveForce...\n"); }
  
  /* Now we get to work. */

  A = buildRigidityMatrix(inLink,tube_radius,lambda,inState);

  // Debugging code. We scan A to see if that's the thing that valgrind's complaining about. 
  //vgscan_taucs_ccs_matrix(A);

  if (A->n == 0) {  /* There are no constraints yet. 
		       Copy dl to dvdt and quit. */		
    
    for(vItr = 0;vItr < inState->totalVerts;vItr++) { dVdt[vItr] = dl[vItr]; }
    if (inState != NULL) {inState->residual = 1.0;}
    
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
    
    minusDL = (double*)calloc((3*plc_num_verts(inLink)),sizeof(double));
    fatalifnull_(minusDL);

    double minusDLnorm = 0;
    
    for( dlItr=0; dlItr<plc_num_verts(inLink); dlItr++ ) {
      
      minusDL[dlItr*3+0] = -dl[dlItr].c[0];
      minusDL[dlItr*3+1] = -dl[dlItr].c[1];
      minusDL[dlItr*3+2] = -dl[dlItr].c[2];
      
      minusDLnorm += plc_M_dot(dl[dlItr],dl[dlItr]);  // Hitch a ride on the loop to compute norm
      
    }
    
    minusDLnorm = sqrt(minusDLnorm);
    
    // solve AX = -dl, x is strut compressions. We set inRelErrTolerance to 2
    // to make sure that we don't try to do a final lsqr step inside t_snnls
    // (this is way too slow to do right now). We also set inPrintErrorWarnings
    // to TRUE (1), since we want to know if the linear algebra goes bad. 

    // Note to self: We really should implement some kind of fail-safe on t_snnls,
    // which causes the code to time out after some number of minutes. 

    if (VERBOSITY >= 10) { logprintf("\tCalling t_snnls..."); }
 
    double l2ResidualNorm;
    compressions = t_snnls_fallback(A, minusDL, &l2ResidualNorm, 2, 1);

    /* We now dump the matrix if state requires us to. */

    if (inState != NULL) { if (inState->dumpAxb) { dumpAxb_sparse(inState,A,compressions,minusDL); }}
       
    char *terr;

    if (tsnnls_error(&terr)) {  /* If tsnnls throws a code, log it. */
	
      logprintf("%s",terr);
      if (inState != NULL) {dumpAxb_sparse(inState,A,NULL,minusDL);}

    }

    if (compressions == NULL) {   /* We really failed. We need to free the memory that we've allocated, and return. */

      logprintf("resolve_force: Linear algebra failure. Returning control to stepper.\n");
      free(dVdt);
      taucs_ccs_free(A);
      free(minusDL);
      return NULL;

    }

    /* We now (being paranoid) check to make sure that the compressions are actually all positive. */

    double compmin = 1000000;
    int compminloc = -1;
    assert(compminloc == -1);

    for(dlItr=0;dlItr<A->n;dlItr++) { 

      if (compressions[dlItr] < compmin) {

	compmin = compressions[dlItr];
	compminloc = dlItr;

      }

      if (compressions[dlItr] < 0) {

	char errmsg[1024];
	
	sprintf(errmsg,
		"resolve_force: tsnnls returned a negative compression on column %d of matrix A.\n"
		"               This indicates a bug in tsnnls which needs to be corrected.\n"
		"               Contact cantarella@math.uga.edu with dump data (A,x,b).\n"
		"               Dumping rigidity matrix (to A), compressions (to x), \n",dlItr);

	if (inState != NULL) {dumpAxb_sparse(inState,A,compressions,minusDL);}

	FatalError(errmsg,__FILE__,__LINE__);
	
      }

    } 

    /* We have gathered compmin, which was essentially free, but we won't usually use it. */

     /*  /\* We failed the linear algebra step. It is likely that our matrix can be saved by deleting  *\/ */
/*       /\* a column. The question is simple: which one? We can only find out by a search algorithm... *\/ */

/*       taucs_ccs_matrix *Ap; */
/*       int col; */
/*       bool success = false; */
/*       double rcond; */

/*       logprintf("Linear algebra failed. Trying to delete a column to find nonsingular matrix.\n"); */

/*       for(col=0;col<A->n && !success;col++) { */

/* 	Ap = taucs_ccs_delete_column(A,col); */
/* 	rcond = taucs_rcond(Ap); */
	
/* 	if (rcond > 1e-12) {      /\* We might be able to solve it! *\/ */

/* 	  printf("\t col deleted: %d rcond: %g, trying again...\n",col,rcond);  */

/* 	  compressions = t_snnls_fallback(Ap, minusDL, &l2ResidualNorm, 2, 1); /\* Try to solve it *\/ */
/* 	  success = (compressions != NULL); */
/* 	  if (success) { logprintf("Succeeded by deleting column %d. Run should continue.\n",col); } */

/* 	} */

/* 	taucs_ccs_free(Ap); */

/*       } */

/*     } */
	  
/*     if (compressions == NULL) { */
  
/*     dumpAxb_sparse(inState,A,NULL,minusDL); */
/*       dumpStruts(inLink,inState,strutdumpname); */
/*       dumpLink(inLink,inState,vectdumpname); */

/*       sprintf(errmsg, */
/* 	      "ridgerunner: Linear algebra failure.\n" */
/* 	      "             rcond of matrix %g.\n" */
/* 	      "             Dumped link to %s.\n" */
/* 	      "             Dumped matrix to A.mat, A.sparse. \n" */
/* 	      "             Dumped minusDL to b.mat.\n" */
/* 	      "             Dumped struts to %s.\n", */
/* 	      taucs_rcond(A),vectdumpname,strutdumpname); */

/*       FatalError(errmsg, __FILE__ , __LINE__ ); */

    if (VERBOSITY >= 10) { logprintf("ok\n"); }

    if (inState != NULL) { 

      inState->residual = l2ResidualNorm/minusDLnorm;  
      /* We now update the strutfree residual if we are keeping track of this... */

      if (gStrutFreeResidual) {

	inState->strutfreeresidual = strut_free_residual(inLink,inState);

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

      if (gNoRcond) {
	inState->rcond = 10; 
      } else {
	inState->rcond = taucs_rcond(A);
      }
      
      inState->tsnnls_evaluations++;

    }
      
    /* We now build dVdt = dl + cleanA*compressions; */

    dVdtflat = (double *)calloc(3*plc_num_verts(inLink),sizeof(double));
    fatalifnull_(dVdtflat);

    ourtaucs_ccs_times_vec(A,compressions,dVdtflat);

    for(dlItr=0; dlItr<plc_num_verts(inLink); dlItr++) {

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
    
    for(dlItr=0;dlItr<plc_num_verts(inLink);dlItr++) {

      plc_M_add_vect(dVdt[dlItr],dl[dlItr]);

    }

    free(minusDL);
    free(compressions);

  }

  /* At this point, we take a snapshot of the computation every state.snapinterval */

  if (inState != NULL) {

    if (inState->steps % inState->snapinterval == 0) {
      
      snapshot(inLink,dVdt,dl,inState);
      
    }

  }

  /* We have now generated dVdt. Go ahead and free memory and then return it. */

  taucs_ccs_free(A);

  /* If we are running with a symmetry, we should now symmetrize the resolved force. */
  /* This will do nothing if inLink->G == NULL. */

  plc_symmetrize_variation(inLink,dVdt);

  return dVdt;

}

void
normalizeStruts( plc_vector* strutDirections, octrope_strut* strutSet, 
		 plCurve* inLink, int strutCount )

// No longer called.
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

 double strut_free_residual( plCurve *L, search_state *state) 

 /* Returns the residual of the curve which results from strut-free portions of the curve only. */

 {
   int nStrutFree;
   bool *freeFlag;
   plc_vector *dLen = NULL;
   double dLenNorm;

   dLen = (plc_vector *)(calloc(plc_num_verts(L), sizeof(plc_vector)));
  /* Note: this must be calloc since the xxxForce procedures add to given buffer. */
  fatalifnull_(dLen);

  freeFlag = malloc_or_die(sizeof(bool)*plc_num_verts(L), __FILE__ , __LINE__ );
  nStrutFree = strut_free_vertices(L,state->tube_radius,freeFlag);
  assert(nStrutFree >= 0);
  
  dlenForce(dLen,L,state);

  int cmpItr,vItr;
  double l2free = 0;

  for(cmpItr=0;cmpItr < L->nc;cmpItr++) {

    for(vItr=0;vItr < L->cp[cmpItr].nv;vItr++) {

      if (freeFlag[vItr]) {

	l2free += pow(plc_norm(dLen[dlenPos(L,cmpItr,vItr)]),2.0);

      }

    } 

  }

  l2free = sqrt(l2free);
  free(freeFlag);

  /* We now compute the l2norm of dLen itself. */

  dLenNorm = plc_l2norm(dLen,state->totalVerts);
   
  free(dLen);

  return l2free/dLenNorm;

}
