/*
 *  stepper.c
 *  ridgerunner
 *
 *  Created by Michael Piatek on Fri May 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "tsnnls.h"

#include "vector.h"

#include "stepper.h"
#include "dlen.h"
#include "eqedge.h"
#include "settings.h"
#include "../errors.h"

octrope_link*		bsearch_step( octrope_link* inLink, search_state* inState );
void				step( octrope_link* inLink, double stepSize, octrope_vector* dVdt );
void				firstVariation( octrope_vector* inOutDvdt, octrope_link* inLink, search_state* inState,
						octrope_strut** outStruts, int* outStrutsCount, int dlenStep);
						
// lesser utility functions
int					equalStruts( const octrope_strut* s1, const octrope_strut* s2 );
void				normalizeStruts( octrope_vector* strutDirections, 
							octrope_strut* strutSet, octrope_link* inLink, int strutCount );

void		placeVertexBars( double* A, octrope_link* inLink, int contactStruts, int totalBarVerts, int totalBars, search_state* inState );
int displayEveryFrame = 0;

static void
export_pushed_edges( octrope_link* L, search_state* inState, double* pushes, char* fname, int colorParam)
{
	int i;
	double  maxPush = 0;
	
	for( i=0; i<inState->totalVerts; i++ )
	{
		if( maxPush < pushes[i] )
			maxPush = pushes[i];
	}		
	
	if( maxPush == 0 )  
		maxPush = 1e-6;
	
	int j;              /* Counter for the for loops */ 
	int nverts = 0;       /* Total number of vertices of all components */
	int colors = inState->totalVerts;       /* Total number of colors of all components */

	  /* Now we begin work. */
  for(i=0;i<L->nc;i++) {
    nverts += L->cp[i].nv;
//    colors += L->cp[i].cc;
  }
	
	FILE* file = fopen(fname, "w");

  /* We are ready to write the link. */
  fprintf(file,"VECT \n");
  fprintf(file,"%d %d %d \n",L->nc,nverts,colors);
  
  for(i=0;i<L->nc;i++) {
    if (L->cp[i].acyclic) {
      fprintf(file,"%d ",L->cp[i].nv); 
    } else {
      fprintf(file,"%d ",-L->cp[i].nv);
    }
  }
  fprintf(file,"\n");

  for(i=0;i<L->nc;i++) {
    fprintf(file,"%d ", L->cp[i].nv);
  }
  fprintf(file,"\n");

  /* Now we write the vertex data . . . */
  for(i=0;i<L->nc;i++) {
    for(j=0;j<L->cp[i].nv;j++) {
      fprintf(file,"%g %g %g \n", L->cp[i].vt[j].c[0], L->cp[i].vt[j].c[1],
                                  L->cp[i].vt[j].c[2]);
    }
  }

  /* . . . and the color data. */
  int colorItr = 0;
  for (i=0; i < L->nc; i++) {
    for (j=0; j < L->cp[i].nv; j++) {
		if( colorParam == 0 )
		  fprintf(file,"%f %f %f %f\n", 0.0, 1.0*(pushes[colorItr]/maxPush), 0.0, 1.0);
		else
		  fprintf(file,"%f %f %f %f\n", 1.0*(pushes[colorItr]/maxPush), 1.0*(pushes[colorItr]/maxPush), 1.0*(pushes[colorItr]/maxPush), 0.0);
		colorItr++;
	}
  }
  
  fclose(file);
}


void
reloadDump( double* A, int rows, int cols, double* x, double* b )
{
	// reload expects things in this form
	static FILE* fp = NULL;
	if( fp == NULL )
	{
		fp = fopen("examples", "w");
	}
	int rItr, cItr;
	
	fprintf( fp, "%d %d\n", rows, cols );
	// A first
	for( rItr=0; rItr<rows; rItr++ )
	{
		for( cItr=0; cItr<cols; cItr++ )
			fprintf( fp, "%10.16lf ", A[rItr*cols + cItr] );
		fprintf( fp, "\n" );
	}
	
	// then x
	for( cItr=0; cItr<cols; cItr++ )
	{
		if( x != NULL )
			fprintf( fp, "%10.16lf\n", x[cItr] );
		else
			fprintf( fp, "0\n" );
	}
	
	// finally b
	for( rItr=0; rItr<rows; rItr++ )
		fprintf( fp, "%lf\n", b[rItr] );
		
	fflush(fp);
//	fclose(fp);
}


int gOutputFlag = 0;
int gConditionCheck = 0;

void 
bsearch_stepper( octrope_link** inLink, unsigned int inMaxSteps, search_state* inState )
{
	unsigned int i, offset = 0;
	int cItr;
	int lastEq = 0;
	unsigned int cSteps = 0;
	
	double maxmaxmin = 0;
	double minthickness = 500;
	double  nextMovieOutput = 0.0;
	
	gConditionCheck = 0;
	int firstRun = 1, secondRun = 0;
	
	FILE*   gnuplotPipes[kTotalGraphTypes];
	FILE*   gnuplotDataFiles[kTotalGraphTypes];
	
	clock_t startTime;
	startTime = clock();
	
	inState->steps = 0;
	
	/* inititalize piping if we are to use it */
	for( i=0; i<kTotalGraphTypes; i++ )
	{
		gnuplotPipes[i] = NULL;
		gnuplotDataFiles[i] = NULL;
		if( inState->graphing[i] != 0 )
		{
			char	fname[255];
			char	cmd[512];
			sprintf( fname, "%dpipe", i );
			sprintf( cmd, "rm -f %s ; mkfifo %s", fname, fname );
			system(cmd);
			
			// get gnuplot going and reading our commands
			sprintf( cmd, "/sw/bin/gnuplot < %s &", fname );
			system(cmd);
			
			// open the pipe for our input
			gnuplotPipes[i] = fopen(fname, "w");
		
			fprintf( gnuplotPipes[i], "set term aqua\n" );
		
			sprintf(fname, "%ddat", i);
			gnuplotDataFiles[i] = fopen(fname, "w");
		} 
	}
	
	int stepItr;
	for( stepItr=0; stepItr<inMaxSteps; stepItr++ )
	{
		int lastSet;
	
		lastSet = inState->lastStepStrutCount;
			
		//if( (i%50)==0 )
			gOutputFlag = 1;
	
		if( inState->shortest < .9999 )
			inState->curvature_step = 0;
		else
		{
			inState->curvature_step = 1;
			cSteps++;
		}
																		
	//	inState->curvature_step = ((inState->curvature_step+1)%2);
	//	inState->curvature_step = 1;
							
		*inLink = bsearch_step(*inLink, inState);
		inState->steps++;
				
	//	if( lastSet != state.lastStepStrutCount || (i%50)==0 )
		if( (i%50) == 0 ) // check things out every 50 steps
		{
	//		refresh_display(*inLink);
	//		printf( "eqing\n" );
		}
		if( inState->shortest > 2*inState->injrad )
			inState->shortest = 2*inState->injrad;
		
		if( inState->shortest < minthickness && inState->shortest != 0 )
			minthickness = inState->shortest;
			
		printf( "s: %d ms: %d len: %lf r: %lf ssize: %lf dcsd: %lf minrad: %lf avgdvdt: %lf residual: %e time: %lf\n", 
					inState->lastStepStrutCount, inState->lastStepMinradStrutCount,
					inState->length, inState->ropelength, inState->stepSize, inState->shortest, inState->minrad, 
					inState->avgDvdtMag, inState->residual, inState->time );
		
		
	/*	if( firstRun )
		{
			firstRun = 0;
			if( inState->shortest != 0 )
			{
				link_scale(*inLink, (2.0*inState->injrad)/inState->shortest);
				*inLink = octrope_fixlength(*inLink);
			}
			secondRun=1;
		}
		/*if( secondRun )
		{
			secondRun = 0;
			if( inState->shortest != 0 )
			{
				link_scale(*inLink, (2.0*inState->injrad)/inState->shortest);
			}
		}*/
				
		if( inState->oldRopeTime == 0 )
		{
			inState->oldRopeTime = inState->time;
			inState->oldRopelength = inState->ropelength;
		}
			
	//	if( inState->residual < 0.1 && inState->residual > 0 /*inState->minrad >= 0.5 /*&& inState->avgDvdtMag < 0.01*/ &&
	//		inState->curvature_step == 1/* && inState->time > 30*/ )
		if( (inState->oldRopeTime + inState->checkDelta < inState->time) && fabs(inState->oldRopelength-inState->ropelength) < .1 )
		{		
			// if we're going to die, it'll be after this, so save our best
			FILE* bestFile = NULL;
			char    fname[512];
			char adjustedName[512];
			
			if( inState->batching == 0 )
			{
				strcpy(adjustedName, inState->fname);
				sprintf( fname, "%s_%d.best", adjustedName, inState->totalVerts );
				bestFile = fopen(fname, "w");
			}
			else
			{
				strcpy(adjustedName, inState->fname);
				strstr(adjustedName,".vect")[0] = '\0';
				sprintf( fname, "%s_%d.best", adjustedName, inState->totalVerts );
				bestFile = fopen(fname, "w");
			}
			
			octrope_link_write(bestFile, *inLink);
			fclose(bestFile);
		
		//	if( inState->totalVerts > 4*inState->ropelength )
			{
				// we are DONE!
				break;
			}
			
			// the strut set can get violent after this change, so let's make sure none exist.
			// it's a lot easier on the condition number if it reemerges quickly than if it's changing
			// rapidly
			link_scale(*inLink, 1.05*(2.0*inState->injrad)/inState->shortest);
		
			octrope_link* oldLink = *inLink;
			*inLink = octrope_double_edges(*inLink);
			octrope_link_free(oldLink);
			
			gConditionCheck = 20; // check condition number for next 10 evals -- this is where we'll fail
			
			// we need to update some state information also
			inState->totalVerts = 0;
			for( cItr=0; cItr<(*inLink)->nc; cItr++ )
			{
				inState->totalVerts += (*inLink)->cp[cItr].nv;
			}
			
			updateSideLengths(*inLink, inState);
			
		//	inState->stepSize = inState->avgDvdtMag*inState->avgDvdtMag;
		//	if( kStepScale*octrope_link_short_edge(*inLink) < inState->stepSize )
		//		inState->stepSize = kStepScale*octrope_link_short_edge(*inLink);
			offset = 0;
			for( i=0; i<(*inLink)->nc; i++ )
			{
				inState->compOffsets[i] = offset;
				offset += (*inLink)->cp[i].nv;
			}
			
			gOutputFlag = 1;
			
			inState->residual = 500; // make this big
		}
		
		if( (inState->oldRopeTime + inState->checkDelta) < inState->time )
		{
			inState->oldRopelength = inState->ropelength;
			inState->oldRopeTime = inState->time;
			printf( "* Checked delta rope and continuing\n" );
		}
		
		double maxmin;
		maxmin = maxovermin(*inLink);
		if( maxmin > maxmaxmin )
			maxmaxmin = maxmin;
		//inState->eqMultiplier = pow(10,maxmin);
		printf( "maxovermin: %3.5lf eqMultiplier: %3.5lf (max max/min: %3.5lf) (min thickness: %3.5lf) last rcond: %lf ssize: %f cstep: %d\n", 
				    maxmin, inState->eqMultiplier, maxmaxmin, minthickness, inState->rcond, inState->stepSize, inState->curvature_step );
		printf( "cstep/step ratio: %lf delta rope: %lf next check: %lf\n", 
				((double)cSteps)/((double)stepItr+1), fabs(inState->oldRopelength-inState->ropelength),
				inState->oldRopeTime+inState->checkDelta );

		if( inState->time >= nextMovieOutput )
		{
			char	fname[384];
			FILE*   frame = NULL;
			
			
			
			if( inState->batching != 0 )
			{
			    sprintf( fname, "restart_%s", inState->fname );
			    (strstr(fname,".vect"))[0] = '\0';
			    printf( "saved restart: %s\n", fname );
			    frame = fopen( fname, "w" );
			    octrope_link_write(frame, *inLink);
			    fclose(frame);
			    
			    nextMovieOutput += 0.05;
			}
			else
			{
			    sprintf( fname, "rmov.%lf_evals-%d_strts-%d_minrads-%d.vect", inState->time, 
							    inState->tsnnls_evaluations, 
							    inState->lastStepStrutCount,
							    inState->lastStepMinradStrutCount );
			  /*  frame = fopen( fname, "w" );
			    octrope_link_write(frame, *inLink);
			    fclose(frame);
			    */
			    nextMovieOutput += 0.05;
			    
			    printf( "movie frame output (tsnnls evals: %d)\n", inState->tsnnls_evaluations );
			    gOutputFlag = 1;
			    
			    /* also do graph outputs */
			    for( i=0; i<kTotalGraphTypes; i++ )
			    {
				    if( inState->graphing[i] != 0 )
				    {
					    char fname[255];
					    sprintf(fname, "%ddat", i);
					    if( i != kConvergence )
						    fprintf( gnuplotDataFiles[i], "%lf\t", inState->time );
					    else
						    fprintf( gnuplotDataFiles[i], "%u\t", inState->steps );
					    switch( i )
					    {
						    case kLength: fprintf( gnuplotDataFiles[i], "%lf", inState->length ); break;
						    case kConvergence: 
						    case kRopelength: fprintf( gnuplotDataFiles[i], "%lf", inState->ropelength ); break;
						    case kStrutCount: fprintf( gnuplotDataFiles[i], "%d", inState->lastStepStrutCount ); break;
						    case kStepSize: fprintf( gnuplotDataFiles[i], "%lf", inState->stepSize ); break;
						    case kThickness: fprintf( gnuplotDataFiles[i], "%lf", inState->shortest ); break;
						    case kMinrad: fprintf( gnuplotDataFiles[i], "%lf", inState->minrad ); break;
						    case kResidual: fprintf( gnuplotDataFiles[i], "%lf", inState->residual ); break;
						    case kMaxOverMin: fprintf( gnuplotDataFiles[i], "%lf", maxmin ); break;
						    case kRcond: fprintf( gnuplotDataFiles[i], "%lf", inState->rcond ); break;
						    case kWallTime: fprintf( gnuplotDataFiles[i], "%lf", (clock() - startTime)/(CLOCKS_PER_SEC / (double) 1000.0) ); break;
						    case kMaxVertexForce: fprintf( gnuplotDataFiles[i], "%lf", (inState->maxPush)/(inState->length/(double)inState->totalVerts) ); break;
					    }
					    fprintf( gnuplotDataFiles[i], "\n" );
					    fflush(gnuplotDataFiles[i]);
					    // show only last 5 seconds
				    //	fprintf( gnuplotPipes[i], "set xrange [%lf:%lf]\n", inState->time-5, inState->time );
					    
					    fprintf( gnuplotPipes[i], "plot \"%s\" using 1:2 w lines title \'", fname );
					    switch( i )
					    {
						    case kLength: fprintf( gnuplotPipes[i], "Length" ); break;
						    case kRopelength: fprintf( gnuplotPipes[i], "Ropelength" ); break;
						    case kStrutCount: fprintf( gnuplotPipes[i], "Strut count" ); break;
						    case kStepSize: fprintf( gnuplotPipes[i], "Step size" ); break;
						    case kThickness: fprintf( gnuplotPipes[i], "Thickness" ); break;
						    case kMinrad: fprintf( gnuplotPipes[i], "Minrad" ); break;
						    case kResidual: fprintf( gnuplotPipes[i], "Residual" ); break;
						    case kMaxOverMin: fprintf( gnuplotPipes[i], "Max/min" ); break;
						    case kRcond: fprintf( gnuplotPipes[i], "reciprocal condition number of A^T A" ); break;
						    case kWallTime: fprintf( gnuplotPipes[i], "process computation time" ); break;
						    case kMaxVertexForce: fprintf( gnuplotPipes[i], "maximum compression sum" ); break;
						    case kConvergence: fprintf( gnuplotPipes[i], "convergence (rope vs steps)" ); break;
					    }
					    fprintf( gnuplotPipes[i], "\'\n" ); 
					    fflush( gnuplotPipes[i] );
				    }
			    }
			} // if not batching
		}
		
		if( gOutputFlag == 1 )
		{
			FILE* verts = fopen("/tmp/verts.vect", "w");
			octrope_link_draw(verts, *inLink);
			fclose(verts);
		}
									
		// eq ourselves
	//	maxovermin(*inLink);
		//lapack_eqedge( *inLink, 4 );
	
		//lapack_eqedge(*inLink, 4);		
		
		if( inState->eqThreshold != 0 && maxmin > inState->eqThreshold )
		{
			int eqsteps;
			octrope_link* oldlink;
			
			lastEq = i;
		
			printf( "**\tEqing\t**\n" );
			
			for( eqsteps=0; eqsteps<4; eqsteps++ )
			{
				oldlink = *inLink;
				*inLink = octrope_fixlength(oldlink);
				octrope_link_free(oldlink);
			}
		}
		
		
	//	maxovermin(*inLink);
	}
	
	/* inititalize piping if we are to use it */
	for( i=0; i<kTotalGraphTypes; i++ )
	{
		if( inState->graphing[i] != 0 )
		{
			char	cmd[512];
			char	fname[255];
			sprintf( fname, "%dpipe", i );
			// this will quit gnuplot
			sprintf( cmd, "echo quit >> %s", fname );
			system(cmd);
			fclose(gnuplotPipes[i]);
			sprintf(cmd, "rm -f %s", fname);
			system(cmd);
			
			sprintf( fname, "%ddat", i );
			fclose(gnuplotDataFiles[i]);
		
			// you might want to keep these, actually
		//	sprintf( cmd, "rm -f %s", fname );
		//	system(cmd);
		} 
	}
	
	// clean up
	free(inState->compOffsets);
	
}

static void
dumpVertsStruts(octrope_link* link, octrope_strut* strutSet, int strutCount)
{
	int cItr, vItr, totalVerts=0;
	FILE* fp = fopen("ridgeverts","w");
	
	for( cItr=0; cItr<link->nc; cItr++ )
	{
		totalVerts += link->cp[cItr].nv;
	}
	
	fprintf( fp, "%d\n", totalVerts );
	
	for( cItr=0; cItr<link->nc; cItr++ )
	{
		for( vItr=0; vItr<link->cp[cItr].nv; vItr++ )
		{
			fprintf(fp, "%lf\n%lf\n%lf\n", link->cp[cItr].vt[vItr].c[0],
				link->cp[cItr].vt[vItr].c[1], link->cp[cItr].vt[vItr].c[2] );
		}
	}
	
	fprintf( fp, "%d\n", strutCount );
	// strut indices and positions. ASSUME ONE COMPONENT
	for( vItr=0; vItr<strutCount; vItr++ )
		fprintf( fp, "%d %d\n", strutSet[vItr].lead_vert[0], strutSet[vItr].lead_vert[1] );
	for( vItr=0; vItr<strutCount; vItr++ )
		fprintf( fp, "%lf %lf\n", strutSet[vItr].position[0], strutSet[vItr].position[1] );
	
	fclose(fp);
}

static void
dumpAxb( taucs_ccs_matrix* A, double* x, double* b )
{
	int i, j;
	double* vals = taucs_convert_ccs_to_doubles(A);
	FILE* fp = fopen("/Users/michaelp/examples", "w");
	
	// out dims
	fprintf( fp, "%d %d\n", A->m, A->n );
	
	// output A
	for( i=0; i<A->m; i++ )
	{
		for( j=0; j<A->n; j++ )
			fprintf( fp, "%lf ", vals[i*A->n + j] );
		fprintf(fp, "\n");
	}
	
	// output x
	for( i=0; i<A->n; i++ )
		fprintf( fp, "%lf\n", x[i] );
	
	// output b
	for( i=0; i<A->m; i++ )
		fprintf( fp, "%lf\n", b[i] );
	
	free(vals);
	fclose(fp);
}

static void
dumpDvdt( octrope_vector* dvdt, int size )
{
	int vItr;
	FILE* fp = fopen("/Users/michaelp/dvdt","a");
	fprintf( fp, "%d\n", size );
	
		for( vItr=0; vItr<size; vItr++ )
		{
			fprintf(fp, "%lf\n%lf\n%lf\n", dvdt[vItr].c[0],
					dvdt[vItr].c[1], dvdt[vItr].c[2] );
		}
		fclose(fp);
	
}

void
step( octrope_link* inLink, double stepSize, octrope_vector* dVdt )
{
	int cItr, vItr, dVdtItr;
	for( cItr=0, dVdtItr=0; cItr<inLink->nc; cItr++ )
	{
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++, dVdtItr++ )
		{
			inLink->cp[cItr].vt[vItr].c[0] += stepSize*dVdt[dVdtItr].c[0];
			inLink->cp[cItr].vt[vItr].c[1] += stepSize*dVdt[dVdtItr].c[1];
			inLink->cp[cItr].vt[vItr].c[2] += stepSize*dVdt[dVdtItr].c[2];
		}
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
barForce( octrope_vector* dVdt, octrope_link* inLink, search_state* inState )
{
	/* solve the unconstrained least squares problem Ax=b where A is the rigidity
	 * matrix formed by bars on the vertices of the link, and b is the change in each
	 * component's x, y, z from its original length
	 */
	double*		A = NULL, *b = NULL;
	int			barVerts=0, cItr, vItr;
	int			bars=0, aIndexer;
	
	octrope_link_fix_wrap(inLink);

	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		if( inState->conserveLength[cItr] != 0 )
		{
			barVerts += inLink->cp[cItr].nv;
			bars += octrope_pline_edges(&inLink->cp[cItr]);
		}
	}
	
	A = (double*)malloc(sizeof(double)*barVerts*3*bars);
	placeVertexBars(A, inLink, 0, barVerts, bars, inState);
	b = (double*)calloc(barVerts*3, sizeof(double));

	/* now compute b, the change in each dimension of each vertex required to get to eq */
	octrope_vector* sides;
	double*			lengths;
	int				edges;
		
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		if( inState->conserveLength[cItr] != 0 )
		{
			double length;
			edges = octrope_pline_edges(&inLink->cp[cItr]);
				
			for( vItr=0; vItr<edges; vItr++ )
			{
				octrope_vector s1, s2, side;
				s1.c[0] = inLink->cp[cItr].vt[vItr].c[0];
				s1.c[1] = inLink->cp[cItr].vt[vItr].c[1];
				s1.c[2] = inLink->cp[cItr].vt[vItr].c[2];
				s2.c[0] = inLink->cp[cItr].vt[vItr+1].c[0];
				s2.c[1] = inLink->cp[cItr].vt[vItr+1].c[1];
				s2.c[2] = inLink->cp[cItr].vt[vItr+1].c[2];
				side.c[0] = s2.c[0] - s1.c[0];
				side.c[1] = s2.c[1] - s1.c[1];
				side.c[2] = s2.c[2] - s1.c[2];
				
				length = octrope_norm(side);
				
				// unit-ize sides since we will use these as the tangential motion basis below
				side.c[0] /= length;
				side.c[1] /= length;
				side.c[2] /= length;
												
				b[inState->compOffsets[cItr]+vItr + 0] = -((inState->sideLengths[inState->compOffsets[cItr]+vItr] - length) * side.c[0]);
				b[inState->compOffsets[cItr]+vItr + 1] = -((inState->sideLengths[inState->compOffsets[cItr]+vItr] - length) * side.c[1]);
				b[inState->compOffsets[cItr]+vItr + 2] = -((inState->sideLengths[inState->compOffsets[cItr]+vItr] - length) * side.c[2]);
			
			//	b[inState->compOffsets[cItr]+vItr + 0] = dVdt[inState->compOffsets[cItr]+vItr + 0].c[0];
			//	b[inState->compOffsets[cItr]+vItr + 1] = dVdt[inState->compOffsets[cItr]+vItr + 0].c[1];
			//	b[inState->compOffsets[cItr]+vItr + 2] = dVdt[inState->compOffsets[cItr]+vItr + 0].c[2];
																		
				printf( "adjusting (%d) %d: %lf\n", cItr, vItr, (inState->sideLengths[inState->compOffsets[cItr]+vItr] - length) );
			}
		}
	}
	
	if( barVerts != 0 )
	{
		double* compressions;
		taucs_ccs_matrix* sparseA;
		
		sparseA = taucs_construct_sorted_ccs_matrix(A, bars, 3*barVerts);
		compressions = t_lsqr( sparseA, b );
	//	taucs_print_ccs_matrix(sparseA);
		taucs_ccs_free(sparseA);
		
		for( cItr=0; cItr<inLink->nc; cItr++ )
		{
			if( inState->conserveLength[cItr] != 0 )
			{
				double partialMult = 0;
				int totalStruts = bars;
				int dlItr, sItr;
				
				for( dlItr=inState->compOffsets[cItr]; dlItr<inState->compOffsets[cItr]+inLink->cp[cItr].nv; dlItr++ )
				{
					int cItr2;
					/* of course, dlItr is not quite right here since we are including ALL vertices in its computation, but all vertices
					 * are not accounted for in A. We need to adjust A's indexer here to account for this, which amounts to subtracting
					 * from dlItr all the components which preceed the one we are currently processing that are NOT bar composed
					 */
					aIndexer = dlItr;
					for( cItr2=0; cItr2<cItr; cItr2++ )
					{
						if( inState->conserveLength[cItr2] == 0 )
							aIndexer -= inLink->cp[cItr2].nv;
					}
				
					for( sItr=0; sItr<totalStruts; sItr++ )
					{
						if( (totalStruts*3*aIndexer)+sItr+0 >= bars*3*barVerts )
							printf( "*** crap!\n" );
						partialMult += compressions[sItr]*A[(totalStruts*3*aIndexer)+sItr+0];
					}
					dVdt[dlItr].c[0] = dVdt[dlItr].c[0] + partialMult;
				
					partialMult=0;
					for( sItr=0; sItr<totalStruts; sItr++ )
					{
						if( (totalStruts*3*aIndexer)+sItr+0 >= bars*3*barVerts )
							printf( "*** crap!\n" );
						partialMult += compressions[sItr]*A[(totalStruts*3*aIndexer)+sItr+totalStruts];
					}
					dVdt[dlItr].c[1] = dVdt[dlItr].c[1] + partialMult;
					
					partialMult = 0;
					for( sItr=0; sItr<totalStruts; sItr++ )
					{
						if( (totalStruts*3*aIndexer)+sItr+0 >= bars*3*barVerts )
							printf( "*** crap!\n" );
						partialMult += compressions[sItr]*A[(totalStruts*3*aIndexer)+sItr+(2*totalStruts)];
					}
					dVdt[dlItr].c[2] = dVdt[dlItr].c[2] + partialMult;
				}
			}
		}
		
		free(compressions);
	}
	
	// debug code!
	exportVect(dVdt, inLink, "/tmp/dVdt.vect");

	free(A);
	free(b);
}

octrope_link*
bsearch_step( octrope_link* inLink, search_state* inState )
{	
	int stepAttempts;
	int dvdtItr;
	double  lastDCSD, eps = inState->stepSize*inState->stepSize;
	octrope_link* workerLink;

	// create initial vector field for which we want to move along in this step
	octrope_vector* dVdt;
	
	dVdt = calloc(inState->totalVerts, sizeof(octrope_vector));
	
	/* 
	 * we need to get the initial strut set to compare our steps to.
	 * this must occur at the start of every step -- we cannot use
	 * the old value since we have changed the curve (and possibly
	 * the strut set) during equilateralization
	 */
	
	// the ordering of these forces is important -- first we introduce the movement force 
	// which we must resolve self contacts over -- but since we don't want to cancel the 
	// force of redistributing the curve (bar force) we place that after strut resolution
	
//	specialForce(dVdt, inLink, inState);
	firstVariation(dVdt, inLink, inState, NULL, &inState->lastStepStrutCount, inState->curvature_step);
	if( inState->shortest > 2*inState->injrad )
		inState->shortest = 2*inState->injrad;
	
	lastDCSD = inState->shortest;

	stepAttempts = 0;
	workerLink = NULL;
	do
	{			
		stepAttempts++;
				
		// move along dVdt
		if( workerLink != NULL )
			octrope_link_free(workerLink);
		workerLink = octrope_link_copy(inLink);
		step(workerLink, inState->stepSize, dVdt);
		if( gOutputFlag == 1 )
		{
			exportVect(dVdt, inLink, "/tmp/dVdt.vect");
			printf( "visualization output\n" );
			gOutputFlag = 0;
		}

		free(dVdt);
		
		// get new dvdt and basically see if we've overrun our error bound for
		// strut length
		dVdt = calloc(inState->totalVerts, sizeof(octrope_vector));
		firstVariation(dVdt, workerLink, inState, NULL, &inState->lastStepStrutCount, inState->curvature_step);
				
		if( inState->shortest < lastDCSD-eps )
		{
			//inState->stepSize = fmax(inState->stepSize*0.5, 0.0039);
			inState->stepSize *= 0.5;
			
			/* If our step size is getting too small, we can surmise that we are at a point of discontinuity and 
			 * we just need to step over it. This will happen whenever a strut is created, since error proportional 
			 * to ssize^2 will always win out with half stepping BEFORE the strut is created (since that's 0 error!)
			 * So it's easy to see that step size spirals down (because tolerated error goes down) and the strut is just 
			 * never created, but this fixes that, which is handy, and the corrector steps will smooth things out for us later
			 * by making sure we're in the 'good zone' of strutness (which is to say struts that are arbitrarily close to thickness)
			 */
			if( inState->stepSize < kMinStepSize ) // notice that this means eps is ss^2, which is close to max double precision!
			{
				inState->stepSize *= 2;
				break;
			}
			continue;
		}
		else
		{
			break;
		}
		
	} while( 1==1 );
	
	inState->time += inState->stepSize;
	
	octrope_link_free(inLink);
	inLink = workerLink;
		
	// we're good, double step size and continue jamming
	if( stepAttempts == 1 && inState->curvature_step != 0 )
	{
		inState->stepSize *= 2;
	}
	
	// grab average dvdt now that we have finished mungering it
	inState->avgDvdtMag=0;
	for( dvdtItr=0; dvdtItr<inState->totalVerts; dvdtItr++ )
	{
		inState->avgDvdtMag += octrope_norm(dVdt[dvdtItr]);
	}
	inState->avgDvdtMag /= inState->totalVerts;

	// we should make sure stepsize isn't > avgDvdt^2 as that's the minrad control bound
	if( inState->stepSize > inState->avgDvdtMag*inState->avgDvdtMag && inState->curvature_step != 0 &&
		inState->avgDvdtMag*inState->avgDvdtMag > kMinStepSize) // this keeps us from zeroing on corrector steps
	{
		inState->stepSize = inState->avgDvdtMag*inState->avgDvdtMag;
	}
	
	// also shouldn't be > 10% of edgelength
	if( inState->stepSize > inState->length/inState->totalVerts*.1 )
		inState->stepSize = inState->length/inState->totalVerts*.1;
	
	free(dVdt);
	
	return inLink;
}

static void
placeContactStruts( double* A, octrope_link* inLink, octrope_strut* strutSet, int strutCount, search_state* inState, int minradStruts )
{
	int sItr, totalStruts;
	
	// in constructing the ridigity matrix, we will need the struts as viewed 
	// as force vectors on the edges, so we create normalized vectors for each strut
	// here
	octrope_vector* strutDirections = (octrope_vector*)calloc(strutCount, sizeof(octrope_vector));
	normalizeStruts( strutDirections, strutSet, inLink, strutCount );

	totalStruts = minradStruts + strutCount;
	
	for( sItr=0; sItr<strutCount; sItr++ )
	{
		int		entry;
		
		// temporarily increment the strut's verts based on their component interactions
		// we undo this change at the end of the for loop in case the user
		// wants to keep the strut set
		strutSet[sItr].lead_vert[0] += inState->compOffsets[strutSet[sItr].component[0]];
		strutSet[sItr].lead_vert[1] += inState->compOffsets[strutSet[sItr].component[1]];
		
		// entry is the offset in A which begin this strut's influce
		// it corresponds to the x influence on the lead_vert[0]th vertex
		// after this line, entry+1 is y, +2: z.
		entry = (totalStruts*3*strutSet[sItr].lead_vert[0])+sItr;
	
		// the strut information includes the position from the strut.lead_vert
		// so we assign "1-position[0]" of the force to the lead vert and "position[0]"
		// of the force to lead_vert+1
		A[entry] = (1-strutSet[sItr].position[0]) * strutDirections[sItr].c[0];
		A[entry+totalStruts] = (1-strutSet[sItr].position[0]) * strutDirections[sItr].c[1];
		A[entry+(2*totalStruts)] = (1-strutSet[sItr].position[0]) * strutDirections[sItr].c[2];
		
		// now for the next vertex, receiving "position[0]" of the force, this is 
		// potential wrapping case
		if( (strutSet[sItr].lead_vert[0]-inState->compOffsets[strutSet[sItr].component[0]]) == (inLink->cp[strutSet[sItr].component[0]].nv-1) &&
			(inLink->cp[strutSet[sItr].component[0]].acyclic == 0) )
		{
			entry = 0;
		}
		else
		{
			entry = (totalStruts*3*(strutSet[sItr].lead_vert[0]+1))+sItr;
		}
		A[entry] = (strutSet[sItr].position[0]) * strutDirections[sItr].c[0];
		A[entry+totalStruts] = (strutSet[sItr].position[0]) * strutDirections[sItr].c[1];
		A[entry+(2*totalStruts)] = (strutSet[sItr].position[0]) * strutDirections[sItr].c[2];
		
		
		// we do the same thing at the opposite end of the strut, except now the 
		// force is negated
		entry = (totalStruts*3*strutSet[sItr].lead_vert[1])+sItr;
					
		A[entry] = (1-strutSet[sItr].position[1]) * -strutDirections[sItr].c[0];
		A[entry+totalStruts] = (1-strutSet[sItr].position[1]) * -strutDirections[sItr].c[1];
		A[entry+(2*totalStruts)] = (1-strutSet[sItr].position[1]) * -strutDirections[sItr].c[2];
		if( (strutSet[sItr].lead_vert[1]-inState->compOffsets[strutSet[sItr].component[1]]) == (inLink->cp[strutSet[sItr].component[1]].nv-1) &&
			(inLink->cp[strutSet[sItr].component[1]].acyclic == 0) )
		{
			entry = 0;
		}
		else
		{
			entry = (totalStruts*3*(strutSet[sItr].lead_vert[1]+1))+sItr;
		}
		A[entry] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[0];
		A[entry+totalStruts] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[1];
		A[entry+(2*totalStruts)] = (strutSet[sItr].position[1]) * -strutDirections[sItr].c[2];
	
		strutSet[sItr].lead_vert[0] -= inState->compOffsets[strutSet[sItr].component[0]];
		strutSet[sItr].lead_vert[1] -= inState->compOffsets[strutSet[sItr].component[1]];
	}
	
	free(strutDirections);
}

void
updateSideLengths( octrope_link* inLink, search_state* inState )
{
	int cItr, vItr;
	
	inState->totalSides = 0;
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		inState->totalSides += inLink->cp[cItr].nv;
	}
	
	if( inState->sideLengths != NULL )
		free(inState->sideLengths);
	
	inState->sideLengths = (double*)malloc(inState->totalSides*sizeof(double));
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
		{
			octrope_vector s1, s2, side;
			s1.c[0] = inLink->cp[cItr].vt[vItr].c[0];
			s1.c[1] = inLink->cp[cItr].vt[vItr].c[1];
			s1.c[2] = inLink->cp[cItr].vt[vItr].c[2];
			s2.c[0] = inLink->cp[cItr].vt[vItr+1].c[0];
			s2.c[1] = inLink->cp[cItr].vt[vItr+1].c[1];
			s2.c[2] = inLink->cp[cItr].vt[vItr+1].c[2];
			side.c[0] = s2.c[0] - s1.c[0];
			side.c[1] = s2.c[1] - s1.c[1];
			side.c[2] = s2.c[2] - s1.c[2];
		
			inState->sideLengths[inState->compOffsets[cItr] + vItr] = octrope_norm(side);
		}
	}
}

void
placeVertexBars( double* A, octrope_link* inLink, int contactStruts, int totalBarVerts, int totalBars, search_state* inState )
{
	int totalStruts = contactStruts + totalBars;
	int cItr, vItr, sItr, next;
	
	octrope_vector  strutDirection;
	octrope_vector  points[2];
	double norm = 0;

	int thresh = ((totalBarVerts+contactStruts)*3)*(totalBars+contactStruts);
	
	octrope_strut*  vertStruts;
	
	vertStruts = (octrope_strut*)malloc(sizeof(octrope_strut)*totalBars);
	
	sItr=0;
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		if( inState->conserveLength[cItr] != 0 )
		{
			for( vItr=0; vItr<octrope_pline_edges(&inLink->cp[cItr]); vItr++ )
			{
				next = ((vItr+1)%inLink->cp[cItr].nv);
			
				vertStruts[sItr].component[0] = cItr;
				vertStruts[sItr].component[1] = cItr;
				vertStruts[sItr].position[0] = 0;
				vertStruts[sItr].position[1] = 0;
				vertStruts[sItr].lead_vert[0] = vItr;
				vertStruts[sItr].lead_vert[1] = next;
				
				sItr++;
			}
		}
	}
	
	// and now we can place these just like minrad struts
	for( sItr=0; sItr<totalBars; sItr++ )
	{
		int entry;
		
	//	printf( "(%d) (%d) %d %d\n", vertStruts[sItr].component[0], vertStruts[sItr].component[1], vertStruts[sItr].lead_vert[0], vertStruts[sItr].lead_vert[1] );
		
		// calc the norm of the converted strut
		octrope_strut_ends( inLink, &vertStruts[sItr], points );
		
		// the normalized difference of pointOne, pointTwo is the strut force vector
		strutDirection.c[0] = points[0].c[0] - points[1].c[0];
		strutDirection.c[1] = points[0].c[1] - points[1].c[1];
		strutDirection.c[2] = points[0].c[2] - points[1].c[2];
		
		norm = octrope_norm(strutDirection);
		strutDirection.c[0] /= norm;
		strutDirection.c[1] /= norm;
		strutDirection.c[2] /= norm;
		
		/* temporarily increment the strut's verts based on their component interactions, 
		 * but in the vertex bar rigidity matrix, only some component's vertices are counted --
		 * namely the ones that are composed of bars. We need to increase our offset, but only by those
		 * guys we have skipped who are bar composed.
		 */
		for( cItr=0; cItr<vertStruts[sItr].component[0]; cItr++ )
		{
			if( inState->conserveLength[cItr] != 0 )
			{
				/* we update both here even though we only check first strut end because vertex bars are always
				 * on the same component */
				vertStruts[sItr].lead_vert[0] += inLink->cp[cItr].nv;
				vertStruts[sItr].lead_vert[1] += inLink->cp[cItr].nv;
			}
		}
		
		// entry is the offset in A which begin this strut's influce
		// it corresponds to the x influence on the lead_vert[0]th vertex
		// after this line, entry+1 is y, +2: z.
		entry = (totalStruts*3*vertStruts[sItr].lead_vert[0])+contactStruts+sItr;
	
		if( entry+(2*totalStruts) >= thresh )
			printf( "****crap!\n" );	
	
		// the strut information includes the position from the strut.lead_vert
		// so we assign "1-position[0]" of the force to the lead vert and "position[0]"
		// of the force to lead_vert+1
		A[entry] = (1-vertStruts[sItr].position[0]) * strutDirection.c[0];
		A[entry+totalStruts] = (1-vertStruts[sItr].position[0]) * strutDirection.c[1];
		A[entry+(2*totalStruts)] = (1-vertStruts[sItr].position[0]) * strutDirection.c[2];
		
		/***************** we don't need to do this since there are no midpoint vertex bars *****************/
		
		// now for the next vertex, receiving "position[0]" of the force, this is 
		// potential wrapping case
	/*	if( (vertStruts[sItr].lead_vert[0]-inState->compOffsets[vertStruts[sItr].component[0]]) == (inLink->cp[vertStruts[sItr].component[0]].nv-1) &&
			(inLink->cp[vertStruts[sItr].component[0]].acyclic == 0) )
		{
			entry = 0;
		}
		else
		{
			entry = (totalStruts*3*(vertStruts[sItr].lead_vert[0]+1))+contactStruts+sItr;
		}
		
		if( entry+(2*totalStruts) >= thresh )
			printf( "****crap!\n" );
		
		A[entry] = (vertStruts[sItr].position[0]) * strutDirection.c[0];
		A[entry+totalStruts] = (vertStruts[sItr].position[0]) * strutDirection.c[1];
		A[entry+(2*totalStruts)] = (vertStruts[sItr].position[0]) * strutDirection.c[2];
		*/
		
		// we do the same thing at the opposite end of the strut, except now the 
		// force is negated
		entry = (totalStruts*3*vertStruts[sItr].lead_vert[1])+contactStruts+sItr;
		
		if( entry+(2*totalStruts) >= thresh )
			printf( "****crap!\n" );			
		
		A[entry] = (1-vertStruts[sItr].position[1]) * -strutDirection.c[0];
		A[entry+totalStruts] = (1-vertStruts[sItr].position[1]) * -strutDirection.c[1];
		A[entry+(2*totalStruts)] = (1-vertStruts[sItr].position[1]) * -strutDirection.c[2];
		
		/*if( (vertStruts[sItr].lead_vert[1]-inState->compOffsets[vertStruts[sItr].component[1]]) == (inLink->cp[vertStruts[sItr].component[1]].nv-1) &&
			(inLink->cp[vertStruts[sItr].component[1]].acyclic == 0) )
		{
			entry = 0;
		}
		else
		{
			entry = (totalStruts*3*(vertStruts[sItr].lead_vert[1]+1))+contactStruts+sItr;
		}
		
		if( entry+(2*totalStruts) >= thresh )
			printf( "****crap!\n" );	
		
		A[entry] = (vertStruts[sItr].position[1]) * -strutDirection.c[0];
		A[entry+totalStruts] = (vertStruts[sItr].position[1]) * -strutDirection.c[1];
		A[entry+(2*totalStruts)] = (vertStruts[sItr].position[1]) * -strutDirection.c[2];
*/
		/* readjusting is problematic for the same reason that adjusting is (see above) and it's not necessary in this 
		 * strut placement, so we just don't do it
		 */
		 
		//vertStruts[sItr].lead_vert[0] -= inState->compOffsets[vertStruts[sItr].component[0]];
		//vertStruts[sItr].lead_vert[1] -= inState->compOffsets[vertStruts[sItr].component[1]];
	}
	
	free(vertStruts);
}

static void
placeMinradStruts( double* A, octrope_link* inLink, octrope_mrloc* minradStruts, 
	int minradLocs, search_state* inState, int contactStruts )
{
	int mItr, i;
	double norm = 0;

	int totalStruts = contactStruts + minradLocs;
		
	for( mItr=0; mItr<minradLocs; mItr++ )
	{
		// we can place these just like regular struts as in placeContactStruts if we just
		// convert the representation
		int		entry;
		double c1, c2;
		
		// view mrloc as center guy, figure out whether ab or bc is longer
		octrope_vector  ab, bc, ba, cb, abperp, bcperp;
		// b-a
		ab.c[0] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[0] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert-1].c[0];
		ab.c[1] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[2] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert-1].c[1];
		ab.c[2] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[1] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert-1].c[2];
		
		// a-b
		ba.c[0] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert-1].c[0] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[0];
		ba.c[1] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert-1].c[2] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[1];
		ba.c[2] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert-1].c[1] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[2];
		
		// c-b
		bc.c[0] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert+1].c[0] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[0];
		bc.c[1] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert+1].c[1] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[1];
		bc.c[2] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert+1].c[2] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[2];
		
		// b-c
		cb.c[0] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[0] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert+1].c[0];
		cb.c[1] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[1] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert+1].c[1];
		cb.c[2] = inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert].c[2] - 
						inLink->cp[minradStruts[mItr].component].vt[minradStruts[mItr].vert+1].c[2];
		
		
		bcperp = octrope_cross(octrope_cross(ba, bc), bc);
		abperp = octrope_cross(octrope_cross(ba, bc), ab);
		
		octrope_vector  As, Bs, Cs;
		octrope_vector  prevSide, thisSide;
		double dot, prevLen, thisLen, angle;
		int vItr = minradStruts[mItr].vert;
		int cItr = minradStruts[mItr].component;

		prevSide.c[0] = inLink->cp[cItr].vt[vItr].c[0] - inLink->cp[cItr].vt[vItr-1].c[0];
		prevSide.c[1] = inLink->cp[cItr].vt[vItr].c[1] - inLink->cp[cItr].vt[vItr-1].c[1];
		prevSide.c[2] = inLink->cp[cItr].vt[vItr].c[2] - inLink->cp[cItr].vt[vItr-1].c[2];
		
		thisSide.c[0] = inLink->cp[cItr].vt[vItr+1].c[0] - inLink->cp[cItr].vt[vItr].c[0];
		thisSide.c[1] = inLink->cp[cItr].vt[vItr+1].c[1] - inLink->cp[cItr].vt[vItr].c[1];
		thisSide.c[2] = inLink->cp[cItr].vt[vItr+1].c[2] - inLink->cp[cItr].vt[vItr].c[2];
		
		dot = octrope_dot(prevSide, thisSide);
		prevLen = octrope_norm(prevSide);
		thisLen = octrope_norm(thisSide);
		
		// this needs to be corrected if we're straight or something
		angle = acos(dot/(prevLen*thisLen));
																																																																																						
		if( octrope_norm(ab) < octrope_norm(bc) )
		{			
			c1 = (1.0/(2*tan(angle/2)*prevLen));
			c2 = prevLen/(8*sin(angle/2)*sin(angle/2));
		
			As.c[0] = c1*ba.c[0];
			As.c[1] = c1*ba.c[1];
			As.c[2] = c1*ba.c[2];
			
			As.c[0] += -c2*(abperp.c[0]/(prevLen*prevLen));
			As.c[1] += -c2*(abperp.c[1]/(prevLen*prevLen));
			As.c[2] += -c2*(abperp.c[2]/(prevLen*prevLen));
		
			Bs.c[0] = c1*ab.c[0];
			Bs.c[1] = c1*ab.c[1];
			Bs.c[2] = c1*ab.c[2];
			
			Bs.c[0] += c2*( (ab.c[0]/(prevLen*prevLen)) + (bcperp.c[0]/(thisLen*thisLen)) );
			Bs.c[1] += c2*( (ab.c[1]/(prevLen*prevLen)) + (bcperp.c[1]/(thisLen*thisLen)) );
			Bs.c[2] += c2*( (ab.c[2]/(prevLen*prevLen)) + (bcperp.c[2]/(thisLen*thisLen)) );
			
			c1 = (-prevLen)/(8*sin(angle/2)*sin(angle/2)*thisLen*thisLen);
			Cs.c[0] = c1*bcperp.c[0];
			Cs.c[1] = c1*bcperp.c[1];
			Cs.c[2] = c1*bcperp.c[2];
		}
		else // bc < ab
		{
			c1 = (1.0/(2*tan(angle/2)*thisLen));
			c2 = prevLen/(8*sin(angle/2)*sin(angle/2));
		
			As.c[0] = -c2*(abperp.c[0]/(prevLen*prevLen));;
			As.c[1] = -c2*(abperp.c[1]/(prevLen*prevLen));;
			As.c[2] = -c2*(abperp.c[2]/(prevLen*prevLen));
					
			Bs.c[0] = c1*cb.c[0];
			Bs.c[1] = c1*cb.c[1];
			Bs.c[2] = c1*cb.c[2];
			
			Bs.c[0] += c2*( (ab.c[0]/(prevLen*prevLen)) + (bcperp.c[0]/(thisLen*thisLen)) );
			Bs.c[1] += c2*( (ab.c[1]/(prevLen*prevLen)) + (bcperp.c[1]/(thisLen*thisLen)) );
			Bs.c[2] += c2*( (ab.c[2]/(prevLen*prevLen)) + (bcperp.c[2]/(thisLen*thisLen)) );
			
			As.c[0] = (1.0/(2*tan(angle/2)*prevLen))*ba.c[0];
			As.c[1] = (1.0/(2*tan(angle/2)*prevLen))*ba.c[1];
			As.c[2] = (1.0/(2*tan(angle/2)*prevLen))*ba.c[2];
			c1 = (-prevLen)/(8*sin(angle/2)*sin(angle/2)*thisLen*thisLen);
			Cs.c[0] += c1*bcperp.c[0];
			Cs.c[1] += c1*bcperp.c[1];
			Cs.c[2] += c1*bcperp.c[2];			
		}
		
		// temporarily increment the strut's verts based on their component interactions
		// we undo this change at the end of the for loop in case the user
		// wants to keep the strut set
		int aVert, bVert, cVert;
		aVert = minradStruts[mItr].vert-1;
		bVert = minradStruts[mItr].vert;
		cVert = minradStruts[mItr].vert+1;

		aVert += inState->compOffsets[minradStruts[mItr].component];
		bVert += inState->compOffsets[minradStruts[mItr].component];
		cVert += inState->compOffsets[minradStruts[mItr].component];
		
		entry = (totalStruts*3*aVert)+contactStruts+mItr;
		if( entry+(2*totalStruts) > (3*inState->totalVerts)*totalStruts )
			printf( "*** crap!\n" );
	
		A[entry] = As.c[0];
		A[entry+totalStruts] = As.c[1];
		A[entry+(2*totalStruts)] = As.c[2];

		entry = (totalStruts*3*bVert)+contactStruts+mItr;
		if( entry+(2*totalStruts) > (3*inState->totalVerts)*totalStruts )
			printf( "*** crap!\n" );
	
		A[entry] = Bs.c[0];
		A[entry+totalStruts] = Bs.c[1];
		A[entry+(2*totalStruts)] = Bs.c[2];

		entry = (totalStruts*3*cVert)+contactStruts+mItr;
		if( entry+(2*totalStruts) > (3*inState->totalVerts)*totalStruts )
			printf( "*** crap!\n" );
	
		A[entry] = Cs.c[0];
		A[entry+totalStruts] = Cs.c[1];
		A[entry+(2*totalStruts)] = Cs.c[2];		
	}
}

static void
checkDuplicates( octrope_strut* struts, int num )
{
	int i, j;
	for( i=0; i<num; i++ )
	{
		// look for the variety of cases that signal a duplicate of strut struts[i] in the
		// rest of the list
		for( j=0; j<num; j++ )
		{
			if( j==i )
				continue;
		
			if( struts[i].component[0] == struts[j].component[0] &&
				struts[i].component[1] == struts[j].component[1] &&
				
				struts[i].lead_vert[0] == struts[j].lead_vert[0] &&
				struts[i].lead_vert[1] == struts[j].lead_vert[1] )
				
			//	fabs(struts[i].position[0] - struts[j].position[0]) < 0.1 &&
			//	fabs(struts[i].position[1] - struts[j].position[1]) < 0.1 )
			{
				printf( "test1: matching strut!\n" );
				exit(-1);
			}
			
			// different sense
			if( struts[i].component[0] == struts[j].component[0] &&
				struts[i].component[1] == struts[j].component[1] &&
				
				struts[i].lead_vert[0] == struts[j].lead_vert[1] &&
				struts[i].lead_vert[1] == struts[j].lead_vert[0] )
				
			//	fabs(struts[i].position[0] - struts[j].position[0]) < 0.1 &&
			//	fabs(struts[i].position[1] - struts[j].position[1]) < 0.1 )
			{
				printf( "test2: matching strut!\n" );
				exit(-1);
			}
			
		}
	}
}

static void
collapseStruts( octrope_strut** struts, int* count )
{
	// dumb method written to test whether almost-useless struts
	// are causing our ill-conditioned matrix problems
	int initialCount = *count;
	int newCount = 0;
	int sItr, compItr, good;
	octrope_strut* newSet = (octrope_strut*)malloc(sizeof(octrope_strut)*initialCount);
	
	for( sItr=0; sItr<initialCount; sItr++ )
	{
		good = 1;
		for( compItr=0; compItr<newCount; compItr++ )
		{
			if( (*struts)[sItr].component[0] == (newSet)[compItr].component[0] &&
				(*struts)[sItr].component[1] == (newSet)[compItr].component[1] &&
				
				(*struts)[sItr].lead_vert[0] == (newSet)[compItr].lead_vert[0] &&
				(*struts)[sItr].lead_vert[1] == (newSet)[compItr].lead_vert[1] )
			{
				good = 0;
				break;
			}
		}
		
		// not a dup, add to list with edge position rounded
		if( good == 1 )
		{
			newSet[newCount].component[0] = (*struts)[sItr].component[0];
			newSet[newCount].component[1] = (*struts)[sItr].component[1];
			
			newSet[newCount].lead_vert[0] = (*struts)[sItr].lead_vert[0];
			newSet[newCount].lead_vert[1] = (*struts)[sItr].lead_vert[1];
			
			newSet[newCount].position[0] = (int)(*struts)[sItr].position[0];
			newSet[newCount].position[1] = (int)(*struts)[sItr].position[1];
			
			newCount++;
		}
	}
	
	free(*struts);
	*struts = newSet;
	*count = newCount;
}

#define kMinradForceScale 0.1
static void
minradForce( octrope_vector* dlen, octrope_link* inLink, search_state* inState )
{
	int cItr, vItr;
	
	// note to self - this can be made much faster by not being dumb
	
	octrope_link_fix_wrap(inLink);
	
	// grab unit edges
    double* norms;
    //vector* units;
    int nextVert, totalVerts=0, dlItr;
    octrope_vector* diffVectors;
    octrope_vector* unitsCopy;
    
    totalVerts = 0;
	for( cItr=0; cItr<inLink->nc; cItr++ )
		totalVerts += inLink->cp[cItr].nv;
    
    norms = (double*)malloc(sizeof(double)*totalVerts);
    //units = (vector*)malloc(sizeof(struct vector_type)*totalVerts);
    diffVectors = (octrope_vector*)calloc(totalVerts, sizeof(struct octrope_vector_type));
    unitsCopy = (octrope_vector*)malloc(sizeof(struct octrope_vector_type)*totalVerts);
    
    fatalifnull_(norms);
    fatalifnull_(diffVectors);
    fatalifnull_(unitsCopy);
	
//	if( (rand() % 2) == 1 )
	{
		int i=0;
		inState->curvature_step = 0;
		for( cItr=0; cItr<inLink->nc; cItr++ )
		{		
			for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
			{
				nextVert = (vItr+1)%inLink->cp[cItr].nv;
				diffVectors[i].c[0] = inLink->cp[cItr].vt[nextVert].c[0] - inLink->cp[cItr].vt[vItr].c[0];
				diffVectors[i].c[1] = inLink->cp[cItr].vt[nextVert].c[1] - inLink->cp[cItr].vt[vItr].c[1];
				diffVectors[i].c[2] = inLink->cp[cItr].vt[nextVert].c[2] - inLink->cp[cItr].vt[vItr].c[2];
				
				// ddot is dot product
				norms[i] = sqrt(cblas_ddot(3, &diffVectors[i].c[0], 1, &diffVectors[i].c[0], 1));
				i++;
			}
		}
		
		// dscal is scale
		i=0;
		for( cItr=0; cItr<inLink->nc; cItr++ )
		{
			// effectively divide by norm
			for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++, i++ )
				cblas_dscal( 3, (1/norms[i]), &diffVectors[i].c[0], 1 );
		}
		
		// now diffVectors are forward units, sum with opposite of previous to get dLen field
		// duplicate the units since daxpy will operate in place
		memcpy( unitsCopy, diffVectors, sizeof(struct octrope_vector_type)*totalVerts );
		i=0;
		int prev;
		for( cItr=0; cItr<inLink->nc; cItr++ )
		{
			for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++)
			{
				// cblas_daxpy - adds a constant times a vector to another vector
				prev = (vItr == 0) ? (inLink->cp[cItr].nv - 1) : (vItr-1);
				prev += i;
				cblas_daxpy( 3, -1, &unitsCopy[prev].c[0], 1, &diffVectors[i+vItr].c[0], 1 );
			}
			i += inLink->cp[cItr].nv;
		
			// if this component is open, zero the end guys so things aren't screwy
			if( inLink->cp[cItr].acyclic != 0  )
			{
				diffVectors[i-inLink->cp[cItr].nv].c[0] = 0;
				diffVectors[i-inLink->cp[cItr].nv].c[1] = 0;
				diffVectors[i-inLink->cp[cItr].nv].c[2] = 0;
				diffVectors[i-1].c[0] = 0;
				diffVectors[i-1].c[1] = 0;
				diffVectors[i-1].c[2] = 0;
			}
		}
    }
//	else
//	{
//		printf( "*" );
//		inState->nocurvature_step = 1;
//	}
	
    free(norms);
    free(unitsCopy);
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		if( inState->conserveLength[cItr] != 0 )
		{
			for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
			{
				diffVectors[vItr+inState->compOffsets[cItr]].c[0] = 0;
				diffVectors[vItr+inState->compOffsets[cItr]].c[1] = 0;
				diffVectors[vItr+inState->compOffsets[cItr]].c[2] = 0;
			}
		}
	}
	
	double scaleFactor = 0.5;
	
	double* pushes = (double*)calloc(inState->totalVerts, sizeof(double));
	int pItr=0;			
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		// the first thing to do is grab edge lengths
		int edges;
		octrope_link* L = inLink;
		
		int i,j;
		double mr = {DBL_MAX}, alpha, rad;
		octrope_vector in,out;
		double normin, normout;
		double rplus, rminus, this_mr;
		int start,end;
		double target = inState->injrad;
		
		edges = octrope_pline_edges(&inLink->cp[cItr]);
								
		i = cItr;
		if (L->cp[i].acyclic) {
			start = 1;
			end = L->cp[i].nv - 1;
			out = octrope_vminus(L->cp[i].vt[1],L->cp[i].vt[0]);      
			normout = octrope_norm(out); 
		} else {
			start = 0;
			end = L->cp[i].nv;
			out = octrope_vminus(L->cp[i].vt[0],L->cp[i].vt[-1]);      
			normout = octrope_norm(out); 
		}

		/* Now we handle the main loop. */
		for (j = start; j < end; j++) 
		{
			in     = out;                
			normin = normout;               

			out = octrope_vminus(L->cp[i].vt[j+1],L->cp[i].vt[j]);
			normout = octrope_norm(out);

			alpha = acos(octrope_dot(in,out)/(normin*normout));
			rad   = 1/(2*tan(alpha/2));

			rminus = normin*rad;
			rplus  = normout*rad;

			this_mr = (rminus < rplus) ? rminus : rplus;
			mr = (mr < this_mr) ? mr : this_mr;
			
			pushes[inState->compOffsets[cItr]+j] = this_mr;

			if (this_mr < target) 
			{
				int prev, next;
				
				prev = j==0 ? end-1 : j-1;
				next = j==end ? 0 : j+1;
				
				// diffVectors is curvature force, so we want to add a little negative curvature force at prev and next
				// and a little positive curvature force at j
			/*	dlen[prev].c[0] += scaleFactor*-1*diffVectors[prev].c[0];
				dlen[prev].c[1] += scaleFactor*-1*diffVectors[prev].c[1];
				dlen[prev].c[2] += scaleFactor*-1*diffVectors[prev].c[2];
				
				dlen[next].c[0] += scaleFactor*-1*diffVectors[next].c[0];
				dlen[next].c[1] += scaleFactor*-1*diffVectors[next].c[1];
				dlen[next].c[2] += scaleFactor*-1*diffVectors[next].c[2];
				
				dlen[j].c[0] += scaleFactor*diffVectors[j].c[0];
				dlen[j].c[1] += scaleFactor*diffVectors[j].c[1];
				dlen[j].c[2] += scaleFactor*diffVectors[j].c[2];
			*/}
		}

	}	
	
	if( gOutputFlag )
	{
		// since this function exportes higher values as 'brighter', we make the minrad values 1/minrad
		for( pItr=0; pItr<inState->totalVerts; pItr++ )
		{
			if( fabs(pushes[pItr]) < 1e-5 )
				pushes[pItr] = 0;
			else
				pushes[pItr] = 1.0/pushes[pItr];
		}	
		export_pushed_edges( inLink, inState, pushes, "/tmp/curvature.vect", 1 );
	}
	free(pushes);
	
	free(diffVectors);
}

static void
spinForce( octrope_vector* dlen, octrope_link* inLink, search_state* inState )
{
	int cItr, vItr;
	
	// note to self - this can be made much faster by not being dumb
	
	octrope_link_fix_wrap(inLink);
		
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		// the first thing to do is grab edge lengths
		int edges;
		double* lengths;
		octrope_vector* sides;
		octrope_vector* adjustments;
		double  averageLength;
		
		edges = octrope_pline_edges(&inLink->cp[cItr]);
		
		lengths = (double*)malloc(sizeof(double)*edges);
		sides = (octrope_vector*)malloc(sizeof(octrope_vector)*edges);
		adjustments = (octrope_vector*)calloc(edges, sizeof(octrope_vector));
		
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
		{
			octrope_vector s1, s2;
			s1.c[0] = inLink->cp[cItr].vt[vItr].c[0];
			s1.c[1] = inLink->cp[cItr].vt[vItr].c[1];
			s1.c[2] = inLink->cp[cItr].vt[vItr].c[2];
			s2.c[0] = inLink->cp[cItr].vt[vItr+1].c[0];
			s2.c[1] = inLink->cp[cItr].vt[vItr+1].c[1];
			s2.c[2] = inLink->cp[cItr].vt[vItr+1].c[2];
			sides[vItr].c[0] = s2.c[0] - s1.c[0];
			sides[vItr].c[1] = s2.c[1] - s1.c[1];
			sides[vItr].c[2] = s2.c[2] - s1.c[2];
			
			lengths[vItr] = octrope_norm(sides[vItr]);
			
			// unit-ize sides since we will use these as the tangential motion basis below
			sides[vItr].c[0] /= lengths[vItr];
			sides[vItr].c[1] /= lengths[vItr];
			sides[vItr].c[2] /= lengths[vItr];
											
			averageLength += lengths[vItr];
		}
		averageLength /= edges;
					
		double spinFactor = 0;
		// fix zero and compute the tangential change necessary for the rest of the vertices
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
		{
			spinFactor = 0.5;
						
			adjustments[(vItr)].c[0] = (spinFactor)*sides[vItr].c[0];
			adjustments[(vItr)].c[1] = (spinFactor)*sides[vItr].c[1];
			adjustments[(vItr)].c[2] = (spinFactor)*sides[vItr].c[2];			
		
			dlen[((vItr)) + inState->compOffsets[cItr]].c[0] += adjustments[(vItr)].c[0];
			dlen[((vItr)) + inState->compOffsets[cItr]].c[1] += adjustments[(vItr)].c[1];
			dlen[((vItr)) + inState->compOffsets[cItr]].c[2] += adjustments[(vItr)].c[2];
		}
		
		free(lengths);
		free(sides);
		free(adjustments);
	}
	exportVect( dlen, inLink, "/tmp/adjustedDL.vect" );
}

void
specialForce( octrope_vector* dlen, octrope_link* inLink, search_state* inState )
{
	int cItr, vItr;
	
	// note to self - this can be made much faster by not being dumb
	
	octrope_link_fix_wrap(inLink);
		
/*	for( cItr=0; cItr<1; cItr++ )
	{
		// the first thing to do is grab edge lengths
		// fix zero and compute the tangential change necessary for the rest of the vertices
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
		{
						
		//	adjustments[(vItr)].c[0] = (spinFactor)*sides[vItr].c[0];
		//	adjustments[(vItr)].c[1] = (spinFactor)*sides[vItr].c[1];
		//	adjustments[(vItr)].c[2] = (spinFactor)*sides[vItr].c[2];			
					
			dlen[((vItr)) + inState->compOffsets[cItr]].c[0] += 0;
			dlen[((vItr)) + inState->compOffsets[cItr]].c[1] += -0.1;
			dlen[((vItr)) + inState->compOffsets[cItr]].c[2] += 0;
		}
	}
*/

	dlen[7].c[1] += -1;
	dlen[8].c[1] += -1;

	if( gOutputFlag )
		exportVect( dlen, inLink, "/tmp/adjustedDL.vect" );
}

static void
eqForce( octrope_vector* dlen, octrope_link* inLink, search_state* inState )
{
	int cItr, vItr;
	
	// note to self - this can be made much faster by not being dumb
	
	octrope_link_fix_wrap(inLink);
	
	if( gOutputFlag )
		exportVect( dlen, inLink, "/tmp/originalDL.vect" );
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		if( inState->conserveLength[cItr] != 0 )
			continue;
	
		// the first thing to do is grab edge lengths
		int edges;
		double* lengths;
		octrope_vector* sides;
		octrope_vector* adjustments;
		double  averageLength;
		double  usedLength;
		double  scaleFactor;
		
		edges = octrope_pline_edges(&inLink->cp[cItr]);
		
	//	printf( "edges: %d\n", edges );
		
		lengths = (double*)malloc(sizeof(double)*edges);
		sides = (octrope_vector*)malloc(sizeof(octrope_vector)*edges);
		adjustments = (octrope_vector*)calloc(edges, sizeof(octrope_vector));
		
		averageLength = 0;
		
		for( vItr=0; vItr<edges; vItr++ )
		{
			octrope_vector s1, s2;
			s1.c[0] = inLink->cp[cItr].vt[vItr].c[0];
			s1.c[1] = inLink->cp[cItr].vt[vItr].c[1];
			s1.c[2] = inLink->cp[cItr].vt[vItr].c[2];
			s2.c[0] = inLink->cp[cItr].vt[vItr+1].c[0];
			s2.c[1] = inLink->cp[cItr].vt[vItr+1].c[1];
			s2.c[2] = inLink->cp[cItr].vt[vItr+1].c[2];
			sides[vItr].c[0] = s2.c[0] - s1.c[0];
			sides[vItr].c[1] = s2.c[1] - s1.c[1];
			sides[vItr].c[2] = s2.c[2] - s1.c[2];
			
			lengths[vItr] = octrope_norm(sides[vItr]);
			
			// unit-ize sides since we will use these as the tangential motion basis below
			sides[vItr].c[0] /= lengths[vItr];
			sides[vItr].c[1] /= lengths[vItr];
			sides[vItr].c[2] /= lengths[vItr];
								
	//		printf( "length(%d): %lf\n", vItr, lengths[vItr] );
			
			averageLength += lengths[vItr];
		}
	//	printf( "total length: %lf\n", averageLength );
		averageLength /= edges;
	//	printf( "target: %lf\n", averageLength );
		
		usedLength = 0;
			
		int procVert;
		// fix zero and compute the tangential change necessary for the rest of the vertices
		for( vItr=0; vItr<edges-1; vItr++ )
		{
			procVert = (vItr+1)%edges;
			usedLength += lengths[vItr];
			scaleFactor = ((vItr+1))*averageLength - usedLength;
			
			// scale the side by the amount we need to move
		//	printf( "need to change %d by: %lf / position: %3.5lf v*avg: %3.5lf\n", procVert, scaleFactor, usedLength, ((vItr+1))*averageLength );
			
			scaleFactor *= inState->eqMultiplier;
			
			//scaleFactor = 1;
			adjustments[(vItr+1)%edges].c[0] = (scaleFactor)*sides[vItr].c[0];
			adjustments[(vItr+1)%edges].c[1] = (scaleFactor)*sides[vItr].c[1];
			adjustments[(vItr+1)%edges].c[2] = (scaleFactor)*sides[vItr].c[2];
			
			// adjustment should have same magnitude as unadjusted dlen vector at this vertex
		/*	adjustments[(vItr+1)%edges].c[0] /= octrope_norm(dlen[(vItr+1)%edges]);
			adjustments[(vItr+1)%edges].c[1] /= octrope_norm(dlen[(vItr+1)%edges]);
			adjustments[(vItr+1)%edges].c[2] /= octrope_norm(dlen[(vItr+1)%edges]);
		*/
		
		/*	adjustments[(vItr+1)%edges].c[0] /= octrope_norm(adjustments[(vItr+1)%edges]);
			adjustments[(vItr+1)%edges].c[1] /= octrope_norm(adjustments[(vItr+1)%edges]);
			adjustments[(vItr+1)%edges].c[2] /= octrope_norm(adjustments[(vItr+1)%edges]);
		*/	
		/*	adjustments[(vItr+1)%edges].c[0] *= octrope_norm(dlen[(vItr+1)%edges]);
			adjustments[(vItr+1)%edges].c[1] *= octrope_norm(dlen[(vItr+1)%edges]);
			adjustments[(vItr+1)%edges].c[2] *= octrope_norm(dlen[(vItr+1)%edges]);
		*/		
		//	printf( "adjustment %d: %lf %lf %lf\n", (vItr+1)%edges, adjustments[(vItr+1)%edges].c[0],
		//						adjustments[(vItr+1)%edges].c[1], adjustments[(vItr+1)%edges].c[2] );
		
			dlen[((vItr+1)%edges) + inState->compOffsets[cItr]].c[0] += adjustments[(vItr+1)%edges].c[0];
			dlen[((vItr+1)%edges) + inState->compOffsets[cItr]].c[1] += adjustments[(vItr+1)%edges].c[1];
			dlen[((vItr+1)%edges) + inState->compOffsets[cItr]].c[2] += adjustments[(vItr+1)%edges].c[2];
		}
		
		free(lengths);
		free(sides);
		free(adjustments);
	}
	
	if( gOutputFlag )
		exportVect( dlen, inLink, "/tmp/adjustedDL.vect" );
}

int gFeasibleThreshold = 0;

void
firstVariation( octrope_vector* dl, octrope_link* inLink, search_state* inState,
		octrope_strut** outStruts, int* outStrutsCount, int dlenStep )
{
	/*
	 * this function creates dVdt, the vector field under which we should actually be 
	 * moving at this point. We resolve dlen curvature force over the strut field as 
	 * returned by liboctrope.
	 */
	
	octrope_vector*		dVdt = NULL;
	int					strutStorageSize = 0;
	int					strutCount = 0;
	octrope_strut*		strutSet = NULL;
	octrope_mrloc*		minradSet = NULL;
	double				thickness = 2*inState->injrad;
	int					cItr, sItr; // loop iterators over components, verticies, and struts
	int					minradLocs;
	int					dlItr;
			
	double*				A = NULL; // the rigidity matrix
	
	if( outStrutsCount != NULL )
		*outStrutsCount = 0;
	if( outStruts != NULL )
		*outStruts = NULL;
	
	// grab dlen
	if( dl == NULL )
		dl = calloc(inState->totalVerts, sizeof(double));
	
	if( dlenStep != 0 )
		dlenForce(dl, inLink, inState);
	
//	dl = (octrope_vector*)calloc(inState->totalVerts, sizeof(struct octrope_vector_type));
		
	// acquire the strut set, we probably haven't increased the strut count much 
	// beyond inState.lastStrutCount, so it's probably safe to assume that we can fit
	// the strut set in 2*inState.lastStrutCount, but maybe not, so we loop over 
	// strut finder until we have stored less than we possibly can.
	strutStorageSize = (inState->lastStepStrutCount != 0) ? (2*inState->lastStepStrutCount) : 10;
	strutSet = NULL;
	do
	{
	
		if( strutSet != NULL )
			free(strutSet);
		if( minradSet != NULL )
			free(minradSet);
			
		strutSet = (octrope_strut*)calloc(strutStorageSize, sizeof(octrope_strut));
		minradSet = (octrope_mrloc*)malloc(strutStorageSize * sizeof(octrope_mrloc));
		octrope(	inLink, 
					DBL_MAX,		// epsilon is inf here
					thickness,		// all struts less than thickness
					1,  // factor 1
					&inState->ropelength,
					// minrad struts
					minradSet, strutStorageSize,
					// minrad info
					&inState->minrad,
					&minradLocs,
					// strut info
					strutSet,
					strutStorageSize,
					&inState->shortest,
					&strutCount,
					&inState->length,
					NULL, 0 );
					
		minradLocs = 0; // these actually aren't helpful
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
			minradLocs += octrope_pline_edges(&inLink->cp[cItr]);
		}
	}

	if( dlenStep != 0 )
	{
		eqForce(dl, inLink, inState);
		spinForce(dl, inLink, inState);
	}
//	specialForce(dl, inLink, inState);
//	//	minradForce(dl, inLink, inState);

	//specialForce(dl, inLink, inState);
	
	// if there are no struts, dVdt is completely controlled by dlen as there 
	// are no self contact constraints to worry about
	if( strutCount+minradLocs == 0 )
	{
		// we don't need this space if there are no struts
		free(strutSet);
		free(minradSet);
	
		// if there are no struts...
		//inState->shortest = 2*inState->injrad;

		//dVdt = dl;
		
		inState->avgDvdtMag = 0;
		for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
		{
			inState->avgDvdtMag += octrope_norm(dl[dlItr]);
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
		Acols = (3*inState->totalVerts)*(strutCount+minradLocs);
		A = (double*)calloc(Acols, sizeof(double)); // calloc zeros A
		// first, place contact struts
		placeContactStruts(A, inLink, strutSet, strutCount, inState, minradLocs);
		
//		placeVertexBars(A, inLink, strutCount, barVerts, minradLocs, inState);
		
		// and then minrad struts
		if( minradLocs > 0 )
		{
			placeMinradStruts(A, inLink, minradSet, minradLocs, inState, strutCount);
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
			
		// construct minusDL -- this is a column vector of size totalVerts
		// we must operate using strictly doubles to interoperate with taucs_snnls
		minusDL = (double*)calloc((3*inState->totalVerts), sizeof(double));
		if( dlenStep != 0 )
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
			for( sItr=0; sItr<strutCount; sItr++ )
			{
				if( strutSet[sItr].length < thickness-0.0001 )
//				if( strutSet[sItr].length < thickness )
				{
					/* get the norm, distribute force over verts */
					octrope_vector  ends[2];
					double		norm;
					int offset;
					octrope_strut_ends(inLink, &strutSet[sItr], ends);
					octrope_vsub(ends[0], ends[1]);
					norm = octrope_norm(ends[0]);
					ends[0].c[0] /= norm;
					ends[0].c[0] /= norm;
					ends[0].c[0] /= norm;
					
					offset = 3*(inState->compOffsets[strutSet[sItr].component[0]] + strutSet[sItr].lead_vert[0]);
					minusDL[offset + 0] += (1-strutSet[sItr].position[0])*(thickness-strutSet[sItr].length)*ends[0].c[0];
					minusDL[offset + 1] += (1-strutSet[sItr].position[0])*(thickness-strutSet[sItr].length)*ends[0].c[1];
					minusDL[offset + 2] += (1-strutSet[sItr].position[0])*(thickness-strutSet[sItr].length)*ends[0].c[2];
					
					if( strutSet[sItr].lead_vert[0] == (inLink->cp[strutSet[sItr].component[0]].nv-1) &&
						(inLink->cp[strutSet[sItr].component[0]].acyclic == 0) )
					{
						offset = (3*inState->compOffsets[strutSet[sItr].component[0]]) - 3;
					}
					minusDL[offset + 3] += (strutSet[sItr].position[1])*(thickness-strutSet[sItr].length)*ends[0].c[0];
					minusDL[offset + 4] += (strutSet[sItr].position[1])*(thickness-strutSet[sItr].length)*ends[0].c[1];
					minusDL[offset + 5] += (strutSet[sItr].position[1])*(thickness-strutSet[sItr].length)*ends[0].c[2];
					
					offset = 3*(inState->compOffsets[strutSet[sItr].component[1]] + strutSet[sItr].lead_vert[1]);
					minusDL[offset + 0] += (1-strutSet[sItr].position[0])*(thickness-strutSet[sItr].length)*-ends[0].c[0];
					minusDL[offset + 1] += (1-strutSet[sItr].position[0])*(thickness-strutSet[sItr].length)*-ends[0].c[1];
					minusDL[offset + 2] += (1-strutSet[sItr].position[0])*(thickness-strutSet[sItr].length)*-ends[0].c[2];
					
					if( strutSet[sItr].lead_vert[1] == (inLink->cp[strutSet[sItr].component[1]].nv-1) &&
						(inLink->cp[strutSet[sItr].component[1]].acyclic == 0) )
					{
						offset = (3*inState->compOffsets[strutSet[sItr].component[1]]) - 3;
					}
					minusDL[offset + 3] += (strutSet[sItr].position[1])*(thickness-strutSet[sItr].length)*-ends[0].c[0];
					minusDL[offset + 4] += (strutSet[sItr].position[1])*(thickness-strutSet[sItr].length)*-ends[0].c[1];
					minusDL[offset + 5] += (strutSet[sItr].position[1])*(thickness-strutSet[sItr].length)*-ends[0].c[2];
				}
			}
			
			for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
			{
				octrope_vector  vert;
				double  norm;
				vert.c[0] = minusDL[3*dlItr+0];
				vert.c[1] = minusDL[3*dlItr+1];
				vert.c[2] = minusDL[3*dlItr+2];
				norm = octrope_norm(vert);
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
			}
		
	/*		if( gOutputFlag == 1 )
			{
				octrope_vector  debug[500];
				for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
				{
					debug[dlItr].c[0] = 0.5*minusDL[3*dlItr+0];
					debug[dlItr].c[1] = 0.5*minusDL[3*dlItr+1];
					debug[dlItr].c[2] = 0.5*minusDL[3*dlItr+2];
				}
				exportVect(debug, inLink, "/tmp/minusDL.vect");
			}
	*/	
		}
							
		// solve AX = -dl, x is strut compressions
		taucs_ccs_matrix* sparseA = taucs_construct_sorted_ccs_matrix(A, strutCount+minradLocs, 3*inState->totalVerts);
	//	taucs_print_ccs_matrix(sparseA);
		//if( strutCount != 0 )
		
	/*	double* foobear = (double*)malloc(sizeof(double)*sparseA->n);
		dumpAxb( sparseA, foobear, minusDL );
		free(foobear);
	*/	
			compressions = t_snnls(sparseA, minusDL, &inState->residual, 2, 0);
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
		
		if( gOutputFlag == 1 && dlenStep != 0 )
			export_struts(inLink, strutSet, strutCount, compressions);

		// if we are graphing rcond, we should record it here
		if( inState->graphing[kRcond] != 0 )
		{
			taucs_ccs_matrix* apda = taucs_ccs_aprime_times_a(sparseA);
			inState->rcond = taucs_rcond(apda);
			taucs_ccs_free(apda);
		}
		
		if( compressions == NULL )
		{
			// we've probably hit the rcond wall, save as best.vect and bail
			exit(0);
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
		
		dVdt = (octrope_vector*)calloc(inState->totalVerts, sizeof(octrope_vector));
		// dVdt = dl + A*compressions is this for loop
		for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
		{
			double partialMult = 0;
			int totalStruts = strutCount + minradLocs;
			
			for( sItr=0; sItr<totalStruts; sItr++ )
				partialMult += compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+0];
			dVdt[dlItr].c[0] = dl[dlItr].c[0] + partialMult;
		
			partialMult=0;
			for( sItr=0; sItr<totalStruts; sItr++ )
				partialMult += compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+totalStruts];
			dVdt[dlItr].c[1] = dl[dlItr].c[1] + partialMult;
			
			partialMult = 0;
			for( sItr=0; sItr<totalStruts; sItr++ )
				partialMult += compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+(2*totalStruts)];
			dVdt[dlItr].c[2] = dl[dlItr].c[2] + partialMult;
		}
		
		if( gOutputFlag )
		{
			int next0, next1, pItr;
		
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
			
			if( gOutputFlag && dlenStep != 0 )
				export_pushed_edges(inLink, inState, totalPushes, "/tmp/compress_push.vect", 0);
			free(totalPushes);
		}
		
		// clean up
		taucs_ccs_free(sparseA);
		free(A);
		free(compressions);
		free(minusDL);
		// we only free this if the user wasn't interested in keeping it
		if( outStruts == NULL )
			free(strutSet);
		free(minradSet);
		// we free this here and not outside the else{} since dVdt _IS_ dl if
		// strutCount == 0, which is the other case.
		if( gOutputFlag == 1 )
			exportVect(dl, inLink, "/tmp/dl.vect");
		//free(dl);
		// throw dVdt back into dl, which on output, now includes dVdt
		for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
		{
			dl[dlItr].c[0] = dVdt[dlItr].c[0];
			dl[dlItr].c[1] = dVdt[dlItr].c[1];
			dl[dlItr].c[2] = dVdt[dlItr].c[2];
		}
		free(dVdt);
	}
	
//	return dVdt;
}

#pragma mark -

void
normalizeStruts( octrope_vector* strutDirections, octrope_strut* strutSet, octrope_link* inLink, int strutCount )
{
	int sItr;
	for( sItr=0; sItr<strutCount; sItr++ )
	{
		octrope_vector  points[2];
		double norm = 0;
		octrope_strut_ends( inLink, &strutSet[sItr], points );
		
		// the normalized difference of pointOne, pointTwo is the strut force vector
		strutDirections[sItr].c[0] = points[0].c[0] - points[1].c[0];
		strutDirections[sItr].c[1] = points[0].c[1] - points[1].c[1];
		strutDirections[sItr].c[2] = points[0].c[2] - points[1].c[2];
		
		norm = octrope_norm(strutDirections[sItr]);
		strutDirections[sItr].c[0] /= norm;
		strutDirections[sItr].c[1] /= norm;
		strutDirections[sItr].c[2] /= norm;
	}

}

double
maxovermin( octrope_link* inLink )
{
	int cItr, vItr;
	double max = 0, min = DBL_MAX;
	octrope_vector s1, s2;
	double max_maxovermin = 0;
	int edges;
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		max = 0;
		min = DBL_MAX;
		
		edges = octrope_pline_edges(&inLink->cp[cItr]);
		
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
			
			octrope_vsub(s1, s2);
			norm = octrope_norm(s1);
		
			if( norm < min )
				min = norm;
			if( norm > max )
				max = norm;
		}
		if( (max/min) > max_maxovermin )
			max_maxovermin = (max/min);
	}
	
	//printf( "max/min: %lf\n", max_maxovermin );
	return max_maxovermin;
}
