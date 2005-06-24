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
#include "lsqr.h"

#include "octrope_vector.h"

#include "stepper.h"
#include "dlen.h"
#include "eqedge.h"
#include "settings.h"
#include "../errors.h"

octrope_link*		bsearch_step( octrope_link* inLink, search_state* inState );
void				step( octrope_link* inLink, double stepSize, octrope_vector* dVdt, search_state* inState );
void				firstVariation( octrope_vector* inOutDvdt, octrope_link* inLink, search_state* inState,
						octrope_strut** outStruts, int* outStrutsCount, int dlenStep);
void				computeCompressPush( octrope_link* inLink, octrope_strut* strutSet,
						octrope_mrloc* minradSet, int strutCount, int minradLocs );

extern void sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod );						

// lesser utility functions
int					equalStruts( const octrope_strut* s1, const octrope_strut* s2 );
void				normalizeStruts( octrope_vector* strutDirections, 
							octrope_strut* strutSet, octrope_link* inLink, int strutCount );

void		placeVertexBars( double* A, octrope_link* inLink, int contactStruts, int totalBarVerts, int totalBars, search_state* inState );
int displayEveryFrame = 0;

extern int gQuiet;
extern int gSurfaceBuilding;
extern int gPaperInfoInTmp;

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
      fprintf(file,"%3.16lf %3.16lf %3.16lf \n", L->cp[i].vt[j].c[0], L->cp[i].vt[j].c[1],
                                  L->cp[i].vt[j].c[2]);
    }
  }

  /* . . . and the color data. */
  int colorItr = 0;
  for (i=0; i < L->nc; i++) {
    for (j=0; j < L->cp[i].nv; j++) {
		if( colorParam == 0 )
		  fprintf(file,"%3.16lf %3.16lf %3.16lf %3.16lf\n", 0.0, 1.0*(pushes[colorItr]/maxPush), 0.0, 1.0);
		else
		  fprintf(file,"%3.16lf %3.16lf %3.16lf %3.16lf\n", 1.0*(pushes[colorItr]/maxPush), 1.0*(pushes[colorItr]/maxPush), 1.0*(pushes[colorItr]/maxPush), 0.0);
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


int gOutputFlag = 1;
int gConditionCheck = 0;

extern int gFastCorrectionSteps;

#define kOutputItrs 1
#define SECS(tv)        (tv.tv_sec + tv.tv_usec / 1000000.0)

int	gEQevents = 0;

static void
placeMinradStruts2( double* rigidityA, octrope_link* inLink, octrope_mrloc* minradStruts, 
	int minradLocs, search_state* inState, int contactStruts )
{
	int mItr;
	int totalStruts = contactStruts + minradLocs;
	
	for( mItr=0; mItr<minradLocs; mItr++ )
	{
		octrope_vector B, A, cross, As, Bs, Cs;
		double	bmag, amag;
		double value, angle;
		double kappa, dot, prevLen, thisLen;
		octrope_vector  prevSide, thisSide, N, fancyL, fancyM, fancyN;
		
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

		// B = b-v
		B.c[0] = inLink->cp[cItr].vt[vItr+1].c[0] - inLink->cp[cItr].vt[vItr].c[0];
		B.c[1] = inLink->cp[cItr].vt[vItr+1].c[1] - inLink->cp[cItr].vt[vItr].c[1];
		B.c[2] = inLink->cp[cItr].vt[vItr+1].c[2] - inLink->cp[cItr].vt[vItr].c[2];
		
		// A = a-v
		A.c[0] = inLink->cp[cItr].vt[vItr-1].c[0] - inLink->cp[cItr].vt[vItr].c[0];
		A.c[1] = inLink->cp[cItr].vt[vItr-1].c[1] - inLink->cp[cItr].vt[vItr].c[1];
		A.c[2] = inLink->cp[cItr].vt[vItr-1].c[2] - inLink->cp[cItr].vt[vItr].c[2];
		
		bmag = octrope_norm(B);
		amag = octrope_norm(A);
		
		value = dot/(prevLen*thisLen);
		if( value >= 1 )
		{
			angle = 0;
		}
		else if( value <= -1 )
		{
			angle = M_PI;
		}
		else
		{
			// this is dangerous for the reasons of Ted's talk. Change me!
			angle = acos(value);
		}
		
		if( thisLen < prevLen )
		{
			// says... maple?
			kappa = -bmag/(2-2*cos(angle));
			
			// BxA
			cross.c[0] = B.c[1] * A.c[2] - B.c[2] * A.c[1];
			cross.c[1] = B.c[2] * A.c[0] - B.c[0] * A.c[2];
			cross.c[2] = B.c[0] * A.c[1] - B.c[1] * A.c[0];
			
			N.c[0] = cross.c[0]/octrope_norm(cross);
			N.c[1] = cross.c[1]/octrope_norm(cross);
			N.c[2] = cross.c[2]/octrope_norm(cross);
			
			double Lconst, Mconst, Nconst;
			
			Lconst = (1/(2*tan(angle/2) * bmag));
			fancyL.c[0] = Lconst * B.c[0];
			fancyL.c[1] = Lconst * B.c[1];
			fancyL.c[2] = Lconst * B.c[2];
			
			Mconst = kappa*(1/(amag*amag));
			// A x N
			cross.c[0] = A.c[1] * N.c[2] - A.c[2] * N.c[1];
			cross.c[1] = A.c[2] * N.c[0] - A.c[0] * N.c[2];
			cross.c[2] = A.c[0] * N.c[1] - A.c[1] * N.c[0];
			fancyM.c[0] = Mconst * cross.c[0];
			fancyM.c[1] = Mconst * cross.c[1];
			fancyM.c[2] = Mconst * cross.c[2];
			
			Nconst = kappa*(1/(bmag*bmag));
			// N x B
			cross.c[0] = N.c[1] * B.c[2] - N.c[2] * B.c[1];
			cross.c[1] = N.c[2] * B.c[0] - N.c[0] * B.c[2];
			cross.c[2] = N.c[0] * B.c[1] - N.c[1] * B.c[0];
			fancyN.c[0] = Nconst * cross.c[0];
			fancyN.c[1] = Nconst * cross.c[1];
			fancyN.c[2] = Nconst * cross.c[2];
			
			As.c[0] = fancyM.c[0];
			As.c[1] = fancyM.c[1];
			As.c[2] = fancyM.c[2];
			
			Bs.c[0] = -fancyM.c[0] - fancyN.c[0] - fancyL.c[0];
			Bs.c[1] = -fancyM.c[1] - fancyN.c[1] - fancyL.c[1];
			Bs.c[2] = -fancyM.c[2] - fancyN.c[2] - fancyL.c[2];
			
			Cs.c[0] = fancyN.c[0] + fancyL.c[0];
			Cs.c[1] = fancyN.c[1] + fancyL.c[1];
			Cs.c[2] = fancyN.c[2] + fancyL.c[2];
		}
		else
		{
			// says... maple?
			kappa = -amag/(2-2*cos(angle));
			
			// BxA
			cross.c[0] = B.c[1] * A.c[2] - B.c[2] * A.c[1];
			cross.c[1] = B.c[2] * A.c[0] - B.c[0] * A.c[2];
			cross.c[2] = B.c[0] * A.c[1] - B.c[1] * A.c[0];
			
			N.c[0] = cross.c[0]/octrope_norm(cross);
			N.c[1] = cross.c[1]/octrope_norm(cross);
			N.c[2] = cross.c[2]/octrope_norm(cross);
			
			double Lconst, Mconst, Nconst;
			
			Lconst = (1/(2*tan(angle/2) * amag));
			fancyL.c[0] = Lconst * A.c[0];
			fancyL.c[1] = Lconst * A.c[1];
			fancyL.c[2] = Lconst * A.c[2];
			
			Mconst = kappa*(1/(amag*amag));
			// A x N
			cross.c[0] = A.c[1] * N.c[2] - A.c[2] * N.c[1];
			cross.c[1] = A.c[2] * N.c[0] - A.c[0] * N.c[2];
			cross.c[2] = A.c[0] * N.c[1] - A.c[1] * N.c[0];
			fancyM.c[0] = Mconst * cross.c[0];
			fancyM.c[1] = Mconst * cross.c[1];
			fancyM.c[2] = Mconst * cross.c[2];
			
			Nconst = kappa*(1/(bmag*bmag));
			// N x B
			cross.c[0] = N.c[1] * B.c[2] - N.c[2] * B.c[1];
			cross.c[1] = N.c[2] * B.c[0] - N.c[0] * B.c[2];
			cross.c[2] = N.c[0] * B.c[1] - N.c[1] * B.c[0];
			fancyN.c[0] = Nconst * cross.c[0];
			fancyN.c[1] = Nconst * cross.c[1];
			fancyN.c[2] = Nconst * cross.c[2];
			
			As.c[0] = fancyM.c[0] + fancyL.c[0];;
			As.c[1] = fancyM.c[1] + fancyL.c[1];;
			As.c[2] = fancyM.c[2] + fancyL.c[2];;
			
			Bs.c[0] = -fancyM.c[0] - fancyN.c[0] - fancyL.c[0];
			Bs.c[1] = -fancyM.c[1] - fancyN.c[1] - fancyL.c[1];
			Bs.c[2] = -fancyM.c[2] - fancyN.c[2] - fancyL.c[2];
			
			Cs.c[0] = fancyN.c[0]; //+ fancyL.c[0];
			Cs.c[1] = fancyN.c[1]; //+ fancyL.c[1];
			Cs.c[2] = fancyN.c[2]; //+ fancyL.c[2];
		}
		
		double norm;
		int nItr;
		norm = 0;
		for( nItr=0; nItr<3; nItr++ )
		{
			norm += As.c[nItr]*As.c[nItr];
			norm += Bs.c[nItr]*Bs.c[nItr];
			norm += Cs.c[nItr]*Cs.c[nItr];
		}
		norm = sqrt(norm);
		
		// if we're doing fast tsnnls based correction steps, the normalization 
		// doesn't really matter
		if( inState->curvature_step != 0 || gFastCorrectionSteps != 0 )
		{
			As.c[0] /= norm;
			As.c[1] /= norm;
			As.c[2] /= norm;
			Bs.c[0] /= norm;
			Bs.c[1] /= norm;
			Bs.c[2] /= norm;
			Cs.c[0] /= norm;
			Cs.c[1] /= norm;
			Cs.c[2] /= norm;
		}
					
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
		
		if( minradStruts[mItr].vert == 0 &&
			(inLink->cp[minradStruts[mItr].component].acyclic == 0) )
		{
			entry = (totalStruts*3)*(inState->compOffsets[minradStruts[mItr].component]);
		}
		else
		{
			entry = (totalStruts*3*aVert)+contactStruts+mItr;
		}

		if( entry+(2*totalStruts) > (3*inState->totalVerts)*totalStruts )
			printf( "*** crap!\n" );
	
		rigidityA[entry] = As.c[0];
		rigidityA[entry+totalStruts] = As.c[1];
		rigidityA[entry+(2*totalStruts)] = As.c[2];

		entry = (totalStruts*3*bVert)+contactStruts+mItr;
		if( entry+(2*totalStruts) > (3*inState->totalVerts)*totalStruts )
			printf( "*** crap!\n" );
	
		rigidityA[entry] = Bs.c[0];
		rigidityA[entry+totalStruts] = Bs.c[1];
		rigidityA[entry+(2*totalStruts)] = Bs.c[2];


		if( minradStruts[mItr].vert+1 == (inLink->cp[minradStruts[mItr].component].nv) &&
			(inLink->cp[minradStruts[mItr].component].acyclic == 0) )
		{
			entry = (totalStruts*3)*(inLink->cp[minradStruts[mItr].component].nv-1 + inState->compOffsets[minradStruts[mItr].component]);
		}
		else
		{
			entry = (totalStruts*3*cVert)+contactStruts+mItr;
		}
		if( entry+(2*totalStruts) > (3*inState->totalVerts)*totalStruts )
			printf( "*** crap!\n" );
	
		rigidityA[entry] = Cs.c[0];
		rigidityA[entry+totalStruts] = Cs.c[1];
		rigidityA[entry+(2*totalStruts)] = Cs.c[2];	
		
	} // for over minrad struts
	
}

static void
glom( octrope_link* inLink, search_state* inState )
{
	octrope_mrloc strut;
	
	strut.component = 0;
	strut.vert = 1;
	strut.mr = octrope_minradval(inLink);

	int Acols = (3*3)*(1);
	double* A = (double*)calloc(Acols, sizeof(double)); // calloc zeros A
	
	double* blah = (double*)calloc(Acols, sizeof(double));
	
	double dpsum = 0, dx;
	double diff;
	int i, mItr;
	
	for( i=0; i<1000; i++ )
	{
		octrope_vector motion[3];
		octrope_vector foo;
		double initialMR;
		
		dpsum = 0;
		
		for( mItr=0; mItr<3; mItr++ )
		{
			motion[mItr].c[0] = (drand48()-0.5)*2;
			motion[mItr].c[1] = (drand48()-0.5)*2;
			motion[mItr].c[2] = (drand48()-0.5)*2;
		}
		
		initialMR = octrope_minradval(inLink);
		
		placeMinradStruts2(A, inLink, &strut, 1, inState, 0);
		
		foo.c[0] = A[0];
		foo.c[1] = A[1];
		foo.c[2] = A[2];
		dpsum += octrope_dot(foo,motion[0]);
		foo.c[0] = A[3];
		foo.c[1] = A[4];
		foo.c[2] = A[5];
		dpsum += octrope_dot(foo,motion[1]);
		foo.c[0] = A[6];
		foo.c[1] = A[7];
		foo.c[2] = A[8];
		dpsum += octrope_dot(foo,motion[2]);
		
		for( dx=1e-1; dx>1e-10; dx *= 0.1 )
		{
			octrope_link* workerLink = octrope_link_copy(inLink);
			octrope_link_fix_wrap(workerLink);
			step( workerLink, dx, motion, inState );
			
			// if everything's right: (minrad of workerLink) - (minrad of inLink) <--- ~= dpsum*dx. we'll see!
			diff = octrope_minradval(workerLink) - octrope_minradval(inLink);
			printf( "diff: %e / dpsum*dx: %e / (diff - dpsum*dx): %e / still!!: %e\n", 
				diff, dpsum*dx, diff-(dpsum*dx), (diff-(dpsum*dx))/dx );
			
			octrope_link_free(workerLink);
		}
	}
	
	free(A);
}

extern int gSuppressOutput;
extern int gVerboseFiling;

int	gCorrectionAttempts = 0;

void 
bsearch_stepper( octrope_link** inLink, search_state* inState )
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
			
			// we only use gnuplot if we have a display, otherwise
			// just create data files
			if( getenv("DISPLAY") != NULL )
			{
				sprintf( fname, "%dpipe", i );
				sprintf( cmd, "rm -f %s ; mkfifo %s", fname, fname );
				system(cmd);
				
				// get gnuplot going and reading our commands
				sprintf( cmd, "/sw/bin/gnuplot < %s &", fname );
				system(cmd);
				
				// open the pipe for our input
				gnuplotPipes[i] = fopen(fname, "w");
			
				fprintf( gnuplotPipes[i], "set term x11\n" );
			}
		
			sprintf(fname, "%ddat", i);
			gnuplotDataFiles[i] = fopen(fname, "w");
		} 
	}
	
	int stepItr;
	struct rusage startStepTime, stopStepTime;
	
	for( stepItr=0; 1==1; stepItr++ )
	{
		int lastSet;
		
		if( stepItr > inState->maxItrs && inState->maxItrs > 0 && inState->curvature_step == 1 )
		{
			FILE* maxFile = NULL;
			char    fname[512];
			char adjustedName[512];
			
			strcpy(adjustedName, inState->fname);
			sprintf( fname, "%s_%d.maxItr", adjustedName, inState->totalVerts );
			maxFile = fopen(fname, "w");
			
			octrope_link_write(maxFile, *inLink);
			fclose(maxFile);

			break;
		}
			
		lastSet = inState->lastStepStrutCount;
			
		//if( (i%50)==0 )
		
	//	inState->curvature_step = 1;
		
		// we need to eq if we're below threshold or if last step hasn't brough us back 
		// within an acceptable range
	/*	if( inState->eqThreshold <= inState->lastMaxMin || 
			(inState->eq_step == 1 && inState->lastMaxMin > (1+((inState->eqThreshold - 1)/2))) )
		{
			inState->eq_step = 1;
			inState->curvature_step = 0;
			
			FILE* verts = fopen("/tmp/verts.vect", "w");
			octrope_link_draw(verts, *inLink);
			fclose(verts);
			
			printf( "" );
		}
		else*/
			inState->eq_step = 0;
			
		if( inState->curvature_step == 0 && gFastCorrectionSteps == 0 )
		{
			if( gQuiet == 0 )
				printf( "last correction step lsqr residual: %e\n", inState->ofvResidual );
		}

		if( (inState->shortest < ((2*inState->injrad)-(inState->overstepTol)*(2*inState->injrad))) ||
			(inState->minrad < inState->minradOverstepTol && !inState->ignore_minrad) )
		{
			// reset counter, this is our first correction attempt
			if( inState->curvature_step != 0 )
			{
				// before reset, output last number if we're recording this
				if( gPaperInfoInTmp != 0 )
				{
					char	fname[512];
					preptmpname(fname,"correction_convergence",inState);
					FILE* tcc = fopen(fname,"a");
					fprintf(tcc, "%d\n", gCorrectionAttempts);
					fclose(tcc);
				}
				
				gCorrectionAttempts = 1;
			}
					
			inState->curvature_step = 0;
		}
		else
		{
			double thickness = 2*inState->injrad;
			double greenZone = thickness - (thickness*inState->overstepTol)*0.5;
			if( (inState->curvature_step == 0 && inState->shortest < greenZone	) ) //||
//				(inState->curvature_step == 0 && inState->minrad < 0.49999 && inState->ignore_minrad==0) )
			{
				// we haven't finished yet, record
				gCorrectionAttempts++;
				
				inState->curvature_step = 0;
			}
			else
				inState->curvature_step = 1;
			
			if( inState->eq_step == 0 )
				cSteps++;
			else
				inState->curvature_step = 0;
		}
		
																		
	//	inState->curvature_step = ((inState->curvature_step+1)%2);
	
//		inState->maxStepSize = 1e-5;
			
	/*		double user;
			
			// grab time for this frame and reset counter
			getrusage(RUSAGE_SELF, &startStepTime);
		*/
		
		if( (stepItr%50)==0 || gVerboseFiling != 0 )
			gOutputFlag = 1;
		
		if( gSuppressOutput == 1 )
			gOutputFlag = 0;
				
		*inLink = bsearch_step(*inLink, inState);
		
/*		if( inState->minrad < 0.49 && inState->ignore_minrad == 0 )
		{
			fprintf(stderr, "minrad was probably going to crumble things\n" );
			error_write(inLink);
			exit(-1);
		}
*/
		/*	getrusage(RUSAGE_SELF, &stopStepTime);
			user = SECS(stopStepTime.ru_utime) - SECS(startStepTime.ru_utime);
			printf( "STEP TIME (user): %f strts: %d\n", user, inState->lastStepStrutCount );
		*/

		inState->steps++;
		
	//	gOutputFlag = 1;
				
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
		
		if( (stepItr%kOutputItrs) == 0 && gQuiet == 0 )
		{
			printf( "s: %d ms: %d len: %lf r: %lf ssize: %e dcsd: %lf minrad: %lf rsdl: %e t: %lf\n", 
						inState->lastStepStrutCount, inState->lastStepMinradStrutCount,
						inState->length, 2*inState->ropelength, inState->stepSize, inState->shortest, inState->minrad, 
						inState->residual, inState->time );
		}
		
				
		if( inState->oldLengthTime == 0 )
		{
			inState->oldLengthTime = inState->cstep_time;
			inState->oldLength = inState->length;
		}
			
		if( (inState->oldLengthTime + inState->checkDelta < inState->cstep_time) && 
			fabs(inState->oldLength-inState->length) < inState->checkThreshold &&
			inState->residualThreshold >= inState->residual
			/*((double)cSteps)/((double)stepItr+1) > 0.05*/ )
		{
			// if we're going to die, it'll be after this, so save our best
			FILE* bestFile = NULL;
			char    fname[512];
			char adjustedName[512];
			
			strcpy(adjustedName, inState->fname);
			sprintf( fname, "%s_%d.best", adjustedName, inState->totalVerts );
			bestFile = fopen(fname, "w");
			
			octrope_link_write(bestFile, *inLink);
			fclose(bestFile);
		
			if( inState->totalVerts > inState->refineUntil*(2*inState->ropelength) )
			{
				// we are DONE!
				break;
			}
			
			// the strut set can get violent after this change, so let's make sure none exist.
			// it's a lot easier on the condition number if it reemerges quickly than if it's changing
			// rapidly
			link_scale(*inLink, /*1.05**/(2.0*inState->injrad)/inState->shortest);
		
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
		
		if( (inState->oldLengthTime + inState->checkDelta) < inState->cstep_time )
		{
			inState->oldLength = inState->length;
			inState->oldLengthTime = inState->cstep_time;
			if( gQuiet == 0 )
				printf( "* Checked delta rope and continuing\n" );
		}
		
		double maxmin;
		maxmin = maxovermin(*inLink, inState);
		inState->lastMaxMin = maxmin;
		if( maxmin > maxmaxmin )
			maxmaxmin = maxmin;
	//	inState->eqMultiplier = 1;
		
/*		if( maxmin > 1.0001 )
			inState->eq_step = 1;
		else
			inState->eq_step = 0;
*/		
		if( (stepItr%kOutputItrs) == 0 && gQuiet == 0 )
		{
			printf( "mm: %3.5lf eqM: %3.5lf (max mm: %3.5lf) (min thick: %3.5lf) lrcond: %e cstep: %d eqAvgDif: %e\n", 
						maxmin, inState->eqMultiplier, maxmaxmin, minthickness, inState->rcond, inState->curvature_step,
						inState->eqAvgDiff );
			printf( "cstep/step ratio: %lf delta length: %lf next check: %lf check threshold: %lf\n", 
					((double)cSteps)/((double)stepItr+1), fabs(inState->oldLength-inState->length),
					inState->oldLengthTime+inState->checkDelta, inState->checkThreshold );
		}

		if( inState->time >= nextMovieOutput )
		{
			char	fname[512];
			FILE*   frame = NULL;
			
			struct rusage stopTime;
			double user;
			
			// grab time for this frame and reset counter
			getrusage(RUSAGE_SELF, &stopTime);
			user = SECS(stopTime.ru_utime) - SECS(inState->frameStart.ru_utime);
			if( gQuiet == 0 )
				printf( "FRAME TIME (user): %f strts: %d\n", user, inState->lastStepStrutCount );
			else
			{
				// same stats when quiet, but not just once / viz output
				printf( "s: %d ms: %d len: %lf r: %lf ssize: %e dcsd: %lf minrad: %lf maxmin: %lf residual: %e time: %lf cstep ratio: %lf\n", 
						inState->lastStepStrutCount, inState->lastStepMinradStrutCount,
						inState->length, 2*inState->ropelength, inState->stepSize, inState->shortest, inState->minrad, 
						maxmin, inState->residual, inState->time, ((double)cSteps)/((double)stepItr+1) );
			}
			getrusage(RUSAGE_SELF, &inState->frameStart);
			
			if( inState->saveConvergence != 0 )
			{
				preptmpname(fname, "rrconvergence.txt", inState);
				FILE* conv = fopen(fname,"a");
				fprintf( conv, "%3.14lf %3.14lf\n", inState->time, max(2*inState->length,2*inState->ropelength) );
				// flushes, most importantly
				fclose(conv);
			}
			
		//	nextMovieOutput += 0.05;
			// make things 24 fps
			nextMovieOutput += 0.041666666667;
					
			sprintf( fname, "restart_%s", inState->fname );
		//	(strstr(fname,".vect"))[0] = '\0';
			if( gQuiet == 0 )
				printf( "saved restart: %s\n", fname );
			frame = fopen( fname, "w" );
			octrope_link_write(frame, *inLink);
			fclose(frame);
			
			if( inState->movie != 0 )
			{
				sprintf( fname, "movie%s/rmov.%lf_evals-%d_strts-%d_length-%lf.vect", 
								inState->fname,
								inState->time, 
								inState->tsnnls_evaluations, 
								inState->lastStepStrutCount,
								inState->length );
				frame = fopen( fname, "w" );
				octrope_link_write(frame, *inLink);
				fclose(frame);
				
				char cmd[1024];
				if( inState->lastStepStrutCount != 0 )
				{
					char	foobear[1024];
					preptmpname(foobear,"struts.vect",inState);
					sprintf(cmd, "cp %s movie%s/struts.%lf.vect", foobear, inState->fname, inState->time);
					system(cmd);
				}
				
				if( gQuiet == 0 )
					printf( "movie frame output (tsnnls evals: %d)\n", inState->tsnnls_evaluations );

			}
			
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
						case kRopelength: fprintf( gnuplotDataFiles[i], "%lf", 2*inState->ropelength ); break;
						case kStrutCount: fprintf( gnuplotDataFiles[i], "%d", inState->lastStepStrutCount ); break;
						case kStepSize: fprintf( gnuplotDataFiles[i], "%lf", inState->stepSize ); break;
						case kThickness: fprintf( gnuplotDataFiles[i], "%lf", inState->shortest ); break;
						case kMinrad: fprintf( gnuplotDataFiles[i], "%lf", inState->minrad ); break;
						case kResidual: fprintf( gnuplotDataFiles[i], "%lf", inState->residual ); break;
						case kMaxOverMin: fprintf( gnuplotDataFiles[i], "%lf", maxmin ); break;
						case kRcond: fprintf( gnuplotDataFiles[i], "%lf", inState->rcond ); break;
						case kWallTime: fprintf( gnuplotDataFiles[i], "%lf", (clock() - startTime)/((double)CLOCKS_PER_SEC) ); break;
						case kMaxVertexForce: fprintf( gnuplotDataFiles[i], "%lf", (inState->maxPush)/(inState->length/(double)inState->totalVerts) ); break;
						case kEQVariance: fprintf( gnuplotDataFiles[i], "%3.15lf %3.15lf", inState->eqVariance, inState->eqAvgDiff ); break;
					}
					fprintf( gnuplotDataFiles[i], "\n" );
					fflush(gnuplotDataFiles[i]);
					// show only last 5 seconds
			//		fprintf( gnuplotPipes[i], "set xrange [%lf:%lf]\n", inState->time-5, inState->time );
					if( i == kLength )
					{
						if( getenv("DISPLAY") != NULL )
							fprintf( gnuplotPipes[i], "set xrange [%lf:%lf]\n", inState->time-1, inState->time );
					}
					
					if( getenv("DISPLAY") != NULL )
					{
						if( i != kEQVariance )
							fprintf( gnuplotPipes[i], "plot \"%s\" using 1:2 w lines title \'", fname );
						else
							fprintf( gnuplotPipes[i], "plot \"%s\" u 1:3 w lines, \"%s\" using 1:2 w lines title \'", fname, fname );
						
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
							case kRcond: fprintf( gnuplotPipes[i], "reciprocal condition number of A (rigidity matrix)" ); break;
							case kWallTime: fprintf( gnuplotPipes[i], "process computation time" ); break;
							case kMaxVertexForce: fprintf( gnuplotPipes[i], "maximum compression sum" ); break;
							case kConvergence: fprintf( gnuplotPipes[i], "convergence (rope vs steps)" ); break;
							case kEQVariance: fprintf( gnuplotPipes[i], "edge length variance from average" ); break;
						}
						fprintf( gnuplotPipes[i], "\'\n" ); 
						fflush( gnuplotPipes[i] );
					} // if $DISPLAY defined
				}
			} // for graph outputs
		}
		
		if( gOutputFlag == 1 && gSuppressOutput == 0 )
		{
			char	fname[1024];
			preptmpname(fname,"verts.vect",inState);
			FILE* verts = fopen(fname, "w");
			octrope_link_draw(verts, *inLink);
			fclose(verts);
		}

																
		// eq ourselves
	//	maxovermin(*inLink);
		//lapack_eqedge( *inLink, 4 );
	
		//lapack_eqedge(*inLink, 4);		
				
	//	maxovermin(*inLink);
	}
	
	/* inititalize piping if we are to use it */
	for( i=0; i<kTotalGraphTypes; i++ )
	{
		if( inState->graphing[i] != 0 )
		{
			char	cmd[512];
			char	fname[255];
			
			if( getenv("DISPLAY") != NULL )
			{
				sprintf( fname, "%dpipe", i );
				// this will quit gnuplot
				sprintf( cmd, "echo quit >> %s", fname );
				system(cmd);
				fclose(gnuplotPipes[i]);
				sprintf(cmd, "rm -f %s", fname);
				system(cmd);
			} 
			
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
step( octrope_link* inLink, double stepSize, octrope_vector* dVdt, search_state* inState )
{
	int cItr, vItr, dVdtItr;
	for( cItr=0, dVdtItr=0; cItr<inLink->nc; cItr++ )
	{
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++, dVdtItr++ )
		{
			/* Since our eq strategy isn't perfect, and we have non-eq theory
			 * we use it here. Scale force at point by the deviation from the 
			 * average edge length of the averge of adjacent edges
			 */
			double avg, scale=1.0;
	/*		avg = (inState->sideLengths[inState->compOffsets[cItr] + vItr] + inState->sideLengths[inState->compOffsets[cItr] + vItr+1])/2.0;
			scale = avg / inState->avgSideLength;
	*/		
			inLink->cp[cItr].vt[vItr].c[0] += scale*stepSize*dVdt[dVdtItr].c[0];
			inLink->cp[cItr].vt[vItr].c[1] += scale*stepSize*dVdt[dVdtItr].c[1];
			inLink->cp[cItr].vt[vItr].c[2] += scale*stepSize*dVdt[dVdtItr].c[2];
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
	char	fname[1024];
	preptmpname(fname,"dVdt.vect",inState);
	exportVect(dVdt, inLink, fname);

	free(A);
	free(b);
}

static void
normalizeVects( octrope_vector* dvdt, int size )
{
	double sum = 0;
	int i;
	for( i=0; i<size; i++ )
	{
		sum += dvdt[i].c[0]*dvdt[i].c[0];
		sum += dvdt[i].c[1]*dvdt[i].c[1];
		sum += dvdt[i].c[2]*dvdt[i].c[2];
	}
	
	sum = sqrt(sum);
	
	for( i=0; i<size; i++ )
	{
		dvdt[i].c[0] = 10*(dvdt[i].c[0]/sum);
		dvdt[i].c[1] = 10*(dvdt[i].c[1]/sum);
		dvdt[i].c[2] = 10*(dvdt[i].c[2]/sum);
	}
}

octrope_link*
bsearch_step( octrope_link* inLink, search_state* inState )
{	
	int stepAttempts;
	int dvdtItr;
	double  lastDCSD, lastMR, eps = inState->stepSize*inState->stepSize;
	octrope_link* workerLink;
	
//	double ERROR_BOUND = ((inState->overstepTol)*(2*inState->injrad))/50;
	double ERROR_BOUND = 1e-5;
//	double MR_ERROR_BOUND = ((inState->overstepTol)*inState->injrad)/20;
	double MR_ERROR_BOUND = 1e-5;

	// create initial vector field for which we want to move along in this step
	octrope_vector* dVdt;
	
	dVdt = calloc(inState->totalVerts, sizeof(octrope_vector));
	
	struct rusage stopTime;
	struct rusage startTime;
	double user;

	if( gPaperInfoInTmp != 0 )
	{
		// keep track of linear algebra runtime
		getrusage(RUSAGE_SELF, &startTime);
	}

	firstVariation(dVdt, inLink, inState, NULL, &inState->lastStepStrutCount, inState->curvature_step);
	
	if( gPaperInfoInTmp != 0 )
	{
		getrusage(RUSAGE_SELF, &stopTime);
		user = SECS(stopTime.ru_utime) - SECS(startTime.ru_utime);
		if( inState->curvature_step == 0 )
			printf( "Correction step firstVartion() time: %lf\n", user );
		else
			printf( "Curvature step firstVartion() time: %lf\n", user );
	}
	// this actually onlt applies with newton stepping without lsqr, which we don't do anymore.
	// not that we ever did, but the code will be left here if i ever come back to it
/*	if( inState->curvature_step == 0 && gFastCorrectionSteps == 0 )
	{
		inState->minrad = octrope_minradval(inLink);
		inState->shortest = octrope_poca(inLink, NULL, 0);
		inState->length = octrope_curvelength(inLink);

		free(dVdt);

		return inLink; // we updated the verts from the lsqr call in firstVariation
	}
*/	
	
	if( inState->shortest > 2*inState->injrad )
		inState->shortest = 2*inState->injrad;
			
//	lastDCSD = octrope_poca(inLink, NULL, 0);
	lastDCSD = inState->shortest;
	lastMR = inState->minrad;

	stepAttempts = 0;
	workerLink = NULL;
	double curr_error = 0, mr_error=0, newpoca = 0, newmr=0, correctionStepSize=inState->correctionStepDefault;
	short	improvedNorm = 0;
	
//	correctionStepSize = 1e-7;
	
	if( inState->curvature_step != 0 )
	{
		normalizeVects(dVdt, inState->totalVerts);
	}
	
//	if( lastMR < 0.5 )
//		correctionStepSize = 0.25;

	if( gSurfaceBuilding )
		correctionStepSize = 0.25;
	
	double oldLength = inState->length;
	double workingLength = 0;
	
	do
	{			
		stepAttempts++;
		
	/*	if( inState->stepSize == 0 )
		{
			inState->stepSize = kMinStepSize;
			break;
		}
	*/	
		// move along dVdt
		if( workerLink != NULL )
			octrope_link_free(workerLink);
		workerLink = octrope_link_copy(inLink);
		octrope_link_fix_wrap(workerLink);
								
		if( inState->curvature_step != 0 )
		{
			step(workerLink, inState->stepSize, dVdt, inState);
			workingLength = octrope_curvelength(workerLink);
		}
		else
		{
			step(workerLink, correctionStepSize, dVdt, inState);
		}

		if( gOutputFlag == 1 )
		{
			char fname[1024];
			preptmpname(fname,"dVdt.vect", inState);
			exportVect(dVdt, inLink, fname);
			
			if( inState->fancyVisualization != 0 )
			{
				//fprintf(inState->fancyPipe, "(load %s)\n",fname);
				preptmpname(fname,"struts.vect",inState);
				fprintf(inState->fancyPipe, "(load %s)\n", fname);
				fflush(inState->fancyPipe);
			}
			if( gQuiet == 0 )
				printf( "visualization output\n" );
			gOutputFlag = 0;
		}

		newmr = octrope_minradval(workerLink);
		newpoca = octrope_poca(workerLink, NULL, 0);
		
		if( isnan(newmr) || isinf(newmr) )
		{
			fprintf(stderr, "error write due to nan minrad!\n");
			error_write(workerLink);
			
			newmr = octrope_minradval(workerLink);
			printf( "%lf\n", newmr );
			
			exit(-1);
		}
		
		inState->minrad = newmr;
		inState->shortest = newpoca;
		
		if( inState->curvature_step != 0 )
		{
			// check our new thickness
			curr_error = max( lastDCSD-newpoca, 0 );
			if( newmr < inState->injrad )
				mr_error = max(lastMR - newmr, 0);
			else
				mr_error = 0;
			
	/*		if( workingLength > oldLength )
				curr_error = 1;
	*/			
			if( curr_error < ERROR_BOUND && (mr_error < MR_ERROR_BOUND || inState->ignore_minrad) )
				inState->stepSize *= 2;
			else
				inState->stepSize /= 2;		
		}
		
		if( mr_error > MR_ERROR_BOUND && !inState->ignore_minrad )
			curr_error = 1;
					
	} while( (curr_error > ERROR_BOUND ) && inState->stepSize < inState->maxStepSize && inState->stepSize > kMinStepSize );
	
	if( inState->curvature_step != 0 && workingLength > oldLength && 1<0)
	{
		int		attempts = 0;
		double	lh=0, rh=inState->stepSize;
		double	ss = (lh+rh)/2.0, newLen;
		int		win = 0;
		// bsearch for minimum length value up to 7 levels
		for( attempts=0; attempts<10; attempts++ )
		{
			if( workerLink != NULL )
				octrope_link_free(workerLink);
			workerLink = octrope_link_copy(inLink);
			octrope_link_fix_wrap(workerLink);

			step(workerLink, ss, dVdt, inState);
			newLen = octrope_curvelength(workerLink);
		//	printf( "%3.16lf (%e) %3.16lf\n", newLen, fabs(newLen-oldLength), oldLength);
			if( newLen > oldLength )
			{
				rh = ss;
			}
			else
			{
				lh = ss;
				win = 1;
			}
			ss = (lh+rh)/2.0;
		}
		
		if( win == 0 )
		{
			fprintf( stderr, "** Never won during length bsearch\n" );
			
		}
	}
	
	if( gPaperInfoInTmp != 0 )
	{
		char fname[512];
		preptmpname(fname,"bsearch_count",inState);
		FILE* tb = fopen(fname,"a");
		fprintf(tb, "%d\n", stepAttempts);
		fclose(tb);
	}
	
	inState->minrad = newmr;
	
	if( inState->stepSize < kMinStepSize )
		inState->stepSize = kMinStepSize;
				
	if( inState->curvature_step != 0 && inState->eq_step == 0 )
	{
		inState->cstep_time += inState->stepSize;
		inState->time += inState->stepSize;
	}
	
	octrope_link_free(inLink);
	inLink = workerLink;
		
	// we're good, double step size and continue jamming
/*	if( stepAttempts == 1 && inState->curvature_step != 0  )
	{
		inState->stepSize *= 2;
	}
*/	
	if( inState->stepSize > inState->maxStepSize )
	{
		inState->stepSize = inState->maxStepSize;
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
	
	inState->length = octrope_curvelength(inLink);
	
	return inLink;
}

static octrope_link*
old_bsearch_step( octrope_link* inLink, search_state* inState )
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
		octrope_link_fix_wrap(workerLink);
		step(workerLink, inState->stepSize, dVdt, inState);
		if( gOutputFlag == 1 )
		{
			char	fname[1024];
			preptmpname(fname,"dVdt.vect",inState);
			exportVect(dVdt, inLink, fname);
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
	if( inState->curvature_step != 0 && inState->eq_step == 0 )
		inState->cstep_time += inState->stepSize;
	
	octrope_link_free(inLink);
	inLink = workerLink;
		
	// we're good, double step size and continue jamming
	if( stepAttempts == 1 && inState->curvature_step != 0  )
	{
		inState->stepSize *= 2;
	}
	
	if( inState->stepSize > inState->maxStepSize )
	{
		inState->stepSize = inState->maxStepSize;
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
	
	if( strutCount == 0 )
		return;
	
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
			entry = (totalStruts*3*inState->compOffsets[strutSet[sItr].component[0]]);
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
			entry = (totalStruts*3*inState->compOffsets[strutSet[sItr].component[1]]);
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
	
	int tot = 0;
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
			inState->avgSideLength += octrope_norm(side);
			tot++;
		}
	}
	inState->avgSideLength /= (double)tot;
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
		
	octrope_vector dlda, dldv, dldb;
	octrope_vector dxda, dxdv, dxdb;
	double x, dtanconst;

	int totalStruts = contactStruts + minradLocs;
		
	for( mItr=0; mItr<minradLocs; mItr++ )
	{
		octrope_vector  a, v, b, As, Bs, Cs;
		octrope_vector  prevSide, thisSide;
		double dot, prevLen, thisLen;
		int vItr = minradStruts[mItr].vert;
		int cItr = minradStruts[mItr].component;
		
		a.c[0] = inLink->cp[cItr].vt[vItr-1].c[0];
		a.c[1] = inLink->cp[cItr].vt[vItr-1].c[1];
		a.c[2] = inLink->cp[cItr].vt[vItr-1].c[2];
		
		v.c[0] = inLink->cp[cItr].vt[vItr].c[0];
		v.c[1] = inLink->cp[cItr].vt[vItr].c[1];
		v.c[2] = inLink->cp[cItr].vt[vItr].c[2];
		
		b.c[0] = inLink->cp[cItr].vt[vItr+1].c[0];
		b.c[1] = inLink->cp[cItr].vt[vItr+1].c[1];
		b.c[2] = inLink->cp[cItr].vt[vItr+1].c[2];

		prevSide.c[0] = inLink->cp[cItr].vt[vItr].c[0] - inLink->cp[cItr].vt[vItr-1].c[0];
		prevSide.c[1] = inLink->cp[cItr].vt[vItr].c[1] - inLink->cp[cItr].vt[vItr-1].c[1];
		prevSide.c[2] = inLink->cp[cItr].vt[vItr].c[2] - inLink->cp[cItr].vt[vItr-1].c[2];
		
		thisSide.c[0] = inLink->cp[cItr].vt[vItr+1].c[0] - inLink->cp[cItr].vt[vItr].c[0];
		thisSide.c[1] = inLink->cp[cItr].vt[vItr+1].c[1] - inLink->cp[cItr].vt[vItr].c[1];
		thisSide.c[2] = inLink->cp[cItr].vt[vItr+1].c[2] - inLink->cp[cItr].vt[vItr].c[2];
		
		dot = octrope_dot(prevSide, thisSide);
		prevLen = octrope_norm(prevSide);
		thisLen = octrope_norm(thisSide);
		
		double value, angle;
		value = dot/(prevLen*thisLen);
		if( value >= 1 )
		{
			angle = 0;
		}
		else if( value <= -1 )
		{
			angle = M_PI;
		}
		else
		{
			// this is dangerous for the reasons of Ted's talk. Change me!
			angle = acos(value);
		}
				
		if( octrope_norm(prevSide) < octrope_norm(thisSide) )
		{
			printf( "prevSide shorter\n" );
		
			dlda.c[0] = (a.c[0] - v.c[0])/prevLen;
			dlda.c[1] = (a.c[1] - v.c[1])/prevLen;
			dlda.c[2] = (a.c[2] - v.c[2])/prevLen;

			dldv.c[0] = -dlda.c[0];
			dldv.c[1] = -dlda.c[1];
			dldv.c[2] = -dlda.c[2];

			dldb.c[0] = dldb.c[1] = dldb.c[2] = 0;

			x = dot/(prevLen*thisLen);
			
			dxda.c[0] = ( (v.c[0]-b.c[0])*prevLen - (dot)*((a.c[0]-v.c[0])/prevLen) )/(prevLen*prevLen*thisLen);
			dxda.c[1] = ( (v.c[1]-b.c[1])*prevLen - (dot)*((a.c[1]-v.c[1])/prevLen) )/(prevLen*prevLen*thisLen);
			dxda.c[2] = ( (v.c[2]-b.c[2])*prevLen - (dot)*((a.c[2]-v.c[2])/prevLen) )/(prevLen*prevLen*thisLen);
			
			dxdv.c[0] = ( (b.c[0] - 2*v.c[0])*prevLen*thisLen - dot*( (prevLen*(v.c[0]-b.c[0])/thisLen) + (thisLen*(v.c[0]-a.c[0])/prevLen) ) )/(prevLen*prevLen*thisLen*thisLen);
			dxdv.c[1] = ( (b.c[1] - 2*v.c[1])*prevLen*thisLen - dot*( (prevLen*(v.c[1]-b.c[1])/thisLen) + (thisLen*(v.c[1]-a.c[1])/prevLen) ) )/(prevLen*prevLen*thisLen*thisLen);
			dxdv.c[2] = ( (b.c[2] - 2*v.c[2])*prevLen*thisLen - dot*( (prevLen*(v.c[2]-b.c[2])/thisLen) + (thisLen*(v.c[2]-a.c[2])/prevLen) ) )/(prevLen*prevLen*thisLen*thisLen);
	
			dxdb.c[0] = ( (v.c[0]-a.c[0])*thisLen - (dot)*((b.c[0]-v.c[0])/thisLen) )/(prevLen*thisLen*thisLen);
			dxdb.c[1] = ( (v.c[1]-a.c[1])*thisLen - (dot)*((b.c[1]-v.c[1])/thisLen) )/(prevLen*thisLen*thisLen);
			dxdb.c[2] = ( (v.c[2]-a.c[2])*thisLen - (dot)*((b.c[2]-v.c[2])/thisLen) )/(prevLen*thisLen*thisLen);			
			
			dtanconst = (x-1)/pow(1-x*x, 3.0/2.0);
			
			//	                   This is the tan(theta/2)' part
			As.c[0] = ( -2*prevLen *   (dxda.c[0] * dtanconst)     + 2*dlda.c[0]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			As.c[1] = ( -2*prevLen *   (dxda.c[1] * dtanconst)     + 2*dlda.c[1]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			As.c[2] = ( -2*prevLen *   (dxda.c[2] * dtanconst)     + 2*dlda.c[2]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			
			Bs.c[0] = ( -2*prevLen *   (dxdv.c[0] * dtanconst)     + 2*dldv.c[0]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			Bs.c[1] = ( -2*prevLen *   (dxdv.c[1] * dtanconst)     + 2*dldv.c[1]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			Bs.c[2] = ( -2*prevLen *   (dxdv.c[2] * dtanconst)     + 2*dldv.c[2]*tan(angle/2) ) / pow(2*tan(angle/2), 2);	
			
			Cs.c[0] = ( -2*prevLen *   (dxdb.c[0] * dtanconst)     + 2*dldb.c[0]*tan(angle/2) ) / pow(2*tan(angle/2), 2);	
			Cs.c[1] = ( -2*prevLen *   (dxdb.c[1] * dtanconst)     + 2*dldb.c[1]*tan(angle/2) ) / pow(2*tan(angle/2), 2);	
			Cs.c[2] = ( -2*prevLen *   (dxdb.c[2] * dtanconst)     + 2*dldb.c[2]*tan(angle/2) ) / pow(2*tan(angle/2), 2);	
		
		}
		else // bc < ab
		{
			printf( "forward shorter\n" );
			
			dlda.c[0] = dlda.c[1] = dlda.c[2] = 0;

			dldv.c[0] = -thisSide.c[0]/thisLen;
			dldv.c[1] = -thisSide.c[1]/thisLen;
			dldv.c[2] = -thisSide.c[2]/thisLen;

			dldb.c[0] = -dldv.c[0];
			dldb.c[1] = -dldv.c[1];
			dldb.c[2] = -dldv.c[2];
			
			x = dot/(prevLen*thisLen);
			
			dxda.c[0] = ( (v.c[0]-b.c[0])*prevLen - (dot)*((a.c[0]-v.c[0])/prevLen) )/(prevLen*prevLen*thisLen);
			dxda.c[1] = ( (v.c[1]-b.c[1])*prevLen - (dot)*((a.c[1]-v.c[1])/prevLen) )/(prevLen*prevLen*thisLen);
			dxda.c[2] = ( (v.c[2]-b.c[2])*prevLen - (dot)*((a.c[2]-v.c[2])/prevLen) )/(prevLen*prevLen*thisLen);
			
			dxdv.c[0] = ( (b.c[0] - 2*v.c[0])*prevLen*thisLen - dot*( (prevLen*(v.c[0]-b.c[0])/thisLen) + (thisLen*(v.c[0]-a.c[0])/prevLen) ) )/(prevLen*prevLen*thisLen*thisLen);
			dxdv.c[1] = ( (b.c[1] - 2*v.c[1])*prevLen*thisLen - dot*( (prevLen*(v.c[1]-b.c[1])/thisLen) + (thisLen*(v.c[1]-a.c[1])/prevLen) ) )/(prevLen*prevLen*thisLen*thisLen);
			dxdv.c[2] = ( (b.c[2] - 2*v.c[2])*prevLen*thisLen - dot*( (prevLen*(v.c[2]-b.c[2])/thisLen) + (thisLen*(v.c[2]-a.c[2])/prevLen) ) )/(prevLen*prevLen*thisLen*thisLen);
	
			dxdb.c[0] = ( (v.c[0]-a.c[0])*thisLen - (dot)*((b.c[0]-v.c[0])/thisLen) )/(prevLen*thisLen*thisLen);
			dxdb.c[1] = ( (v.c[1]-a.c[1])*thisLen - (dot)*((b.c[1]-v.c[1])/thisLen) )/(prevLen*thisLen*thisLen);
			dxdb.c[2] = ( (v.c[2]-a.c[2])*thisLen - (dot)*((b.c[2]-v.c[2])/thisLen) )/(prevLen*thisLen*thisLen);

			dtanconst = (x-1)/pow(1-x*x, 3.0/2.0);
			
			//	                   This is the tan(theta/2)' part
			As.c[0] = ( -2*thisLen *   (dxda.c[0] * dtanconst)     + 2*dlda.c[0]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			As.c[1] = ( -2*thisLen *   (dxda.c[1] * dtanconst)     + 2*dlda.c[1]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			As.c[2] = ( -2*thisLen *   (dxda.c[2] * dtanconst)     + 2*dlda.c[2]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			
			Bs.c[0] = ( -2*thisLen *   (dxdv.c[0] * dtanconst)     + 2*dldv.c[0]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			Bs.c[1] = ( -2*thisLen *   (dxdv.c[1] * dtanconst)     + 2*dldv.c[1]*tan(angle/2) ) / pow(2*tan(angle/2), 2);
			Bs.c[2] = ( -2*thisLen *   (dxdv.c[2] * dtanconst)     + 2*dldv.c[2]*tan(angle/2) ) / pow(2*tan(angle/2), 2);	
			
			Cs.c[0] = ( -2*thisLen *   (dxdb.c[0] * dtanconst)     + 2*dldb.c[0]*tan(angle/2) ) / pow(2*tan(angle/2), 2);	
			Cs.c[1] = ( -2*thisLen *   (dxdb.c[1] * dtanconst)     + 2*dldb.c[1]*tan(angle/2) ) / pow(2*tan(angle/2), 2);	
			Cs.c[2] = ( -2*thisLen *   (dxdb.c[2] * dtanconst)     + 2*dldb.c[2]*tan(angle/2) ) / pow(2*tan(angle/2), 2);	

		}
		
	/*	As.c[0] /= octrope_norm(As);
		As.c[1] /= octrope_norm(As);
		As.c[2] /= octrope_norm(As);
		
		Bs.c[0] /= octrope_norm(Bs);
		Bs.c[1] /= octrope_norm(Bs);
		Bs.c[2] /= octrope_norm(Bs);
		
		Cs.c[0] /= octrope_norm(Cs);
		Cs.c[1] /= octrope_norm(Cs);
		Cs.c[2] /= octrope_norm(Cs);
	*/	
		
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
		
		if( minradStruts[mItr].vert == 0 &&
			(inLink->cp[minradStruts[mItr].component].acyclic == 0) )
		{
			entry = (totalStruts*3)*(inState->compOffsets[minradStruts[mItr].component]);
		}
		else
		{
			entry = (totalStruts*3*aVert)+contactStruts+mItr;
		}

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


		if( minradStruts[mItr].vert+1 == (inLink->cp[minradStruts[mItr].component].nv) &&
			(inLink->cp[minradStruts[mItr].component].acyclic == 0) )
		{
			entry = (totalStruts*3)*(inLink->cp[minradStruts[mItr].component].nv-1 + inState->compOffsets[minradStruts[mItr].component]);
		}
		else
		{
			entry = (totalStruts*3*cVert)+contactStruts+mItr;
		}
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
	
	if( gOutputFlag == 0 )
		return;
	
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
		char fname[1024];
		preptmpname(fname,"curvature.vect",inState);
		export_pushed_edges( inLink, inState, pushes, fname, 1 );
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
	char fname[1024];
	preptmpname(fname,"adjustedDL.vect",inState);
	exportVect( dlen, inLink, fname );
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
	{
		char fname[1024];
		preptmpname(fname,"adjustedDL.vect",inState);
		exportVect( dlen, inLink, fname );
	}
}

static void
eqForce( octrope_vector* dlen, octrope_link* inLink, search_state* inState )
{
	int cItr, vItr;
	
	// note to self - this can be made much faster by not being dumb
	
	octrope_link_fix_wrap(inLink);
	
	if( gOutputFlag )
	{
		char fname[1024];
		preptmpname(fname,"originalDL.vect",inState);
		exportVect( dlen, inLink, fname );
	}
	
	// if we're graphing kEqVariance
	double* diffFromAvg;
	double* averages;
	if( inState->graphing[kEQVariance] != 0 )
	{
		diffFromAvg = (double*)malloc(sizeof(double)*inState->totalSides);
		averages = (double*)malloc(sizeof(double)*inLink->nc);
	}
	
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
		
		if( inState->graphing[kEQVariance] != 0 )
		{
			for( vItr=0; vItr<edges; vItr++ )
				diffFromAvg[vItr+inState->compOffsets[cItr]] = fabs(lengths[vItr] - averageLength)/averageLength;
			
			averages[cItr] = averageLength;
		}
		
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
			
	//		scaleFactor = 1;


//			** OLD WAY **
/*			adjustments[(vItr+1)%edges].c[0] = (scaleFactor)*sides[vItr].c[0];
			adjustments[(vItr+1)%edges].c[1] = (scaleFactor)*sides[vItr].c[1];
			adjustments[(vItr+1)%edges].c[2] = (scaleFactor)*sides[vItr].c[2];
*/
			/* "we change the eq code so that it moves in the direction perpendicular to 
			* the gradient of length (that is, the average of the two tangent vectors) rather 
			* than in the direction of the tangent vector to the "left" edge at each vertex.
			* This might make a small difference in eq performance (it certainly shouldn't hurt)
			* and it makes the theory nicer." 
			*/
			adjustments[(vItr+1)%edges].c[0] = (scaleFactor)*(sides[vItr].c[0]+sides[vItr+1].c[0])/2.0;
			adjustments[(vItr+1)%edges].c[1] = (scaleFactor)*(sides[vItr].c[1]+sides[vItr+1].c[1])/2.0;
			adjustments[(vItr+1)%edges].c[2] = (scaleFactor)*(sides[vItr].c[2]+sides[vItr+1].c[2])/2.0;
			
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
	
	if( inState->graphing[kEQVariance] != 0 )
	{
		double sum = 0, isum, sampleAvg=0;
		int N = 0;

		// remember you made these all percent diffs 
		
		for( cItr=0; cItr<inLink->nc; cItr++ )
		{
			for( vItr=0; vItr<inState->totalSides; vItr++ )
			{
				sampleAvg += diffFromAvg[vItr+inState->compOffsets[cItr]];
				N++;
			}
		}
		sampleAvg /= (double)N;

		for( cItr=0; cItr<inLink->nc; cItr++ )
		{
			for( vItr=0; vItr<inState->totalSides; vItr++ )
			{
				isum = diffFromAvg[vItr+inState->compOffsets[cItr]] - sampleAvg;
				sum += isum*isum;
			}
		}
		
		inState->eqAvgDiff = sampleAvg;
		inState->eqVariance = (1.0/(double)N)*sum;
		
		free(diffFromAvg);
		free(averages);
	}
	
	if( gOutputFlag )
	{
		char	fname[1024];
		preptmpname(fname,"adjustedDL.vect",inState);
		exportVect( dlen, inLink, fname );
	}
}

static double*
stanford_lsqr( taucs_ccs_matrix* sparseA, double* minusDL, double* residual )
{
	/* if there are no constrained struts, t_snnls won't actually work, so use SOL LSQR */
	lsqr_input   *lsqr_in;
	lsqr_output  *lsqr_out;
	lsqr_work    *lsqr_work;
	lsqr_func    *lsqr_func;
	int bItr;
	double*		result;
						
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
	lsqr_func->mat_vec_prod = sparse_lsqr_mult;
	
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
	double				dummyThick;
			
	double*				A = NULL; // the rigidity matrix
	
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
	
//	dl = (octrope_vector*)calloc(inState->totalVerts, sizeof(struct octrope_vector_type));
		
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
						0.0005,
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
			minradLocs += octrope_pline_edges(&inLink->cp[cItr]);
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
			placeMinradStruts2(A, inLink, minradSet, minradLocs, inState, strutCount);
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
			
		taucs_ccs_matrix* sparseA = NULL;
		
		// fill this in now if we're not correcting, otherwise we might 
		// have to flip strut gradients for those already in the green zone
		if( dlenStep != 0 || inState->eq_step != 0 )
			sparseA = taucs_construct_sorted_ccs_matrix(A, strutCount+minradLocs, 3*inState->totalVerts);
			
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
				
				int sIndex, mItr;
				for( sIndex=0; sIndex<greenZoneCount; sIndex++ )
				{
					// flip the entries in this guy's column
					int totalStruts = strutCount+minradLocs;
					sItr = greenZoneStruts[sIndex];
					int entry = (totalStruts*3*strutSet[sItr].lead_vert[0])+sItr;
					A[entry] *= -1;
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
					A[entry+(2*totalStruts)] *= -1;
				}
							
				// done with these
				free(greenZoneStruts);
				
				// now we can build since we've flipped gradients in full A.
				sparseA = taucs_construct_sorted_ccs_matrix(A, strutCount+minradLocs, 3*inState->totalVerts);
				sparseAT = taucs_ccs_transpose(sparseA);
		
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

				
				/*octrope_vector* ofvVec;
				ofvVec = (octrope_vector*)malloc(sizeof(octrope_vector)*inState->totalVerts);
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
					octrope_vector  vert;
					double  norm;
					vert.c[0] = ofv[3*dlItr+0];
					vert.c[1] = ofv[3*dlItr+1];
					vert.c[2] = ofv[3*dlItr+2];
					norm = octrope_norm(vert);
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
				}*/
			
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
			} // gFastCorrectionSteps
			else
			{
				int rItr;
				// do this. follow grads.
				double	scaleFactor;
				double secondGreenZone = thickness - (thickness*inState->overstepTol)*.25;
				double greenZone = thickness - (thickness*inState->overstepTol)*0.5;
				
				sparseA = taucs_construct_sorted_ccs_matrix(A, strutCount+minradLocs, 3*inState->totalVerts);	
				
				// struts are rows
				for( sItr=0; sItr<strutCount; sItr++ )
				{
					if( strutSet[sItr].length < greenZone )
					{
						scaleFactor = thickness-strutSet[sItr].length;
						// all we need to do is copy over the column in an appropriately 
						// scaled way.
						for( rItr=sparseA->colptr[sItr]; rItr<sparseA->colptr[sItr+1]; rItr++ )
						{
							minusDL[sparseA->rowind[rItr]] = scaleFactor*sparseA->values.d[rItr];
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
						
						for( rItr=sparseA->colptr[sItr]; rItr<sparseA->colptr[sItr+1]; rItr++ )
						{
							minusDL[sparseA->rowind[rItr]] += scaleFactor*sparseA->values.d[rItr];
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
				taucs_print_ccs_matrix(sparseA, matrixF);
				fclose(matrixF);
			}
										
		// solve AX = -dl, x is strut compressions
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
			inState->rcond = taucs_rcond(sparseA);
			
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
		
		if( compressions == NULL )
		{
			// we've probably hit the rcond wall, scaling should smooth transition
			//exit(-1);
			printf( "****** NULL compressions!\n" );
			fprintf( stderr, "****** NULL compressions!\n" );
//			link_scale(inLink, 1.001);
			exit(-1);
			return;
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
			{
				if( sItr < strutCount )
					partialMult += compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+0];
				else // we're in minrad land (-- actually this means nothing now)
					partialMult += 1*compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+0];
			}
			dVdt[dlItr].c[0] = dl[dlItr].c[0] + partialMult;
		
			partialMult=0;
			for( sItr=0; sItr<totalStruts; sItr++ )
			{
				if( sItr < strutCount )
					partialMult += compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+totalStruts];
				else
					partialMult += 1*compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+totalStruts];
			}
			dVdt[dlItr].c[1] = dl[dlItr].c[1] + partialMult;
			
			partialMult = 0;
			for( sItr=0; sItr<totalStruts; sItr++ )
			{
				if( sItr < strutCount )
					partialMult += compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+(2*totalStruts)];
				else
					partialMult += 1*compressions[sItr]*A[(totalStruts*3*dlItr)+sItr+(2*totalStruts)];
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
		taucs_ccs_free(sparseA);
		free(A);
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
maxovermin( octrope_link* inLink, search_state* inState )
{
	int cItr, vItr;
	double max = 0, min = DBL_MAX;
	octrope_vector s1, s2;
	double max_maxovermin = 0, len=0;
	int edges, totalEdges=0;
	int maxVert, maxComp, minVert, minComp;
	
	inState->avgSideLength = 0;
	int tot = 0;
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		max = 0;
		min = DBL_MAX;
		
		edges = octrope_pline_edges(&inLink->cp[cItr]);
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
			
			octrope_vsub(s1, s2);
			norm = octrope_norm(s1);
			
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
computeCompressPush( octrope_link* inLink, octrope_strut* strutSet,
				octrope_mrloc* minradSet, int strutCount, int minradLocs )
{
	
}
