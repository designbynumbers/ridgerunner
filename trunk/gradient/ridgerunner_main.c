/*
 *  ridgerunner_main.c
 *  ridgerunner
 *
 *  Created by Michael Piatek on Fri May 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
 
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/dir.h>

#include "errors.h"
#include "stepper.h"
#include "eqedge.h"
#include "settings.h"

void usage();
void reload_handler( int sig );
void initializeState( search_state* state, octrope_link** inLink, const char* fname );

int gVerboseFiling = 0;
int gSuppressOutput = 0;
int gQuiet = 0;
int gSurfaceBuilding = 0;
int gAvoidTmpConflicts = 0;
int gPaperInfoInTmp = 0;
int gFastCorrectionSteps = 0;

#define min(a,b) ((a)<(b)) ? (a) : (b)

/*	
	struct options longopts[] = {
		{"avoidTmpConflicts",	no_argument,		&gAvoidTmpConflicts,	1},
		{"saveConvergenceData", no_argument,		&saveConvergence,		1},
		{"quiet",				no_argument,		&gQuiet,				1},
		{"buildSurface",		no_argument,		&gSurfaceBuilding,		1},
		{"verbose",				no_argument,		&gVerboseFiling,		1},
		{"suppressFiles",		no_argument,		&gSuppressOutput,		1},
		{"movie",				no_argument,		&movie,					1},
		{"ignoreCurvature",		no_argument,		&ignorecurvature,		1},
		{"autoscale",			no_argument,		&autoscale,				1},
		{"fancyViz",			no_argument,		&fancyViz,				1},
		{"infile",				required_argument,	NULL,					'a'},
		{"scale",				required_argument,	NULL,					'b'},
		{"checkDelta",			required_argument,	NULL,					'c'},
		{"injrad",				required_argument,	NULL,					'd'},
		{"eqMultiplier",		required_argument,	NULL,					'e'},
		{"residualThreshold",	required_argument,	NULL,					'f'},
		{"refineUntil",			required_argument,	NULL,					'g'},
		{"checkThreshold",		required_argument,	NULL,					'h'},
		{"maxSteps",			required_argument,	NULL,					'i'},
		{"overstepTolerance",	required_argument,	NULL,					'j'},
		{"maxStepSize",			required_argument,	NULL,					'k'},
		{"double",				no_argument,		NULL,					'l'},
		{0, 0, 0, 0}
		};
*/
	
int
main( int argc, char* argv[] )
{
	octrope_link*	link = NULL;
	FILE*			linkFile = NULL;
	int				cItr=0;
	search_state	state;
	char			opt;
	char			cmd[1024];
	short			autoscale = 0, fixlengths = 0, movie = 0;
	int				refineUntil = 0, ignorecurvature=0;
	double			checkDelta = 1;
	double			injrad = 0.5, overstepTol=0.0001;
	double			residualThreshold = 65000;
	double			checkThreshold = 0.05;
	double			eqMult = 2.0;
	int				doubleCount = 0;
	char			fname[1024];
	short			fancyViz = 0, saveConvergence = 0;
	double			scaleAmt = 1;
	double			maxStep = -1;
	long			maxItrs = -1;
	int				graphIt[kTotalGraphTypes];
	int				gItr;
	double			correctionStepSize = 0.25;
	double			minradOverstepTol = 0.499975;
	
	bzero(&graphIt[0], sizeof(int)*kTotalGraphTypes);
		
	printf( "cvs client build: %s (%s)\n", __DATE__, __TIME__ );
	
	srand(time(NULL));
	
	while( (opt = getopt(argc, argv, "vlf:mnuap:t:d:k:g:b:c:r:i:o:e:yqjx:sh:w:z")) != -1 )
	{
		switch(opt)
		{
			case 'g':
				switch((int)(optarg[0]))
				{
					case 'a':
						for( gItr=0; gItr<kTotalGraphTypes; gItr++ )
							graphIt[gItr] = 1;
						break;
					
					case 'l': graphIt[kLength] = 1; break;
					case 'c': graphIt[kWallTime] = 1; break;
					case 'r': graphIt[kRcond] = 1; break;
					case 'm': graphIt[kMinrad] = 1; break;
					case 's': graphIt[kStrutCount] = 1; break;
					case 'v': graphIt[kEQVariance] = 1; break;
						
					default: 
						fprintf(stderr,"Unknown graphing option\n");
						break;
				}
				break;
			
			case 'k':
				correctionStepSize = atof(optarg);
				break;
		
			case 'z':
				gAvoidTmpConflicts = 1;
				break;
		
			case 'w':
				//saveConvergence = 1;
				switch( (int)(optarg[0]) )
				{
					case 'f':
						gFastCorrectionSteps = 1;
						break;
					 
					default: 
						printf( "unknown width option\n" ); 
						break;
				}
				break;
		
			case 'q':
				gQuiet = 1;
				break;
				
			case 'x':
				sscanf(optarg, "%lf", &maxStep);
				break;
			
			case 'u':
				gSurfaceBuilding = 1;
				break;
				
			case 'v':
				gVerboseFiling = 1;
				break;
								
			case 'd':
				minradOverstepTol = atof(optarg);
				break;
				
			case 's':
				gSuppressOutput = 1;
				break;
				
			case 'o': // min overstep, percentage of thickness
				if( optarg == NULL )
					usage();
				else
					overstepTol = atof(optarg);
				break;
			
			case 'p':
				maxItrs = atoi(optarg);
				break;
				
			case 'h': // check delta threshold
				if( optarg == NULL )
					usage();
				else
					checkThreshold = atof(optarg);
				break;
			
			case 'r': // refine until discretization at x(times)rope
				if( optarg == NULL )
					usage();
				refineUntil = atoi(optarg);
				if( refineUntil < 0 )
					usage();
				break;
			
			case 't':
				if( optarg == NULL )
					usage();
				residualThreshold = atof(optarg);
				if( residualThreshold < 0 )
					usage();
				break;
			
			case 'e':
				if( optarg == NULL )
					usage();
				eqMult = atof(optarg);
				if( eqMult < 0 )
					usage();
				break;
			
			case 'i':
				if( optarg == NULL )
					usage();
				injrad = atof(optarg);
				if( injrad < 0 )
					usage();
				break;
		
			case 'c': // check delta
				if( optarg == NULL )
					usage();
				checkDelta = atof(optarg);
				if( checkDelta < 0 )
					usage();
				break;
				
			case 'b':
				scaleAmt = atof(optarg);
				break;
				
			case 'l': // fixlengths
				fixlengths = 1;
				break;
		
			case 'f': // file
				linkFile = fopen(optarg, "r");
				if( linkFile == NULL )
				{
					fprintf(stderr, "Specified file: %s doesn't exist\n", optarg);
					return -1;
				}
				link = octrope_link_read(linkFile);
				strcpy(fname, optarg);
				fclose(linkFile);
				break;
				
			case 'm': // movie
				movie = 1;
				break;
			
			case 'n':
				ignorecurvature = 1;
				break;
				
			case 'j':
				gPaperInfoInTmp = 1;
				break;
			
			case 'a':
				autoscale = 1;
				break;
			
			case 'y':
				fancyViz = 1;
				break;
			
			default:
				usage();
				break;
		}
	}
	
	if( link == NULL )
		usage();
		
	if( fixlengths == 1 )
	{
		octrope_link* tempLink = link;
		link = octrope_fixlength(tempLink);
		octrope_link_free(tempLink);
		tempLink = link;
		link = octrope_fixlength(tempLink);
		octrope_link_free(tempLink);
		tempLink = link;
		link = octrope_fixlength(tempLink);
		octrope_link_free(tempLink);
	}
	
	int dItr;
	for( dItr=0; dItr<doubleCount; dItr++ )
	{
		octrope_link* tempLink = link;
		link = octrope_double_edges(tempLink);
		octrope_link_free(tempLink);
	}
	
	initializeState( &state, &link, fname);
	updateSideLengths(link, &state);
	
	state.refineUntil = refineUntil;
	state.checkDelta = checkDelta;
	state.checkThreshold = checkThreshold;
	state.movie = movie;
	
	state.injrad = injrad;
	state.overstepTol = overstepTol;
	state.minradOverstepTol = minradOverstepTol;
	state.minminrad = minradOverstepTol;
	
	state.eqMultiplier = eqMult;
	
	state.residualThreshold = residualThreshold;

	printf( "Overstep tolerance: %f of thickness %f (%f) minminrad: %3.8lf\n", state.overstepTol, state.injrad*2, 
			state.injrad*2 - state.overstepTol*state.injrad*2, 
			state.minminrad );
	
	// scale to thickness 1
	if( autoscale == 1 )
	{
		printf( "autoscaling with injrad %f, scale factor: %e\n", state.injrad, ((2*state.injrad)/octrope_thickness(link, 1, NULL, 0))-(2*state.injrad) );
		link_scale(link, (2*state.injrad)/octrope_thickness(link, 1, NULL, 0) );
	}
	else
	{
		link_scale(link, scaleAmt );
	}
		
	// Create directory to store movie frames if we're making a movie
	if( state.movie != 0 )
	{
		sprintf( cmd, "mkdir movie%s", state.fname );
		system(cmd);
	}
	
	// use amd column ordering if the user hasn't specified something else
	setenv("COL_ORDERING", "amd", 0);
				
	// install signal handler
	signal(SIGUSR1, reload_handler);
	
	state.fancyVisualization = fancyViz;
	
	if( gSuppressOutput == 0 )
	{
		preptmpname(fname, "orig.vect", &state);
		FILE* verts = fopen(fname, "w");
		octrope_link_draw(verts, link);
		fclose(verts);
		
		if( state.fancyVisualization != 0 )
		{
			// create fancypipe
			preptmpname(fname, "fancypipe", NULL );
			char cmd[1024];
			sprintf(cmd, "togeomview %s geomview -nopanels < /dev/null", fname );
			system(cmd);
			state.fancyPipe = fopen(fname,"w");
			fprintf(state.fancyPipe, "(normalization g0 none)\n");
			fprintf(state.fancyPipe, "(bbox-draw g0 no)\n");
			fprintf(state.fancyPipe, "(merge-ap g0 appearance {*shading smooth} )\n");
			fflush(state.fancyPipe);
		}
	}

	// realtime graphing! -- it is not advisable to do more than a few of these at once
	// (or any, if you're not on a mac using fink installed gnuplot)
	if( fancyViz != 0 || graphIt[kLength] )
		state.graphing[kLength] = 1;
	state.graphing[kRopelength] = graphIt[kRopelength];
	state.graphing[kStrutCount] = graphIt[kStrutCount];
	state.graphing[kStepSize] = graphIt[kStepSize];
	state.graphing[kThickness] = graphIt[kThickness];
	state.graphing[kMinrad] = graphIt[kMinrad];
	state.graphing[kResidual] = graphIt[kResidual]; 
	state.graphing[kMaxOverMin] = graphIt[kMaxOverMin];
	state.graphing[kRcond] = graphIt[kRcond];
	state.graphing[kWallTime] = graphIt[kWallTime];
	state.graphing[kMaxVertexForce] = graphIt[kMaxVertexForce];
	state.graphing[kConvergence] = graphIt[kConvergence];
	state.graphing[kEQVariance] = graphIt[kEQVariance];
	
	if( ignorecurvature != 0 )
		state.ignore_minrad = 1;
		
	if( maxStep > 0 )
	{
		printf("max step size: %e\n", maxStep);
		state.maxStepSize = maxStep;
	}
	state.stepSize = min(state.stepSize, state.maxStepSize);
	
	state.maxItrs = maxItrs;
	state.saveConvergence = saveConvergence;
	if( saveConvergence != 0 )
	{
		char fname[1024];
		preptmpname(fname, "rrconvergence.txt", &state);
		// empty file and touch, essentially
		FILE* conv = fopen(fname,"w");
		fclose(conv);
	}
	
	state.correctionStepDefault = correctionStepSize;
	
	if( gPaperInfoInTmp != 0 )
	{
		char	fname[512];
		preptmpname(fname, "bsearch_count", &state);
		remove(fname);
		preptmpname(fname, "correction_convergence",&state);
		remove(fname);
		preptmpname(fname, "rcond", &state);
		remove(fname);
	}
							
	bsearch_stepper(&link, &state);
	octrope_link_free(link);
		
	return kNoErr;
}

static short
rgetline( char* outLine, int inSize, FILE* fp )
{
	outLine[0] ='\0';
	int pos = 0;
	
	do
	{
		outLine[pos] = fgetc(fp);
	} while( outLine[pos++] != '\n' && !feof(fp) );
	outLine[pos] = '\0';
	
	if( strlen(outLine) == 0 )
		return 0;
	return 1;
}

void
initializeState( search_state* state, octrope_link** inLink, const char* fname )
{
	int cItr, i, offset=0;
	char	line[1024];
	FILE* fp = fopen(fname, "r"); // we would have failed earlier if this didn't exist
	
	// zero all those pointers
	bzero(state, sizeof(search_state));
	
	strncpy(state->fname, fname, sizeof(state->fname)-1);
	
	state->injrad = 0.5; // fix this for now
		
	state->maxStepSize = 0.1*octrope_link_short_edge(*inLink);
		
	if( state->maxStepSize > 1e-3 )
		state->maxStepSize = 1e-3;
//	state->maxStepSize = 1e-4;
	state->shortest = 2*state->injrad;
	state->totalVerts = 0;
	state->eqThreshold = 1.05;
	state->eqMultiplier = 1;
	state->factor = 1;
	
	state->minminrad = 0.499975;
	
	state->minrad = octrope_minradval(*inLink);
	for( cItr=0; cItr<(*inLink)->nc; cItr++ )
	{
		state->totalVerts += (*inLink)->cp[cItr].nv;
	}
	
	state->stepSize = 0.5*state->maxStepSize;
	if( ((double)1)/(double)state->totalVerts < state->stepSize )
	{
		state->stepSize = ((double)1/(double)state->totalVerts);
	}
	state->compOffsets = malloc(sizeof(int)*(*inLink)->nc);
	for( i=0; i<(*inLink)->nc; i++ )
	{
		state->compOffsets[i] = offset;
		offset += (*inLink)->cp[i].nv;
	}
	state->conserveLength = (int*)calloc((*inLink)->nc, sizeof(int));
	for( i=0; i<(*inLink)->nc; i++ )
	{
		state->conserveLength[i] = 0;
	}
	
	state->residual = 0;
	 
	state->checkDelta = 1.0;	// default check delta
	
	fclose(fp);
	
	getrusage(RUSAGE_SELF, &state->frameStart);
}

void
usage()
{
	fprintf( stdout, 
"ridgerunner \n\n \
-h threshold for check delta, will stop if length improvement < threshld in checkDelta time\n \
-d minrad control -- specify min minrad (should be less than thickness/2)\n \
-o tolerance\t Overstep tolerance (maybe not broken)\n \
-r mult\t Runs until stopping criteron satisfied, then divides, until verts < mult*ropelength\n \
-e mult\t Sets the eq-force multiplier\n \
-b amt\tscales curve by amt\n \
-i injrad\t Sets injrad constraint (probably broken)\n \
-c delta\t Stop if in delta timestep time, we haven't improved length more than 0.05\n \
-l\t\t Forces edge lengths to be approximately equal via Rawdon's tangential runaround\n \
-f file\t Specifies runtile knot (required!)\n \
-m\t\t Will create a directory with timestep frames and strut sets to facilitate movie creation\n \
-x\tmax step size\tspecifies the maxiumum size of integration step\n \
-n\t\t Run will ignore curvature constraint (useful if curvature breaks things)\n \
-a\t\t Autoscale specified knot to thickness 1.0. Useful to save scaling runtime\n \
-t threshold\t Sets additional residual stopping requirement that residual < threshold\n \
-s suppress out Suppresses the progress files in /tmp\n \
-p max itrs\tspecifies the maximum number of integration steps to be taken. default is inf\n \
-v\t\textra verbose filing, will output files _every_ step \n \
-y\t\tEnables geomview visualization fanciness, req. tube, gnuplot, geomview on path\n \
-g\tgraphing with arg: (l)ength (c)locktime (r)cond (m)inrad (r)esidual (s)truts (v)ariance of edge diff from avg (a)ll\n \
-u\tsurface strut gen\tgenerates multiple files in /tmp for use with surfaceBuilder, use with -v\n \
-w\twith: (f)ast correction stepper\n \
-z\tavoid conflicts in /tmp/ between multiple users by stamping with pid and file name\n \
-j\tpaper info, print avg times through bsearch, correction steps to green in /tmp/ and rcond every step if its graphing on\n \
-k\tcorrection step size, default 0.25\n"
);
	exit(kNoErr);
}

void
reload_handler( int sig )
{
	if( sig == SIGUSR1 )
	{
		// immediately reload config
		
		// reinstall handler
		signal(SIGUSR1, reload_handler);
	}
}
