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
#include <signal.h>

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

int
main( int argc, char* argv[] )
{
	octrope_link*	link = NULL;
	FILE*			linkFile = NULL;
	int				cItr=0;
	search_state	state;
	char			opt;
	char			cmd[1024];
	short			autoscale = 0, fixlengths = 0;
	int				checkDelta = 1, refineUntil = 0;
	
	printf( "cvs client build: %s (%s)\n", __DATE__, __TIME__ );
	
	while( (opt = getopt(argc, argv, "lf:mac:r:")) != -1 )
	{
		switch(opt)
		{
			case 'r': // refine until discretization at x(times)rope
				if( optarg == NULL )
					usage();
				refineUntil = atoi(optarg);
				if( refineUntil < 0 )
					usage();
				break;
		
			case 'c': // check delta
				if( optarg == NULL )
					usage();
				checkDelta = atoi(optarg);
				if( checkDelta < 0 )
					usage();
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
				initializeState( &state, &link, optarg);
				updateSideLengths(link, &state);
				fclose(linkFile);
				break;
				
			case 'm': // movie
				state.movie = 1;
				break;
			
			case 'a':
				autoscale = 1;
				break;
				
			default:
				usage();
				break;
		}
	}
	
	if( link == NULL )
		usage();
	
	state.refineUntil = refineUntil;
	state.checkDelta = checkDelta;
	
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
	
	// scale to thickness 1
	if( autoscale == 1 )
	{
		link_scale(link, 1.0/octrope_thickness(link, 1, NULL, 0) );
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
			
	FILE* verts = fopen("/tmp/orig.vect", "w");
	octrope_link_draw(verts, link);
	fclose(verts);

	// realtime graphing! -- it is not advisable to do more than a few of these at once
	// (or any, if you're not on a mac using fink installed gnuplot)
	state.graphing[kLength] = 0;
	state.graphing[kRopelength] = 0;
	state.graphing[kStrutCount] = 0;
	state.graphing[kStepSize] = 0;
	state.graphing[kThickness] = 0;
	state.graphing[kMinrad] = 0;
	state.graphing[kResidual] = 0; 
	state.graphing[kMaxOverMin] = 0;
	state.graphing[kRcond] = 0;
	state.graphing[kWallTime] = 0;
	state.graphing[kMaxVertexForce] = 0;
	state.graphing[kConvergence] = 0;
		
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
	state->maxStepSize = .01;
	state->stepSize = state->maxStepSize;
	state->shortest = 2*state->injrad;
	state->totalVerts = 0;
	state->eqThreshold = 0;
	state->eqMultiplier = 2;
	state->factor = 1;
	for( cItr=0; cItr<(*inLink)->nc; cItr++ )
	{
		state->totalVerts += (*inLink)->cp[cItr].nv;
	}
	
	state->stepSize = kStepScale*octrope_link_short_edge(*inLink);
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
}

void
usage()
{
	fprintf( stdout, 
			"Options information is as follows\n"
            "ridgerunner [runfile|vectfile] options\n\n"
            "    -h, -help			       Print this message\n");
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
