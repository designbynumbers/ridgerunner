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

	octrope_link* link = NULL;
	FILE* linkFile = NULL;
	int aItr=1, cItr=0;
	search_state   state;

	// no real arg processing for now, just assume grad evolve on specified link
	if( argc < 2 )
		usage();
	
	printf( "cvs client\n" );
		
	// use amd column ordering if the user hasn't specified something else
	setenv("COL_ORDERING", "amd", 0);
				
	// install signal handler
	signal(SIGUSR1, reload_handler);
	
	linkFile = fopen(argv[1],"r");
	fatalifnull_(linkFile);
	
	system( "rm -f examples");
	system( "rm -f rmov.*");
	
	if( fgetc(linkFile) == 'V' )
	{
		ungetc('V', linkFile);
	
		link = octrope_link_read(linkFile);
		
	//	link->cp[0].acyclic = 1;

//		for( cItr=0; cItr<link->nc; cItr++ )
//			link->cp[cItr].acyclic = 1;
		
		fatalifnull_(link);
		fclose(linkFile);
		
	//	link = octrope_double_edges(link);
	//	link = octrope_double_edges(link);
				
		initializeState( &state, &link, argv[1]);
				
	//		link->cp[0].acyclic=1;
	//	link->cp[1].acyclic=1;
		
	//	state.conserveLength[1] = 1;
	//	state.conserveLength[0] = 1;
				
		updateSideLengths(link, &state);

		FILE* verts = fopen("/tmp/orig.vect", "w");
		octrope_link_draw(verts, link);
		fclose(verts);

	//	link = octrope_fixlength(link);

//		state.graphing[0] = 1;

		// realtime graphing! -- it is not advisable to do more than a few of these at once
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
		
	//	state.conserveLength[1] = 1;
		
	}
	else
	{
		fclose(linkFile);
		link = NULL;
		// try to parse runfile
		initializeState(&state, &link, argv[1]);
		
		if( link == NULL )
		{
			printf( "no knots to process!\n" );
			return 0;
		}
		
		updateSideLengths(link, &state);
	}
		
/*	printf( "%f\n", maxovermin(link) );
	link = octrope_fixlength(link);
/*	printf( "%f\n", maxovermin(link) );
	link = octrope_fixlength(link);
	printf( "%f\n", maxovermin(link) );
*/	

	FILE* verts = fopen("/tmp/verts.vect", "w");
	octrope_link_draw(verts, link);
	fclose(verts);
	
	// 10000000
	
	octrope_link* store = link;
	link = octrope_fixlength(link);
	octrope_link_free(store);
	link_scale(link,1.0001);
		
	bsearch_stepper(&link, 10000000, &state);
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
	
	if( *inLink == NULL )
	{
		// for now just assume dirscanning and get first *.vect in the dir
		// and call ourselves again
		struct dirent* dp;
		DIR* workingDir = opendir(".");
		if( workingDir != NULL ) 
		{
			while( (dp = readdir(workingDir)) != NULL )
			{
				char	path[512];
				sprintf( path, "%s", dp->d_name );

				struct stat finfo;
				if( stat( path, &finfo ) == 0 )
				{
					if( !(finfo.st_mode & S_IFDIR) )	// if this isn't a directory...
					{
						char* vectPortion = strstr(path, ".vect");
						if( vectPortion != NULL )
						{
							char cmd[1024];
							FILE* linkFile = fopen(path,"r");
							if( linkFile == NULL )
								continue; // someone else got it first
							*inLink = octrope_link_read(linkFile);
							initializeState(state,inLink, path);
							fclose(linkFile);
							
							state->batching = 1;

							printf( "operating on %s\n", path );
							
							strcpy(cmd, "mv ");
							strcat(cmd, path);
							strcat(cmd, " ");
							vectPortion[0] = '\0';
							strcat(cmd, path);
							strcat(cmd, ".working");
							system(cmd);
							
							fclose(fp);
							
							return;
						}
					}
				}
			}
		}
	}
	else
	{
		// zero all those pointers
		bzero(state, sizeof(search_state));
		
		strncpy(state->fname, fname, sizeof(state->fname)-1);
		
		state->injrad = 0.5; // fix this for now
		state->maxStepSize = .01;
		state->stepSize = state->maxStepSize;
		state->shortest = 2*state->injrad;
		state->totalVerts = 0;
		state->eqThreshold = 0;
		state->eqMultiplier = 1;
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
	}
	
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
