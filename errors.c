/*
 *  errors.c
 *  ridgerunner3
 *
 *  Created by Michael Piatek on Fri Jan 16 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "errors.h"
#include "octrope.h"

void 
DebugThrow( int inErr, const char* inFile, long inLine )
{
    printf( "fatal error: %d file: %s line: %ld \n", inErr, inFile, inLine );
    exit(inErr);
}

void
DebugWarning( int inErr, const char* inFile, long inLine )
{
    printf( "warning: %d file: %s line: %ld\n", inErr, inFile, inLine );
    fflush(stdout);
}

void
error_write( octrope_link* inLink )
{
	FILE* fp = fopen("/tmp/err.vect","w");
	octrope_link_write(fp,inLink);
	fclose(fp);
}

