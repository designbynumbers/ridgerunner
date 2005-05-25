/*
 *  settings.c
 *  ridgerunner
 *
 *  Created by Michael Piatek on Thu Jun 03 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "settings.h"

#include <sys/types.h>
#include <unistd.h>

RSettings   gRidgeSettings;

void
preptmpname( char* outName, const char* inName, search_state* inState )
{
	sprintf( outName, "/tmp/%d.%s-%s", getpid(), inState->fname, inName );
}

