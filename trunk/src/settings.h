/*
 *  settings.h
 *  ridgerunner
 *
 *  Created by Michael Piatek on Thu Jun 03 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _H_settings
#define _H_settings

#include "gradient/stepper.h"

typedef struct
{
	double  rescaleThreshold;
} RSettings;

extern RSettings gRidgeSettings;

void preptmpname( char* outName, const char* inName, search_state* inState );

#endif // _H_settings
