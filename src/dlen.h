/*
 *  dlen.h
 *  ridgerunner3
 *
 *  Created by Michael Piatek on Sun Jan 18 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _H_dlen
#define _H_dlen

#include "octrope.h"
#include "stepper.h"

// note that we return dlen as a concatenated list of dlen for 
// each component in ascending order
void dlenForce( plc_vector* ioDL, plCurve* inLink, search_state* inState );

#endif // _H_dlen
