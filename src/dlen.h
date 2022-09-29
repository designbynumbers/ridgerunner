/*
 *  dlen.h
 *  ridgerunner3
 *
 *  Created by Michael Piatek on Sun Jan 18 2004.

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
