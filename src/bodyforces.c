/*
 *  bodyforces.c
 *  ridgerunner
 *
 *  Created by Jason Cantarella on 1/1/08.
 *  Copyright Jason Cantarella. 

 *  Loads a file of "bodyforces" to apply to a given link.

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

#include "ridgerunner.h"
#include "plCurve_readoogl.h"

bodyforce *read_bf_file( FILE *infile, plCurve *inLink );

/* Procedure reads a bodyforce file and checks it against the curve inLink for sanity. */

/* NOT IMPLEMENTED YET! */

