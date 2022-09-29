/*

   mangle.h : Simple interface for mangling filenames by removing an extension
              (if present) and adding a new one. 

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

*/

#ifndef MANGLE_H__
#define MANGLE_H__ 1

#include"config.h"

#ifdef HAVE_STRING_H
#include<string.h>
#endif

#ifdef HAVE_STDLIB_H
#include<stdlib.h>
#endif

#ifdef HAVE_STDIO_H
#include<stdio.h>
#endif

FILE *fmangle(const char *filename,const char *oldextension,const char *newextension);
/* Create an output file by replacing "oldextension" with "newextension" in "filename". */
char *mangle(const char *filename,const char *oldextension,const char *newextension);
/* Create a filename by replacing "oldextension" by "newextension" in
   "filename". Caller must free the returned string. */
void  nmangle(char *newname,int newname_size,
	      const char *filename,const char *oldextension,const char *newextension);
/* Creates a filename in caller-supplied buffer "newname" of size
   "newname_size" by replacing "oldextension" by "newextension" in "filename". */

#endif
