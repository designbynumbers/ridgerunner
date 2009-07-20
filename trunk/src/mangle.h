/*

   mangle.h : Simple interface for mangling filenames by removing an extension
              (if present) and adding a new one. 

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
