/* Ordie duplicates some of the errors.c functionality in a way that 
   doesn't require the whole apparatus of a ridgerunner run. 

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

#include "portability.h"

void FatalError(char *debugmsg,const char *file,int line) 

{
  fprintf(stderr,"%s at %s:%d.\n",debugmsg,file,line);
  exit(1);
}
  

void *malloc_or_die(size_t size, const char *file, const int line)
     
     /* Allocates memory or dies trying. */

{
  void *mret;

  mret = malloc(size);

  if (mret == NULL) {

    fprintf(stderr,"Could not allocate block of size %d in %s:%d\n",
	    (int)(size),file,line);
    exit(1);

  }

  return mret;
}
