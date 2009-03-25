/* Ordie duplicates some of the errors.c functionality in a way that 
   doesn't require the whole apparatus of a ridgerunner run. */

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
