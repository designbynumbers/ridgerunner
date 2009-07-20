#include "mangle.h"

FILE *fmangle(const char *filename,const char *oldextension,const char *newextension)
{
  FILE *outfile;
  char *name;

  name = mangle(filename,oldextension,newextension);
 
  outfile = fopen(name,"w");

  if (outfile == NULL) {

    fprintf(stderr,"fmangle: Could not open filename %s.\n",
	    name);

    exit(1);

  }

  free(name);
  return outfile;
}

char *mangle(const char *filename,const char *oldextension,const char *newextension)
/* Allocate space, then perform an "nmangle". */
{
  char *newname;
  int nnsize;

  nnsize = strlen(filename) + strlen(newextension) + 10;
  newname = calloc(nnsize,sizeof(char));
  
  if (newname == NULL) {

    fprintf(stderr,"mangle: Couldn't allocate string of size %d.\n",nnsize);
    exit(1);

  }

  nmangle(newname,nnsize,filename,oldextension,newextension);

  return newname;
}

void  nmangle(char *newname,int nnsize,
	      const char *filename,const char *oldextension,const char *newextension)

/* Replace the (terminating) string "oldextension" with "newextension" if present in 
"filename". Return the results in "newname". */
{
  
  if (nnsize < strlen(filename) + strlen(newextension) + 2) {

    fprintf(stderr,"nmangle: Buffer newname is not long enough to add extension %s to filename %s.\n",
	    newextension,filename);
    exit(1);

  }

  /* We have already checked that we have enough space, but we're cautious. */

  strncpy(newname,filename,nnsize);

  if (strstr(newname,oldextension) != NULL) {
	
    strcpy(strstr(newname,oldextension),newextension);
    
  } else {
       
    strncat(newname,newextension,nnsize);
       
  }

}
