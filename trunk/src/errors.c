/*
 *  errors.c
 *  ridgerunner3
 *
 *  Created by Michael Piatek on Fri Jan 16 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "errors.h"
#include "octrope.h"

void 
DebugThrow( int inErr, const char* inFile, long inLine )
{
    printf( "fatal error: %d file: %s line: %ld \n", inErr, inFile, inLine );
    exit(inErr);
}

void
DebugWarning( int inErr, const char* inFile, long inLine )
{
    printf( "warning: %d file: %s line: %ld\n", inErr, inFile, inLine );
    fflush(stdout);
}

void
error_write( plCurve* inLink )
{
	FILE* fp = fopen("/tmp/err.vect","w");
	plc_write(fp,inLink);
	fclose(fp);
}

void FatalError(char *debugmsg,const char *file,int line)

     /* Write the error message to the global logfile
	and to stderr and quit. */

{
  
  fprintf(gLogfile,"ridgerunner: Fatal error in file %s, line %d.\n",file,line);
  fprintf(gLogfile,"%s",debugmsg);
  
  fprintf(stderr,"ridgerunner: Fatal error in file %s, line %d.\n",file,line);
  fprintf(stderr,"%s",debugmsg);

#ifdef HAVE_ASCTIME
#ifdef HAVE_LOCALTIME
#ifdef HAVE_TIME
  
  time_t end_time;
  end_time = time(NULL);
  
  fprintf(gLogfile,
	  "Run ended: %s.\n",
	  asctime(localtime(&end_time)));
  
  fprintf(stderr,
	  "Run ended: %s.\n",
	  asctime(localtime(&end_time)));
  
  
#endif
#endif
#endif
  
  exit(1);
  
}

FILE *fopen_or_die(const char *filename,const char *mode,const char *file,const int line) 

     /* Opens the file or dies trying. */

{

  FILE *outfile;
  char  errmsg[1024];

  outfile = fopen(filename,mode);

  if (outfile == NULL) {

    sprintf(errmsg,"ridgerunner: Could not open file %s.\n",filename);
    FatalError(errmsg,file,line);

  }

  return outfile;

}


void
dumpAxb_full( search_state *inState, 
	      double* A, int rows, int cols, 
	      double* x, double* b )
{
  /* Saves a problem Ax = b in MATLAB or Octave format */
  
  static FILE* fp;
  int rItr, cItr;
  char filename[1024],basename[1024];
  char errmsg[1024];
  
  /* Construct filename for A */

  sprintf(filename,"%sA.dat",inState->fprefix);
  fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
  
  /* Now construct the file. */
  
  for( rItr=0; rItr<rows; rItr++ ) {

    for( cItr=0; cItr<cols; cItr++ )
      fprintf( fp, "%10.16lf ", A[rItr*cols + cItr] );
      
    fprintf( fp, "\n" );
  }
  
  fclose(fp);

  /* Construct filename for x. */

  sprintf(basename,"x.dat");
  strcpy(filename,inState->fprefix);
  strcat(filename,basename);
  fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
  
  for( cItr=0; cItr<cols; cItr++ ) {

    if( x != NULL )
      fprintf( fp, "%10.16lf\n", x[cItr] );
    else
      fprintf( fp, "0\n" );
    
  }
  
  fclose(fp);

  /* Construct filename for b. */

  sprintf(basename,"b.dat");
  strcpy(filename,inState->fprefix);
  strcat(filename,basename);
  fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 

  for( rItr=0; rItr<rows; rItr++ )
    fprintf( fp, "%lf\n", b[rItr] );
  
  fclose(fp);
}

static void
dumpAxb_sparse( search_state *inState, taucs_ccs_matrix* A, double* x, double* b )
{
  int i, j;
  double* vals = taucs_convert_ccs_to_doubles(A);

  dumpAxb_full(inState,vals,A->m,A->n,x,b);

  free(vals);

}  
 
static void
dumpVertsStruts(plCurve* link, octrope_strut* strutSet, int strutCount)

     /* This is debugging code which shouldn't be called in a normal
	build of the software. It's not immediately clear what it's
	good for, but I'm leaving it archived here in case I figure
	out a need for it later. */
     
{
  int cItr, vItr, totalVerts=0;
  FILE* fp = fopen("ridgeverts","w");
  
  for( cItr=0; cItr<link->nc; cItr++ )
    {
      totalVerts += link->cp[cItr].nv;
    }
  
  fprintf( fp, "%d\n", totalVerts );
  
  for( cItr=0; cItr<link->nc; cItr++ )
    {
      for( vItr=0; vItr<link->cp[cItr].nv; vItr++ )
	{
	  fprintf(fp, "%lf\n%lf\n%lf\n", link->cp[cItr].vt[vItr].c[0],
		  link->cp[cItr].vt[vItr].c[1], link->cp[cItr].vt[vItr].c[2] );
	}
    }
  
  fprintf( fp, "%d\n", strutCount );
  // strut indices and positions. ASSUME ONE COMPONENT
  for( vItr=0; vItr<strutCount; vItr++ )
    fprintf( fp, "%d %d\n", strutSet[vItr].lead_vert[0], strutSet[vItr].lead_vert[1] );
  for( vItr=0; vItr<strutCount; vItr++ )
    fprintf( fp, "%lf %lf\n", strutSet[vItr].position[0], strutSet[vItr].position[1] );
  
  fclose(fp);
}


static void
dumpDvdt( search_state *inState, plc_vector* dvdt, int size )

     /* Again, this is debugging code which shouldn't usually be called. */

{
  int vItr;
  FILE* fp;
  char filename[1024];

  sprintf(filename,"%sDvDt.dat",inState->fprefix);
  fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 

  fprintf( fp, "%d\n", size );
  
  for( vItr=0; vItr<size; vItr++ ) {

    fprintf(fp, "%g %g %g\n", plc_M_clist(dvdt[vItr]));

  }

  fclose(fp);
	
}

void dumpLink( plCurve *inLink, search_state *inState, char *dumpname) 

     /* Writes a copy of the link to the current error directory. */
     /* Sets dumpname to the (full) filename of the dumped curve. */
     /* We expect dumpname to be large enough to hold the filename. */

{
  FILE *dumpfile;
  
  sprintf(dumpname,"./%s.rr/dumpLink.vect",inState->fname);
  dumpfile = fopen_or_die(dumpname,"w", __FILE__ , __LINE__ );
  plc_write(dumpfile,inLink);
  fclose(dumpfile);

}

static void
checkDuplicates( octrope_strut* struts, int num )

     /* This is debugging code that isn't called in a standard build. */
     /* Used when one is suspicious about a strut set. */

{
  int i, j;
  for( i=0; i<num; i++ ) {

    // look for the variety of cases that signal a duplicate of strut struts[i] in the
    // rest of the list
    for( j=0; j<num; j++ ) {

      if( j==i )
	continue;
      
      if( struts[i].component[0] == struts[j].component[0] &&
	  struts[i].component[1] == struts[j].component[1] &&
	  
	  struts[i].lead_vert[0] == struts[j].lead_vert[0] &&
	  struts[i].lead_vert[1] == struts[j].lead_vert[1] )
	
	//	fabs(struts[i].position[0] - struts[j].position[0]) < 0.1 &&
	//	fabs(struts[i].position[1] - struts[j].position[1]) < 0.1 )
	{
	  printf( "test1: matching strut!\n" );
	  exit(-1);
	}
      
      // different sense
      if( struts[i].component[0] == struts[j].component[0] &&
	  struts[i].component[1] == struts[j].component[1] &&
	  
	  struts[i].lead_vert[0] == struts[j].lead_vert[1] &&
	  struts[i].lead_vert[1] == struts[j].lead_vert[0] )
	
	//	fabs(struts[i].position[0] - struts[j].position[0]) < 0.1 &&
	//	fabs(struts[i].position[1] - struts[j].position[1]) < 0.1 )
	{
	  printf( "test2: matching strut!\n" );
	  exit(-1);
	}
      
    }
  }
}

static void
collapseStruts( octrope_strut** struts, int* count )

     /* Not called in a standard build of the software, this converts
	all struts to vertex-vertex struts in order to see whether
	struts which are very close to (but not on) vertices cause a
	major change in condition number of matrices. */

{
  // dumb method written to test whether almost-useless struts
  // are causing our ill-conditioned matrix problems

  int initialCount = *count;
  int newCount = 0;
  int sItr, compItr, good;
  octrope_strut* newSet = (octrope_strut*)malloc(sizeof(octrope_strut)*initialCount);
  
  for( sItr=0; sItr<initialCount; sItr++ )
    {
      good = 1;
      for( compItr=0; compItr<newCount; compItr++ )
	{
	  if( (*struts)[sItr].component[0] == (newSet)[compItr].component[0] &&
	      (*struts)[sItr].component[1] == (newSet)[compItr].component[1] &&
	      
	      (*struts)[sItr].lead_vert[0] == (newSet)[compItr].lead_vert[0] &&
	      (*struts)[sItr].lead_vert[1] == (newSet)[compItr].lead_vert[1] )
	    {
	      good = 0;
	      break;
	    }
	}
      
      // not a dup, add to list with edge position rounded
      if( good == 1 )
	{
	  newSet[newCount].component[0] = (*struts)[sItr].component[0];
	  newSet[newCount].component[1] = (*struts)[sItr].component[1];
	  
	  newSet[newCount].lead_vert[0] = (*struts)[sItr].lead_vert[0];
	  newSet[newCount].lead_vert[1] = (*struts)[sItr].lead_vert[1];
	  
	  newSet[newCount].position[0] = (int)(*struts)[sItr].position[0];
	  newSet[newCount].position[1] = (int)(*struts)[sItr].position[1];
	  
	  newCount++;
	}
    }
  
  free(*struts);
  *struts = newSet;
  *count = newCount;
}

