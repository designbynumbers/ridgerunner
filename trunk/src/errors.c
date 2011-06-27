/*
 *  errors.c
 *  ridgerunner3
 *
 *  Created by Michael Piatek on Fri Jan 16 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "errors.h"
#include "ridgerunner.h"


void 
DebugThrow( int inErr, const char* inFile, long inLine )
{
  char errmsg[1024];

  if (inErr == kNULLPointer) {

    sprintf(errmsg,"ridgerunner: Null pointer error.\n");

  } else {
    
    sprintf(errmsg,"ridgerunner: Error number %d.\n",inErr);

  }

  FatalError(errmsg, inFile, inLine);

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
  /* We may not have lived long enough to open gLogfile. */

  if (gLogfile != NULL) {

    fprintf(gLogfile,"ridgerunner: Fatal error in file %s, line %d.\n",file,line);
    fprintf(gLogfile,"%s",debugmsg);

  }
  
  fprintf(stderr,"ridgerunner: Fatal error in file %s, line %d.\n",file,line);
  fprintf(stderr,"%s",debugmsg);

#ifdef HAVE_ASCTIME
#ifdef HAVE_LOCALTIME
#ifdef HAVE_TIME
  
  time_t end_time;
  end_time = time(NULL);
  
  if (gLogfile != NULL) {

    fprintf(gLogfile,
	    "Run ended: %s.\n",
	    asctime(localtime(&end_time)));
    
  }

  fprintf(stderr,
	  "Run ended: %s.\n",
	  asctime(localtime(&end_time)));
  
  
#endif
#endif
#endif

  fflush(stdout);
  fflush(stderr);

  if (gLogfile != NULL) {

    fflush(gLogfile);
    fclose(gLogfile);

  }

  exit(1);
  
}

void NonFatalError(char *debugmsg,const char *file,int line)

     /* Write the error message to the global logfile
	and to stderr, then return. */

{
  /* We may not have lived long enough to open gLogfile. */

  if (gLogfile != NULL) {

    fprintf(gLogfile,"ridgerunner: Error in file %s, line %d.\n",file,line);
    fprintf(gLogfile,"%s",debugmsg);

  }
  
  fprintf(stderr,"ridgerunner: Error in file %s, line %d.\n",file,line);
  fprintf(stderr,"%s",debugmsg);

#ifdef HAVE_ASCTIME
#ifdef HAVE_LOCALTIME
#ifdef HAVE_TIME
  
  time_t end_time;
  end_time = time(NULL);
  
  if (gLogfile != NULL) {

    fprintf(gLogfile,
	    "Error logged: %s.\n",
	    asctime(localtime(&end_time)));
    
  }

  fprintf(stderr,
	  "Error logged: %s.\n",
	  asctime(localtime(&end_time)));
  
  
#endif
#endif
#endif
  
}

FILE *fopen_or_die(const char *filename,const char *mode,
		   const char *file,const int line) 

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

int system_or_die(char *cmdline,const char *file,int line)

{

  int sysresult=0;
  char errmsg[1024];

  sysresult = system(cmdline);

  if (sysresult != 0) {

    sprintf(errmsg,
	    "ridgerunner: Command line \n"
	    "               %s \n"
	    "             failed with result %d.\n",
	    cmdline,sysresult);

    FatalError(errmsg,file,line);
    
  }

  return sysresult;
}

void *malloc_or_die(size_t size, const char *file, const int line)
     
     /* Allocates memory or dies trying. */

{
  void *mret;
  char  errmsg[1024];

  mret = malloc(size);

  if (mret == NULL) {

    sprintf(errmsg,"ridgerunner: Could not allocate block of size %d.\n",
	    (int)(size));
    FatalError(errmsg,file,line);

  }

  return mret;
}

void remove_or_die(char *filename,const char *file, const int line)

/* Deletes the file "filename" or dies trying. */

{
  
  int retcode;
  char errmsg[1024];

  retcode = remove(filename);

  if (retcode != 0) {

    sprintf(errmsg,
	    "ridgerunner: Couldn't remove file %s due to %s\n",
	    filename,strerror(retcode));
    
    FatalError(errmsg,file,line);

  }
  
}

void rename_or_die(char *oldname,char *newname,const char *file, const int line)

/* Deletes the file "filename" or dies trying. */

{
  
  int retcode;
  char errmsg[1024];

  retcode = rename(oldname,newname);

  if (retcode == -1) {

    sprintf(errmsg,
	    "ridgerunner: Couldn't rename file %s to %s due to %s\n",
	    oldname,newname,strerror(retcode));
    
    FatalError(errmsg,file,line);

  }
  
}

DIR *opendir_or_die(char *dirname,const char *file, const int line)

/* Opens the directory "dirname" or dies trying. */

{
  
  DIR *ret;
  char errmsg[1024];

  ret = opendir(dirname);

  if (ret == NULL) {

    sprintf(errmsg,
	    "ridgerunner: Couldn't open directory %s due to %s\n",
	    dirname,strerror(errno));
    
    FatalError(errmsg,file,line);

  }
  
  return ret;

}

int   mkstemp_or_die(char *template,const char *file, const int line) 

/* Construct a temporary file from the given template or die trying. */

{

  int retcode;
  char errmsg[1024];

  retcode = mkstemp(template);

  if (retcode == -1) {

    sprintf(errmsg,
	    "ridgerunner: Couldn't mkstemp file based on template %s\n",
	    template);
    
    FatalError(errmsg,file,line);

  }

  return retcode;
  
}

FILE *fdopen_or_die(int fd, const char *opentype, const char *file, const int line)

/* Open a stream from the given file descriptor or die trying. */

{
  FILE *retfile;
  char errmsg[1024];

  retfile = fdopen(fd,opentype);

  if (retfile == NULL) {

    sprintf(errmsg,
	    "ridgerunner: Couldn't associate stream to fd %d.\n",
	    fd);
    
    FatalError(errmsg,file,line);

  }

  return retfile;
  
}  

void logprintf(char *format, ... )

/* Function prints a message both to the screen and the logfile. */
/* If the curses interface is running, will (eventually) write   */
/* to an appropriate area of the screen. */
  
/* For this function we flush as soon as we write if VERBOSITY >= 10
   since we're trying to detect crashes of the system then. */

/* In any case, we flush after every 5 writes. */

{
  va_list args;
  char msgbuf[2048];
  static int nwrites=0;
  
  va_start(args,format);
  vsprintf(msgbuf,format,args);
  va_end(args);

#ifdef CURSES_DISPLAY

  if (gLogwin != NULL) { wprintw(gLogwin,"%s",msgbuf); }
  else { printf("%s",msgbuf); }	/* The error might be in the windowing system! */

#else

  printf("%s",msgbuf);
  if (VERBOSITY >= 10) {fflush(stdout);}

#endif
  
  if (gLogfile != NULL) {

    fprintf(gLogfile,"%s",msgbuf);
    if (VERBOSITY >= 10 || nwrites%5 == 0) {fflush(gLogfile);}
  }

  nwrites++;
}

void
dumpAxb_full( search_state *inState, 
	      double* A, int rows, int cols, 
	      double* x, double* b )
{
  /* Saves a problem Ax = b in MATLAB or Octave format */
  
  FILE* fp;
  int rItr, cItr;
  char filename[1024];
  
  /* Construct filename for A */

  if (A != NULL) {

    sprintf(filename,"%sA.mat",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    
    /* Now construct the file. */
    
    fprintf(fp,
	    "# Created by ridgerunner (dumpAxb_full) \n"
	    "# name: A\n"
	    "# type: matrix\n"
	    "# rows: %d\n"
	    "# columns: %d\n",rows,cols);
    
    for( rItr=0; rItr<rows; rItr++ ) {
      
      for( cItr=0; cItr<cols; cItr++ )
	fprintf( fp, "%10.16lf ", A[rItr*cols + cItr] );
      
      fprintf( fp, "\n" );
    }
    
    fclose(fp);

  }

  /* Construct filename for x. */

  if (x != NULL) {
    
    sprintf(filename,"%sx.mat",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    colvector_write_mat(fp,x,cols,"x");
    fclose(fp);

  }

  /* Construct filename for b. */

  if (b != NULL) {

    sprintf(filename,"%sb.mat",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    colvector_write_mat(fp,b,cols,"b");    
    fclose(fp);

  }

}

void
dumpAxb_sparse( search_state *inState, taucs_ccs_matrix* A, double* x, double* b )

/* This is debugging code which is only called as part of an exit from error. */

{
  FILE* fp;
  char filename[1024];
  
  /* Construct filename for A */

  if (A != NULL) {

    sprintf(filename,"%sA.sparse",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    taucs_ccs_write_sparse(fp,A);
    fclose(fp);

    sprintf(filename,"%sA.mat",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    taucs_ccs_write_mat(fp,A);
    fclose(fp);

    sprintf(filename,"%sA.dat",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    taucs_ccs_write_dat(fp,A);
    fclose(fp);

  }

  /* Construct filename for x. */

  if (x != NULL) {
    
    sprintf(filename,"%sx.mat",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    colvector_write_mat(fp,x,A->n,"x");
    fclose(fp);

    sprintf(filename,"%sx.dat",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    colvector_write_dat(fp,x,A->n,"x");
    fclose(fp);


  }

  /* Construct filename for b. */

  if (b != NULL) {

    sprintf(filename,"%sb.mat",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    colvector_write_mat(fp,b,A->m,"b");    
    fclose(fp);

    sprintf(filename,"%sb.dat",inState->fprefix);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    colvector_write_dat(fp,b,A->m,"b");    
    fclose(fp);

  }

}

 
void
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

void 
dumpStruts( plCurve *inLink, search_state *inState, char *dumpname)

     /* Writes the current strut set to a file in the run directory. */
     /* Assumes that the strut set is given inside inState. */
{
  FILE *strutsVect, *strutsTed;
  char tedname[1024];

  sprintf(dumpname,"%s%s.struts.vect",inState->fprefix,inState->basename);
  strutsVect = fopen_or_die(dumpname,"w", __FILE__ , __LINE__ );   
  strut_vectfile_write(inLink,
		       inState->lastStepStruts,inState->lastStepStrutCount,
		       strutsVect);
  fclose(strutsVect);

  sprintf(tedname,"%s%s.struts",inState->fprefix,inState->basename);
  strutsTed = fopen_or_die(tedname,"w", __FILE__ , __LINE__ );   
  octrope_strutfile_write(inState->lastStepStrutCount,
			  inState->lastStepStruts,
			  inState->lastStepMinradStrutCount,
			  inState->lastStepMRlist,
			  strutsTed);
  fclose(strutsTed);
}

void
dumpDvdt( plc_vector* dvdt, plCurve *inLink, search_state *inState )

     /* Again, this is debugging code which shouldn't usually be called. */

{
  FILE* fp;
  char filename[1024];
  plCurve *dVdTcurve;
  plc_color blue;

  sprintf(filename,"%sDvDt.vect",inState->fprefix);
  fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 

  dVdTcurve = vectorfield_to_plCurve(dvdt,inLink);
  blue = plc_build_color(0,0,1,1);
  plc_set_color(dVdTcurve,blue);

  plc_write(fp,dVdTcurve);
  plc_free(dVdTcurve);

}

void
dumpdLen( plc_vector* dLen, plCurve *inLink, search_state *inState )

     /* Again, this is debugging code which shouldn't usually be called. */

{
  FILE* fp;
  char filename[1024];
  plCurve *dLencurve;
  plc_color red;

  sprintf(filename,"%sdLen.vect",inState->fprefix);
  fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 

  dLencurve = vectorfield_to_plCurve(dLen,inLink);
  red = plc_build_color(1,0,0,1);
  plc_set_color(dLencurve,red);

  plc_write(fp,dLencurve);
  plc_free(dLencurve);

}

void snapshot( plCurve *inLink,
	       plc_vector *dVdt,plc_vector *dlen,
	       search_state *inState )

{
  FILE *fp,*gfp;
  char filename[1024];
  plCurve *VFcurve;
  plc_color col;
  taucs_ccs_matrix*  A = NULL;
 
  /* We drop a complete snapshot of the link, together with a geomview
     file tying everything together, to the snapshot directory. */

  if (VERBOSITY >= 10) { logprintf("Snapshot at step %d.\n",inState->steps); }

  sprintf(filename,"%s.%d.geom",inState->snapprefix,inState->steps);
  gfp = fopen_or_die(filename,"w", __FILE__ , __LINE__ );
  fprintf(gfp,"LIST {\n");

  sprintf(filename,"%s.%d.vect",inState->snapprefix,inState->steps);
  fprintf(gfp,"\t{ < %s.%d.vect }\n",inState->basename,inState->steps);
  fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ );
  plc_write(fp,inLink);
  fclose(fp);

  if (inState->lastStepStrutCount > 0) {

    sprintf(filename,"%s.%d.struts.vect",inState->snapprefix,inState->steps);
    fprintf(gfp,"\t{ < %s.%d.struts.vect }\n",inState->basename,inState->steps);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ );
    strut_vectfile_write(inLink,inState->lastStepStruts,
			 inState->lastStepStrutCount,fp);
    fclose(fp);

  }

  sprintf(filename,"%s.%d.dlen.vect",inState->snapprefix,inState->steps);
  fprintf(gfp,"\t{ < %s.%d.dlen.vect }\n",inState->basename,inState->steps);
  fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ );
  VFcurve = vectorfield_to_plCurve(dlen,inLink);
  col = plc_build_color(0,0,1,1);
  plc_set_color(VFcurve,col);
  plc_write(fp,VFcurve);
  plc_free(VFcurve);
  fclose(fp);

  sprintf(filename,"%s.%d.dVdt.vect",inState->snapprefix,inState->steps);
  fprintf(gfp,"\t{ < %s.%d.dVdt.vect }\n",inState->basename,inState->steps);
  fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ );
  VFcurve = vectorfield_to_plCurve(dVdt,inLink);
  col = plc_build_color(1,0,0,1);
  plc_set_color(VFcurve,col);
  plc_write(fp,VFcurve);
  plc_free(VFcurve);
  fclose(fp);

  fprintf(gfp,"}\n");
  fclose(gfp);

  /* Part of the snapshot is a drop of the current rigidity matrix A */

  A = buildRigidityMatrix(inLink,inState->tube_radius,gLambda,inState);

  if (A != NULL && A->n > 0 && A->m > 0) {  /* There has to BE a matrix to drop */

    sprintf(filename,"%s.%d.A.sparse",inState->snapprefix,inState->steps);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    taucs_ccs_write_sparse(fp,A);
    fclose(fp);

    sprintf(filename,"%s.%d.A.mat",inState->snapprefix,inState->steps);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    taucs_ccs_write_mat(fp,A);
    fclose(fp);

    sprintf(filename,"%s.%d.A.dat",inState->snapprefix,inState->steps);
    fp = fopen_or_die(filename,"w", __FILE__ , __LINE__ ); 
    taucs_ccs_write_dat(fp,A);
    fclose(fp);

  }

  taucs_ccs_free(A);
  
}
  
void dumpLink( plCurve *inLink, search_state *inState, char *dumpname) 

     /* Writes a copy of the link to the current error directory. */
     /* Sets dumpname to the (full) filename of the dumped curve. */
     /* We expect dumpname to be large enough to hold the filename. */

{
  FILE *dumpfile;
  
  sprintf(dumpname,"%s%s.dump.vect",inState->fprefix,inState->basename);
  dumpfile = fopen_or_die(dumpname,"w", __FILE__ , __LINE__ );
  plc_write(dumpfile,inLink);
  fclose(dumpfile);

}

void dumpNamedLink( plCurve *inLink, search_state *inState, char *tag, char *dumpname)

/* Writes a copy of the link to the current error directory. */
/* Sets dumpname to the (full) filename of the dumped curve. */
/* We expect dumpname to be large enough to hold the filename. */
  
{
  FILE *dumpfile;
  
  sprintf(dumpname,"%s%s.%s.vect",inState->fprefix,inState->basename,tag);
  dumpfile = fopen_or_die(dumpname,"w", __FILE__ , __LINE__ );
  plc_write(dumpfile,inLink);
  fclose(dumpfile);

}


// this is debugging code which is not used in a standard build
double
rigidityEntry( taucs_ccs_matrix* A, int strutCount, int minradLocs, int row, int col )
{
  int vItr;
  
  // the most common case will be a 0, so get that out of the way quick --
  // movd this check to firstVariation for speed purposes
  
  // also note that this exploits the t_snnls enforced ordering of row
  // entries and doesn't work for generic ccs matrices
  /*	if( row < A->rowind[A->colptr[col]] || row > A->rowind[A->colptr[col+1]-1] )
	return 0;
  */
  // now we look.
  if( col < strutCount ) {

    // we can manually unroll this loop -- data dependencies in here
    // result in large slowdown according to shark.
    /*		for( vItr=A->colptr[col]; vItr<A->colptr[col+1]; vItr++ )
		{
		if( A->rowind[vItr] == row )
		return A->values.d[vItr];
		}
    */
    /* Remember that if this column represents a strut, there at most
       12 rows in the column (less only if the strut is vertex-vertex,
       or vertex-edge).  So we need only search through 12 possible
       rowind values looking for the desired row. 

       If there are too few rows in this column, the search could in 
       principle wrap around and return a row from the next column.*/

    vItr=A->colptr[col];

    if( A->rowind[vItr+0] == row ) return A->values.d[vItr+0];
    if( A->rowind[vItr+1] == row ) return A->values.d[vItr+1];
    if( A->rowind[vItr+2] == row ) return A->values.d[vItr+2];
    if( A->rowind[vItr+3] == row ) return A->values.d[vItr+3];
    
    if( A->rowind[vItr+4] == row ) return A->values.d[vItr+4];
    if( A->rowind[vItr+5] == row ) return A->values.d[vItr+5];
    if( A->rowind[vItr+6] == row ) return A->values.d[vItr+6];
    if( A->rowind[vItr+7] == row ) return A->values.d[vItr+7];
    
    if( A->rowind[vItr+8] == row ) return A->values.d[vItr+8];
    if( A->rowind[vItr+9] == row ) return A->values.d[vItr+9];
    if( A->rowind[vItr+10] == row ) return A->values.d[vItr+10];
    if( A->rowind[vItr+11] == row ) return A->values.d[vItr+11];
 
  } else {

    /* This is not a strut constraint, so we're not sure how many 
       columns to search. Check them all just to be sure. */
    
    for( vItr=A->colptr[col]; vItr<A->colptr[col+1]; vItr++ ) {

      if( A->rowind[vItr] == row )
	return A->values.d[vItr];
    }
  }

  /* We should never reach this point in the code. */

  char errmsg[1024];

  sprintf(errmsg,
	  "ridgerunner: illegal call to rigidityEntry\n"
	  "             tried for entry (%d,%d) of %d x %d matrix A\n"
	  "             strutcount %d\n"
	  "             minradLocs %d\n",
	  row,col,A->n,A->m,strutCount,minradLocs);

  FatalError(errmsg, __FILE__ , __LINE__ );

  return 0;
}


void
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

void
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

