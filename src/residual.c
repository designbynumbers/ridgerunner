/* 
 * residual.c : This short program computes the residual of a given vect file.
                We link with a number of the ridgerunner object files to do this.
 *
 */

#include "portability.h"
#include "octrope.h"
#include "ridgerunner.h"
#include <argtable2.h>

int main(int argc,char *argv[]) {

  FILE *infile_fptr = {NULL};
  double lambda = 1.0,tuberad = 0.5;
  plCurve *L;

  int nerrors;
    
  struct arg_file *infile = arg_filen(NULL,NULL,"<VECT file>", 1, 10000, "input VECT file(s)");
  struct arg_dbl  *tube_radius = arg_dbl0("r","radius","<x>","radius of tube");
  struct arg_dbl  *arg_lambda = arg_dbl0("l","lambda","<x>","minimum diameter of curvature");
  struct arg_lit  *help        = arg_lit0("h","help","display help message");
  struct arg_end  *end = arg_end(20);

  void *argtable[] = {infile,tube_radius,arg_lambda,help,end}; 
  
  struct arg_end  *helpend = arg_end(20);

  void *helptable[] = {help, helpend};
  int filenum;
  char revision[20] = "$Revision: 1.20 $";
  char *dollar;

  char winning_file[1024] = "No file";
  double winning_residual = 2, winning_ropelength = -1, winning_thickness = -1;
  int    winning_strutcount = -1, winning_mrstruts  = -1;

  dollar = strchr(&revision[1],'$');
  dollar[0] = '\0';
  printf("residual v%s, ridgerunner v%s\n",&revision[11],PACKAGE_VERSION);
 /* The PACKAGE_VERSION preprocessor symbol is defined in "config.h" by autoconf. */

  printf("  computes the ridgerunner residual of a polygonal knot.\n");
  
  /* We start by parsing the command-line arguments with argtable. */
  
  if (arg_nullcheck(argtable) != 0 || arg_nullcheck(helptable) != 0)
    printf("error: insufficient memory\n");

  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */

    nerrors = arg_parse(argc,argv,helptable);

    if (nerrors > 0) {  /* The help table didn't match either-- the first set of */
                        /* errors was probably more helpful, so we display it. */

      arg_print_errors(stdout,helpend,"helptable");
      arg_print_errors(stdout,end,"struts");
      exit(1);

    } else {  /* The help table matched, which means we asked for help or gave nothing */
  
      printf("residual computes the fraction of the gradient of length for a polygonal knot \n"
	     "which is unresolved by strut and kink compression forces. This provides a measure\n"
	     "of how close the configuration is to being (polygonal) ropelength-critical.\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }
  
  /* Convert command-line arguments into proper global variables */
  
  if (arg_lambda->count > 0) {

    lambda = arg_lambda->dval[0];
    gLambda = lambda;

  }

  if (tube_radius->count > 0) {

    tuberad = tube_radius->dval[0];

  }
  
  search_state inState;
  inState.dumpAxb = 0;
  inState.snapinterval = 1000;
  inState.steps = 1; /* Make sure that we don't try to create a snapshot! */
  inState.lastStepStruts = NULL;
  inState.lastStepMRlist = NULL;
  inState.tube_radius = 0.5;

  const char *FILENAME;

  for (filenum = 0;filenum < infile->count;filenum++) {
    
    FILENAME = infile->filename[filenum];
    
    infile_fptr = fopen(infile->filename[filenum],"r");
    
    if (infile_fptr == NULL) {
      
      fprintf(stderr,"residual: Couldn't open file %s.\n",infile->filename[filenum]);
      continue;
      
    }
    
    octrope_error_num = 0;
    L = plc_read(infile_fptr,&octrope_error_num,octrope_error_str,80);
    
    /* We now demonstrate the octrope library's error handling protocol: */
    
    if (octrope_error_num > 0) {   /* This is the signal for an error. */
      
      fprintf(stderr,"residual: link reading error\n%s\n",octrope_error_str);
      continue;
      
    }
    
    fclose(infile_fptr);
  
    /* The actual computation is performed as a side effect of stepDirection. */

    stepDirection(L,tuberad,0,lambda,&inState);

    /* We now have the data inside inState. */

    printf("%s Rop: %g Thi: %g Struts: %d MrStruts: %d \nResidual: %g \n",
	   infile->basename[filenum],inState.ropelength,inState.thickness,
	   inState.lastStepStrutCount,inState.lastStepMinradStrutCount,inState.residual);

    plc_free(L);

    if (infile->count > 1) { /* Keep a record of the minimum residual file among those processed */

      if (inState.residual < winning_residual) { 

	strncpy(winning_file,infile->basename[filenum],sizeof(winning_file));
	winning_residual = inState.residual;
	winning_ropelength = inState.ropelength;
	winning_thickness = inState.thickness;
	winning_strutcount = inState.lastStepStrutCount;
	winning_mrstruts = inState.lastStepMinradStrutCount;

      }

    }

  }

  if (infile->count > 1) {

    printf("----------------------------------------------------------\n");
    printf("         Least Residual among %d files checked            \n",infile->count);
    printf("----------------------------------------------------------\n");

    printf("%s Rop: %g Thi: %g Struts: %d MrStruts: %d \nResidual: %g \n",
	   winning_file,winning_ropelength,winning_thickness,winning_strutcount,winning_mrstruts,winning_residual);

    printf("\n\n");

  }

  exit(0);
  
}
	   
