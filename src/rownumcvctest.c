/*******************************************************************

   rownumcvctest.c : 

   A big part of the ridgerunner algorithm involves converting knots
   and links from the structured format of plCurve to a flat format
   consisting of a single buffer of plc_vectors or doubles.

   Ridgerunner provides several utility functions to accomplish this
   conversion:

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include"ridgerunner.h"
#include"mangle.h"
#include<argtable2.h>

int rownum(plCurve *inLink, int cmp, int vert, int coord);

/* Computes the row number in the rigidity matrix corresponding
   to the given data. Handles wraparound correctly, and does
   error checking. We should use this to get stuff into rows of
   vectors and the matrix A exclusively, rather than try to
   recode the computation at various places in the program. */

void cvc(plCurve *inLink,int rownum,int *cmp,int *vert,int *coord);

/* Computes the component, vector and coordinate corresponding to a given row
   in the rigidity matrix, using the componentOffsets in inState. cvc should be 
   the inverse of rownum. */

/* This test function creates some plCurves with random numbers of components and vertices
   and iterates over all the coordinates of all the vertices, checking that rownum and cvc
   are indeed inverses of one another. */


bool VERBOSE = false;
bool gFail = false;

bool test_cvc_rownum(plCurve *L)
{
  plCurve *link;
  int ntests = 0;

  if (L == NULL) { 

    /* Make up a link from scratch. */

    int nv[10],cc[10];
    bool open[10];

    int nc,cp;

    nc = rand() % 9;

    for(cp=0;cp < nc;cp++) {

      nv[cp] = rand() % 213 + 1;
      cc[cp] = 0;
      open[cp] = (rand() % 2 == 0);

    }

    link = plc_new(cp,nv,open,cc);

    assert(link != NULL);

  } else {

    link = L;

  }

  /* Now we have a link. We go ahead and test. */

  int cp,vt,c;
  int rcp,rvt,rc,rnum;

  for(cp=0;cp<link->nc;cp++) {

    for(vt=0;vt<link->cp[cp].nv;vt++) {

      for(c=0;c<3;c++) {

	rnum = rownum(link,cp,vt,c);
	cvc(link,rnum,&rcp,&rvt,&rc);

	if (rcp != cp || rvt != vt || rc != c) {

	  gFail = true;
	  printf("Rownumber %d for (cmp,vert,coord) = (%d,%d,%d) fails to match reported (cmp,vert,coord) of (%d,%d,%d).  FAIL.\n",
		 rnum,cp,vt,c,rcp,rvt,rc);
	  return false;

	}

	ntests++;

      }

    }

  }
  
  printf("Passed %d rownum/cvc consistency checks.\n",ntests);
  return true;

}
	  	

int main(int argc,char *argv[]) 
{
  int nerrors;
  
  struct arg_file *arg_infile =  arg_filen(NULL,NULL,"<file>",0,100000,"input files");
  struct arg_lit  *verbose = arg_lit0("v","verbose","display lots of human-readable debugging information");
  struct arg_lit  *help    = arg_lit0("h","help","display help message");
  struct arg_end  *end     = arg_end(20);
  
  void *argtable[] = {arg_infile,help,verbose,end};
  
  struct arg_end  *helpend = arg_end(20);
  void *helptable[] = {help, helpend};
    
  printf("rownumcvctest, ridgerunner v%s\n\n",PACKAGE_VERSION);
  /* The PACKAGE_VERSION preprocessor symbol is defined in "config.h" by autoconf. */
  
  printf("Test the code in ridgerunner for converting to and from component, vertex, coordinate notation to row numbers in\n"
	 "the rigidity matrix. Can be run either on an input file or (if no input file is given) on randomly generated\n"
	 "plCurves.\n\n"
	 "Running this program (rownumcvctest) from the command line with --verbose will give more details on tests.\n\n");
  
  /* We start by parsing the command-line arguments with argtable. */
  
  if (arg_nullcheck(argtable) != 0 || arg_nullcheck(helptable) != 0)
    printf("error: insufficient memory\n");
  
  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the first set of */
                        /* errors was probably more helpful, so we display it. */
      
      arg_print_errors(stdout,helpend,"helptable");
      arg_print_errors(stdout,end,"rownumcvctest");
      exit(1);
      
    } else {  /* The help table matched, which means we asked for help or gave nothing */
      
      printf("rownumcvctest performs self-tests on the rigidity matrix building code in ridgerunner.\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);
      
    }
    
  }
  
  if (verbose->count > 0) {

    VERBOSE = true;

  }
  
  gLambda = 1.0;

  if (arg_infile->count == 0) {

    test_cvc_rownum(NULL);

  } else {

    printf("Filename                  Result\n");
    printf("--------------------------------\n");

    int infilenum;
    
    for(infilenum=0;infilenum < arg_infile->count;infilenum++) {
      
      plCurve *link;
      FILE *infile_fptr;
      
      infile_fptr = fopen(arg_infile->filename[infilenum],"r");
      
      if (infile_fptr == NULL) {
	
	fprintf(stderr,"rownumcvctest: Couldn't open file %s.\n",
		arg_infile->filename[infilenum]);
	continue;  /* Try the next file */
	
      }
      
      int plr_error_num;
      char plr_error_str[1024];
      
      link = plc_read(infile_fptr,
		      &plr_error_num,plr_error_str,sizeof(plr_error_str));
      
      /* We now demonstrate the octrope library's error handling protocol: */
      
      if (plr_error_num > 0) {   /* This is the signal for an error. */
	
	fprintf(stderr,"rownumcvctest: link reading error\n%s\n",plr_error_str);
	continue;  /* Try the next file */
	
      }
      
      fclose(infile_fptr);

      printf("%-22s   ",arg_infile->basename[infilenum]);
      
      if (test_cvc_rownum(link)) { printf("   pass\n"); }
      else { printf("   FAIL.\n"); }
      
    }

  }
    
  if (gFail) { exit(1); } 
  else {exit(0);}
  
}
