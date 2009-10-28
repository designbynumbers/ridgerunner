/*******************************************************************

   strutgradtest.c : Perform some consistency checks on the minrad gradient 
                  code to make sure that it behaves like the actual gradient
                  of minrad. 

                  We do this by coming up with some 3 vertex vect files, 
                  adjusting the angle between the edges to suit our purposes
                  varying by the minrad gradient and checking minrad before 
                  and after. 

                  If the minrad gradient is correct, then 

                  lim_{h -> 0} (h <var, gradMR> - (MR(x + h grad MR) - MR(x))/h = 0.

                  This is our first test. 

	   ************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<ridgerunner.h>
#include<mangle.h>
#include<argtable2.h>
#include<gsl/gsl_fit.h>

bool VERBOSE = false;
bool gFail = false;

search_state *gState;

void placeContactStruts ( taucs_ccs_matrix* A, plCurve* inLink, 
			  octrope_strut* strutSet, int strutCount,  
			  search_state* inState );

void cvc(search_state *inState,plCurve *inLink,int rownum,int *cmp,int *vert,int *coord);

/* These live in stepper.c and aren't public, so we include a prototype directly. */

/* int octrope_struts(plCurve *L,  */
/*                    const double strut_cutoff, const double epsilon,  */
/*                    octrope_strut *strutlist, int sl_size, */
/*                    double *shortest, void *mem, const int memsize); */

octrope_strut poca(plCurve *L)
/* Computes the poca strut (only) using a clever call to octrope_struts. The magic here*/
/* is that by setting epsilon to 0 we store only the shortest of all the struts. */
{
  octrope_strut poca;
  double pval;
  
  octrope_struts(L,100000,0,&poca,1,&pval,NULL,0);

  return poca;
}

plc_vector *compute_poca_gradient(plCurve *inLink,octrope_strut strut)
/* We do a bit of work to create a fake rigidity matrix, place the strut, and then extract the data. */
/* Uses the global gState */
{
  taucs_ccs_matrix *cleanA;

  cleanA = taucs_ccs_new(12,3*plc_num_verts(inLink),1); /* Creates an nstrut column ccs matrix */
  placeContactStruts(cleanA,inLink,&strut,1,gState);    /* This is the code we're trying to debug, so it's important to call it. */

  /* We now want to extract the data from cleanA into the correct components of first and second. */

  plc_vector *grad;
  grad = calloc(plc_num_verts(inLink),sizeof(plc_vector));

  int j,cptr;
  double val;
  int cp,vt,c;

  for(j=0;j<cleanA->n;j++) { /* Loop over columns (presumably only one) */

    for(cptr=cleanA->colptr[j];cptr<cleanA->colptr[j+1];cptr++) {

      j= cleanA->rowind[cptr];
      val= cleanA->values.d[cptr];
      grad[(int)(floor(j/3.0))].c[j % 3] = val;

      /* We now sanity check the j value. */
      
      cvc(gState,inLink,j,&cp,&vt,&c); 

      switch (cptr % 12) {

      case 0 : 

	assert(cp == strut.component[0]);
	assert(vt == strut.lead_vert[0]);
	assert(c  == 0);
	break; 

      case 1 : 

	assert(cp == strut.component[0]);
	assert(vt == strut.lead_vert[0]);
	assert(c  == 1);
	break; 

      case 2 : 

	assert(cp == strut.component[0]);
	assert(vt == strut.lead_vert[0]);
	assert(c  == 2);
	break; 

      case 3 : 

	assert(cp == strut.component[0]);
	assert(vt == strut.lead_vert[0]+1);
	assert(c  == 0);
	break; 

      case 4 : 

	assert(cp == strut.component[0]);
	assert(vt == strut.lead_vert[0]+1);
	assert(c  == 1);
	break; 

      case 5 : 

	assert(cp == strut.component[0]);
	assert(vt == strut.lead_vert[0]+1);
	assert(c  == 2);
	break; 
      
      case 6 : 

	assert(cp == strut.component[1]);
	assert(vt == strut.lead_vert[1]);
	assert(c  == 0);
	break; 

      case 7 : 

	assert(cp == strut.component[1]);
	assert(vt == strut.lead_vert[1]);
	assert(c  == 1);
	break; 

      case 8 : 

	assert(cp == strut.component[1]);
	assert(vt == strut.lead_vert[1]);
	assert(c  == 2);
	break; 

      case 9 : 

	assert(cp == strut.component[1]);
	assert(vt == strut.lead_vert[1]+1);
	assert(c  == 0);
	break; 

      case 10 : 

	assert(cp == strut.component[1]);
	assert(vt == strut.lead_vert[1]+1);
	assert(c  == 1);
	break; 

      case 12 : 

	assert(cp == strut.component[1]);
	assert(vt == strut.lead_vert[1]+1);
	assert(c  == 2);
	break; 

      }

    }

  }

  /* We have now transcribed the entire column into the grad array. */

  taucs_ccs_free(cleanA);
  return grad;

}


void step( plCurve* inLink, double stepSize, plc_vector* dVdt ); // This is in stepper.c, but not exposed to public.


plCurve *simplecross(double len1, double len2, double angle)

/* generates a test curve consisting of segments of variable lengths in the planes z = 0 and z = 1 */
{
  plCurve *R;
  int nv[2] = {2,2}, cc[2] = {0,0};
  bool open[2] = {true,true};

  R = plc_new(2,nv,open,cc);

  R->cp[0].vt[0] = plc_build_vect(-len1/2.0,0,0);
  R->cp[0].vt[1] = plc_build_vect(len1/2.0,0,0);
  R->cp[1].vt[0] = plc_build_vect(-len2*cos(angle),-len2*sin(angle),1);
  R->cp[1].vt[1] = plc_build_vect(len2*cos(angle),len2*sin(angle),1);

  /* Now we fix gState. */

  if (gState->compOffsets != NULL) { free(gState->compOffsets); gState->compOffsets = NULL; }

  gState->compOffsets = calloc(2,sizeof(int));
  gState->compOffsets[0] = 0;
  gState->compOffsets[1] = 2;

  return R;
}

void test_pocagrad_var(int steps)

/* Perform poca gradient test at angles between 0 and pi. We don't include 0 or pi. */

{
  double PI = 3.14159265358979323846;
  double theta,tstep;
  int i;

  plCurve *L;
  double poca_before,poca_after;
  double hdata[256],quotdata[256];
  int ndata;
  bool newTheta = true;

  printf("Angle  Edgelengths              Poca                       quot fit                                 Result\n");
  printf("----------------------------------------------------------------------------------------------------------\n");

  for(tstep = PI/(steps+2),i=1,theta=tstep;
      i<steps+1;i++,theta+=tstep,newTheta = true) {

    double len1,len2=2.0;

    for(len1=2-1e-6;len1<2+1e-6;len1+=1e-7) {

      L = simplecross(len1,len2,theta);
      poca_before = octrope_poca(L,NULL,0);

      plc_vector *grad;
      grad = compute_poca_gradient(L,poca(L));

      if (newTheta) {

	printf("%6.2g ",theta);
	newTheta = false;
	
      } else {

	printf("       ");

      }

      printf("(%2.8f,%2.8f)  %2.6f  ",
	     len1,len2,poca_before);
      
      if (VERBOSE) {
	
	printf("\n");
	
	octrope_strut Lpoca;

	Lpoca = poca(L);

	printf("\t poca from component %d edge (%d-%d) at %g to component %d edge (%d-%d) at %g, with length %g.\n",
	       Lpoca.component[0],Lpoca.lead_vert[0],Lpoca.lead_vert[0]+1,Lpoca.position[0],
	       Lpoca.component[1],Lpoca.lead_vert[1],Lpoca.lead_vert[1]+1,Lpoca.position[1],
	       Lpoca.length);

	printf("\t poca_gradient = (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f)\n",
	       plc_M_clist(grad[gState->compOffsets[Lpoca.component[0]] + Lpoca.lead_vert[0]]) ,
	       plc_M_clist(grad[gState->compOffsets[Lpoca.component[0]] + (Lpoca.lead_vert[0]+1)]) ,
	       plc_M_clist(grad[gState->compOffsets[Lpoca.component[1]] + Lpoca.lead_vert[1]]) ,
	       plc_M_clist(grad[gState->compOffsets[Lpoca.component[1]] + (Lpoca.lead_vert[1]+1)]));     
      }

      double h,quot;

      plc_vector *stepDir;
      stepDir = calloc(plc_num_verts(L),sizeof(plc_vector));

      double stepNorm = 0;
      int i;
      
      for(i=0;i<plc_num_verts(L);i++) { stepDir[i] = plc_random_vect(); stepNorm += plc_dot_prod(stepDir[i],stepDir[i]);}
      for(i=0;i<plc_num_verts(L);i++) { stepDir[i] = plc_scale_vect(1/stepNorm,stepDir[i]); }

      if (VERBOSE)  
	{ 
	  
	  printf("\t Variation       = (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f),  (%+2.8f,%+2.8f,%+2.8f)\n",
		 plc_M_clist(stepDir[0]),plc_M_clist(stepDir[1]),plc_M_clist(stepDir[2]),plc_M_clist(stepDir[3]));
	  
	}

      double dder = 0;
      
      for(i=0;i<plc_num_verts(L);i++) { dder += plc_dot_prod(stepDir[i],grad[i]); }

      /* We have now computed the directional derivative of poca in the step direction. It remains to actually vary the curve in this 
         direction and see to what extent this predicts the actual change in poca. */
      
      for(h=1e-5,ndata=0;h>1e-9;h/=2.0,ndata++) {
	
	plCurve *workerLink;

	workerLink = plc_copy(L);
	step(workerLink,h,stepDir);
	poca_after = octrope_poca(workerLink,NULL,0);
	
	quot = fabs((h*dder - (poca_after - poca_before))/h);

	if (VERBOSE) {
	  printf("\t\t h = %2.6f delta poca = %1.10f quot = %1.10f \n",h,
		 poca_after - poca_before,quot);
	}

	hdata[ndata] = h;
	quotdata[ndata] = quot;
	
      }
           
      double c0,c1,cov00,cov01,cov11,sumsq;  /* We are going to fit to the linear model y = c0 + c1 h. */
      bool tFail = false;

      gsl_fit_linear(hdata,1,quotdata,1,ndata,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
      printf("%10.4g h + %10.4g (%10.4g) ",c1,c0,sumsq);

      if (fabs(c0) > 1e-4 || fabs(sumsq) > 1e-7) {

	tFail = true;
	gFail = true;
       
      }

      if (tFail) {

	printf("FAIL\n");

      } else {

	printf("pass.\n");

      }
	
    }

    printf("\n");

  }

}

int main(int argc,char *argv[]) 
{
  int nerrors;
  
  struct arg_lit  *verbose = arg_lit0("v","verbose","display lots of human-readable debugging information");
  struct arg_lit  *help    = arg_lit0("h","help","display help message");
  struct arg_end  *end     = arg_end(20);
  
  void *argtable[] = {help,verbose,end};
  
  struct arg_end  *helpend = arg_end(20);
  void *helptable[] = {help, helpend};
    
  printf("strutgradtest, ridgerunner v%s\n\n",PACKAGE_VERSION);
  /* The PACKAGE_VERSION preprocessor symbol is defined in "config.h" by autoconf. */
  
  printf("Test the point of closest approach gradient code in ridgerunner by taking small steps in a randomly chosen direction.\n"
	 "We compute the expected directional derivative and take the difference between the predicted and actual\n"
	 "values of poca, then divide by stepsize (h) to get a quotient. We compute\n"
	 "these values for a range of h from 1e-2 to 1e-8. The resulting (h,quot) data should \n"
	 "fit very well to a line through the origin (the slope should be random). We display the best least-squares\n"
	 "fit below (along with the sum of the squares of the residuals). If the constant term and residuals are both\n"
	 "low, the algorithm will score the test as PASS.\n\n"
	 "Running this program (strutgradtest) from the command line with --verbose will give more details on tests.\n\n");
  
  /* We start by parsing the command-line arguments with argtable. */
  
  if (arg_nullcheck(argtable) != 0 || arg_nullcheck(helptable) != 0)
    printf("error: insufficient memory\n");
  
  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the first set of */
                        /* errors was probably more helpful, so we display it. */
      
      arg_print_errors(stdout,helpend,"helptable");
      arg_print_errors(stdout,end,"strutgradtest");
      exit(1);
      
    } else {  /* The help table matched, which means we asked for help or gave nothing */
      
      printf("strutgradtest performs self-tests on the poca gradient code in ridgerunner.\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);
      
    }
    
  }
  
  if (verbose->count > 0) {

    VERBOSE = true;

  }
  
  gState = calloc(1,sizeof(search_state));
  gLambda = 1.0;


  test_pocagrad_var(100);

  if (gFail) { exit(1); } 
  else {exit(0);}
  
}
