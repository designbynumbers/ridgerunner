/*******************************************************************

   mrgradtest.c : Perform some consistency checks on the minrad gradient 
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

#include"ridgerunner.h"
#include"mangle.h"
#include<argtable2.h>
#include<gsl/gsl_fit.h>

bool VERBOSE = false;
bool gFail = false;

void compute_minrad_gradient(plCurve *inLink,octrope_mrloc mrloc,plc_vector *AAs,plc_vector *BBs,plc_vector *CCs);
/* This lives in stepper.c and isn't a public function, so we don't have a prototype in ridgerunner.h */

void compute_pm_minrad_gradient(plCurve *inLink,int cItr,int vItr,
				plc_vector pGrad[3],plc_vector mGrad[3]) 
				 
/* Returns the gradient of plus and minus minrad at inLink->cp[cItr].vt[vItr]. Each gradient has three vectors which 
   describe the motion of the vertex before the target [0], the target [1], and the vertex after the
   target [2]. */

{
  plc_vector B, A, cross;
  double bmag, amag;
  double angle;
  double kappa, prevLen, thisLen;
  plc_vector  prevSide, thisSide, N, fancyL, fancyM, fancyN;
 		
  prevSide = plc_vect_diff(inLink->cp[cItr].vt[vItr],inLink->cp[cItr].vt[vItr-1]);
  thisSide = plc_vect_diff(inLink->cp[cItr].vt[vItr+1],inLink->cp[cItr].vt[vItr]);
  
  /* dot = plc_M_dot(prevSide, thisSide); */
  
  prevLen = plc_M_norm(prevSide);
  thisLen = plc_M_norm(thisSide);
  
  // B = b-v = thisSide. 
  
  B = thisSide;
  A = plc_scale_vect(-1,prevSide);
  
  bmag = plc_M_norm(B);
  amag = plc_M_norm(A);
  
  /* value = dot/(prevLen*thisLen); */
  bool ok = true;
  angle = plc_angle(prevSide,thisSide,&ok);
  
  if (!ok) {
    
    char errmsg[1024];
    sprintf(errmsg,
	    "ridgerunner: Couldn't compute turning angle at vertex %d of cmp %d\n"
	    "             of polyline. \n",
	    vItr,cItr);
    FatalError(errmsg, __FILE__ , __LINE__ );
    
  }
  
  /* We now check that angle is high enough for the following 
     stuff to work. */
  
  double PI = 3.14159265358979323846;
  
  if (angle < 1e-12 || angle > PI - 1e-12) {
    
    char errmsg[1024];
    sprintf(errmsg,
	    "ridgerunner: Can't compute minrad gradient when edges are\n"
	    "             almost colinear. Angle between edges is %g.\n",
	    angle);
    
    FatalError(errmsg, __FILE__ , __LINE__ );
    
  }
  
  /* First, we compute the "plus" gradient, depending on the length of the B or second edge */
    
  // says... maple?
  kappa = -bmag/(2-2*cos(angle));
  
  // BxA
  
  plc_M_cross(N,B,A);
  N = plc_normalize_vect(N,&ok);
  assert(ok);
  
  double Lconst, Mconst, Nconst;
  
  Lconst = (1/(2*tan(angle/2) * bmag));	
  fancyL = plc_scale_vect(Lconst,B);
  
  Mconst = kappa*(1/(amag*amag));
  // A x N
  
  plc_M_cross(cross,A,N);
  fancyM = plc_scale_vect(Mconst,cross);
  
  Nconst = kappa*(1/(bmag*bmag));
  // N x B
  
  plc_M_cross(cross,N,B);
  fancyN = plc_scale_vect(Nconst,cross);
  
  pGrad[0] = fancyM;
  
  pGrad[1].c[0] = -fancyM.c[0] - fancyN.c[0] - fancyL.c[0];
  pGrad[1].c[1] = -fancyM.c[1] - fancyN.c[1] - fancyL.c[1];
  pGrad[1].c[2] = -fancyM.c[2] - fancyN.c[2] - fancyL.c[2];
  
  /* This is really the fastest way to accomplish this operation. The plCurve */
  /* option would be the code: */
  
  /* Bs = plc_scale_vect(-1,plc_vect_sum(plc_vect_sum(fancyM,fancyN),fancyL)); */
  
  /* which is really pretty awful. */
  
  pGrad[2] = plc_vect_sum(fancyN,fancyL);
    
  /* We now compute the "minus" gradient, where the length is taken from the first (A) edge. */
    
  // says... maple?
  kappa = -amag/(2-2*cos(angle));
  
  // BxA
  
  plc_M_cross(N,B,A);
  N = plc_normalize_vect(N,&ok);
  assert(ok);
  
  //double Lconst, Mconst, Nconst;
  
  Lconst = (1/(2*tan(angle/2) * amag));
  
  fancyL = plc_scale_vect(Lconst,A);
  Mconst = kappa*(1/(amag*amag));
  
  // A x N
  
  plc_M_cross(cross,A,N);
  fancyM = plc_scale_vect(Mconst,cross);
  
  Nconst = kappa*(1/(bmag*bmag));
  // N x B
  plc_M_cross(cross,N,B);
  fancyN = plc_scale_vect(Nconst,cross);
  
  mGrad[0] = plc_vect_sum(fancyM,fancyL);
  
  mGrad[1].c[0] = -fancyM.c[0] - fancyN.c[0] - fancyL.c[0];
  mGrad[1].c[1] = -fancyM.c[1] - fancyN.c[1] - fancyL.c[1];
  mGrad[1].c[2] = -fancyM.c[2] - fancyN.c[2] - fancyL.c[2];
  
  /* See above. This is better than using plCurve. */
  
  mGrad[2] = fancyN;
  
}

double pmMinrad(plCurve *L, int cp, int vt,double *pminrad,double *mminrad) 

/* Compute the plus and minus minrads for a given vertex on a plCurve. */
{
  double mr = {DBL_MAX}, rad;
  plc_vector in,out;
  double normin, normout;
  double dot_prod,cross_prod_norm;
  
  /* We start with some initializations. */
  
  in  = plc_vect_diff(L->cp[cp].vt[vt],L->cp[cp].vt[vt-1]);
  out = plc_vect_diff(L->cp[cp].vt[vt+1],L->cp[cp].vt[vt]);

  normin  = plc_norm(in);
  normout = plc_norm(out);
      
  dot_prod = plc_M_dot(in,out);
  cross_prod_norm = plc_norm(plc_cross_prod(in,out));
  
  rad = (normin*normout + dot_prod)/(2*cross_prod_norm);
  
  if (cross_prod_norm < 1e-12) { /* Presumably, cross_prod_norm == 0 */
    
    *pminrad = DBL_MAX;  // The minrad is actually infinity.
    *mminrad = DBL_MAX;

  }
  
  *mminrad = normin*rad;
  *pminrad = normout*rad;
  
  mr = ((*mminrad < *pminrad) ? *mminrad : *pminrad);
  
  return mr;
  
}


void step( plCurve* inLink, double stepSize, plc_vector* dVdt ); // This is in stepper.c, but not exposed to public.


plCurve *anglecurve(double len1, double len2, double angle)
{
  plCurve *R;
  int nv = {3}, cc = {0};
  bool open = {true};

  R = plc_new(1,&nv,&open,&cc);

  R->cp[0].vt[0] = plc_build_vect(-len1,0,0);
  R->cp[0].vt[1] = plc_build_vect(0,0,0);
  R->cp[0].vt[2] = plc_build_vect(len2*cos(angle),len2*sin(angle),0);
 
  return R;
}

void test_mrgrad_var(int steps)

/* Perform mr gradient test at steps angles between 0 and pi. We don't include 0 or pi. */

{
  double PI = 3.14159265358979323846;
  double theta,tstep;
  int i;

  plCurve *L;
  plc_vector pGrad[3], mGrad[3];
  double mr_before,mr_after,pMinrad,mMinrad;
  double hdata[256],pquotdata[256],mquotdata[256];
  int ndata;
  bool newTheta = true;

  printf("Angle  Edgelengths              Minrad                       pquot fit                             mquot fit                               Result\n");
  printf("-------------------------------------------------------------------------------------------------------------------------------------------------\n");

  for(tstep = PI/(steps+2),i=1,theta=tstep;
      i<steps+1;i++,theta+=tstep,newTheta = true) {

    double len1,len2=2.0;

    for(len1=2-1e-6;len1<2+1e-6;len1+=1e-7) {

      L = anglecurve(len1,len2,theta);
      mr_before = pmMinrad(L,0,1,&pMinrad,&mMinrad);

      /* To compute the minrad gradients using Ridgerunner's code, we need to 
	 fake an octrope_mrloc. */

      octrope_mrloc thisloc;

      thisloc.component = 0;
      thisloc.vert = 1;
      thisloc.svert = 2;
      thisloc.mr = pMinrad;
      
      compute_minrad_gradient(L,thisloc,&pGrad[0],&pGrad[1],&pGrad[2]);

      thisloc.svert = 0;
      compute_minrad_gradient(L,thisloc,&mGrad[0],&mGrad[1],&mGrad[2]);

      // compute_pm_minrad_gradient(L,0,1,pGrad,mGrad);  This reference implementation is now used only for debugging. 

      if (newTheta) {

	printf("%6.2g ",theta);
	newTheta = false;
	
      } else {

	printf("       ");

      }

      printf("(%2.8f,%2.8f)  %2.6f (%2.6f,%2.6f) ",
	     len1,len2,mr_before,pMinrad,mMinrad);
      
      if (VERBOSE) {
	
	printf("\n");
	printf("\t pMinradGradient = (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f)\n",
	       plc_M_clist(pGrad[0]),plc_M_clist(pGrad[1]),plc_M_clist(pGrad[2]));
	printf("\t mMinradGradient = (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f)\n",
	       plc_M_clist(mGrad[0]),plc_M_clist(mGrad[1]),plc_M_clist(mGrad[2]));

      }
	
      plc_vector cmpP[3],cmpM[3];

      compute_pm_minrad_gradient(L,0,1,cmpP,cmpM);

      double pdist=0,mdist=0;
      pdist = plc_distance(cmpP[0],pGrad[0]) + plc_distance(cmpP[1],pGrad[1]) + plc_distance(cmpP[2],pGrad[2]);
      mdist = plc_distance(cmpM[0],mGrad[0]) + plc_distance(cmpM[1],mGrad[1]) + plc_distance(cmpM[2],mGrad[2]);
	    
      if (VERBOSE || pdist > 1e-10 || mdist > 1e-10) {
	
	printf("\t Comparison pMinradGradient  = (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f)\n",
	       plc_M_clist(cmpP[0]),plc_M_clist(cmpP[1]),plc_M_clist(cmpP[2]));
	printf("\t Comparison mMinradGradient  = (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f)\n",
	       plc_M_clist(cmpM[0]),plc_M_clist(cmpM[1]),plc_M_clist(cmpM[2]));
      }
      
      if (pdist > 1e-10 || mdist > 1e-10) {

	printf("Comparison Minrad gradients don't match the ridgerunner version. Test FAILS.\n");
	gFail = true;

      }
     
      

      double h,pquot,mquot;

      plc_vector stepDir[3];

      stepDir[0] = plc_random_vect();
      stepDir[1] = plc_random_vect();
      stepDir[2] = plc_random_vect();

      double stepNorm; 

      stepNorm = sqrt(plc_dot_prod(stepDir[0],stepDir[0]) + 
		      plc_dot_prod(stepDir[1],stepDir[1]) +
		      plc_dot_prod(stepDir[2],stepDir[2]));

      stepDir[0] = plc_scale_vect(1/stepNorm,stepDir[0]);
      stepDir[1] = plc_scale_vect(1/stepNorm,stepDir[1]);
      stepDir[2] = plc_scale_vect(1/stepNorm,stepDir[2]);

      if (VERBOSE)  
	{ 
	  
	  printf("\t Variation       = (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f)\n",
		 plc_M_clist(stepDir[0]),plc_M_clist(stepDir[1]),plc_M_clist(stepDir[2]));
	  
	}

      double mdder,pdder;

      mdder = plc_dot_prod(stepDir[0],mGrad[0]) + 
	plc_dot_prod(stepDir[1],mGrad[1]) + plc_dot_prod(stepDir[2],mGrad[2]);
      
      pdder = plc_dot_prod(stepDir[0],pGrad[0]) + 
	plc_dot_prod(stepDir[1],pGrad[1]) + plc_dot_prod(stepDir[2],pGrad[2]);
      
      for(h=1e-5,ndata=0;h>1e-9;h/=2.0,ndata++) {
	
	plCurve *workerLink;
	double pMinrad_after,mMinrad_after;

	workerLink = plc_copy(L);
	step(workerLink,h,stepDir);
	mr_after = pmMinrad(workerLink,0,1,&pMinrad_after,&mMinrad_after);
	
	pquot = fabs((h*pdder - (pMinrad_after - pMinrad))/h);
	mquot = fabs((h*mdder - (mMinrad_after - mMinrad))/h);

	if (VERBOSE) {
	  printf("\t\t h = %2.6f delta pMinrad = %1.10f pquot = %1.10f \t\t delta mMinrad = %1.10f mquot = %1.10f \n",h,
		 pMinrad_after - pMinrad,pquot,
		 mMinrad_after - mMinrad,mquot);
	}

	hdata[ndata] = h;
	pquotdata[ndata] = pquot;
	mquotdata[ndata] = mquot;
	
      }
           
      double c0,c1,cov00,cov01,cov11,sumsq;  /* We are going to fit to the linear model y = c0 + c1 h. */
      bool tFail = false;

      gsl_fit_linear(hdata,1,pquotdata,1,ndata,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
      printf("%10.4g h + %10.4g (%10.4g) ",c1,c0,sumsq);

      if (fabs(c0) > 1e-4 || fabs(sumsq) > 1e-7) {

	tFail = true;
	gFail = true;
       
      }

      gsl_fit_linear(hdata,1,mquotdata,1,ndata,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
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
    
  printf("mrgradtest, ridgerunner v%s\n\n",PACKAGE_VERSION);
  /* The PACKAGE_VERSION preprocessor symbol is defined in "config.h" by autoconf. */
  
  printf("Test the minrad gradient code in ridgerunner by taking small steps in a randomly chosen direction.\n"
	 "We compute the expected directional derivative and take the difference between the predicted and actual\n"
	 "values of plus and minus minrad, then divide by stepsize (h) to get the plus and minus quotients. We compute\n"
	 "these values for a range of h from 1e-2 to 1e-8. The resulting (h,pquot) and (h,mquot) data should \n"
	 "fit very well to lines through the origin (the slopes should be random). We display the best least-squares\n"
	 "fit below (along with the sum of the squares of the residuals). If the constant term and residuals are both\n"
	 "low, the algorithm will score the test as PASS.\n\n"
	 "Running this program (mrgradtest) from the command line with --verbose will give more details on tests.\n\n");
  
  /* We start by parsing the command-line arguments with argtable. */
  
  if (arg_nullcheck(argtable) != 0 || arg_nullcheck(helptable) != 0)
    printf("error: insufficient memory\n");
  
  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the first set of */
                        /* errors was probably more helpful, so we display it. */
      
      arg_print_errors(stdout,helpend,"helptable");
      arg_print_errors(stdout,end,"mrgradtest");
      exit(1);
      
    } else {  /* The help table matched, which means we asked for help or gave nothing */
      
      printf("mrgradtest performs self-tests on the minrad gradient code in ridgerunner.\n"
	     "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);
      
    }
    
  }
  
  if (verbose->count > 0) {

    VERBOSE = true;

  }
  
  test_mrgrad_var(100);

  if (gFail) { exit(1); } 
  else {exit(0);}
  
}
