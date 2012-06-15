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

#include"ridgerunner.h"
#include"mangle.h"
#include<argtable2.h>
#include<gsl/gsl_fit.h>

bool VERBOSE = false;
bool gFail = false;

search_state *gState;

void placeContactStruts ( taucs_ccs_matrix* A, plCurve* inLink, 
			  octrope_strut* strutSet, int strutCount );

void cvc(plCurve *inLink,int rownum,int *cmp,int *vert,int *coord);

/* These live in stepper.c and aren't public, so we include a prototype directly. */

/* int octrope_struts(plCurve *L,  */
/*                    const double strut_cutoff, const double epsilon,  */
/*                    octrope_strut *strutlist, int sl_size, */
/*                    double *shortest, void *mem, const int memsize); */

octrope_strut poca(plCurve *L)
/* Computes the poca strut (only) using a clever call to octrope_struts. The magic here*/
/* is that by setting epsilon to 0 we store only the shortest of all the struts. */
{
  octrope_strut poca = {{0,0},{0,0},{0,0},0,0};
  double pval;
  
  octrope_struts(L,0,0,&poca,1,&pval,NULL,0);

  return poca;
}

plc_vector *compute_poca_gradient(plCurve *inLink,octrope_strut strut)
/* We do a bit of work to create a fake rigidity matrix, place the strut, and then extract the data. */
/* Uses the global gState */
{
  taucs_ccs_matrix *cleanA;

  cleanA = taucs_ccs_new(3*plc_num_verts(inLink),1,12); /* (rows,columns,nnz) ccs matrix */
  cleanA->colptr[0] = 0;
  cleanA->colptr[1] = 12;

  placeContactStruts(cleanA,inLink,&strut,1);    /* This is the code we're trying to debug, so it's important to call it. */

  /* We now want to extract the data from cleanA into the correct components of first and second. */

  plc_vector *grad;
  grad = calloc(plc_num_verts(inLink),sizeof(plc_vector));

  int i,j,cptr;
  double val;
  int cp,vt,c;

  for(i=0;i<cleanA->n;i++) { /* Loop over columns (presumably only one) */

    for(cptr=cleanA->colptr[i];cptr<cleanA->colptr[i+1];cptr++) {

      j= cleanA->rowind[cptr];
      val= cleanA->values.d[cptr];
      grad[(int)(floor(j/3.0))].c[j % 3] = val;

      /* We now sanity check the j value. */
      
      cvc(inLink,j,&cp,&vt,&c); 

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


plCurve *simplecross(double len1, double pos1, double len2, double angle)

/* generates a test curve consisting of segments of variable lengths in the planes z = 0 and z = 1 */
{
  plCurve *R;
  int nv[2] = {2,2}, cc[2] = {0,0};
  bool open[2] = {true,true};

  R = plc_new(2,nv,open,cc);

  R->cp[0].vt[0] = plc_build_vect((pos1)*(-len1)/2.0,0,0);
  R->cp[0].vt[1] = plc_build_vect((1-pos1)*len1/2.0,0,0);
  R->cp[1].vt[0] = plc_build_vect(-(len2/2.0)*cos(angle),-(len2/2.0)*sin(angle),1);
  R->cp[1].vt[1] = plc_build_vect((len2/2.0)*cos(angle),(len2/2.0)*sin(angle),1);

  /* Now we fix gState. */

  if (gState->compOffsets != NULL) { free(gState->compOffsets); gState->compOffsets = NULL; }

  gState->compOffsets = calloc(2,sizeof(int));
  gState->compOffsets[0] = 0;
  gState->compOffsets[1] = 2;

  return R;
}

plCurve *angryguy(double len1, double ht1, double len2, double angle)

/* generates a test curve consisting of a v-shaped guy in the xy plane with vertex at 0,0,1 and a straight segment in the z = 0 plane */

{
  plCurve *R;
  int nv[2] = {3,2}, cc[2] = {0,0};
  bool open[2] = {true,true};

  R = plc_new(2,nv,open,cc);

  R->cp[0].vt[0] = plc_build_vect((-len1/2.0),0,0);
  R->cp[0].vt[1] = plc_build_vect(0,0,ht1);
  R->cp[0].vt[2] = plc_build_vect((len1/2.0),0,0);

  R->cp[1].vt[0] = plc_build_vect(-(len2/2.0)*cos(angle),-(len2/2.0)*sin(angle),1);
  R->cp[1].vt[1] = plc_build_vect((len2/2.0)*cos(angle),(len2/2.0)*sin(angle),1);

  /* Now we fix gState. */

  if (gState->compOffsets != NULL) { free(gState->compOffsets); gState->compOffsets = NULL; }

  gState->compOffsets = calloc(2,sizeof(int));
  gState->compOffsets[0] = 0;
  gState->compOffsets[1] = 3;

  return R;
}

plCurve *guy_from_10_154() 
{
  plCurve *R;
  int nv[2] = {3,2}, cc[2] = {0,0};
  bool open[2] = {true,true};

  R = plc_new(2,nv,open,cc);

  R->cp[0].vt[0] = plc_build_vect(2.0295126472416132, 0.43795076600880289, -1.1589291197160101);
  R->cp[0].vt[1] = plc_build_vect(2.0413493167425178, 0.3666563635753321, -1.136636493529819);
  R->cp[0].vt[2] = plc_build_vect(2.0533828160698868, 0.29984901861703778, -1.110372871991034);

  R->cp[1].vt[0] = plc_build_vect(2.9294924562856028, 0.40246302218472219, -1.596714018083581);
  R->cp[1].vt[1] = plc_build_vect(2.9378304833641922, 0.36187142029496788, -1.5797140655272379);

  /* Now we fix gState. */

  if (gState->compOffsets != NULL) { free(gState->compOffsets); gState->compOffsets = NULL; }

  gState->compOffsets = calloc(2,sizeof(int));
  gState->compOffsets[0] = 0;
  gState->compOffsets[1] = 3;

  return R;
}


void test_10_154_grad_var(int steps)

/* Perform poca gradient test for v-shaped curve with angles between 0 and pi. We don't include 0 or pi. */

{
  int i;

  plCurve *L;
  double poca_before,poca_after;
  double hdata[256],quotdata[256];
  int ndata;
  //bool newTheta = true;

  printf("Guy from 10_154                 quot fit                           Result\n");
  printf("-------------------------------------------------------------------------\n");


  L = guy_from_10_154();
  poca_before = octrope_poca(L,NULL,0);
  
  plc_vector *grad;
  grad = compute_poca_gradient(L,poca(L));
  
  printf("(10_154_dataset)  %2.6f  ",
	 poca_before);
  
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
    
  for(i=0;i<plc_num_verts(L);i++) { stepDir[i] = plc_random_vect(); stepNorm += plc_dot_prod(stepDir[i],stepDir[i]);}
  stepNorm = sqrt(stepNorm);
  for(i=0;i<plc_num_verts(L);i++) { stepDir[i] = plc_scale_vect(1/stepNorm,stepDir[i]); }
  
  //stepDir[1] = plc_build_vect(0,0,-1);
  
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
    
    quot = fabs((h*dder - ((poca_after/2.0) - (poca_before/2.0)))/h);
    
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

  /* Now we're going to do the same drill with dLen. */

  printf("\t (dlen force)       ");

  stepDir = calloc(plc_num_verts(L),sizeof(plc_vector));
  dlenForce(stepDir,L,NULL);
  stepNorm = sqrt(stepNorm);
  for(i=0;i<plc_num_verts(L);i++) { stepDir[i] = plc_scale_vect(1/stepNorm,stepDir[i]); }
  
  //stepDir[1] = plc_build_vect(0,0,-1);
  
  if (VERBOSE)  
    { 
      
      printf("\t Variation       = (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f),  (%+2.8f,%+2.8f,%+2.8f)\n",
	     plc_M_clist(stepDir[0]),plc_M_clist(stepDir[1]),plc_M_clist(stepDir[2]),plc_M_clist(stepDir[3]));
      
    }
  
  for(dder=0,i=0;i<plc_num_verts(L);i++) { dder += plc_dot_prod(stepDir[i],grad[i]); }
  
  /* We have now computed the directional derivative of poca in the step direction. It remains to actually vary the curve in this 
     direction and see to what extent this predicts the actual change in poca. */
  
  for(h=1e-5,ndata=0;h>1e-9;h/=2.0,ndata++) {
    
    plCurve *workerLink;
    
    workerLink = plc_copy(L);
    step(workerLink,h,stepDir);
    poca_after = octrope_poca(workerLink,NULL,0);
    
    quot = fabs((h*dder - ((poca_after/2.0) - (poca_before/2.0)))/h);
    
    if (VERBOSE) {
      printf("\t\t h = %2.6f delta poca = %1.10f quot = %1.10f \n",h,
	     poca_after - poca_before,quot);
    }
    
    hdata[ndata] = h;
    quotdata[ndata] = quot;
    
  }
  
  //double c0,c1,cov00,cov01,cov11,sumsq;  /* We are going to fit to the linear model y = c0 + c1 h. */
  //bool tFail = false;
  
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

  printf("\n");

}


void test_vgrad_var(int steps)

/* Perform poca gradient test for v-shaped curve with angles between 0 and pi. We don't include 0 or pi. */

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

    double ht1,len1=2.0,len2=2.0;

    for(ht1=1e-2;ht1>1e-8;ht1/=2) {

      L = angryguy(len1,ht1,len2,theta);
      poca_before = octrope_poca(L,NULL,0);

      plc_vector *grad;
      grad = compute_poca_gradient(L,poca(L));

      if (newTheta) {

	printf("%6.2g ",theta);
	newTheta = false;
	
      } else {

	printf("       ");

      }

      printf("(%2.8f,%2.8f,%2.8f)  %2.6f  ",
	     len1,ht1,len2,poca_before);
      
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
      stepNorm = sqrt(stepNorm);
      for(i=0;i<plc_num_verts(L);i++) { stepDir[i] = plc_scale_vect(1/stepNorm,stepDir[i]); }

      //stepDir[1] = plc_build_vect(0,0,-1);

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
	
	quot = fabs((h*dder - ((poca_after/2.0) - (poca_before/2.0)))/h);

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

    double pos1,len1=2.0,len2=2.0;

    for(pos1=0.1;pos1<=0.9;pos1+=0.1) {

      L = simplecross(len1,pos1,len2,theta);
      poca_before = octrope_poca(L,NULL,0);

      plc_vector *grad;
      grad = compute_poca_gradient(L,poca(L));

      if (newTheta) {

	printf("%6.2g ",theta);
	newTheta = false;
	
      } else {

	printf("       ");

      }

      printf("(%2.8f,%2.8f,%2.8f)  %2.6f  ",
	     len1,pos1,len2,poca_before);
      
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
      stepNorm = sqrt(stepNorm);
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
	
	quot = fabs((h*dder - ((poca_after/2.0) - (poca_before/2.0)))/h);

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

void test_link(plCurve *L,int ndirs)

/* Test variations of a particular file in a random direction.  */

{
  int i,j;
  double poca_before,poca_after;
  double hdata[256],quotdata[256];
  int ndata;
  
  //printf("Angle  Edgelengths              Poca                       quot fit                                 Result\n");
  //printf("----------------------------------------------------------------------------------------------------------\n");

  poca_before = octrope_poca(L,NULL,0);

  plc_vector *grad;
  octrope_strut Lpoca;
    
  Lpoca = poca(L);
  grad = compute_poca_gradient(L,Lpoca);

  printf("%+2.8f\n",poca_before);

  if (VERBOSE) {
    
    printf("\n");
    
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

  for(j=0;j<ndirs;j++) {
  
    plc_vector *stepDir;
    stepDir = calloc(plc_num_verts(L),sizeof(plc_vector));

    // Experiment. Let's only move the POCA vertices. 
    
    double stepNorm = 0;
    
    for(i=0;i<plc_num_verts(L);i++) { stepDir[i] = plc_random_vect(); stepNorm += plc_dot_prod(stepDir[i],stepDir[i]);}
    stepNorm = sqrt(stepNorm);
    for(i=0;i<plc_num_verts(L);i++) { stepDir[i] = plc_scale_vect(1/stepNorm,stepDir[i]); }
    
    if (VERBOSE)  
      { 
	
	printf("\t Variation       = (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f), (%+2.8f,%+2.8f,%+2.8f),  (%+2.8f,%+2.8f,%+2.8f)\n",
	       plc_M_clist(stepDir[gState->compOffsets[Lpoca.component[0]] + Lpoca.lead_vert[0]]) ,
	       plc_M_clist(stepDir[gState->compOffsets[Lpoca.component[0]] + (Lpoca.lead_vert[0]+1)]) ,
	       plc_M_clist(stepDir[gState->compOffsets[Lpoca.component[1]] + Lpoca.lead_vert[1]]) ,
	       plc_M_clist(stepDir[gState->compOffsets[Lpoca.component[1]] + (Lpoca.lead_vert[1]+1)]));     
	
      }
    
    double dder = 0;
    
    for(i=0;i<plc_num_verts(L);i++) { dder += plc_dot_prod(stepDir[i],grad[i]); }

    if (dder > 0) { /* We don't expect this to work when the derivative is positive, as another poca should take over */

      for(i=0;i<plc_num_verts(L);i++) { stepDir[i] = plc_scale_vect(-1.0,stepDir[i]); }
      dder = 0;
      for(i=0;i<plc_num_verts(L);i++) { dder += plc_dot_prod(stepDir[i],grad[i]); }

    }
    
    /* We have now computed the directional derivative of poca in the step direction. It remains to actually vary the curve in this 
       direction and see to what extent this predicts the actual change in poca. */
    
    for(h=1e-9,ndata=0;h>1e-11;h/=2.0,ndata++) {
      
      plCurve *workerLink;
      octrope_strut afterPoca;
      
      workerLink = plc_copy(L);
      step(workerLink,h,stepDir);
      poca_after = octrope_poca(workerLink,NULL,0);
      afterPoca = poca(workerLink);
      
      quot = fabs((h*dder - ((poca_after/2.0) - (poca_before/2.0)))/h);
      
      if (VERBOSE) {
	printf("\n\t\t\t h = %2.6f delta poca = %1.10f quot = %1.10f \n",h,
	       poca_after - poca_before,quot);
      }
      
      hdata[ndata] = h;
      quotdata[ndata] = quot;
      
    }
    
    double c0,c1,cov00,cov01,cov11,sumsq;  /* We are going to fit to the linear model y = c0 + c1 h. */
    bool tFail = false;
    
    gsl_fit_linear(hdata,1,quotdata,1,ndata,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
    printf("                                                 %10.4g h + %10.4g (%10.4g) ",c1,c0,sumsq);
    
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

  if (arg_infile->count == 0) {

    test_pocagrad_var(100);
    test_vgrad_var(100);
    test_10_154_grad_var(100);

  } else {

    printf("Filename                  Poca                       quot fit                                 Result\n");
    printf("----------------------------------------------------------------------------------------------------\n");

    int infilenum;
    
    for(infilenum=0;infilenum < arg_infile->count;infilenum++) {
      
      plCurve *link;
      FILE *infile_fptr;
      
      infile_fptr = fopen(arg_infile->filename[infilenum],"r");
      
      if (infile_fptr == NULL) {
	
	fprintf(stderr,"strutgradtest: Couldn't open file %s.\n",
		arg_infile->filename[infilenum]);
	continue;  /* Try the next file */
	
      }
      
      int plr_error_num;
      char plr_error_str[1024];
      
      link = plc_read(infile_fptr,
		      &plr_error_num,plr_error_str,sizeof(plr_error_str));
      
      /* We now demonstrate the octrope library's error handling protocol: */
      
      if (plr_error_num > 0) {   /* This is the signal for an error. */
	
	fprintf(stderr,"strutgradtest: link reading error\n%s\n",plr_error_str);
	continue;  /* Try the next file */
	
      }
      
      fclose(infile_fptr);

      printf("%-22s   ",arg_infile->basename[infilenum]);
      
      /* We have now loaded the file. We need to change gState. */
      
      if (gState->compOffsets != NULL) { free(gState->compOffsets); gState->compOffsets = NULL; }
      gState->compOffsets = malloc(sizeof(int)*(link->nc));
      fatalifnull_(gState->compOffsets);
      
      int offset = 0,i;
      
      for( i=0; i<link->nc; i++ )    {
	
	gState->compOffsets[i] = offset;
	offset += link->cp[i].nv;
	
      }
      
      /* Now we can go ahead and run the test. */
      
      test_link(link,10);
      
    }

  }
    
  if (gFail) { exit(1); } 
  else {exit(0);}
  
}
