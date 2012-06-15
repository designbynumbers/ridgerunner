/*******************************************************************

   lookahead.c : Test program for ridgerunner. Given a file, constructs
                 the ridgerunner step direction and samples the ropelength
                 of the curve at various stepsizes between the minimum and 0.01.

	   ************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<stdio.h>
#include<stdlib.h>

#include"ridgerunner.h"
#include"mangle.h"
#include<argtable2.h>

/* Global variables live here. */

void step( plCurve* inLink, double stepSize, plc_vector* dVdt ); // This is in stepper.c, but not exposed to public.

struct arg_dbl  *arg_minstep;
struct arg_dbl  *arg_scale;
struct arg_dbl  *arg_maxstep;
struct arg_int  *arg_nsamples;
struct arg_lit  *arg_logsampling;
struct arg_dbl  *arg_lambda;
struct arg_dbl  *target_thickness; 
struct arg_lit  *arg_dumpAxb;
struct arg_lit  *arg_timewarp;

struct arg_lit  *verbose;
struct arg_file *arg_infile;
struct arg_lit  *help;

struct arg_lit  *quiet;

struct arg_end *end;
struct arg_end *helpend;

FILE    *infile_fptr,*outfile_fptr;

int    QUIET=0;

struct dpoint {
  double x;
  double y;
};

double scale = 0.1;

int compare_data(const void *a,const void *b) {

  const struct dpoint *A = a;
  const struct dpoint *B = b;
  if (A->x > B->x) return 1;
  if (A->x < B->x) return -1;
  return 0;
}

/*plc_vector *inputForce( plCurve *inLink, double tube_radius, double eqMultiplier, double lambda,search_state *inState ); */

plc_vector *dLenDirection(plCurve *inLink,search_state *inState) {

  return inputForce(inLink,inState->tube_radius,inState->eqMultiplier,gLambda,inState);

}


int main(int argc,char *argv[])
{
  int       infilenum,i,nerrors;
  plCurve   *link;
  
  double    overstepTol=0.0001; /* Should this only be for tube_radius = 1? */
  double    tube_radius = 0.5;
  
  void *argtable[] = 
    {

     target_thickness  = arg_dbl0("r","radius","<x>","target tube radius"),
     arg_scale  = arg_dbl0("s","scale","<double>","maximum stepsize to plot"),
     arg_lambda = arg_dbl0("l","lambda","<double>",
			   "minimum radius of curvature for unit rope"),

     arg_timewarp = arg_lit0(NULL,"Timewarp","turn on timewarp to accelerate convergence of strut-free sections of tube"),
     arg_dumpAxb = arg_lit0(NULL,"dumpAxb","dump the rigidity matrix, right-hand side, and solution"),
     verbose = arg_lit0("v","verbose","print debugging information"),
     arg_infile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
    
     arg_minstep = arg_dbl0(NULL,"MinStep","<x>",
			         "minimum step size to check (can be negative)"),
     arg_maxstep = arg_dbl0(NULL,"MaxStep","<x>",
				 "maximum step size to check (should be > MinStep)"),

     arg_nsamples = arg_int0("s","Samples","<n>","plot ropelength at <n> sample locations"),

     /*arg_logsampling = arg_lit0("l","LogSampling","distribute samples logarithmically"),*/

     quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  fprintf(stderr,"lookahead (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");

  if (arg_nullcheck(argtable) != 0)
    printf("lookahead: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"lookahead");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("lookahead plots the ropelength of a VECT file along the standard search direction.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  if (arg_timewarp->count > 0) { gNoTimeWarp = false; } else { gNoTimeWarp = true; }

  QUIET = quiet->count > 0; /* Register if we're in batch mode */

  if (arg_scale->count > 0) { scale = arg_scale->dval[0]; }

  double correctionStepSize = 0.25;
  double minradOverstepTol = 0.00005; /* Used to be abs val of 0.499975 */

  /* We now need to fake an inState to pretend that we're in the middle
     of a ridgerunnner run. */

  search_state state;

  gLambda = 1.0;
  if (arg_lambda->count > 0) { gLambda = arg_lambda->dval[0]; }

  state.maxItrs = 999999;
  state.stop20  = -20;
  state.residualThreshold = -1;
  state.stopTime = 999999;

  if (arg_dumpAxb->count > 0) {
    state.dumpAxb = true;
  } else {
    state.dumpAxb = false;
  } 

  state.moviefactor = 0;    // Nobody knows whether these still work
  state.maxmovieframes = 0;

  state.saveConvergence = 0;
  state.movie = 0;
  state.fancyVisualization = 0;
  state.fancyPipe = 0;

  state.correctionStepDefault = correctionStepSize;
  state.movie = 0;
  state.moviefactor = 1.0;

  state.tube_radius = 0.5;

  if (target_thickness->count > 0) 
    { state.tube_radius = target_thickness->dval[0]; }
  
  state.overstepTol = overstepTol;
  state.minradOverstepTol = minradOverstepTol;
  state.minminrad = gLambda*state.tube_radius*(1 - minradOverstepTol);

  state.lastStepStrutCount = 0;
  state.lastStepStruts = NULL;
  
  state.lastStepMinradStrutCount = 0;
  state.lastStepMRlist = NULL;

  state.last_cstep_attempts = 0;
  state.cstep_count = 0;

  state.steps = 1;
  state.maxStepSize = 0.01;
  state.snapinterval = 10;

#ifdef HAVE_ASCTIME
#ifdef HAVE_LOCALTIME
#ifdef HAVE_TIME
  
  state.start_time = time(NULL);

#endif
#endif
#endif

  gLsqrLogging = false;

  printf("Filename           Data                                                                   \n");
  printf("------------------------------------------------------------------------------------------\n");

  for(infilenum = 0;infilenum < arg_infile->count;infilenum++) {

    printf("%-32s ",arg_infile->basename[infilenum]);

    /* We start by loading the tube from core. */

    FILE *infile_fptr;

    infile_fptr = fopen(arg_infile->filename[infilenum],"r");
  
    if (infile_fptr == NULL) {
      
      fprintf(stderr,"lookahead: Couldn't open file %s.\n",
	      arg_infile->filename[infilenum]);
       continue;  /* Try the next file */
       
    }

    int plr_error_num;
    char plr_error_str[1024];
    
    link = plc_read(infile_fptr,
		    &plr_error_num,plr_error_str,sizeof(plr_error_str));
    
    /* We now demonstrate the octrope library's error handling protocol: */
  
    if (plr_error_num > 0) {   /* This is the signal for an error. */
      
      fprintf(stderr,"lookahead: link reading error\n%s\n",plr_error_str);
      continue;  /* Try the next file */
      
    }
    
    fclose(infile_fptr);

    /* Now we set the state information as needed. */

    char *extptr;
    
    strncpy(state.basename,arg_infile->basename[0],sizeof(state.basename));
    extptr = strstr(state.basename,arg_infile->extension[0]);
    *extptr = 0; /* This chops the extension off of basename. */

    strncpy(state.fname,arg_infile->filename[0],sizeof(state.fname));
    sprintf(state.workingfilename,"./%s.lookahead/%s.working.vect",
	    state.basename,state.basename);
    sprintf(state.workingstrutname,"./%s.lookahead/%s.working.struts",
	    state.basename,state.basename);
    sprintf(state.vectprefix,"./%s.lookahead/vectfiles/%s",
	    state.basename,state.basename);
    sprintf(state.snapprefix,"./%s.lookahead/snapshots/%s",
	    state.basename,state.basename);
    sprintf(state.logfilename,"./%s.lookahead/%s.log",
	    state.basename,state.basename);
    sprintf(state.fprefix,"./%s.lookahead/",state.basename);

    /* We now open the local directory for storing output. */
    
    char cmdline[1024];
    
    sprintf(cmdline,"rm -fr %s.lookahead",state.basename);
    system_or_die(cmdline, __FILE__ , __LINE__ );
    
    sprintf(cmdline,"mkdir %s.lookahead",state.basename);
    system_or_die(cmdline, __FILE__ , __LINE__ );
    
    sprintf(cmdline,"mkdir %s.lookahead/logfiles",state.basename);
    system_or_die(cmdline, __FILE__ , __LINE__ );
    
    //sprintf(cmdline,"mkdir %s.cstest/vectfiles",state.basename);
    //system_or_die(cmdline, __FILE__ , __LINE__ );
    
    sprintf(cmdline,"mkdir %s.lookahead/snapshots",state.basename);
    system_or_die(cmdline, __FILE__ , __LINE__ );

    state.maxStepSize = 0.1*plCurve_short_edge(link);
    state.minrad = octrope_minradval(link);
    state.totalVerts = plc_num_verts(link);
    state.stepSize = 0.5*state.maxStepSize;
    
    if( ((double)1)/(double)state.totalVerts < state.stepSize )
      {
	state.stepSize = ((double)1/(double)state.totalVerts);
      }
    
    if( state.maxStepSize > 1e-3 ) { state.maxStepSize = 1e-3; }
    
    if( arg_maxstep->count > 0 ) { state.maxStepSize = arg_maxstep->dval[0]; }  
    state.stepSize = min(state.stepSize, state.maxStepSize);
    
    state.compOffsets = malloc(sizeof(int)*(link->nc));
    fatalifnull_(state.compOffsets);
    
    int offset = 0;
    
    for( i=0; i<link->nc; i++ )    {
      
      state.compOffsets[i] = offset;
      offset += link->cp[i].nv;
      
    }
    
    state.conserveLength = NULL; /* We aren't using this. */
    
    if (state.maxmovieframes % 2 == 1) { state.maxmovieframes--; } 
    /* Make sure this is even. */
    
    state.minrad = octrope_minradval(link);
    state.thickness = octrope_thickness(link,NULL,0,gLambda);
    state.ropelength = octrope_ropelength(link,NULL,0,gLambda);
    state.residual = 1.0;
    state.shortest = octrope_poca(link,NULL,0);
    state.eqMultiplier = 0;
    
    state.totalVerts = plc_num_verts(link);
    
    /*printf("Overstep tolerance: %f of thickness %f (%f) minminrad: %3.8lf\n", 
      state.overstepTol, state.tube_radius, 
      state.tube_radius*(1 - state.overstepTol), 
      state.minminrad ); */
    
    /* We now initialize the log with a lot of (hopefully) helpful information
       about the run. */
    
    gLogfile = fopen(state.logfilename,"w");
    if (gLogfile == NULL) { 
      fprintf(stderr,"lookahead: Couldn't open logfile %s.\n",state.logfilename);
      exit(1);
    }
    
    fprintf(gLogfile,"lookahead logfile.\n");
    
    fprintf(gLogfile,"lookahead %s\n",PACKAGE_VERSION);
    
    char svntag[1024];

    sprintf(svntag,"%s",SVNVERSION);
    if (!strstr("exported",svntag)) {  /* We were built from svn */
      fprintf(gLogfile,"svn version %s\n",SVNVERSION);
    }
    
    fprintf(gLogfile,"Built %s, %s.\n", __DATE__ , __TIME__ );
    
#ifdef HAVE_ASCTIME
#ifdef HAVE_LOCALTIME
#ifdef HAVE_TIME
    
    time_t start_time;
    state.start_time = start_time = time(NULL);
    
    fprintf(gLogfile,"Run began: %s",asctime(localtime(&start_time)));
    
#else
    
    fprintf(gLogfile,
	    "Warning: system doesn't have 'time' function.\n"
	    "         data logs won't include time information.\n");
    
#endif
#else
    
    fprintf(gLogfile,
	    "Warning: system doesn't have 'localtime' function.\n"
	    "         data logs won't include time information.\n");
    
#endif
#else
  
    fprintf(gLogfile,
	    "Warning: system doesn't have 'asctime' function.\n"
	    "         data logs won't include time information");
#endif
    
    fprintf(gLogfile,"Initial filename: %s.\n",
	    state.fname);
    
    fprintf(gLogfile,"Command line: ");
    
    for(i=0;i<argc;i++) {
      
      fprintf(gLogfile,"%s ",argv[i]);
      
    } 
    
    fprintf(gLogfile,"\n");
    
#ifdef HAVE_GETPID
    
    fprintf(gLogfile,"Process ID: %d.\n",(int)(getpid()));
    
#endif
    
    /* We are now prepared to get the step direction. */

    plc_vector *stepDir,*dLen;

    stepDir = stepDirection(link,tube_radius,0,gLambda,&state); // tube_radius, no Eq force, lambda = gLambda
    //stepDir = calloc(plc_num_verts(link),sizeof(plc_vector));
    //stepDir[0] = plc_build_vect(1,1,1);

    dLen = dLenDirection(link,&state);

    /* Now it would be completely unfair to compare these vectorfields, because they have different norms. */
    
    int k;
    double scaleFactor;
    scaleFactor = plc_l2norm(stepDir,state.totalVerts)/plc_l2norm(dLen,state.totalVerts);
    for(k=0;k<state.totalVerts;k++) { dLen[k] = plc_scale_vect(scaleFactor,dLen[k]); }
    scaleFactor = plc_l2norm(stepDir,state.totalVerts)/plc_l2norm(dLen,state.totalVerts);

    /* Now they have the same scale */

    if (stepDir == NULL || dLen == NULL) {

      char errMsg[1024],dumpname[1024];

      dumpLink(link,&state,dumpname);
      sprintf(errMsg,
	      "Could not compute a step direction at the start of a steepest descent step.\n"
	      "Either this link cannot be run with -c or correction stepping has resulted\n"
	      "in a degenerate configuration.\n"
	      "Dumped link to %s.\n"
	      "Terminating run.\n",dumpname);
      FatalError(errMsg, __FILE__, __LINE__ );

    }

    /* If we've made it here, we're ready to get some data. */

    double *sample_x;
    int samps = 100;
    
    if (arg_nsamples->count > 0) {
      
      samps = arg_nsamples->ival[0];

    }

    sample_x = calloc(samps,sizeof(double));
    assert(sample_x != NULL);

    double a,b;

    a=-scale; 
    b=scale;

    if (arg_minstep->count > 0) { 

      a = arg_minstep->dval[0];

    } 

    if (arg_maxstep->count > 0) {

      b = arg_maxstep->dval[0];

    }

    int i;
    double x;   

    for(i=0,x=a;i<samps;i++,x+=(b-a)/(samps-1)) {

      sample_x[i] = x;

    }

    /* We now gather the actual data. */

    double *scores,*dlscores, *prescores, *predlscores;

    scores = calloc(samps,sizeof(double));
    assert(scores != NULL);

    dlscores = calloc(samps,sizeof(double));
    assert(dlscores != NULL);

    prescores = calloc(samps,sizeof(double));
    assert(prescores != NULL);

    predlscores = calloc(samps,sizeof(double));
    assert(predlscores != NULL);

   /*  /\* DEBUGGING CODE for a particular failure. *\/ */

/*     scores[1] = predict_deltarop(link,dLen,1e-5,tube_radius,gLambda); */

/*     plCurve *wLink; */
/*     wLink = plc_copy(link); */

/*     step(wLink,1e-5,dLen); */

/*     double newrop,newthi,newmr,newpoca,newlen; */
/*     octrope_strut struts[3000]; */
/*     octrope_mrloc mrlocs[3000]; */
/*     int nmr,nstrut; */

/*     octrope(wLink,&newrop,&newthi,&newlen,&newmr,&newpoca, */
/* 	    gLambda*state.tube_radius,0,mrlocs,3000,&nmr, */
/* 	    2*state.tube_radius,0,struts,3000,&nstrut, */
/* 	    NULL,0,gLambda); */
    
/*     scores[0] = stepScore(link,&state,dLen,1e-5) - state.ropelength; */

/*     printf("predicted %g, got %g.\n",scores[1],scores[0]); */

/*     exit(1); */



    for(i=0;i<samps;i++) {

      scores[i] = stepScore(link,&state,stepDir,sample_x[i]) - stepScore(link,&state,stepDir,0);
      dlscores[i] = stepScore(link,&state,dLen,sample_x[i]) - stepScore(link,&state,dLen,0);
      prescores[i] = predict_deltarop(link,stepDir,sample_x[i],tube_radius,gLambda);
      predlscores[i] = predict_deltarop(link,dLen,sample_x[i],tube_radius,gLambda);

    }

    /* We now output a summary for the user. */

    for(i=0;i<8;i++) {

      printf("%13.8f ",sample_x[(int)(floor(i*samps/8))]);

    } 

    printf("\n     scores                      ");

    
    for(i=0;i<8;i++) {

      printf("%13.8f ",scores[(int)(floor(i*samps/8))]);

    } 

    printf("\n     predicted                   ");

    
    for(i=0;i<8;i++) {

      printf("%13.8f ",prescores[(int)(floor(i*samps/8))]);

    } 
    

    printf("\n                 ");


    /* We now output the data at the largest scale. */
    /* We now collect data. */

    FILE *datfile;
    char datname[1024];
    double stepMax[4] = {-100000,-100000,-100000,-100000}, stepMin[4] = {100000000,1000000,10000000,1000000}; 
    /* We will set the plot range to the max of stepDir */
    
    struct dpoint *data, *dldata, *predata, *predldata;

    data = calloc(samps*5,sizeof(struct dpoint));
    dldata = calloc(samps*5,sizeof(struct dpoint));
    predata = calloc(samps*5,sizeof(struct dpoint));        /* Predicted stepDir data */
    predldata = calloc(samps*5,sizeof(struct dpoint));      /* Predicted dLen data */
    
    int dItr = 0;

    for(i=0;i<samps;i++,dItr++) {

      data[dItr].x = sample_x[i];
      data[dItr].y = scores[i];

      dldata[dItr].x = sample_x[i];
      dldata[dItr].y = dlscores[i];

      predata[dItr].x = sample_x[i];
      predata[dItr].y = prescores[i];
       
      predldata[dItr].x = sample_x[i];
      predldata[dItr].y = predlscores[i];

      stepMax[0] = (data[dItr].y > stepMax[0]) ? data[dItr].y : stepMax[0];
      stepMin[0] = (data[dItr].y < stepMin[0]) ? data[dItr].y : stepMin[0];

    }
    
    /* We now repeat at other scales; */

    for(k=0;k<3;k++) {

      a /= 10.0;
      b /= 10.0;
      
      for(i=0,x=a;i<samps;i++,x+=(b-a)/(samps-1),dItr++) {

	data[dItr].x = x;
	data[dItr].y = stepScore(link,&state,stepDir,x) - stepScore(link,&state,stepDir,0);

	dldata[dItr].x = x;
	dldata[dItr].y = stepScore(link,&state,dLen,x) - stepScore(link,&state,dLen,0);

	predata[dItr].x = x;
	predata[dItr].y = predict_deltarop(link,stepDir,x,tube_radius,gLambda);

	predldata[dItr].x = x;
	predldata[dItr].y = predict_deltarop(link,dLen,x,tube_radius,gLambda);

	stepMax[k+1] = (data[dItr].y > stepMax[k+1]) ? data[dItr].y : stepMax[k+1];
	stepMin[k+1] = (data[dItr].y < stepMin[k+1]) ? data[dItr].y : stepMin[k+1];
	
      }

    }

    a *= 10*10*10;
    b *= 10*10*10;

    dldata[dItr].x = 0;
    dldata[dItr].y = 0;

    data[dItr].x = 0;
    data[dItr].y = 0;

    predldata[dItr].x = 0;
    predldata[dItr].y = 0;

    predata[dItr].x = 0;
    predata[dItr].y = 0;

    /* We now sort the data by x. */

    qsort(data,4*samps+1,sizeof(struct dpoint),compare_data);
    qsort(dldata,4*samps+1,sizeof(struct dpoint),compare_data);
    qsort(predata,4*samps+1,sizeof(struct dpoint),compare_data);
    qsort(predldata,4*samps+1,sizeof(struct dpoint),compare_data);

    /* And now we write it to the file. */

    sprintf(datname,"./%s.lookahead/%s.la.dat",
	    state.basename,state.basename);

    datfile = fopen(datname,"w");
    assert(datfile != NULL);

    for(k=0;k<4*samps;k++) {

      fprintf(datfile,"%15.15f %15.15f \n",data[k].x,data[k].y);

    }

    fclose(datfile);
    free(data);

    sprintf(datname,"./%s.lookahead/%s_d.la.dat",
	    state.basename,state.basename);

    datfile = fopen(datname,"w");
    assert(datfile != NULL);

    for(k=0;k<4*samps;k++) {

      fprintf(datfile,"%15.15f %15.15f \n",dldata[k].x,dldata[k].y);

    }

    fclose(datfile);
    free(dldata);

    sprintf(datname,"./%s.lookahead/%s_p.la.dat",
	    state.basename,state.basename);
    
    datfile = fopen(datname,"w");
    assert(datfile != NULL);

    for(k=0;k<4*samps;k++) {

      fprintf(datfile,"%15.15f %15.15f \n",predata[k].x,predata[k].y);

    }

    fclose(datfile);
    free(predata);

    sprintf(datname,"./%s.lookahead/%s_pd.la.dat",
	    state.basename,state.basename);

    datfile = fopen(datname,"w");
    assert(datfile != NULL);

    for(k=0;k<4*samps;k++) {

      fprintf(datfile,"%15.15f %15.15f \n",predldata[k].x,predldata[k].y);

    }

    fclose(datfile);
    free(predldata);

    double margin;

    for(k=0;k<4;k++) { 

      margin = 0.1 * (stepMax[k] - stepMin[k]);
      stepMax[k] += margin;
      stepMin[k] -= margin;

    }
   
    /* We now try to generate and xv a gnuplot graph of the results. */

    FILE *gpfile;
    char gpfilename[1024];
  
    sprintf(gpfilename,"./%s.lookahead/genplot.gnuplot",
	    state.basename);

    gpfile = fopen(gpfilename,"w");
    assert(gpfile != NULL);
    
    fprintf(gpfile,

	    "# \n"
	    "set terminal pdf\n"
	    "set output \"%s.la.pdf\"\n"
	    "set autoscale y2\n"
	    
	    "set style data lines \n"
	    "#set size 1.0, 1.0 \n"
	    "#set origin 0.0, 0.0 \n"
	    "set multiplot layout 2, 2 title \"Delta Ropelength score vs stepsize for %s.vect \" \n"

	    "#set size 0.5, 0.5 \n"
	    "#set origin 0.0, 0.5 \n"
	    "set xrange [%g:%g] \n"
	    "set yrange [%13.21g:%13.21g] \n"
	    "plot '%s.la.dat' title 'stepDir', '%s_d.la.dat' title 'dLen', '%s_p.la.dat' title 'p(stepDir)', '%s_pd.la.dat' title 'p(dLen)' \n"

	    "#set size 0.5, 0.5 \n"
	    "#set origin 0.5, 0.5 \n"
	    "set xrange [%g:%g] \n"
	    "set yrange [%13.21g:%13.21g] \n"
	    "plot '%s.la.dat' title 'stepDir', '%s_d.la.dat' title 'dLen',  '%s_p.la.dat' title 'p(stepDir)', '%s_pd.la.dat' title 'p(dLen)' \n"

	    "#set size 0.5, 0.5 \n"
	    "#set origin 0.0, 0.0 \n"
	    "set xrange [%g:%g] \n"
	    "set yrange [%13.21g:%13.21g] \n"
	    "plot '%s.la.dat' title 'stepDir', '%s_d.la.dat' title 'dLen',  '%s_p.la.dat' title 'p(stepDir)', '%s_pd.la.dat' title 'p(dLen)'  \n"

	    "#set size 0.5, 0.5 \n"
	    "#set origin 0.5, 0.0 \n"
	    "set xrange [%g:%g] \n"
	    "set yrange [%13.21g:%13.21g] \n"
	    "plot '%s.la.dat' title 'stepDir', '%s_d.la.dat' title 'dLen',  '%s_p.la.dat' title 'p(stepDir)', '%s_pd.la.dat' title 'p(dLen)'  \n"
	    "unset multiplot \n"
	    
	    "#clear \n"
	    "# \n"
	    "#set terminal pdf\n"
	    "#set output \"%s.la.pdf\"\n"
	    "set autoscale y\n"
	    "set autoscale y2\n"
	    
	    "set style data lines \n"
	    "#set size 1.0, 1.0 \n"
	    "#set origin 0.0, 0.0 \n"
	    "set multiplot layout 2, 2 title \"Delta Ropelength score vs stepsize for %s.vect \" \n"

	    "#set size 0.5, 0.5 \n"
	    "#set origin 0.0, 0.5 \n"
	    "set xrange [%g:%g] \n"
	    "#set yrange [%13.21g:%13.21g] \n"
	    "plot '%s.la.dat' title 'stepDir', '%s_d.la.dat' title 'dLen',  '%s_p.la.dat' title 'p(stepDir)', '%s_pd.la.dat' title 'p(dLen)' \n"

	    "#set size 0.5, 0.5 \n"
	    "#set origin 0.5, 0.5 \n"
	    "set xrange [%g:%g] \n"
	    "#set yrange [%13.21g:%13.21g] \n"
	    "plot '%s.la.dat' title 'stepDir', '%s_d.la.dat' title 'dLen',  '%s_p.la.dat' title 'p(stepDir)', '%s_pd.la.dat' title 'p(dLen)' \n"

	    "#set size 0.5, 0.5 \n"
	    "#set origin 0.0, 0.0 \n"
	    "set xrange [%g:%g] \n"
	    "#set yrange [%13.21g:%13.21g] \n"
	    "plot '%s.la.dat' title 'stepDir', '%s_d.la.dat' title 'dLen',  '%s_p.la.dat' title 'p(stepDir)', '%s_pd.la.dat' title 'p(dLen)'  \n"

	    "#set size 0.5, 0.5 \n"
	    "#set origin 0.5, 0.0 \n"
	    "set xrange [%g:%g] \n"
	    "#set yrange [%13.21g:%13.21g] \n"
	    "plot '%s.la.dat' title 'stepDir', '%s_d.la.dat' title 'dLen',  '%s_p.la.dat' title 'p(stepDir)', '%s_pd.la.dat' title 'p(dLen)'  \n"
	    "unset multiplot \n"
	    ,
	    state.basename,state.basename,
	    a,b,stepMin[0],stepMax[0],state.basename,state.basename,state.basename,state.basename,
	    a/10.0,b/10.0,stepMin[1],stepMax[1],state.basename,state.basename,state.basename,state.basename,
	    a/100.0,b/100.0,stepMin[2],stepMax[2],state.basename,state.basename,state.basename,state.basename,
	    a/1000.0,b/1000.0,stepMin[3],stepMax[3],state.basename,state.basename,state.basename,state.basename,

	    state.basename,state.basename,
	    a,b,stepMin[0],stepMax[0],state.basename,state.basename,state.basename,state.basename,
	    a/10.0,b/10.0,stepMin[1],stepMax[1],state.basename,state.basename,state.basename,state.basename,
	    a/100.0,b/100.0,stepMin[2],stepMax[2],state.basename,state.basename,state.basename,state.basename,
	    a/1000.0,b/1000.0,stepMin[3],stepMax[3],state.basename,state.basename,state.basename,state.basename
	    );

    fclose(gpfile);

    char cmd[1024];
    
    sprintf(cmd,"bash -c \"cd %s.lookahead; gnuplot genplot.gnuplot; gv --scale=4 --spartan *.pdf &\"",state.basename);
    system(cmd);
    
    /* We now clean up after ourselves */

    free(stepDir);
    free(scores);
    free(sample_x);

    plc_free(link);

  }
  
 if (state.lastStepStruts != NULL) { free(state.lastStepStruts); state.lastStepStruts = NULL;}
 if (state.lastStepMRlist != NULL) { free(state.lastStepMRlist); state.lastStepMRlist = NULL;}

 printf("\n");

 exit(0);

}
    


