/*
 *  ridgerunner_main.c
 *  ridgerunner
 *
 *  Created by Michael Piatek on Fri May 28 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 */
 
#include"ridgerunner.h"

void usage();
void reload_handler( int sig );
void initializeState( search_state* state, plCurve** inLink, const char* fname );

/* globals moved to ridgerunner.h */

#ifndef min
  #define min(a,b) ((a)<(b)) ? (a) : (b)
#endif

/*	
	struct options longopts[] = { {"avoidTmpConflicts",
		no_argument, &gAvoidTmpConflicts, 1},
		{"saveConvergenceData", no_argument, &saveConvergence,
		1}, {"quiet", no_argument, &gQuiet, 1},
		{"buildSurface", no_argument, &gSurfaceBuilding, 1},
		{"verbose", no_argument, &gVerboseFiling, 1},
		{"suppressFiles", no_argument, &gSuppressOutput, 1},
		{"movie", no_argument, &movie, 1}, {"ignoreCurvature",
		no_argument, &ignorecurvature, 1}, {"autoscale",
		no_argument, &autoscale, 1}, {"fancyViz", no_argument,
		&fancyViz, 1}, {"infile", required_argument, NULL,
		'a'}, {"scale", required_argument, NULL, 'b'},
		{"checkDelta", required_argument, NULL, 'c'},
		{"injrad", required_argument, NULL, 'd'},
		{"eqMultiplier", required_argument, NULL, 'e'},
		{"residualThreshold", required_argument, NULL, 'f'},
		{"refineUntil", required_argument, NULL, 'g'},
		{"checkThreshold", required_argument, NULL, 'h'},
		{"maxSteps", required_argument, NULL, 'i'},
		{"overstepTolerance", required_argument, NULL, 'j'},
		{"maxStepSize", required_argument, NULL, 'k'},
		{"double", no_argument, NULL, 'l'}, {0, 0, 0, 0} };
*/
	
int
main( int argc, char* argv[] )
{

  struct arg_file *arg_infile = arg_file1(NULL,NULL,"<VECT file>","input file");

  struct arg_dbl  *arg_lambda = arg_dbl0("l","lambda","<double>",
					 "minimum radius of curvature for unit rope");

  struct arg_rem  *arg_curveopts = arg_rem("","Operations on input curve before run");

  struct arg_lit  *arg_autoscale = arg_lit0("a","autoscale","automatically scale curve "
					    "to thickness 2.001");

  struct arg_dbl  *arg_resolution = arg_dbl0("r","res","<# verts/unit ropelength>",
					     "spline the curve to this resolution");

  struct arg_lit  *arg_eqit = arg_lit0("e","eq","equilateralize curve"); 

  struct arg_rem  *arg_stopopts = arg_rem("","Stopping Criteria (can use more than one)");
  struct arg_dbl  *arg_stop20 = arg_dbl0(NULL,"Stop20","<change in Rop>","stop when "
					 "average change in ropelength over last 20 "
					 "steps is less than this value");

  struct arg_dbl  *arg_stopRes = arg_dbl0(NULL,"StopResidual","<fraction>",
					  "stop when less than this fraction of "
					  "curvature force remains after projection");

  struct arg_int  *arg_stopSteps = arg_int0("s","StopSteps","<#steps>","stop after "
					    "this number of steps");

  struct arg_rem  *arg_fileopts = arg_rem("","File Output Options");

  struct arg_lit  *arg_suppressfiles = arg_lit0(NULL,"NoOutputFiles",
		      "don't save intermediate files during run");
  struct arg_rex  *arg_outpath = arg_rex0(NULL,"OutPath","/*/",
					  "</home/../outdir/>",
					  0,"path for output files");

  struct arg_rem  *arg_progopts = arg_rem("","Program and Algorithm Options");


  struct arg_lit  *arg_quiet  = arg_lit0("q","quiet","quiet");
  struct arg_lit  *arg_verbose = arg_lit0("v","verbose","verbose");
  struct arg_lit  *arg_help = arg_lit0("h","help","print this help and exit");

  struct arg_dbl  *arg_cstep_size = arg_dbl0("k","CorrectionStepSize","<fraction>",
					     "initial size of Newton correction step");
  
  struct arg_dbl  *arg_eqmult = arg_dbl0(NULL,"EqMultiplier","<scalar>","increase to "
					 "make 'equilateralization force' stronger");
  
  struct arg_dbl  *arg_overstep = arg_dbl0("o","OverstepTol","<x>",
					   "start correction step if "
					   "thickness < 2.0 - x");

  struct arg_dbl  *arg_mroverstep = arg_dbl0(NULL,"MinRadOverstepTol","<x>",
					     "start correction step if "
					     "minrad < lambda - x");

  struct arg_dbl  *arg_maxstep = arg_dbl0(NULL,"MaxStep","<scalar>","maximum stepsize "
					  "for a length-reduction step");

  struct arg_end *end = arg_end(20);
  
  void *argtable[] = {arg_infile,arg_lambda,
		      arg_curveopts,arg_autoscale,arg_resolution,arg_eqit,
		      arg_stopopts,arg_stop20,arg_stopRes,arg_stopSteps,
		      arg_fileopts,arg_suppressfiles,arg_outpath,
		      arg_progopts,arg_quiet,arg_verbose,arg_help,
		        arg_cstep_size,arg_eqmult,arg_overstep,arg_mroverstep,
		      arg_maxstep,
		      end};
  int nerrors;


  plCurve*	link = NULL;
  FILE*		linkFile = NULL;
  int		cItr=0;
  search_state	state;
  char		opt;
  char		cmd[1024];
  short		autoscale = 0, fixlengths = 0, movie = 0;
  double	refineUntil = -1;  /* Used to be an int */
  int           ignorecurvature=0;
  
  double	injrad = 0.5, overstepTol=0.0001; /* Should this be injrad = 1? */
  double	residualThreshold = 65000;
  
  /*double	checkThreshold = 0.05;
    double	checkDelta = 1; */

  double        stop20 = 0.0;

  double        eqMult = 2.0;
  int		doubleCount = 0;

  char		fname[1024];
  short		fancyViz = 0, saveConvergence = 0;
  double	scaleAmt = 1;
  double	maxStep = -1;
  long		maxItrs = -1;
  int		graphIt[kTotalGraphTypes];
  int		gItr;
  double	correctionStepSize = 0.25;
  double	minradOverstepTol = 0.499975;
  
  for(i=0;i<kTotalGraphTypes;i++) {graphIt[i] = 0;}
  srand(time(NULL));

  /* Display opening message. */

  printf("Ridgerunner %s (cvs build %s %s)\n",PACKAGE_VERSION, __DATE__ , __TIME__ );
  plc_version(NULL,0);
  octrope_version(NULL,0);
  tsnnls_version(NULL,0); 
  
  /* Parse command-line arguments */

  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries detected, allocations must have failed. */
    fprintf(stderr,"%s: insufficient memory\n",argv[0]);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  nerrors = arg_parse(argc,argv,argtable);

  /* special case: '--help' takes preceence over error reporting */
  if (arg_help->count > 0) {
    fprintf(stderr,"Usage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 0;
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    fprintf(stderr,"\n");
    arg_print_errors(stderr,end,argv[0]);
    fprintf(stderr,"\nUsage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }
   
  if (arg_cstep_size->count > 0) { 

    correctionStepSize = arg_cstep_size->dval[0];

  }

  /* Note: there used to be a way to turn on "gFastCorrectionSteps" from the cmd line */
  
  if (arg_quiet->count > 0) { gQuiet = 1; }

  if (arg_maxstep->count > 0) { maxStep = arg_maxstep->dval[0]; }

  /* Note: there used to be a way to turn on "gSurfaceBuilding" and "gVerboseFiling" */

  if (arg_mroverstep->count > 0) {minradOverstepTol = lambda - arg_mroverstep->dval[0];}

  if (arg_suppressfiles->count > 0) { gSuppressOutput = 1; }

  if (arg_overstep->count > 0) { overstepTol = arg_overstep->dval[0]; }

  if (arg_stopSteps->count > 0) { maxItrs = arg_stopSteps->ival[0]; }

  if (arg_stop20->count > 0) { stop20 = arg_stop20->dval[0]; }

  if (arg_resolution->count > 0) { refineUntil = arg_resolution->dval[0]; }
 
  if (arg_stopRes->count > 0) { residualThreshold = arg_stopRes->dval[0]; }

  if (arg_eqmult->count > 0) { eqMult = arg_eqmult->dval[0]; }

  /* Note: There used to be a way to change "injrad", "scaleamt", and "fixlengths" 
     from the cmdline. */
  
  if (arg_autoscale->count > 0) { autoscale = 1; }
  
  /* Note: There used to be a way to set "movie", "gPaperInfoinTmp", "ignorecurvature",
     and "fancyviz" from the cmdline. */

  /* We now open the given link file. */

  linkFile = fopen(arg_infile->filename[0],"r");

  if( linkFile == NULL ) {
    
    fprintf(stderr, "ridgerunner: Couldn't open file %s\n", arg_infile->filename);
    exit(1);
  }

  int plc_read_error_num;
  char plc_error_string[1024];
  size_t plc_error_size = {1024};

  link = plc_read(linkFile,&plc_read_error_num,plc_error_string,plc_error_size); 

  if (plc_read_error_num != 0) {

    fprintf(stderr,"ridgerunner: file %s had \n"
	    "             plc_read error %s.\n",plc_error_string);
    exit(1);
    
  }

  strncpy(fname,arg_infile->basename,sizeof(fname));
  fclose(linkFile);

  printf("Loaded %d component, %d vertex plCurve from %s.\n",
	 link->nc,plc_num_verts(link),fname);

  /* We now open the local directory for storing output. */

  char cmdline[1024];

  sprintf(cmdline,"rm -fr %s.rr",arg_infile->basename);
  system(cmdline);
  sprintf(cmdline,"mkdir %s.rr",arg_infile->basename);

  if (system(cmdline) > 0) {

    printf("Created directory %s.rr to hold output of this run.\n",
	   arg_infile->basename);

  } else {

    printf("WARNING. Creation of output directory %s.rr may have failed.\n",
	   arg_infile->basename);

  }

  sprintf(cmdline,"mkdir %s.rr/logfiles",arg_infile->basename);

  if (system(cmdline) <= 0) {

    printf("WARNING. Creation of logfile directory %s.rr/logfiles may have failed.\n");

  }

  if (arg_suppressfiles->count == 0) {

    sprintf(cmdline,"mkdir %s.rr/vectfiles",arg_infile->basename);

  }

  /* We archive a version of the initial input file in the run directory. */

  char filename[1024];
  FILE *savefile;

  sprintf(filename,"./%s.rr/%s.vect",arg_infile->basename,arg_infile->basename);
  
  savefile = fopen(filename,"w");
  
  if (savefile == NULL) {

    fprintf(stderr,"ridgerunner: Couldn't open archive copy %s for writing.\n",
	    filename);
    exit(1);

  }

  plc_write(savefile,link);
  fclose(savefile);

  printf("Saved copy of %s to %s.\n",arg_infile->filename,filename);

  /* Now perform initial operations. */

  if (arg_resolution->count > 0) {

    tempLink = plCurve_fixresolution(link,arg_resolution->dval[0]);
    plc_free(link);
    link = tempLink;

    printf("Splined to resolution of %g verts/rop. New curve has %d verts.\n",
	   arg_resolution->dval[0],plc_num_verts(link));

  }

  if( arg_eqit->count > 0 ) {

    for(i=0;i<3;i++) {

      tempLink = link;
      link = octrope_fixlength(tempLink);
      plc_free(tempLink);

    }

    printf("Equilateralized curve has max edgelength/min edgelength %g.\n",
	   plCurve_long_edge(link)/plCurve_short_edge(link));

  }

  if( arg_autoscale->count > 0 ) {

    double thickness;

    thickness = octrope_thickness(link,NULL,0,gLambda);
    printf("Curve has thickness %g. Autoscaling to thickness 2.001.\n",
	   thickness);
    plc_scale(link,2.001/thickness);
    thickness = octrope_thickness(link,NULL,0,gLambda);
    printf("Scaled curve has thickness %g.",thickness);

    if (fabs(thickness - 2.001) > 1e-12) {

      printf("\n");
      fprintf(stderr,"ridgerunner: Failed to scale %s to thickness 2.001."
	      "             Aborting run.\n",fname);
      exit(1);

    } else {

      printf(" Autoscale selftest ok.\n");

    }

  }

  /* Now save the modified file so that we record what we actually ran with. */

  sprintf(filename,"./%s.rr/%s.atstart.vect",arg_infile->basename,arg_infile->basename);
  savefile = fopen(filename,"w");
  
  if (savefile == NULL) {

    fprintf(stderr,"ridgerunner: Couldn't open atstart save version %s for writing.\n",
	    filename);
    exit(1);

  }

  plc_write(savefile,link);
  fclose(savefile);

  printf("Rerez'd, autoscaled, eq'd file written to %s.\n",filename);
    
  /* Now get the rest of the program ready to go. */
  
  initializeState( &state, &link, fname);
  updateSideLengths(link, &state);

  state.stop20 = stop20;
  
  state.refineUntil = refineUntil;
  state.checkDelta = checkDelta;
  state.checkThreshold = checkThreshold;
  state.movie = movie;
  
  state.injrad = injrad;
  state.overstepTol = overstepTol;
  state.minradOverstepTol = minradOverstepTol;
  state.minminrad = minradOverstepTol;
  
  state.eqMultiplier = eqMult;
  state.residualThreshold = residualThreshold;

  sprintf(state.finalfilename,"./%s.rr/%s.final.vect",
	  arg_infile->basename,arg_infile->basename);
  sprintf(state.finalstrutname,"./%s.rr/%s.final.struts",
	  arg_infile->basename,arg_infile->basename);
  sprintf(state.logfilename,"./%s.rr/%s.rrlog.dat",
	  arg_infile->basename,arg_infile->basename);
  sprintf(state.vectprefix,"./%s.rr/vectfiles/%s",
	  arg_infile->basename,arg_infile->basename); 

    /* Piatek's realtime graphing fanciness has been replaced by
     comprehensive data logging, so a bunch of code was deleted
     here. */
    
  if( maxStep > 0 ) {

      printf("max step size: %e\n", maxStep);
      state.maxStepSize = maxStep;
  
  }  
  state.stepSize = min(state.stepSize, state.maxStepSize);
  state.maxItrs = maxItrs;  
  state.correctionStepDefault = correctionStepSize;
  
  printf( "Overstep tolerance: %f of thickness %f (%f) minminrad: %3.8lf\n", 
	  state.overstepTol, state.injrad*2, 
	  state.injrad*2 - state.overstepTol*state.injrad*2, 
	  state.minminrad );
  
  // use amd column ordering if the user hasn't specified something else
  setenv("COL_ORDERING", "amd", 0);
  
  /* We now initialize the log with a lot of (hopefully) helpful information
     about the run. */

  FILE *logfile;
  
  logfile = fopen(state.logfilename,"w");
  if (logfile == NULL) { 
    fprintf(stderr,"Ridgerunner: Couldn't open logfile %s.\n",state.logfilename);
    exit(1);
  }

  fprintf(logfile,"Ridgerunner logfile.\n");

  #ifdef HAVE_ASCTIME
  #ifdef HAVE_LOCALTIME
  #ifdef HAVE_TIME
  
  fprintf(logfile,"Run began: %s.\n",asctime(localtime(time(NULL))));

  #else

  fprintf(logfile,
	  "Warning: system doesn't have 'time' function.\n"
	  "         data logs won't include time information.\n");

  #endif
  #else

  fprintf(logfile,
	  "Warning: system doesn't have 'localtime' function.\n"
	  "         data logs won't include time information.\n");
  
  #endif
  #else

  fprintf(logfile,
	  "Warning: system doesn't have 'asctime' function.\n"
	  "         data logs won't include time information");
  #endif

  fprintf(logfile,"Initial filename: %s.\n",
	  arg_infile->filename);

  fprintf(logfile,"Command line: ");

  for(i=0;i<argc;i++) {

    fprintf(logfile,"%s ",argv[i]);

  } 
  
  fprintf(logfile,"\n");
  fprintf(logfile,"------------------------------------------------\n");
  
  fclose(logfile);

  /* Now we open the data logfiles in state. */

  char tmpfilename[1024];
  
  for(i=0;i<state.nlogs;i++) {

    sprintf(tmpfilename,"./%s.rr/logfiles/%s.dat",
	    arg_infile->basename,
	    state.logfilenames[i]);

    state.logfiles[i] = fopen(tmpfilename,"w");

    if (state.logfiles[i] == NULL) {

      fprintf(stderr,"ridgerunner: Couldn't open logfile %s.\n",
	      tmpfilename);
      exit(1);

    }

  }

  /* Now we actually run the stepper. */
  
  bsearch_stepper(&link, &state);

  /* We're done. Print some concluding text to the main logfile.*/
 
  logfile = fopen(state.logfilename,"a");
  if (logfile == NULL) { 
    fprintf(stderr,"Ridgerunner: Couldn't open logfile %s.\n",state.logfilename);
    exit(1);
  }

  fprintf(logfile,"--------------------------------------------------\n");
  
  #ifdef HAVE_ASCTIME
  #ifdef HAVE_LOCALTIME
  #ifdef HAVE_TIME
  
  fprintf(logfile,
	  "Run ended: %s.\n",
	  asctime(localtime(time(NULL))),
	  arg_infile->filename);

  #endif
  #endif
  #endif

  /* Now close all of the datalogs. */

  for(i=0;i<state.nlogs;i++) {

    fclose(state.logfiles[i]);

  }

  plc_free(link);
  return kNoErr;
}

static short
rgetline( char* outLine, int inSize, FILE* fp )
{
  outLine[0] ='\0';
  int pos = 0;
  
  do
    {
      outLine[pos] = fgetc(fp);
    } while( outLine[pos++] != '\n' && !feof(fp) );
  outLine[pos] = '\0';
  
  if( strlen(outLine) == 0 )
    return 0;
  return 1;
}

void
initializeState( search_state* state, plCurve** inLink, const char* fname )
{
  int cItr, i, offset=0;
  char	line[1024];
  FILE* fp = fopen(fname, "r"); // we would have failed earlier if this didn't exist
  static char *log_fnames[] = {
    "length","ropelength","strutcount","stepsize","thickness","minrad","residual",
    "maxovermin","rcond","walltime","maxvertforce","convergence","edgelenvariance"
  }; 

  // zero all those pointers
  bzero(state, sizeof(search_state));
  
  strncpy(state->fname, fname, sizeof(state->fname)-1);
  
  state->injrad = 0.5; // fix this for now
  
  state->maxStepSize = 0.1*plCurve_short_edge(*inLink);
  
  if( state->maxStepSize > 1e-3 )
    state->maxStepSize = 1e-3;
  //	state->maxStepSize = 1e-4;
  state->shortest = 2*state->injrad;
  state->totalVerts = 0;
  state->eqThreshold = 1.05;
  state->eqMultiplier = 1;
  state->factor = 1;
  
  state->minminrad = 0.499975;
  
  state->minrad = octrope_minradval(*inLink);
  for( cItr=0; cItr<(*inLink)->nc; cItr++ )
    {
      state->totalVerts += (*inLink)->cp[cItr].nv;
    }
  
  state->stepSize = 0.5*state->maxStepSize;
  if( ((double)1)/(double)state->totalVerts < state->stepSize )
    {
      state->stepSize = ((double)1/(double)state->totalVerts);
    }
  state->compOffsets = malloc(sizeof(int)*(*inLink)->nc);
  for( i=0; i<(*inLink)->nc; i++ )
    {
      state->compOffsets[i] = offset;
      offset += (*inLink)->cp[i].nv;
    }
  state->conserveLength = (int*)calloc((*inLink)->nc, sizeof(int));
  for( i=0; i<(*inLink)->nc; i++ )
    {
      state->conserveLength[i] = 0;
    }
	
  state->residual = 0;
  
  state->checkDelta = 1.0;	// default check delta
  
  fclose(fp);
  
  getrusage(RUSAGE_SELF, &state->frameStart);

  state->nlogs = sizeof(log_fnames);
  for(i=0;i<state->nlogs;i++) {state->logfilenames[i] = log_fnames[i];}

}

void
usage()
{
	fprintf( stdout, 
"ridgerunner \n\n \
-h threshold for check delta, will stop if length improvement < threshld in checkDelta time\n \
-d minrad control -- specify min minrad (should be less than thickness/2)\n \
-o tolerance\t Overstep tolerance (maybe not broken)\n \
-r mult\t Runs until stopping criteron satisfied, then divides, until verts < mult*ropelength\n \
-e mult\t Sets the eq-force multiplier\n \
-b amt\tscales curve by amt\n \
-i injrad\t Sets injrad constraint (probably broken)\n \
-c delta\t Stop if in delta timestep time, we haven't improved length more than 0.05\n \
-l\t\t Forces edge lengths to be approximately equal via Rawdon's tangential runaround\n \
-f file\t Specifies runtile knot (required!)\n \
-m\t\t Will create a directory with timestep frames and strut sets to facilitate movie creation\n \
-x\tmax step size\tspecifies the maxiumum size of integration step\n \
-n\t\t Run will ignore curvature constraint (useful if curvature breaks things)\n \
-a\t\t Autoscale specified knot to thickness 1.0. Useful to save scaling runtime\n \
-t threshold\t Sets additional residual stopping requirement that residual < threshold\n \
-s suppress out Suppresses the progress files in /tmp\n \
-p max itrs\tspecifies the maximum number of integration steps to be taken. default is inf\n \
-v\t\textra verbose filing, will output files _every_ step \n \
-y\t\tEnables geomview visualization fanciness, req. tube, gnuplot, geomview on path\n \
-g\tgraphing with arg: (l)ength (c)locktime (r)cond (m)inrad (r)esidual (s)truts (v)ariance of edge diff from avg (a)ll\n \
-u\tsurface strut gen\tgenerates multiple files in /tmp for use with surfaceBuilder, use with -v\n \
-w\twith: (f)ast correction stepper\n \
-z\tavoid conflicts in /tmp/ between multiple users by stamping with pid and file name\n \
-j\tpaper info, print avg times through bsearch, correction steps to green in /tmp/ and rcond every step if its graphing on\n \
-k\tcorrection step size, default 0.25\n"
);
	exit(kNoErr);
}

void
reload_handler( int sig )
{
	if( sig == SIGUSR1 )
	{
		// immediately reload config
		
		// reinstall handler
		signal(SIGUSR1, reload_handler);
	}
}
