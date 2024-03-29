/*
 *  ridgerunner_main.c
 *  ridgerunner
 *
 *  Created by Michael Piatek on Fri May 28 2004.

 Copyright Jason Cantarella.

 This file is part of ridgerunner. ridgerunner is free software: you can
 redistribute it and/or modify it under the terms of the GNU General
 Public License as published by the Free Software Foundation, either
 version 3 of the License, or (at your option) any later version.

 ridgerunner is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 for more details.  You should have received a copy of the GNU General
 Public License along with ridgerunner. If not, see
 <https://www.gnu.org/licenses/>.

*/
 
#include"ridgerunner.h"

/* This "main" file is where all of the global variables live. */

void usage();
void reload_handler( int sig );
void initializeState( search_state* state );

/* globals moved to ridgerunner.h */

#ifndef min
  #define min(a,b) ((a)<(b)) ? (a) : (b)
#endif

double parse_size(const char *sizestring);
void parse_display_arg(search_state *inState, struct arg_str *display);
	
int
main( int argc, char* argv[] )
{

  struct arg_file *arg_infile = arg_file1(NULL,NULL,"<VECT file>","input file");

  struct arg_dbl  *arg_lambda = arg_dbl0("l","Lambda","<double>",
					 "minimum radius of curvature for unit rope");
  struct arg_dbl  *arg_tuberadius = arg_dbl0("t","TubeRadius","<double>","radius of tube around core curve");
  struct arg_rem  *arg_bl0 = arg_rem("","");
  struct arg_rem  *arg_curveopts = arg_rem("","Operations on input curve before run");
  struct arg_rem  *arg_bl1 = arg_rem("","");

  struct arg_lit  *arg_autoscale = arg_lit0("a","Autoscale","always scale curve "
					    "to thickness .501");

  struct arg_lit  *arg_continue = arg_lit0("c","Continue","continue mode -- no initial rescaling of curve");

  struct arg_lit  *arg_timewarp = arg_lit0(NULL,"Timewarp","accelerate shrinking of strut-free sections of curve");

  struct arg_lit  *arg_mangler = arg_lit0(NULL,"MangleMode","instead of trying to reduce length, try to change configuration");

  struct arg_int  *arg_manglesteps = arg_int0(NULL,"MangleSteps","<#steps>","switch mangling algorithm after this many steps");

  struct arg_lit  *arg_cg = arg_lit0(NULL,"ConjugateGradient","turn on conjugate gradient mode (experimental)");

  struct arg_dbl  *arg_minstep = arg_dbl0(NULL,"MinStep","<double>","minimum size for step");

  /*  struct arg_dbl  *arg_resolution = arg_dbl0("r","res","<verts/rop>",
      "spline curve to this resolution"); */

  struct arg_lit  *arg_eqit = arg_lit0(NULL,"EqOn","automatically equilateralizes when max/min edge > 3"); 

  struct arg_rem  *arg_bl2 = arg_rem("","");
  struct arg_rem  *arg_stopopts = arg_rem("","Stopping Criteria (stop when)");
  struct arg_rem  *arg_bl3 = arg_rem("","");

  struct arg_dbl  *arg_stop20 = arg_dbl0(NULL,"Stop20","<change in Rop>",
					 "rop change over last 20 "
					 "steps < value");

  struct arg_dbl  *arg_stopRes = arg_dbl0(NULL,"StopResidual","<fraction>",
					  "fraction of "
					  "curvature force after projection < value");

  struct arg_int  *arg_stopSteps = arg_int0("s","StopSteps","<#steps>",
					    "steps >= value");

  struct arg_int  *arg_stopTime = arg_int0(NULL,"StopTime","<minutes>","stop after this many (wall clock) minutes");

  struct arg_rem  *arg_bl4 = arg_rem("","");
  struct arg_rem  *arg_fileopts = arg_rem("","File Output Options");
  struct arg_rem  *arg_bl5 = arg_rem("","");
  
  struct arg_lit  *arg_suppressfiles = arg_lit0(NULL,"NoOutputFiles",
						"don't save intermediate VECT files for animation");

  struct arg_lit  *arg_nolsqrlog = arg_lit0(NULL,"NoLsqrLog","don't log output from lsqr");

  struct arg_lit  *arg_nopng = arg_lit0(NULL,"NoPNGOutput","don't tube or draw final snapshot with POVRAY");

  struct arg_str  *arg_maxlogsize = arg_str0(NULL,"MaxLogSize","<100K>",
					     "maximum log file size in bytes, K, or M");

  struct arg_str  *arg_maxvectdirsize = arg_str0(NULL,"MaxVectDirSize","<100M>",
					      "maximum size of movie dir in bytes, K, or M");

  struct arg_lit  *arg_nocolor = arg_lit0(NULL,"NoColor", "don't color curves at all");

  struct arg_lit  *arg_nohighlight = arg_lit0(NULL,"NoHighlight", "don't highlight straight or kinked portions of curve"); 
  
  // struct arg_rex  *arg_outpath = arg_rex0(NULL,"OutPath","/*/",
  //					  "</home/../outdir/>",
  //					  0,"path for output files");

  struct arg_rem  *arg_bl6 = arg_rem("","");
  struct arg_rem  *arg_progopts = arg_rem("","Program and Algorithm Options");
  struct arg_rem  *arg_bl7 = arg_rem("","");

  struct arg_lit  *arg_quiet  = arg_lit0("q","Quiet","quiet");
  struct arg_lit  *arg_verbose = arg_lit0("v","Verbose","verbose");
  struct arg_lit  *arg_vverbose = arg_lit0(NULL,"VeryVerbose","very verbose (debugging)");

  struct arg_lit  *arg_help = arg_lit0("h","help","print this help and exit");

  struct arg_lit  *arg_animation = arg_lit0(NULL,"AnimationStepper","produce smooth, but slow, movie. won't converge to low residual");

  struct arg_dbl  *arg_cstep_size = arg_dbl0("k","CorrectionStepSize","<fraction>",
					     "initial size of Newton correction step");
  
  struct arg_dbl  *arg_eqmult = arg_dbl0(NULL,"EqMultiplier","<scalar>","increase to "
					 "make 'equilateralization force' stronger");

  struct arg_lit  *arg_eq = arg_lit0(NULL,"EqForceOn","turn on equilateralization force during the run (can interfere with stepper)");
  
  struct arg_dbl  *arg_overstep = arg_dbl0("o","OverstepTol","<x>",
					   "start correction step if "
					   "thickness < tuberad*(1 - x)");

  struct arg_dbl  *arg_mroverstep = arg_dbl0(NULL,"MinRadOverstepTol","<x>",
					     "start correction step if "
					     "minrad < lambda*(1 - x)");

  struct arg_dbl  *arg_maxstep = arg_dbl0(NULL,"MaxStep","<scalar>","maximum stepsize "
					  "for a length-reduction step");

  struct arg_int  *arg_maxcorr = arg_int0(NULL,"MaxCorrectionAttempts","<n>",
					  "maximum # of Newton steps in error correction");

  struct arg_lit  *arg_trynewton = arg_lit0(NULL,"TryNewtonCorrection","always attempt to correct with Newton steps");

  struct arg_lit  *arg_sono = arg_lit0(NULL,"SONO","emulate Pieranski's SONO");

  struct arg_lit  *arg_spin = arg_lit0(NULL,"spinForce","add spin force");

  struct arg_int  *arg_snapinterval = arg_int0(NULL,"SnapshotInterval","<n>",
					       "save a complete snapshot "
					       "of computation every <n> steps");
  
  struct arg_lit *arg_rcond = arg_lit0(NULL,"Rcond","log condition number for matrices");

  struct arg_lit *arg_sfr = arg_lit0(NULL,"StrutFreeResidual","log the portion of residual on strut-free sections of curve");

  struct arg_str  *arg_symmetry = arg_strn(NULL,"Symmetry","Z/pZ,D2,cplanes",
					   0,1,"rotation/reflection symmetry group (around z-axis), or across all coordinate planes to enforce on curve");

  struct arg_rem  *arg_bl8 = arg_rem("","");
  struct arg_rem  *arg_dispopts = arg_rem("","Display Options");
  struct arg_rem  *arg_bl9 = arg_rem("","");

  struct arg_str  *arg_display = arg_strn(NULL,"Display","ropelength/strutcount/thickness",0,128,"name of logfile to display on screen during run");
  
  struct arg_end *end = arg_end(20);
  
  void *argtable[] = {arg_infile,arg_lambda,arg_tuberadius,
		      arg_bl0,arg_curveopts,arg_bl1,
		      arg_autoscale,/* arg_resolution, */arg_continue,

		      arg_bl6,arg_progopts,arg_bl7,
		      arg_quiet,arg_verbose,arg_vverbose,arg_help,
		      arg_animation, arg_timewarp,arg_mangler,arg_manglesteps,arg_cg,arg_eqit,/*arg_cstep_size,arg_maxcorr,*/
		      arg_eqmult,arg_eq,arg_overstep,arg_mroverstep,
		      arg_maxstep,arg_minstep,arg_snapinterval,arg_trynewton,
		      arg_sono, arg_spin, arg_rcond, arg_sfr, arg_symmetry,

		      arg_bl8,arg_dispopts,arg_bl9,
		      arg_display,

		      arg_bl2,arg_stopopts,arg_bl3,
		      arg_stop20,arg_stopRes,arg_stopSteps,arg_stopTime,

		      arg_bl4,arg_fileopts,arg_bl5,
		      arg_suppressfiles, arg_nopng,
		      arg_nolsqrlog, arg_maxlogsize, 
		      arg_maxvectdirsize, arg_nocolor, arg_nohighlight, /* arg_outpath, */
		      
		      end};
  int nerrors;
  int snapinterval = 10000;

  plCurve*	link = NULL;
  FILE*		linkFile = NULL;
  search_state	state;
  short		movie = 0; // autoscale
  /*  double	refineUntil = -1; *//*  Used to be an int */
  
  double	overstepTol=0.0001; /* Should this only be for tube_radius = 1? */
  double	residualThreshold = 0; /* Don't stop on residual by default. */
  
  /*double	checkThreshold = 0.05;
    double	checkDelta = 1; */

  double        stop20 = -1000;        /* For this to fail, ropelength would have to climb! */
  double        eqMult = 0.0;

  double	maxStep = -1,minStep = 1e-6;
  long		maxItrs = 10000000;
  double	correctionStepSize = 0.25;
  double	minradOverstepTol = 0.00005; /* Used to be abs val of 0.499975 */
  int           i,j;

  /* Display opening message. */

  printf("Ridgerunner %s\n",PACKAGE_VERSION);  
  printf("Built %s, %s.\n", __DATE__ , __TIME__ );

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

  /* Now we move the arguments into local variables. */
   
  if (arg_cstep_size->count > 0) { correctionStepSize = arg_cstep_size->dval[0]; }

  if (arg_maxcorr->count > 0) { gMaxCorrectionAttempts = arg_maxcorr->ival[0]; }

  /* Note: there used to be a way to turn on "gFastCorrectionSteps" from the cmd line */
  
  if (arg_quiet->count > 0) { gQuiet = 1; }

  if (arg_maxstep->count > 0) { maxStep = arg_maxstep->dval[0]; }

  if (arg_minstep->count > 0) { minStep = arg_minstep->dval[0]; }

  /* Note: there used to be a way to turn on "gSurfaceBuilding" 
     and "gVerboseFiling" */

  if (arg_snapinterval->count > 0) {snapinterval = arg_snapinterval->ival[0];}

  if (arg_rcond->count > 0) { gNoRcond = 0; } else {gNoRcond = 1;}

  if (arg_sfr->count > 0) { gStrutFreeResidual = 1; } else { gStrutFreeResidual = 0;}

  if (arg_cg->count > 0) {gConjugateGradient = 1;} else {gConjugateGradient = 0;}

  if (arg_mroverstep->count > 0) {minradOverstepTol = arg_mroverstep->dval[0];}

  if (arg_trynewton->count > 0) {gTryNewton = 1;} 

  if (arg_suppressfiles->count > 0) { gSuppressOutput = 1; }

  if (arg_nolsqrlog->count > 0) { gLsqrLogging = 0; }

  if (arg_overstep->count > 0) { overstepTol = arg_overstep->dval[0]; }

  if (arg_stopSteps->count > 0) { maxItrs = arg_stopSteps->ival[0]; }

  if (arg_stop20->count > 0) { stop20 = arg_stop20->dval[0]; }

  if (arg_mangler->count > 0) { 

    if (arg_animation->count == 0) {

      printf("Warning: --MangleMode requires --AnimationStepper. Faking --AnimationStepper.\n"); 
      gAnimationStepper = 1; 
      eqMult = 1.0;

    }

    gMangleMode = 1;

  }

  /*   if (arg_resolution->count > 0) { refineUntil = arg_resolution->dval[0]; } */
 
  if (arg_stopRes->count > 0) { residualThreshold = arg_stopRes->dval[0]; }

  if (arg_eqmult->count > 0) { eqMult = arg_eqmult->dval[0]; }

  if (arg_eq->count > 0) { eqMult = 1.0; } else {eqMult = 0.0; }

  /* Note: There used to be a way to change "tube_radius", "scaleamt", and "fixlengths" 
     from the cmdline. */
  
  if (arg_quiet->count > 0) { VERBOSITY = 0; }

  if (arg_verbose->count > 0) { VERBOSITY = 5; }

  if (arg_vverbose->count > 0) { VERBOSITY = 10; }

  if (arg_lambda->count > 0) { gLambda = arg_lambda->dval[0]; }

  if (arg_timewarp->count > 0) { gNoTimeWarp = 0; }

  if (arg_animation->count > 0) { gAnimationStepper = 1; eqMult = 1.0; }
  
  if (arg_eqit->count > 0) { gEqIt = 1; }

  if (arg_sono->count > 0) { gSONO = 1; }

  if (arg_spin->count > 0) { gSpinForce = 1; }

  if (arg_symmetry->count > 0) { 

    if (arg_animation->count == 0) {

      printf("Warning: --Symmetry requires --AnimationStepper. Faking --AnimationStepper.\n"); 
      gAnimationStepper = 1; 
      eqMult = 1.0;

    }

    if (gEqIt == 1) {

      printf("Warning: --EqOn is incompatible with --Symmetry. Turning off --EqOn.\n");
      gEqIt = 0;

    }

    /* We will revisit arg_symmmetry later AFTER loading the curve. */

  }

  /* Note: There used to be a way to set "movie", "gPaperInfoinTmp", "ignorecurvature",
     and "fancyviz" from the cmdline. */

  /* We now set as much of "state" as we can from this amount of data. */
  
  initializeState( &state );  

#define BYTES_PER_ENTRY 32.0

  if (arg_maxlogsize->count > 0) { 
    
    state.maxlogsize = (int)(parse_size(arg_maxlogsize->sval[0])/BYTES_PER_ENTRY);

  } else {

    state.maxlogsize = (int)(10*1024*1024.0/BYTES_PER_ENTRY);  /* 10 Mb */

  }

  state.stop20 = stop20;
  state.maxItrs = maxItrs;  
  state.residualThreshold = residualThreshold;
  
  if (arg_stopTime->count > 0) {

    state.stopTime = 60*arg_stopTime->ival[0]; // The argument is in minutes, the var in secs

  } else {

    state.stopTime = 60*60*24*7; // One week.

  }
  
  if (arg_tuberadius->count > 0) {

    state.tube_radius = arg_tuberadius->dval[0];

  } else {

    state.tube_radius = 0.5;

  }

  state.newDir = NULL;
  state.lastGrad = NULL;
  state.correctionStepDefault = correctionStepSize;
  state.movie = movie;
  state.moviefactor = 1.0;
  
  state.overstepTol = overstepTol;
  state.minradOverstepTol = minradOverstepTol;
  state.minminrad = gLambda*state.tube_radius*(1 - minradOverstepTol);
  
  state.eqMultiplier = eqMult;
  state.snapinterval = snapinterval;

  state.manglemode = torusrotate;
  state.this_manglemode_startstep = 0;

  if (arg_manglesteps->count > 0) {

    state.steps_in_mode = arg_manglesteps->ival[0];

  } else {

    state.steps_in_mode = 30000;

  }

  /* We now do a little bit of filename mangling. */

  char *extptr;

  strncpy(state.basename,arg_infile->basename[0],sizeof(state.basename));
  extptr = strstr(state.basename,arg_infile->extension[0]);
  *extptr = 0; /* This chops the extension off of basename. */
  strncpy(state.fname,arg_infile->filename[0],sizeof(state.fname));

  sprintf(state.workingfilename,"./%s.rr/%s.working.vect",
	  state.basename,state.basename);
  sprintf(state.workingstrutname,"./%s.rr/%s.working.struts",
	  state.basename,state.basename);
  sprintf(state.vectprefix,"./%s.rr/vectfiles/%s",
	  state.basename,state.basename);
  sprintf(state.snapprefix,"./%s.rr/snapshots/%s",
	  state.basename,state.basename);
  sprintf(state.logfilename,"./%s.rr/%s.log",
	  state.basename,state.basename);
  sprintf(state.fprefix,"./%s.rr/",state.basename);

  /* We now open the local directory for storing output. */
  
  char cmdline[1024];
  
  sprintf(cmdline,"rm -fr %s.rr",state.basename);
  system_or_die(cmdline, __FILE__ , __LINE__ );

  sprintf(cmdline,"mkdir %s.rr",state.basename);
  system_or_die(cmdline, __FILE__ , __LINE__ );
  
  sprintf(cmdline,"mkdir %s.rr/logfiles",state.basename);
  system_or_die(cmdline, __FILE__ , __LINE__ );
  
  if (arg_suppressfiles->count == 0) {
    
    sprintf(cmdline,"mkdir %s.rr/vectfiles",state.basename);
    system_or_die(cmdline, __FILE__ , __LINE__ );
    
  }
  
  sprintf(cmdline,"mkdir %s.rr/snapshots",state.basename);
  system_or_die(cmdline, __FILE__ , __LINE__ );

  /* We now initialize the log with a lot of (hopefully) helpful information
     about the run. */
  
  gLogfile = fopen(state.logfilename,"w");
  if (gLogfile == NULL) { 
    fprintf(stderr,"Ridgerunner: Couldn't open logfile %s.\n",state.logfilename);
    exit(1);
  }
  
  fprintf(gLogfile,"Ridgerunner logfile.\n");

  fprintf(gLogfile,"Ridgerunner %s\n",PACKAGE_VERSION);
   
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

  fprintf(gLogfile,"Stopping criteria: Max # of steps %ld, Min residual %g, Stop20 %g, stopTime %d (sec).\n",
	  state.maxItrs,state.residualThreshold,state.stop20,state.stopTime);

#ifdef HAVE_GETPID

  fprintf(gLogfile,"Process ID: %d.\n",(int)(getpid()));

#endif
  
  fprintf(gLogfile,"------------------------------------------------\n");
  
  /* We now open the given link file. */
  
  linkFile = fopen_or_die(state.fname,"r", __FILE__ , __LINE__ );
  
  int plc_read_error_num;
  char plc_error_string[1024];
  size_t plc_error_size = {1024};
  char  errmsg[1024];
  
  link = plc_read(linkFile,&plc_read_error_num,plc_error_string,plc_error_size); 

  if (plc_read_error_num != 0) {
    
    sprintf(errmsg,"ridgerunner: file %s had \n"
	    "             plc_read error %s.\n",
	    state.fname,plc_error_string);
    FatalError(errmsg, __FILE__, __LINE__ );
    
  }
  
  fclose(linkFile);

  /* We have some fairly serious randomness needs when running manglemode. Try to get a seed from dev/random. */
  
  FILE *devrand;
  unsigned int seed;

  devrand = fopen("/dev/random","r");
  if (devrand != NULL) {
  
    fread(&seed,sizeof(unsigned int),1,devrand);
    fclose(devrand);
    logprintf("Seeded from dev/random with seed %d.\n",seed);

  } else {
 
    seed = time(NULL);
    logprintf("Seeded from time with seed %d.\n",seed);

  }

  srand(seed);


  if (gMangleMode) {

    logprintf("Running in MangleMode.\n");

  }

  if (gAnimationStepper) { 

    logprintf("Running with Animation Stepper.\n");
    
  }

  if (arg_symmetry->count > 0) {

    logprintf("Running with symmetry %s.\n",arg_symmetry->sval[0]);

  }

  logprintf("Loaded %d component, %d vertex plCurve from %s.\n",
	    link->nc,plc_num_verts(link),state.fname);
  
  /* We archive a version of the initial input file in the run directory. */
  
  char filename[1024];
  FILE *savefile;
  
  sprintf(filename,"./%s.rr/%s.vect",state.basename,state.basename);
  savefile = fopen_or_die(filename,"w", __FILE__, __LINE__ );
  plc_write(savefile,link);
  fclose(savefile);
  
  logprintf("Saved copy of %s to %s.\n",state.fname,filename);

  /* We now handle adding color to the curve if needed. */

  if (arg_nocolor->count == 0) {  /* The nocolor flag is not on, so we will color the curve, stomping existing colors. */

    for(i=0;i<link->nc;i++) {

      plc_resize_colorbuf(link,i,link->cp[i].nv);

      for(j=0;j<link->cp[i].nv;j++) {

	link->cp[i].clr[j] = gTubeColors[i % gNumTubeColors ];

      }

    }

  }
 
  /* Now complete state initializations which depend on the curve. */

  state.maxStepSize = 100; // 0.1*plCurve_short_edge(link); /* We are experimenting with large steps. */
  state.minStep = minStep;
  state.minrad = octrope_minradval(link);
  state.totalVerts = plc_num_verts(link);
  state.stepSize = 0.01;

  if( ((double)1)/(double)state.totalVerts < state.stepSize )
    {
      state.stepSize = ((double)1/(double)state.totalVerts);
    }

  if( state.maxStepSize > 1e3 ) { state.maxStepSize = 1e3; }

  if( maxStep > 0 ) {
    
    printf("max step size: %e\n", maxStep);
    state.maxStepSize = maxStep;
    
  } 
  
  state.stepSize = min(state.stepSize, state.maxStepSize);

  state.compOffsets = malloc(sizeof(int)*(link->nc));
  fatalifnull_(state.compOffsets);
  
  int offset = 0;

  for( i=0; i<link->nc; i++ )    {
    
    state.compOffsets[i] = offset;
    offset += link->cp[i].nv;
    
  }

  state.conserveLength = NULL; /* We aren't using this. */
  
  /* Now perform initial operations. */
  
  
/*   if (arg_resolution->count > 0) { */
    
/*     logprintf("Original curve has %d verts, %g verts/rop.\n",plc_num_verts(link), */
/* 	      plc_num_verts(link)*octrope_thickness(link,NULL,0,gLambda)/plc_arclength(link,NULL)); */

/*     tempLink = plCurve_fixresolution(link,arg_resolution->dval[0]); */
/*     plc_free(link); */
/*     link = tempLink; */
    
/*     logprintf("Splined to resolution at least %g verts/rop. \nNew curve has %d verts, %g verts/rop.\n", */
/* 	   arg_resolution->dval[0],plc_num_verts(link),plc_num_verts(link)*octrope_thickness(link,NULL,0,gLambda)/plc_arclength(link,NULL)); */
        
/* 	   }  */
  
  double thickness;
  thickness = octrope_thickness(link,NULL,0,gLambda);
  double t_margin = 0.0001;

  if (arg_continue->count == 0) {

    if(( arg_autoscale->count > 0 || thickness < state.tube_radius + t_margin)) {  
      
      logprintf("Curve has thickness %g. Scaling to thickness %g.\n",
		thickness,state.tube_radius + 10*t_margin);
      
      plc_scale(link,(state.tube_radius + 10*t_margin)/thickness);
      thickness = octrope_thickness(link,NULL,0,gLambda);
      
      logprintf("Scaled curve has thickness %g.\n",thickness);
      
      if (thickness < state.tube_radius + t_margin) {

	/* Make sure that our thickness is close to the desired value */
	
	sprintf(errmsg,"ridgerunner: Failed to scale %s to thickness %g."
		"             Aborting run.\n",state.basename,state.tube_radius + t_margin);
	FatalError(errmsg, __FILE__ , __LINE__ );
	
      } else {
	
	printf(" Autoscale selftest ok.\n");
	
      }
      
    }
    
  }

  state.oktoscale = true;
  
  plc_constraint *thisCst;
  
  for(thisCst = link->cst; 
      thisCst != NULL; 
      thisCst = thisCst->next) {

    /* Look for "line" constraints. For various technical reasons,
       these are best handled by user replacement with a pair of
       "plane" constraints, so we refuse to run on anything with these
       in place. */

    if (thisCst->kind == line) {

      sprintf(errmsg,
	      "ridgerunner: This version does not implement 'line' constraints\n"
	      "             Replace them with a pair of 'plane' constraints.\n");
      FatalError(errmsg, __FILE__ , __LINE__ );

    }

    if (thisCst->kind == fixed) {

      state.oktoscale = false;

    }

  }

      
  /* Now save the modified file so that we record what we actually ran with. */

  sprintf(filename,"./%s.rr/%s.atstart.vect",
	  state.basename,state.basename);
  savefile = fopen_or_die(filename,"w", __FILE__ , __LINE__ );
  plc_write(savefile,link);
  fclose(savefile);

  logprintf("Rerez'd, autoscaled, eq'd file written to %s.\n",filename);
     
  /* We now initialize the link-dependent parts of "state". */

  /* Piatek's realtime graphing fanciness has been replaced by
     comprehensive data logging, so a bunch of code was deleted
     here. */

  sprintf(filename,"./%s.rr/%s.atstart.vect",
	  state.basename,state.basename);
  
  struct stat savefilebuf;

  if (stat(filename,&savefilebuf)) {

    sprintf(errmsg,"ridgerunner: Failed to stat the file %s.\n",filename);
    FatalError(errmsg, __FILE__ , __LINE__ );

  }

  if (arg_maxvectdirsize->count > 0) { 
    
    state.maxmovieframes = (int)
      (parse_size(arg_maxvectdirsize->sval[0]))/((float)(savefilebuf.st_size));
       
  } else {
    
    state.maxmovieframes = (int)(50*1024*1024.0/((float)(savefilebuf.st_size)));  /* 50 Mb */

  }

  if (state.maxmovieframes % 2 == 1) { state.maxmovieframes--; } /* Make sure this is even. */

  state.dumpAxb = false;
  
  state.minrad = octrope_minradval(link);
  state.thickness = octrope_thickness(link,NULL,0,gLambda);
  state.ropelength = octrope_ropelength(link,NULL,0,gLambda);
  state.shortest = octrope_poca(link,NULL,0);
  state.residual = DBL_MAX;
  state.strutfreeresidual = DBL_MAX; // This is a nonsense value to be logged if we don't turn this computation on
  state.score = stepScore(link,&state,NULL,0);

  state.totalVerts = plc_num_verts(link);

  printf( "Overstep tolerance: %f of thickness %f (%f) minminrad: %3.8lf\n", 
	  state.overstepTol, state.tube_radius, 
	  state.tube_radius*(1 - state.overstepTol), 
	  state.minminrad );
  
  // use amd column ordering if the user hasn't specified something else
  setenv("COL_ORDERING", "amd", 0);
    
  open_runtime_logs(&state,'w');
  parse_display_arg(&state,arg_display);
  init_runtime_display(&state);

  logprintf("ridgerunner: Starting run. Will stop if \n"
	 "               step number >= %ld \n"
	 "               residual < %g \n"
	 "               ropelength decrease over last 20 steps < %g.\n\n",
	 state.maxItrs,state.residualThreshold,state.stop20);

  /* We are now going to parse the symmetry argument, if present, 
     and attempt to associate the corresponding symmetry group 
     with the actual link plCurve. */

  if (arg_symmetry->count > 0) {

    /* The first thing we have to do is try to scanf the string and recognize it. */

    int p;
    if (sscanf(arg_symmetry->sval[0],"Z/%dZ",&p) == 1) {

      plc_vector zaxis = {{0,0,1}};
      /* We recognize the symmetry as rotational. */

      logprintf("Recognized symmetry as rotational (Z/%dZ) around z-axis.\n",p);
      link->G = plc_rotation_group(link,zaxis,p);
      
      if (link->G == NULL) {

	char errmsg[1024];
	sprintf(errmsg,"Couldn't build Z/%dZ symmetry on input link. Aborting run.\n",p);
	FatalError(errmsg, __FILE__, __LINE__);

      }

      plc_symmetrize(link);
      logprintf("Built Z/%dZ symmetry on link and symmetrized to initial error %g.\n",p,plc_symmetry_group_check(link));
      
    } else if (strcmp(arg_symmetry->sval[0],"D2") == 0) {

      plc_vector zaxis = {{0,0,1}};
      logprintf("Recognized symmetry as reflection (D2) over z-axis.\n");
      link->G = plc_reflection_group(link,zaxis);

      if (link->G == NULL) {

	char errmsg[1024];
	sprintf(errmsg,"Couldn't build D2 reflection symmetry on input link. Aborting run.\n");
	FatalError(errmsg, __FILE__, __LINE__);

      }

      plc_symmetrize(link);
      logprintf("Built D2 reflection symmetry on link and symmetrized to initial error %g.\n",p,plc_symmetry_group_check(link));
    
    } else if (strcmp(arg_symmetry->sval[0],"cplanes") == 0) {

      logprintf("Recognized symmetry as reflections over coordinate planes (cplanes).\n");
      link->G = plc_coordplanes_reflection_group(link);

      if (link->G == NULL) {

	char errmsg[1024];
	sprintf(errmsg,"Couldn't build cplanes reflection symmetry on input link. Aborting run.\n");
	FatalError(errmsg, __FILE__, __LINE__);

      }

      plc_symmetrize(link);
      logprintf("Built cplanes reflection symmetry on link and symmetrized to initial error %g.\n",
		p,plc_symmetry_group_check(link));

    } else {

      	char errmsg[1024];
	sprintf(errmsg,"Couldn't parse symmetry group %s. Aborting run.\n",arg_symmetry->sval[0]);
	FatalError(errmsg, __FILE__, __LINE__);

    }

  }
      
  /* Now we actually run the stepper. */
  
  bsearch_stepper(&link, &state);

  /* We're done. Print some concluding text to the main logfile.*/
 
  fprintf(gLogfile,"--------------------------------------------------\n");
  
  #ifdef HAVE_ASCTIME
  #ifdef HAVE_LOCALTIME
  #ifdef HAVE_TIME
  
  time_t end_time;
  end_time = time(NULL);

  fprintf(gLogfile,
	  "Run ended: %s.\n",
	  asctime(localtime(&end_time)));


  #endif
  #endif
  #endif

  /* We now color the straight segments, kinks, helices and so forth for final visualization purposes. */

  if (arg_nohighlight->count == 0) { 

    highlight_curve(link,&state);

  }

  /* Now write the concluding file to the appropriate directory. */

  char tmpfilename[2048];

  sprintf(tmpfilename,"./%s.rr/%s.final.vect",
	  state.basename,state.basename);  
  savefile = fopen_or_die(tmpfilename,"w", __FILE__ , __LINE__ );
  plc_write(savefile,link);
  fclose(savefile);

  /* Now write the final strut set (in octrope format). */

  sprintf(tmpfilename,"./%s.rr/%s.final.struts",
	  state.basename,state.basename);
  
  savefile = fopen_or_die(tmpfilename,"w", __FILE__ , __LINE__ );

  octrope_strutfile_write(state.lastStepStrutCount,
			  state.lastStepStruts,
			  state.lastStepMinradStrutCount,
			  state.lastStepMRlist,
			  savefile);

  fclose(savefile);

  /* We can provide some convenient post-processing of the results if
     we have various utilities on this system. */

  if (arg_nopng->count == 0) {
    
#ifdef HAVE_TUBE

    char tmpcommand[8096];
    
    sprintf(tmpcommand,"cd ./%s.rr; tube -r %g %s.final.vect",
	    state.basename,state.tube_radius,state.basename);  
    
    system(tmpcommand);
    
#ifdef HAVE_POVRAY 

#ifdef HAVE_POVSNAP 

    printf("Running povsnap to generate image of final configuration.\n");
    sprintf(tmpcommand,"cd ./%s.rr; orient -a 2 %s.final.tube.off; povsnap -q -s %s.final.tube.off; rm -fr %s.final.tube",state.basename,state.basename,state.basename,state.basename);
    system(tmpcommand);
    
#endif

#endif

#endif

  }

  /* Now clean up allocated memory and close files. */

  close_runtime_logs(&state);
  close_runtime_display();
  free_search_state(&state);
  plc_free(link);

  fclose(gLogfile);

  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

  /* Finally, signal the user that we terminated normally. */

  printf("ridgerunner: Run complete.\n");

  return kNoErr;
}

void
initializeState( search_state* state )

     /* Performs all initializations which _don't_ depend on inLink. */
{
  int  i;

  static char *log_fnames[] = {
    "length","ropelength","strutcount","stepsize","thickness","minrad","residual",
    "maxovermin","rcond","walltime","maxvertforce","csteps_to_converge","edgelenvariance",
    "lsqroutput","memused","effectiveness","score","strutfreeresidual","symmetryerror"
  }; 

  search_state *zerostate;
  
  // zero everything
  zerostate = calloc(1,sizeof(search_state));
  fatalifnull_(zerostate);
  *state = *zerostate;
  free(zerostate);

  state->eqMultiplier = 1;
  state->factor = 1;
  
  state->minminrad = 0.499975;
  

  state->tsnnls_evaluations = 0;
  state->octrope_calls = 0;
  state->cstep_count = 0;	
  state->residual = 0;

  assert(kTotalLogTypes == sizeof(log_fnames)/sizeof(char *));
  for(i=0;i<kTotalLogTypes;i++) {state->logfilenames[i] = log_fnames[i];}
  state->loginterval = 1;
}

void free_search_state(search_state *inState)

     /* Free memory. */

{

  free(inState->compOffsets);  
  free(inState->lastStepStruts);
  free(inState->lastStepMRlist);
  
}

double parse_size(const char *argstring)

/* Function parses a human-readable size and returns a value in bytes. */

{

  char sizestring[1024];
  double unitcvrt = 1;
  int len;
  double sizef;
  double result;

  strncpy(sizestring,argstring,sizeof(sizestring));

  /* First we parse and remove any letter representing units (M or K). */
  
  len = strlen(sizestring);
  
  if (sizestring[len-1] == 'G' || sizestring[len-1] == 'g') { /* Gigabytes */

    unitcvrt = 1024*1024*1024;
    sizestring[len-1] = 0;
    
  } else if (sizestring[len-1] == 'M' || sizestring[len-1] == 'm') { /* Megabytes */

    unitcvrt = 1024*1024;
    sizestring[len-1] = 0;

  } else if (sizestring[len-1] == 'K' || sizestring[len-1] == 'k') { /* kilobytes */

    unitcvrt = 1024;
    sizestring[len-1] = 0;

  } 

  /* Now we attempt to parse the number. */
  
  if (sscanf(sizestring,"%lg",&sizef) != 1) {

    logprintf("parse_size: Could not parse size argument %s.\n",sizestring);
    exit(0);

  }

  if (2.0*1024.0*1024.0*1024.0 < unitcvrt*sizef) {

    logprintf("parse_size: Warning! Size greater than 2GB will be truncated to 2GB.\n");    
    result = 2*1024*1024*1024.0;

  } else {

    result = floor(unitcvrt*sizef);
    
  }

  return result;

}


