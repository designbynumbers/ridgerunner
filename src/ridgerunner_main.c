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
		{"tube_radius", required_argument, NULL, 'd'},
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
					    "to thickness 1.001");

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
  search_state	state;
  short		autoscale = 0, movie = 0;
  double	refineUntil = -1;  /* Used to be an int */
  
  double	tube_radius = 0.5, overstepTol=0.0001; /* Should this be tube_radius = 1? */
  double	residualThreshold = 65000;
  
  /*double	checkThreshold = 0.05;
    double	checkDelta = 1; */

  double        stop20 = 0.0;

  double        eqMult = 2.0;


  char		fname[1024];
  double	maxStep = -1;
  long		maxItrs = -1;
  int		graphIt[kTotalGraphTypes];
  double	correctionStepSize = 0.25;
  double	minradOverstepTol = 0.499975;
  int           i;
  
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

  /* Now we move the arguments into local variables. */
   
  if (arg_cstep_size->count > 0) { 

    correctionStepSize = arg_cstep_size->dval[0];

  }

  /* Note: there used to be a way to turn on "gFastCorrectionSteps" from the cmd line */
  
  if (arg_quiet->count > 0) { gQuiet = 1; }

  if (arg_maxstep->count > 0) { maxStep = arg_maxstep->dval[0]; }

  /* Note: there used to be a way to turn on "gSurfaceBuilding" and "gVerboseFiling" */

  if (arg_mroverstep->count > 0) {minradOverstepTol = gLambda - arg_mroverstep->dval[0];}

  if (arg_suppressfiles->count > 0) { gSuppressOutput = 1; }

  if (arg_overstep->count > 0) { overstepTol = arg_overstep->dval[0]; }

  if (arg_stopSteps->count > 0) { maxItrs = arg_stopSteps->ival[0]; }

  if (arg_stop20->count > 0) { stop20 = arg_stop20->dval[0]; }

  if (arg_resolution->count > 0) { refineUntil = arg_resolution->dval[0]; }
 
  if (arg_stopRes->count > 0) { residualThreshold = arg_stopRes->dval[0]; }

  if (arg_eqmult->count > 0) { eqMult = arg_eqmult->dval[0]; }

  /* Note: There used to be a way to change "tube_radius", "scaleamt", and "fixlengths" 
     from the cmdline. */
  
  if (arg_autoscale->count > 0) { autoscale = 1; }
  
  /* Note: There used to be a way to set "movie", "gPaperInfoinTmp", "ignorecurvature",
     and "fancyviz" from the cmdline. */

  /* We now set as much of "state" as we can from this amount of data. */
  
  initializeState( &state, &link, fname);
  updateSideLengths(link, &state);

  state.stop20 = stop20;
  
  state.refineUntil = refineUntil;
  state.movie = movie;
  
  state.overstepTol = overstepTol;
  state.minradOverstepTol = minradOverstepTol;
  state.minminrad = minradOverstepTol;
  
  state.eqMultiplier = eqMult;
  state.residualThreshold = residualThreshold;

  sprintf(state.finalfilename,"./%s.rr/%s.final.vect",
	  arg_infile->basename[0],arg_infile->basename[0]);
  sprintf(state.finalstrutname,"./%s.rr/%s.final.struts",
	  arg_infile->basename[0],arg_infile->basename[0]);
  sprintf(state.vectprefix,"./%s.rr/vectfiles/%s",
	  arg_infile->basename[0],arg_infile->basename[0]);
  sprintf(state.logfilename,"./%s.rr/%s.log",
	  arg_infile->basename[0],arg_infile->basename[0]);
  sprintf(state.fprefix,"./%s.rr/",arg_infile->basename[0]);

  /* We now open the local directory for storing output. */
  
  char cmdline[1024];
  
  sprintf(cmdline,"rm -fr %s.rr",arg_infile->basename[0]);
  system(cmdline);
  sprintf(cmdline,"mkdir %s.rr",arg_infile->basename[0]);

  if (system(cmdline) > 0) {

    printf("Created directory %s.rr to hold output of this run.\n",
	   arg_infile->basename[0]);

  } else {

    printf("WARNING. Creation of output directory %s.rr may have failed.\n",
	   arg_infile->basename[0]);

  }

  sprintf(cmdline,"mkdir %s.rr/logfiles",arg_infile->basename[0]);

  if (system(cmdline) <= 0) {

    printf("WARNING. Creation of logfile directory %s.rr/logfiles may have failed.\n",
	   arg_infile->basename[0]);

  }

  if (arg_suppressfiles->count == 0) {

    sprintf(cmdline,"mkdir %s.rr/vectfiles",arg_infile->basename[0]);

  }

  /* We now initialize the log with a lot of (hopefully) helpful information
     about the run. */

  gLogfile = fopen(state.logfilename,"w");
  if (gLogfile == NULL) { 
    fprintf(stderr,"Ridgerunner: Couldn't open logfile %s.\n",state.logfilename);
    exit(1);
  }

  fprintf(gLogfile,"Ridgerunner logfile.\n");

  #ifdef HAVE_ASCTIME
  #ifdef HAVE_LOCALTIME
  #ifdef HAVE_TIME
  
  time_t start_time;
  start_time = time(NULL);

  fprintf(gLogfile,"Run began: %s.\n",asctime(localtime(&start_time)));

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
	  arg_infile->filename[0]);

  fprintf(gLogfile,"Command line: ");

  for(i=0;i<argc;i++) {

    fprintf(gLogfile,"%s ",argv[i]);

  } 
  
  fprintf(gLogfile,"\n");
  fprintf(gLogfile,"------------------------------------------------\n");

  /* We now open the given link file. */

  linkFile = fopen_or_die(arg_infile->filename[0],"r", __FILE__ , __LINE__ );

  int plc_read_error_num;
  char plc_error_string[1024];
  size_t plc_error_size = {1024};

  link = plc_read(linkFile,&plc_read_error_num,plc_error_string,plc_error_size); 

  if (plc_read_error_num != 0) {

    sprintf(errmsg,"ridgerunner: file %s had \n"
	    "             plc_read error %s.\n",
	    arg_infile->filename[0],plc_error_string);
    FatalError(errmsg, __FILE__, __LINE__ );
    
  }

  strncpy(fname,arg_infile->basename[0],sizeof(fname));
  fclose(linkFile);

  printf("Loaded %d component, %d vertex plCurve from %s.\n",
	 link->nc,plc_num_verts(link),fname);
  fprintf(gLogfile,"Loaded %d component, %d vertex plCurve from %s.\n",
	 link->nc,plc_num_verts(link),fname);

  /* We archive a version of the initial input file in the run directory. */

  char filename[1024];
  FILE *savefile;

  sprintf(filename,"./%s.rr/%s.vect",arg_infile->basename[0],arg_infile->basename[0]);
  savefile = fopen_or_die(filename,"w", __FILE__, __LINE__ );
  plc_write(savefile,link);
  fclose(savefile);

  printf("Saved copy of %s to %s.\n",arg_infile->filename[0],filename);
  fprintf(gLogfile,"Saved copy of %s to %s.\n",arg_infile->filename[0],filename);

  /* Now perform initial operations. */

  plCurve *tempLink;

  if (arg_resolution->count > 0) {

    tempLink = plCurve_fixresolution(link,arg_resolution->dval[0]);
    plc_free(link);
    link = tempLink;

    printf("Splined to resolution of %g verts/rop. New curve has %d verts.\n",
	   arg_resolution->dval[0],plc_num_verts(link));
    fprintf(gLogfile,"Splined to resolution of %g verts/rop. New curve has %d verts.\n",
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
    fprintf(gLogfile,"Equilateralized curve has max edgelength/min edgelength %g.\n",
	   plCurve_long_edge(link)/plCurve_short_edge(link));

  }

  double thickness;
  thickness = octrope_thickness(link,NULL,0,gLambda);

  if( arg_autoscale->count > 0 || thickness < state.tube_radius + 0.001)) {

    printf("Curve has thickness %g. Scaling to thickness %g.\n",
	   thickness,state.tube_radius + 0.001);
    fprintf(gLogfile,"Curve has thickness %g. Scaling to thickness %g.\n",
	    thickness,state.tube_radius + 0.001);

    plc_scale(link,(state.tube_radius + 0.001)/thickness);
    thickness = octrope_thickness(link,NULL,0,gLambda);
    
    printf("Scaled curve has thickness %g.",thickness);
    fprintf(gLogfile,"Scaled curve has thickness %g.",thickness);
    
    if (fabs(thickness - (state.tube_radius + 0.001) > 1e-12) {

      sprintf(errmsg,"ridgerunner: Failed to scale %s to thickness 1.001."
	      "             Aborting run.\n",fname);
      FatalError(errmsg, __FILE__ , __LINE__ );

    } else {

      printf(" Autoscale selftest ok.\n");

    }

  }

  /* Now save the modified file so that we record what we actually ran with. */

  sprintf(filename,"./%s.rr/%s.atstart.vect",
	  arg_infile->basename[0],arg_infile->basename[0]);
  savefile = fopen_or_die(filename,"w", __FILE__ , __LINE__ );
  plc_write(savefile,link);
  fclose(savefile);

  printf("ridgerunner: Rerez'd, autoscaled, eq'd file written to %s.\n",filename);
  fprintf(gLogfile,"ridgerunner: Rerez'd, autoscaled, eq'd file written to %s.\n",
	  filename);
    
  /* We now initialize the link-dependent parts of "state". */

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

  state.minrad = octrope_minradval(link);
  state.thickness = octrope_thickness(link,NULL,0,gLambda);
  state.ropelength = plc_arclength(link,NULL)/state.thickness;
  state.residual = DBL_MAX;
  state.shortest = state.thickness+1.0;

  printf( "Overstep tolerance: %f of thickness %f (%f) minminrad: %3.8lf\n", 
	  state.overstepTol, state.tube_radius, 
	  state.tube_radius - state.overstepTol*state.tube_radius, 
	  state.minminrad );
  
  // use amd column ordering if the user hasn't specified something else
  setenv("COL_ORDERING", "amd", 0);
    
  /* Now we open the data logfiles in state. */

  char tmpfilename[1024];
  
  for(i=0;i<state.nlogs;i++) {

    sprintf(tmpfilename,"./%s.rr/logfiles/%s.dat",
	    arg_infile->basename[0],
	    state.logfilenames[i]);

    state.logfiles[i] = fopen_or_die(tmpfilename,"w", __FILE__ , __LINE__ );

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

  /* Now write the concluding file to the appropriate directory. */

  sprintf(tmpfilename,"./%s.rr/%s.final.vect",
	  arg_infile->basename[0],arg_infile->basename[0]);  
  savefile = fopen_or_die(tmpfilename,"w", __FILE__ , __LINE__ );
  plc_write(savefile,link);
  fclose(savefile);

  /* Now write the final strut set (in octrope format). */

  sprintf(tmpfilename,"./%s.rr/%s.final.struts",
	  arg_infile->basename[0],arg_infile->basename[0]);
  
  savefile = fopen_or_die(tmpfilename,"w", __FILE__ , __LINE__ );

  octrope_strutfile_write(state.lastStepStrutCount,
			  state.lastStepStruts,
			  state.lastStepMinradStrutCount,
			  state.lastStepMRlist,
			  savefile);

  fclose(savefile);

  /* Now clean up allocated memory and close files. */
  
  free_search_state(state);
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

  static char *log_fnames[] = {
    "length","ropelength","strutcount","stepsize","thickness","minrad","residual",
    "maxovermin","rcond","walltime","maxvertforce","csteps_to_converge","edgelenvariance"
  }; 

  // zero all those pointers
  bzero(state, sizeof(search_state));
  
  strncpy(state->fname, fname, sizeof(state->fname)-1);
  
  state->tube_radius = 0.5; // fix this for now
  
  state->maxStepSize = 0.1*plCurve_short_edge(*inLink);
  
  if( state->maxStepSize > 1e-3 )
    state->maxStepSize = 1e-3;
  //	state->maxStepSize = 1e-4;
  state->shortest = 2*state->tube_radius;
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
  
  state->nlogs = sizeof(log_fnames);
  for(i=0;i<state->nlogs;i++) {state->logfilenames[i] = log_fnames[i];}

}

void free_search_state(search_state *inState)

     /* Free memory and close filehandles. */

{

  int i;

  for(i=0;i<kTotalGraphTypes;i++) { fclose(inState->logfiles[i]); }

  free(inState->compOffsets);  
  free(inState->lastStepStruts);
  free(inState->lastStepMRlist);
  
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
