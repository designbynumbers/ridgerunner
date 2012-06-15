/*******************************************************************

   csteptest.c : Test program for the correction stepper in ridgerunner.
                 Inputs a vectfile and attempts to correct it to thickness 0.5 
                 using the ridgerunner stepper. The results are saved to the 
                 file file.corrected.vect. 

	   ************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include"ridgerunner.h"
#include"mangle.h"
#include<argtable2.h>

/* Global variables live here. */

struct arg_dbl  *target_thickness;
struct arg_dbl  *arg_lambda;
struct arg_dbl  *arg_cstep_size;
struct arg_dbl  *arg_maxstep;
struct arg_dbl  *arg_overstep;
struct arg_dbl  *arg_mroverstep;
struct arg_int  *arg_maxcorr;


struct arg_lit  *verbose;
struct arg_file *arg_infile;
struct arg_lit  *help;

struct arg_lit  *quiet;

struct arg_end *end;
struct arg_end *helpend;

FILE    *infile_fptr,*outfile_fptr;

int    QUIET=0;

int main(int argc,char *argv[])
{
  int       infilenum,i,nerrors;
  plCurve   *link;
  bool      allpass = true;
  
  double    overstepTol=0.0001; /* Should this only be for tube_radius = 1? */
  
  void *argtable[] = 
    {

     target_thickness  = arg_dbl0("r","radius","<x>","target tube radius"),
     arg_lambda = arg_dbl0("l","lambda","<double>",
			   "minimum radius of curvature for unit rope"),
     verbose = arg_lit0("v","verbose","print debugging information"),
     arg_infile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     arg_cstep_size = arg_dbl0("k","CorrectionStepSize","<fraction>",
			       "initial size of Newton correction step"),
     arg_mroverstep   = arg_dbl0(NULL,"MinRadOverstepTol","<x>",
				 "start correction step if "
				 "minrad < lambda*(1 - x)"),
     arg_overstep = arg_dbl0("o","OverstepTol","<x>",
			     "start correction step if "
			     "thickness < tuberad*(1 - x)"),
     arg_maxstep = arg_dbl0(NULL,"MaxStep","<scalar>","maximum stepsize "
			    "for a length-reduction step"),
     arg_maxcorr = arg_int0(NULL,"MaxCorrectionAttempts","<n>",
			    "maximum # of Newton steps in error correction"),
     quiet = arg_lit0("q","quiet","suppress almost all output (for scripting)"), 
     help = arg_lit0(NULL,"help","display help message"),
     end = arg_end(20)};
  
  void *helptable[] = {help,helpend = arg_end(20)};
  
  /* First, we parse the command-line arguments using argtable. */

  fprintf(stderr,"csteptest (" PACKAGE_STRING ") compiled " __DATE__ " " __TIME__ "\n");

  if (arg_nullcheck(argtable) != 0)
    printf("csteptest: Insufficient memory to allocate argument table.\n");

  nerrors = arg_parse(argc,argv,argtable);

  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */
    
    nerrors = arg_parse(argc,argv,helptable);
    
    if (nerrors > 0) {  /* The help table didn't match either-- the
                         first set of errors was probably more
                         helpful, so we display it. */

      arg_print_errors(stdout,end,"csteptest");
      exit(1);

    } else {  /* The help table matched, which means we 
		 asked for help or gave nothing */
  
      printf("csteptest attempts to correct a Geomview VECT file to given thickness.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }

  QUIET = quiet->count > 0; /* Register if we're in batch mode */

  double correctionStepSize = 0.25;
  double minradOverstepTol = 0.00005; /* Used to be abs val of 0.499975 */

  if (arg_cstep_size->count > 0) 
    { correctionStepSize = arg_cstep_size->dval[0]; }

  if (arg_maxcorr->count > 0) 
    { gMaxCorrectionAttempts = arg_maxcorr->ival[0]; } // In globals.c
  
  if (arg_mroverstep->count > 0) {minradOverstepTol = arg_mroverstep->dval[0];}

  if (arg_overstep->count > 0) { overstepTol = arg_overstep->dval[0]; }

  /* We now need to fake an inState to pretend that we're in the middle
     of a ridgerunnner run. */

  search_state state;

  state.maxItrs = 999999;
  state.stop20  = -20;
  state.residualThreshold = -1;
  state.stopTime = 999999;

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

  state.steps = 0;
  state.maxStepSize = 0.01;

#ifdef HAVE_ASCTIME
#ifdef HAVE_LOCALTIME
#ifdef HAVE_TIME
  
  state.start_time = time(NULL);

#endif
#endif
#endif

  gLsqrLogging = false;

  printf("Filename                         Verts Initial Thi Iter Final Thi (Poca,MR)         Result\n");
  printf("------------------------------------------------------------------------------------------\n");

  for(infilenum = 0;infilenum < arg_infile->count;infilenum++) {

    printf("%-32s ",arg_infile->basename[infilenum]);

    /* We start by loading the tube from core. */

    FILE *infile_fptr;

    infile_fptr = fopen(arg_infile->filename[infilenum],"r");
  
    if (infile_fptr == NULL) {
      
      fprintf(stderr,"csteptest: Couldn't open file %s.\n",
	      arg_infile->filename[infilenum]);
       continue;  /* Try the next file */
       
    }

    int plr_error_num;
    char plr_error_str[1024];
    
    link = plc_read(infile_fptr,
		    &plr_error_num,plr_error_str,sizeof(plr_error_str));
    
    /* We now demonstrate the octrope library's error handling protocol: */
  
    if (plr_error_num > 0) {   /* This is the signal for an error. */
      
      fprintf(stderr,"csteptest: link reading error\n%s\n",plr_error_str);
      continue;  /* Try the next file */
      
    }
    
    fclose(infile_fptr);

    /* Report to user */

    printf("%5d %11g",plc_num_verts(link),octrope_thickness(link,NULL,0,gLambda)); 

    /* Now we set the state information as needed. */

    char *extptr;
    
    strncpy(state.basename,arg_infile->basename[0],sizeof(state.basename));
    extptr = strstr(state.basename,arg_infile->extension[0]);
    *extptr = 0; /* This chops the extension off of basename. */

    strncpy(state.fname,arg_infile->filename[0],sizeof(state.fname));
    sprintf(state.workingfilename,"./%s.cstest/%s.working.vect",
	    state.basename,state.basename);
    sprintf(state.workingstrutname,"./%s.cstest/%s.working.struts",
	    state.basename,state.basename);
    sprintf(state.vectprefix,"./%s.cstest/vectfiles/%s",
	    state.basename,state.basename);
    sprintf(state.snapprefix,"./%s.cstest/snapshots/%s",
	    state.basename,state.basename);
    sprintf(state.logfilename,"./%s.cstest/%s.log",
	    state.basename,state.basename);
    sprintf(state.fprefix,"./%s.cstest/",state.basename);

    /* We now open the local directory for storing output. */
    
    char cmdline[1024];
    
    sprintf(cmdline,"rm -fr %s.cstest",state.basename);
    system_or_die(cmdline, __FILE__ , __LINE__ );
    
    sprintf(cmdline,"mkdir %s.cstest",state.basename);
    system_or_die(cmdline, __FILE__ , __LINE__ );
    
    sprintf(cmdline,"mkdir %s.cstest/logfiles",state.basename);
    system_or_die(cmdline, __FILE__ , __LINE__ );
    
    //sprintf(cmdline,"mkdir %s.cstest/vectfiles",state.basename);
    //system_or_die(cmdline, __FILE__ , __LINE__ );
    
    sprintf(cmdline,"mkdir %s.cstest/snapshots",state.basename);
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
    
    state.dumpAxb = false;
    state.conserveLength = NULL; /* We aren't using this. */
    
    if (state.maxmovieframes % 2 == 1) { state.maxmovieframes--; } 
    /* Make sure this is even. */
    
    state.minrad = octrope_minradval(link);
    state.thickness = octrope_thickness(link,NULL,0,gLambda);
    state.ropelength = octrope_ropelength(link,NULL,0,gLambda);
    state.residual = 1.0;
    state.shortest = octrope_poca(link,NULL,0);
    
    state.totalVerts = plc_num_verts(link);
    
    /*printf("Overstep tolerance: %f of thickness %f (%f) minminrad: %3.8lf\n", 
      state.overstepTol, state.tube_radius, 
      state.tube_radius*(1 - state.overstepTol), 
      state.minminrad ); */
    
    
    /* We now initialize the log with a lot of (hopefully) helpful information
       about the run. */
    
    gLogfile = fopen(state.logfilename,"w");
    if (gLogfile == NULL) { 
      fprintf(stderr,"csteptest: Couldn't open logfile %s.\n",state.logfilename);
      exit(1);
    }
    
    fprintf(gLogfile,"csteptest logfile.\n");
    
    fprintf(gLogfile,"csteptest %s\n",PACKAGE_VERSION);
    
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
    
    /* We are now prepared to run the correction stepper. */

    bool pass = true;

    if (!correct_thickness(link,&state)) {

      pass = false;

    }

    /* Now we report the results. */

    double newThi,newPoca,newMR;
   
    newThi = octrope_thickness(link,NULL,0,gLambda);
    newPoca = octrope_poca(link,NULL,0);
    newMR = octrope_minradval(link)/gLambda;

    if (newPoca < state.tube_radius - state.overstepTol) {

      pass = false;
    
    }

    if (newMR < state.tube_radius - state.minradOverstepTol) {

      pass = false;

    }

    printf("%5d %8g(%7g,%7g)  ",state.last_cstep_attempts,newThi,newPoca,newMR);
    if (pass) { printf("pass\n"); } else { printf("FAIL\n"); allpass = false; }

    outfile_fptr = fmangle(state.basename,".vect",".corrected.vect");

    plc_write(outfile_fptr,link);
    plc_free(link);

  }

  printf("\nFinal Result: ");
  if (allpass) { printf("PASS\n"); exit(0); } else { printf("FAIL\n"); exit(1); }
  
}
    


