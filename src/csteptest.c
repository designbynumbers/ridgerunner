/*******************************************************************

   csteptest.c : Test program for the correction stepper in ridgerunner.
                 Inputs a vectfile and attempts to correct it to thickness 0.5 
                 using the ridgerunner stepper. The results are saved to the 
                 file file.corrected.vect. 

	   ************************************************************/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<ridgerunner.h>
#include<mangle.h>
#include<argtable2.h>

/* Global variables live here. */

struct arg_dbl  *target_thickness;
struct arg_dbl  *arg_lambda;
struct arg_dbl  *arg_cstep_size;
struct arg_dbl  *arg_maxstep;
struct arg_dbl  *arg_mroverstep;
struct arg_int  *arg_maxcorr;


struct arg_lit  *verbose;
struct arg_file *corefile;
struct arg_lit  *help;

struct arg_lit  *quiet;

struct arg_end *end;
struct arg_end *helpend;

plCurve *core;
FILE    *infile_fptr,*outfile_fptr;

int    QUIET=0;

int main(int argc,char *argv[])
{
  int       infilenum,i,nerrors;
  plCurve   *link;

  void *argtable[] = 
    {

     target_thickness  = arg_dbl0("r","radius","<x>","target tube radius"),
     arg_lambda = arg_dbl0("l","lambda","<double>",
			   "minimum radius of curvature for unit rope"),
     verbose = arg_lit0("v","verbose","print debugging information"),
     corefile  = arg_filen(NULL,NULL,"<file>",1,100000,"input files"),
     arg_cstep_size = arg_dbl0("k","CorrectionStepSize","<fraction>",
			       "initial size of Newton correction step"),
     arg_mroverstep   = arg_dbl0(NULL,"MinRadOverstepTol","<x>",
				 "start correction step if "
				 "minrad < lambda*(1 - x)"),
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
    printf("csteptestt: Insufficient memory to allocate argument table.\n");

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
  state.movie = movie;
  state.moviefactor = 1.0;

  state.tube_radius = 0.5;

  if (target_thickness->count > 0) 
    { state.tube_radius = target_thickness->dval[0]; }
  
  state.overstepTol = overstepTol;
  state.minradOverstepTol = minradOverstepTol;
  state.minminrad = gLambda*state.tube_radius*(1 - minradOverstepTol);

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

  int infilenum;

  for(infilenum = 0;infilenum < corefile->count;infilenum++) {

    /* We start by loading the tube from core. */

    FILE *infile_fptr;

    infile_fptr = fopen(corefile->filename[infilenum],"r");
  
    if (infile_fptr == NULL) {
      
      fprintf(stderr,"csteptest: Couldn't open file %s.\n",
	      corefile->filename[infilenum]);
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

    /* Now we set the state information as needed. */
    
    state.maxStepSize = 0.1*plCurve_short_edge(link);
    state.minrad = octrope_minradval(link);
    state.totalVerts = plc_num_verts(link);
    state.stepSize = 0.5*state.maxStepSize;
    
    if( ((double)1)/(double)state.totalVerts < state.stepSize )
      {
	state.stepSize = ((double)1/(double)state.totalVerts);
      }
    
    if( state.maxStepSize > 1e-3 ) { state.maxStepSize = 1e-3; }
    
    if( maxStep > 0 ) { state.maxStepSize = maxStep; }  
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
    state.residual = DBL_MAX;
    state.shortest = state.thickness+1.0;
    
    state.totalVerts = plc_num_verts(link);
    
    /*printf("Overstep tolerance: %f of thickness %f (%f) minminrad: %3.8lf\n", 
      state.overstepTol, state.tube_radius, 
      state.tube_radius*(1 - state.overstepTol), 
      state.minminrad ); */
    
    
    
