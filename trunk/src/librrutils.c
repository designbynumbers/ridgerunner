/*

librrutils : Contains utility functions for establishing a state and
making other global initializations required to run code from the
ridgerunner core as if it were "library" code. 

*/

#include<ridgerunner.h>

void initialize_rr_globals()

/* Initialize the globals for ridgerunner to default values. */

{
  gLambda = 1.0;
  gLsqrLogging = false;
  gNoTimeWarp = true;
}

search_state initialize_rr_state(plCurve *link)

/* Initializes search_state acceptably, given a link. We don't fill in filenames, because
   they probably depend on the particular debugging program being run. 
   
   We are guaranteed that every other value will be initialized, with an acceptable default 
   value. We assume that the globals are initialized with initialize_rr_globals.  */

{
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

  state.correctionStepDefault = 0.025;
  state.movie = 0;
  state.moviefactor = 1.0;
  state.tube_radius = 0.5;
  
  state.overstepTol = 0.0001;
  state.minradOverstepTol = 0.00005;
  state.minminrad = gLambda*state.tube_radius*(1 - state.minradOverstepTol);

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

  state.maxStepSize = 0.1*plCurve_short_edge(link);
  state.minrad = octrope_minradval(link);
  state.totalVerts = plc_num_verts(link);
  state.stepSize = 0.5*state.maxStepSize;
  
  if( ((double)1)/(double)state.totalVerts < state.stepSize )
    {
      state.stepSize = ((double)1/(double)state.totalVerts);
    }
  
  if( state.maxStepSize > 1e-3 ) { state.maxStepSize = 1e-3; }
  
  state.stepSize = state.stepSize < state.maxStepSize ? state.stepSize : state.maxStepSize;
  
  state.compOffsets = calloc(link->nc,sizeof(int));
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

#ifdef HAVE_TIME
    
    time_t start_time;
    state.start_time = start_time = time(NULL);
    
#else
    
}
