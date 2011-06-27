/*

   Display.c : code to update runtime displays as the program runs. 

*/

#include "ridgerunner.h"

FILE *gclpipe = NULL;

#define RUNSTART_TEXT_YX   0,52
#define RUNSTART_FIELD_YX  0,64

#define RUNCLOCK_TEXT_YX   1,52
#define RUNCLOCK_FIELD_YX  1,64

#define ITER_TEXT_YX       5,0
#define ITER_FIELD_YX      5,18

#define TSCALLS_TEXT_YX    6,2
#define TSCALLS_FIELD_YX   6,18

#define ORCALLS_TEXT_YX    7,2
#define ORCALLS_FIELD_YX   7,18

#define MAXITER_TEXT_YX    10,0
#define MAXITER_FIELD_YX   10,18

#define CSATT_TEXT_YX      9,0
#define CSATT_FIELD_YX     9,20

#define ATT_TEXT_YX        10,0
#define ATT_FIELD_YX       10,18

#define FILEPATH_TEXT_YX   20,0
#define FILEPATH_FIELD_YX  20,10

#define THI_TEXT_YX  5,40
#define THI_FIELD_YX 5,55

#define ROP_TEXT_YX  6,40
#define ROP_FIELD_YX 6,55

#define MR_TEXT_YX   7,40
#define MR_FIELD_YX  7,55

#define STRUTS_TEXT_YX 9,40
#define STRUTS_FIELD_YX 9,55

#define MRSTRUTS_TEXT_YX 10,40
#define MRSTRUTS_FIELD_YX 10,55

#define CONSTRAINT_TEXT_YX 11,0
#define CONSTRAINT_FIELD_YX 11,15

#define LAMBDA_TEXT_YX 12,0
#define LAMBDA_FIELD_YX 12,15

#define TUBERAD_TEXT_YX 13,0
#define TUBERAD_FIELD_YX 13,15

#define MAXMIN_TEXT_YX 12,40
#define MAXMIN_FIELD_YX 12,55

void parse_display_arg(search_state *inState, struct arg_str *display)
{
  inState->ndisplay = 3;
  inState->runtime_display[0] = kRopelength;
  inState->runtime_display[1] = kStrutCount;
  inState->runtime_display[2] = kThickness;
  
  int n;

  for(n=0;n<display->count;n++) {
    
    if (!strcmp(display->sval[n],"length")) {
      
      inState->runtime_display[inState->ndisplay++] = kLength;
      
    } else if (!strcmp(display->sval[n],"ropelength")) {
      
      inState->runtime_display[inState->ndisplay++] = kRopelength;
      
    } else if (!strcmp(display->sval[n],"strutcount")) {
      
      inState->runtime_display[inState->ndisplay++] = kStrutCount;
      
    } else if (!strcmp(display->sval[n],"stepsize")) {
      
      inState->runtime_display[inState->ndisplay++] = kStepSize;
      
    } else if (!strcmp(display->sval[n],"thickness")) {
      
      inState->runtime_display[inState->ndisplay++] = kThickness;
      
    } else if (!strcmp(display->sval[n],"minrad")) {
      
      inState->runtime_display[inState->ndisplay++] = kMinrad;
      
    } else if (!strcmp(display->sval[n],"residual")) {
      
      inState->runtime_display[inState->ndisplay++] = kResidual;
      
    } else if (!strcmp(display->sval[n],"maxovermin")) {
      
      inState->runtime_display[inState->ndisplay++] = kMaxOverMin;
      
    } else if (!strcmp(display->sval[n],"rcond")) {
      
      inState->runtime_display[inState->ndisplay++] = kRcond;
      
    } else if (!strcmp(display->sval[n],"walltime")) {
      
      inState->runtime_display[inState->ndisplay++] = kWallTime;
      
    } else if (!strcmp(display->sval[n],"maxvertforce")) {
      
      inState->runtime_display[inState->ndisplay++] = kMaxVertexForce;
      
    } else if (!strcmp(display->sval[n],"csteps_to_converge")) {
      
      inState->runtime_display[inState->ndisplay++] = kCorrectionStepsNeeded;
      
    } else if (!strcmp(display->sval[n],"edgelenvariance")) {
      
      inState->runtime_display[inState->ndisplay++] = kEQVariance;
      
    } else if (!strcmp(display->sval[n],"memused")) {
      
      inState->runtime_display[inState->ndisplay++] = kMemused;
      
    } else if (!strcmp(display->sval[n],"effectiveness")) {
      
      inState->runtime_display[inState->ndisplay++] = kEffectiveness;
      
    } else if (!strcmp(display->sval[n],"score")) {

      inState->runtime_display[inState->ndisplay++] = kScore;

    } else if (!strcmp(display->sval[n],"strutfreeresidual")) {

      inState->runtime_display[inState->ndisplay++] = kStrutFreeResidual;

    } else if (!strcmp(display->sval[n],"symmetryerror")) {

      inState->runtime_display[inState->ndisplay++] = kSymmetryError;

    }
    
  }

}

void init_runtime_display(search_state *inState)
{

#ifdef CURSES_DISPLAY

  char errmsg[1024];
  char plc_ver[1024],octrope_ver[1024],tsnnls_ver[1024];
  char disp[1024];

  plc_version(plc_ver,sizeof(plc_ver));
  octrope_version(octrope_ver,sizeof(octrope_ver));
  tsnnls_version(tsnnls_ver,sizeof(tsnnls_ver));

  initscr();
  sprintf(disp,"Ridgerunner %s (cvs build %s %s)",PACKAGE_VERSION, __DATE__ , __TIME__ );
  printw(disp);

  sprintf(disp,"plCurve %s",plc_ver);
  mvprintw(1,0,disp);

  sprintf(disp,"octrope %s",octrope_ver);
  mvprintw(2,0,disp);

  sprintf(disp,"tsnnls %s",tsnnls_ver);
  mvprintw(3,0,disp);

#ifdef HAVE_ASCTIME
#ifdef HAVE_LOCALTIME
#ifdef HAVE_TIME
  
  mvprintw(RUNSTART_FIELD_YX,"%s",asctime(localtime(&(inState->start_time))));

  mvprintw(RUNSTART_TEXT_YX,"runstart:");
  mvprintw(RUNCLOCK_TEXT_YX,"runclock:");

#endif
#endif
#endif

  mvprintw(ITER_TEXT_YX,   "Step number    :"); 
  mvprintw(TSCALLS_TEXT_YX,"tsnnls calls :");
  mvprintw(ORCALLS_TEXT_YX,"octrope calls:");

  mvprintw(CSATT_TEXT_YX,  "correction steps:");

  mvprintw(THI_TEXT_YX,       "thickness   :");
  mvprintw(STRUTS_TEXT_YX,    "struts      :");
  mvprintw(MRSTRUTS_TEXT_YX,  "mr struts   :");
  mvprintw(CONSTRAINT_TEXT_YX,"constraints :");

  attron(A_BOLD);
  mvprintw(ROP_TEXT_YX,       "ropelength  :");
  attroff(A_BOLD);

  mvprintw(MR_TEXT_YX,        "minrad      :");

  mvprintw(LAMBDA_TEXT_YX,    "lambda      :");
  mvprintw(TUBERAD_TEXT_YX,   "tube radius :");
  mvprintw(MAXMIN_TEXT_YX,    "max/min side:");


  mvprintw(15,0,"Ridgerunner log messages:");
  gLogwin = newwin(0, 0, 16, 1);
  
  if (gLogwin == NULL) {

    printf(errmsg,"ridgerunner: Couldn't open curses window gLogwin.\n");
    FatalError(errmsg, __FILE__ , __LINE__ );

  }
  
  wrefresh(gLogwin);
  idlok(gLogwin,TRUE);
  scrollok(gLogwin,TRUE);
  curs_set(0);  /* Set cursor to invisible during run */

  refresh();

#else 

  printf("ridgerunner: Curses screen display disabled, using stdout.\n" 
	 "             Define CURSES_DISPLAY to enable better screen display.\n");

#endif

}  
  
void update_runtime_display(plCurve *inLink,search_state *inState) 
{

#ifdef CURSES_DISPLAY

#ifdef HAVE_DIFFTIME
#ifdef HAVE_TIME

  time_t now;
  double rclock;
  int hrs,mins,secs;

  now = time(NULL);
  rclock = difftime(now,inState->start_time);

  /* rclock is the runclock in seconds. Convert to hrs:mins:secs. */

  hrs = floor(rclock/(60*60));
  rclock -= (60*60*hrs);

  mins = floor(rclock/(60));
  rclock -= (60*mins);

  secs = floor(rclock);

  mvprintw(RUNCLOCK_FIELD_YX,"%d:%d:%02d",hrs,mins,secs);

#endif
#endif
  
  mvprintw(ITER_FIELD_YX,"%d (%d)",inState->steps,inState->maxItrs);
  mvprintw(TSCALLS_FIELD_YX,"%d",inState->tsnnls_evaluations);
  mvprintw(ORCALLS_FIELD_YX,"%d",inState->octrope_calls);
  mvprintw(CSATT_FIELD_YX,"%d (%d)",inState->cstep_count,inState->last_cstep_attempts);

  mvprintw(THI_FIELD_YX,"%3.10g",inState->thickness);
  mvprintw(ROP_FIELD_YX,"%3.10g",inState->ropelength);
  mvprintw(MR_FIELD_YX,"%3.10g",inState->minrad);
  mvprintw(STRUTS_FIELD_YX,"%d",inState->lastStepStrutCount);
  mvprintw(MRSTRUTS_FIELD_YX,"%d",inState->lastStepMinradStrutCount);
 
  mvprintw(CONSTRAINT_FIELD_YX,"%d",plCurve_score_constraints(inLink));
  mvprintw(LAMBDA_FIELD_YX,"%g",gLambda);
  mvprintw(TUBERAD_FIELD_YX,"%g",inState->tube_radius);
  mvprintw(MAXMIN_FIELD_YX,"%.8g",inState->lastMaxMin);

  refresh();
  wrefresh(gLogwin);

#else 

  printf("%4d ",inState->steps);

  int i;

  for(i=0;i<inState->ndisplay;i++) {

    switch (inState->runtime_display[i]) {

    case kLength : printf("Len: %3.7f ",inState->length); break;
    case kRopelength : printf("Rop: %3.11f ",inState->ropelength); break;
    case kStrutCount : printf("Str: %3d MRstruts: %3d ",inState->lastStepStrutCount,inState->lastStepMinradStrutCount); break;
    case kStepSize : printf("Step size: %1.7f ",inState->stepSize); break;
    case kThickness : printf("Thi: %1.7f ",inState->thickness); break;
    case kMinrad : printf("Minrad: %1.7f ", inState->minrad); break;
    case kResidual : printf("Residual: %1.8f ",inState->residual); break;
    case kMaxOverMin : printf("Max/min edge: %1.7f ",inState->lastMaxMin); break;
    case kRcond : printf("Rcond: %g ",inState->rcond); break;
    case kEffectiveness : printf("Eff: %1.4f MrEff: %1.4f ",inState->lastStepPocaEffectiveness,inState->lastStepMREffectiveness); break;
    case kScore : printf("Sco: %3.11f ",inState->score); break;
    case kStrutFreeResidual : printf("Sfr: %1.8f ",inState->strutfreeresidual); break;
    case kSymmetryError : printf("SymErr: %1.8f ",plc_symmetry_group_check(inLink)); break;

    default: break;
      
    }

  }

  printf("\n");

#endif

}
  
void close_runtime_display() {

#ifdef CURSES_DISPLAY

  endwin();

#endif

}

  
 

void init_gv_display()
{	
}

void
shutdown_gv_display()
{
  fclose(gclpipe);
}

void refresh_gv_display(plCurve *L) 

     /* Procedure displays the curve on the linked copy of Geomview. */
     /* We note that OOGLPIPE must be running in Geomview, and we must */
     /* be able to write to /tmp. */

{

  FILE *curvepipe;

  /* Open the curve file and update it. */

  curvepipe = fopen("/tmp/curvepipe.oogl","w");
  plCurve_draw(curvepipe,L);
  fclose(curvepipe);

  /* Now tell Geomview about the situation. */

  fprintf(gclpipe,"(geometry curve { < /tmp/curvepipe.oogl})\n");
  fprintf(gclpipe,"(geometry struts { < /tmp/struts.vect})\n");
  fprintf(gclpipe,"(geometry dl { < /tmp/dl.vect})\n");
  fprintf(gclpipe,"(geometry dVdt { < /tmp/dVdt.vect})\n");
    
  fprintf(gclpipe,"(look g0)\n");

  fflush(gclpipe);

}

void strut_vectfile_write(plCurve *inLink, octrope_strut *strutlist, 
			  int strutCount, FILE *fp)

     /* Write struts to a VECT file for later display. */

{
  int i;
  double maxCompression;
  
  fprintf(fp, "VECT\n");
  fprintf(fp, "%d %d %d\n", strutCount, 2*strutCount, strutCount);
  fprintf(fp, "2");
  for( i=1; i<strutCount; i++ )
    fprintf(fp, " 2");
  fprintf(fp, "\n");
  
  fprintf(fp, "1");
  for( i=1; i<strutCount; i++ )
    fprintf(fp, " 1");
  fprintf(fp, "\n");
  
  maxCompression = 0;
  
  for( i=0; i<strutCount; i++ ) {
    
    maxCompression = (strutlist[i].compression > maxCompression) ? 
      strutlist[i].compression : maxCompression;
    
  }
  
  if( maxCompression == 0 )
    maxCompression = 1e-6;
  
  for( i=0; i<strutCount; i++ ) {

    plc_vector  points[2];
    
    octrope_strut_ends( inLink, &strutlist[i], points );
    
    fprintf(fp, "%lf %lf %lf\n", plc_M_clist(points[0]));
    fprintf(fp, "%lf %lf %lf\n", plc_M_clist(points[1]));
    
  }
  
  for( i=0; i<strutCount; i++ ) {

    fprintf(fp, "%f,%f,%f,0\n", strutlist[i].compression/maxCompression, 0.0, 0.0);
  }
  
}


plCurve *vectorfield_to_plCurve(plc_vector *vf, plCurve *inLink) 

     /* Create a plCurve representing vf. */

{
  int vItr, vfItr, cItr;
  
  plCurve *outlink;
  int *nv;
  bool *open;
  int *cc;
  int verts;

  /* First, a bit of error checking */

  fatalifnull_(vf);
  fatalifnull_(inLink);

  /* Now allocate the plCurve... */

  verts = plc_num_verts(inLink);

  nv = (int *)malloc_or_die(verts*sizeof(int), __FILE__ , __LINE__ ); 
  open = (bool *)malloc_or_die(verts*sizeof(bool), __FILE__ , __LINE__ );
  cc = (int *)malloc_or_die(verts*sizeof(int), __FILE__ , __LINE__ );

  for(cItr = 0;cItr < verts;cItr++) {

    nv[cItr] = 2;
    open[cItr] = true;
    cc[cItr] = 1;

  }

  outlink = plc_new(verts,nv,open,cc);
  fatalifnull_(outlink);

  /* Now we enter the appropriate locations in outlink. */

  for(vfItr=0,cItr=0;cItr < inLink->nc;cItr++) {

    for(vItr=0;vItr < inLink->cp[cItr].nv;vItr++,vfItr++) {

      outlink->cp[vfItr].vt[0] = inLink->cp[cItr].vt[vItr];
      outlink->cp[vfItr].vt[1] = plc_vect_sum(outlink->cp[vfItr].vt[0],
					      vf[vfItr]);

    }
  }

  /* We can free memory, then we're done. */

  free(nv); free(open); free(cc);

  return outlink;

}
      

