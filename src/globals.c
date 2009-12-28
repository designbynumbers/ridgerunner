/* This file contains definitions of the various global variables for rr */

#include"ridgerunner.h"

int VERBOSITY;
FILE *gLogfile;

int gSuppressOutput = 0;
int gQuiet = 0;
double gLambda = 1.0;   /* lambda-stiffness of rope */
int gMaxCorrectionAttempts = 25;
int gNoRcond = 0;
int gLsqrLogging = 1;
int gNoTimeWarp = 1;  /* By Default, we turn timewarp OFF */
int gEqIt = 0; /* By default, we never equilateralize */
int gAnimationStepper = 0; /* By default, we use the steepest descent stepper. */
int gConjugateGradient = 0;
int gTryNewton = 0;

int gNumTubeColors = 5;   /* Note: We must keep this in sync with the gTubeColors array below */

plc_color gTubeColors[5] = { \
  {0.9647,0.9098,0.7647},  /* cream */ \
  {0.8470,0.7019,0.3960},  /* ltbrown */ \
  {0.78039,0.9176,0.8980},  /* ltgreen */	\
  {0.00392,0.4,0.3686},    /* dkbrown */ \
  {0.5490,0.3176,0.0392},  /* dkgreen */ \

};

plc_color gStraightSegColor = 
    {128.0/255.0,177.0/255.0,211.0/255.0,1.0}; /* Blue */

plc_color gKinkColor = 
    {251.0/255.0,128.0/255.0,114.0/255.0,1.0}; /* Red */

plc_color gHelixColor = 
    {141.0/255.0,211.0/255.0,199.0/255.0,1.0}; /* Gray-greeny blue */

plc_color gStrutColor =                        /* Kind of a dark green */
  {35/255.0,132/255.0,67/255.0,1.0};
