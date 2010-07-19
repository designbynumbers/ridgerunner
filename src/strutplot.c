/*
 * strutplot.c : This test program reads link information from a VECT file and
 *            outputs the strut information in various (hopefully useful)
 *            formats. 
 *
 *  $Id: strutplot.c,v 1.20 2006-04-18 19:14:37 ashted Exp $
 */

/* Copyright 2004 The University of Georgia. */

/* This file is part of ridgerunner. */

#include "portability.h"
#include "octrope.h"
#include "ridgerunner.h"
#include <argtable2.h>

#ifndef min 
  #define min(A,B) ((A < B) ? A : B)
#endif

/* Helpful data types for sorting and processing strutset data. */

typedef struct augmented_strut_type {
  const octrope_strut *strut;
  double s,t;
  double end;
} augmented_strut;

typedef struct tag_type {

  double plotx,ploty;
  char   text[100];
  char   code;      /* L -> left tag, R -> right tag, C -> ctr tag, S -> skipped tag */
  double tagx,tagy;

} tag;

#define pline_edges(P) (((P).open) ? (P).nv-1 : (P).nv)


/* These globals are used to set scale in the PostScript output of 
   strutset and curvature graph plots. */

int    NBOXES,NTICKS = {10};
double boxwidth,plotscale,tickplotstep;
double PLOTWIDTH = {5};  /* Measured in inches, the size of the entire plot */
int    NSTRUTCOLORS = {5};
double TEXTGRAY = {0}; /* Black */
double BOXMULT = {1.0};  /* Size of strut boxes (in edgelengths) */
double sbox[2] = {0,1}, tbox[2] = {0,1};  /* Clipping box for s,t coords. */
int    DATA_IN_COMMENTS = {0};

/* These globals are used in the key */

const char   *FILENAME;
double  shortest,max_curv,min_curv,curv_length;

/* These are the main colors for the plot. */

plc_color checkercolors[2] = {

   {1.02 * 0.9647,1.02 * 0.9098,1.02  * 0.7647},  //cream 
   // {0.78039,0.9176,0.8980}, //ltgreenblue 
   {0.6784,0.8,0.79215} // ltblue
  //{1.0,1.0,179.0/255.0,1.0},	                /* Cream */
  //{190.0/255.0,186.0/255.0,218.0/255.0,1.0}   /* Purple */
    
  };

plc_color sciencechecker[2] = {

  {79.0/255.0,83.0/255.0,105.0/255.0},
  {86.0/255.0,89.0/255.0,99.0/255.0}

};

plc_color strutcolors[5] = 
    
  /* These are from ColorBrewer 7 class sequentil YlGn (last 4 classes) */
    {
      {173.0/255.0,221.0/255.0,142.0/255.0,1.0},
      {120/255.0,198/255.0,121/255.0,1.0}, 
      {65/255.0,171/255.0,93/255.0,1.0},
      {35/255.0,132/255.0,67/255.0,1.0},
      {0/255.0,90/255.0,50/255.0,1.0} 

    };

const plc_color sciencestrutcolor = {216.0/255.0,222.0/255.0,45.0/255.0};

const plc_color compressioncolor = {12.0/255,44.0/255,132.0/255,1.0};

/* These globals are used to control aspects of the plot from the command line */

int SS_HIGHLIGHT = {true},
  KAPPAPLOT = {true},
  DRAWGRID  = {true},
  KINK_HIGHLIGHT = {true},
  EPSPLOT = {false},
  PLOT_COMPRESSIONS = {false},
  PLOT_KEY = {true},
  DEBUGLEVEL = {0},
  PLOT_BACKGROUND = {true},
  USE_DSC = {true};

double gLambda;

/**********************************************************************/
/*  Horizontal tag plotting (crude)                                   */
/**********************************************************************/

int compare_htags(const void *a, const void *b) 
{
  tag *A, *B;

  A = (tag *)(a);
  B = (tag *)(b);

  return A->plotx - B->plotx;
}

void display_htags(double tagwidth,int ntags,tag *tags,FILE *outfile)

/* Procedure cleverly displays as many of the horizontal tags
   in the "htags" array as it can without overlap. Tags that 
   are "squeezed out" are simply dropped, with a warning 
   printed to the user.
*/

{

  int i;

  qsort(tags,ntags,sizeof(tag),compare_htags); /* First, sort by x */

  for(i=0;i<ntags-1;i++) {

    if (tags[i+1].plotx - tags[i].plotx < tagwidth) { /* Too close */

      if (tags[i].code == 'R') { /* This is already squeezed */

	tags[i].code = 'S'; /* So drop it */

      } else { /* We can move this left in several cases*/

	if (i == 0) {

	  tags[i].code = 'L';

	} else if (tags[i].plotx - tags[i-1].plotx > 1.5*tagwidth) {
	  
	  tags[i].code = 'L';
	  
	} else {
	  
	  tags[i].code = 'S';
	  
	}

      }
      
      tags[i+1].code = 'R';
      
    }

  }

  /* We have now been through the list of tags and are prepared to write them */

  for(i=0;i<ntags;i++) {

    fprintf(outfile,"gsave %g %g translate %g setgray %s ",
	    tags[i].plotx,tags[i].ploty,TEXTGRAY,tags[i].text);

    switch (tags[i].code) {
      
    case 'L' :      
      fprintf(outfile,"dlt ");
      break;

    case 'R' :
      fprintf(outfile,"drt ");
      break;
      
    case 'C' :
      fprintf(outfile,"ht ");
      break;
      
    }
      
    fprintf(outfile,"grestore\n");
    
  }

}

/*******************************************************************/
/* Vertical tag placement (relative sophisticated gradient method) */
/*******************************************************************/  

int compare_vtags(const void *a, const void *b) 
{
  tag *A, *B;

  A = (tag *)(a);
  B = (tag *)(b);

  return A->ploty - B->ploty;

}

double tag_grad_norm(int ntags,double *tag_grad) 

{
  int i;
  double sum = 0;

  for(i=0;i<ntags;i++) { sum += tag_grad[i]*tag_grad[i]; }
  return sqrt(sum);
}

void compute_tag_grad(double taght, int ntags, tag *tags, double *tag_grad)

/* Computes an "adjustment vector" designed to keep tags from overlap
   while keeping them as close to their indicated points as possible. */

{
  int i;
  double POS = {0.5}, OVERLAP = {15};

  /* The variables "POS" and "OVERLAP" determine the relative strength 
     of the algorithm's preference for insisting on labels close to the 
     positions they represent (POS) and avoidance of overlaps (OVERLAP).

     You can tune these at will, but remember that setting OVERLAP too
     low may result in overlapping tags which will then be cut by the 
     layout algorithm later. */

  /* First compute "divergence from desired position gradient" */

  for (i=0;i<ntags;i++) {

    tag_grad[i] = -POS*(tags[i].tagy - tags[i].ploty)/taght;

  }

  /* Now compute "overlap gradient" */

  for (i=1;i<ntags;i++) {

    if (tags[i].tagy - tags[i-1].tagy < 1.2*taght) { 
      
      if (tags[i].tagy - tags[i-1].tagy > 0.5*taght) {

	tag_grad[i]   += OVERLAP * pow((tags[i].tagy - tags[i-1].tagy)/taght,-3.0) 
	  - OVERLAP*pow(1.2,-3.0);
	tag_grad[i-1] += OVERLAP * pow((tags[i-1].tagy - tags[i].tagy)/taght,-3.0)
	  - OVERLAP*pow(-1.2,-3.0);
  
      } else { /* We are really close, so change the function to avoid huge steps */
	
	tag_grad[i]   += OVERLAP * ( 8 - pow(1.2,-3.0) 
				     + (0.5 - (tags[i].tagy - tags[i-1].tagy)/taght) * 25);
      
	tag_grad[i-1] -= tag_grad[i];
	
      }

    }

  }

}

void display_vtags(double taght,int ntags,tag *tags,FILE *outfile)

/* Procedure displays vertical tags as best it can, keeping the tags
   as close to their points of origin as possible. */

{
  int    i,j;
  double *tag_grad;
  double alpha = {0.01};

  qsort(tags,ntags,sizeof(tag),compare_vtags); /* First, sort by y. */

  /* Now compute an initial gradient */

  tag_grad = calloc(ntags,sizeof(double));
  compute_tag_grad(taght,ntags,tags,tag_grad);

  printf("\tLaying out vertical tags with energy method...\n");
  printf("\t\tInitial gradient norm = %g.\n",tag_grad_norm(ntags,tag_grad));

  for(i=0;i<100000 && tag_grad_norm(ntags,tag_grad) > 0.001;i++) {

    for(j=0;j<ntags;j++) {
      
      tags[j].tagy += alpha*tag_grad[j];

    }
    
    compute_tag_grad(taght,ntags,tags,tag_grad);

  }

  printf("\t\t%d steps taken, final gradient norm = %g.\n",
	 i,tag_grad_norm(ntags,tag_grad));

  free(tag_grad);

  /* Now check for overlaps */

  for(i=1;i<ntags;i++) {

    if (tags[i].tagy - tags[i-1].tagy < taght) {

      printf("\t Warning! Tags %d and %d conflict (and were deleted).\n",
	     i,i-1);

      tags[i].code = tags[i-1].code = 'S';

    }

  }

  /* Now actually display the tags */

  fprintf(outfile,"tls\n");

  for (i=0;i<ntags;i++) {

    if (tags[i].code != 'S') {

      fprintf(outfile,"%s %g %g %g %g rttag\n",
	      tags[i].text,
	      tags[i].tagx,tags[i].tagy,
	      tags[i].plotx,tags[i].ploty);

    }
    
  }

}


int compare_augmented_struts(const void *a, const void *b)
{
  augmented_strut *A,*B;

  A = (augmented_strut *)(a);
  B = (augmented_strut *)(b);

  if (A->s < B->s) {
    
    return -1; 

  } else {

    return 1;

  }

}

int compare_mrlocs(const void *a, const void *b)
{
  octrope_mrloc *A,*B;

  A = (octrope_mrloc *)(a);
  B = (octrope_mrloc *)(b);

  if (A->component != B->component) {

    return A->component - B->component;

  } 

  if (A->vert != B->vert) {
    
    return A->vert - B->vert; 

  } 

  return 0;

}

void set_ps_color(plc_color color,FILE *outfile)

/* Procedure changes colors in PostScript */

{
  fprintf(outfile,"%g %g %g setrgbcolor\n",color.r,color.g,color.b);
}
 

 


void create_st_plot(plCurve *L, 
		    int n_struts, 
		    const octrope_strut *strutlist, 
		    int num_min_rad_locs,
		    octrope_mrloc *min_rad_locs,
                    int whole,
                    int pscomments,
		    FILE *outfile)

    /* Create a very informative plot showing curvature, struts, 
       and regions of interest such as straight segments for L 
       in arclength form. */

{ 
  double tot_length,edgelength;
  int cmp,nv,vert,strut,i,j,k;
  int end;
  plc_strand *cp;
  const octrope_strut *st;
  
  plc_vector temp_vect;
  augmented_strut *as_list;
  
  double **s_positions;  
  double avg_edgelength;

  double *compressions;
 
  int slo, shi;
  double maxskip, s_skip=0.0;
  double straightlength=0.0;
  double tickplotstep,tick_s_step,tickpos;

  int intpos;
  double curv,curvpixel;
  int printedpt,pt_ok,printed_total;

  double yi,yj,mij,predicted_curvk;
  double cmpwidth,cmpstart,cmpheight;

  tag *kinktags;
  int   nkinks = 0;

  tag  vtags[1000];
  int  nvtags = 0;

  int  struts_output = 0;
  
  double maxcompression,mincompression,comp_level[100],comp;
  int  *comp_bins[100],binsize[100];
  
  if (L == NULL) {
    return;
  }
 
  /******************************************************************/
  /*                                                                */
  /* 1. Link preprocessing.                                         */
  /*                                                                */
  /******************************************************************/
  
  NBOXES = plc_num_edges(L);
  
  if (DEBUGLEVEL > 0) {
    printf("Building s,t plot for %d edge link with %d struts...\n\n",
	   NBOXES,n_struts);  
  }
  
  /* First, compute length of link and s positions of verts. */  
  
  tot_length = 0;
  
  s_positions = (double **)(calloc(L->nc,sizeof(double *)));
  for (cmp = 0; cmp < L->nc; cmp++) {
    
    s_positions[cmp] = calloc(L->cp[cmp].nv+1,sizeof(double));   
    cp = &L->cp[cmp];
    nv = (cp->open) ? cp->nv-1 : cp->nv;
    for (vert = 0; vert < nv; vert++) {
      s_positions[cmp][vert] = tot_length;  /* record s position of vertex */ 
      
      temp_vect = cp->vt[vert+1];
      plc_M_sub_vect(temp_vect,cp->vt[vert]);
      tot_length += plc_M_norm(temp_vect);
      
    }

    s_positions[cmp][nv] = tot_length;
    s_positions[cmp][cp->nv] = tot_length; /* for open components */

  }
  
  boxwidth = tot_length/(double)(NBOXES); // Changing to 1/100 from (double)(NBOXES)
  avg_edgelength = tot_length/(double)(plc_num_edges(L));
  

  if (DEBUGLEVEL > 0) {
    
    printf("Found arclength positions for vertices...\n");

  }
  
  /* Now compute augmented strut list. We list every strut twice,
     switching endpoints to ensure */
  /* that we have both s and t coordinates in the augmented list. */
  
  as_list = malloc(2 *n_struts * sizeof(augmented_strut));
  
  for (strut=0;strut < n_struts;strut++) {
    
    st = &(strutlist[strut]);
    
    for (end=0;end<2;end++) {
      
      temp_vect = L->cp[st->component[end]].vt[st->lead_vert[end]+1];
      plc_M_sub_vect(temp_vect,
          L->cp[st->component[end]].vt[st->lead_vert[end]]);
      edgelength = plc_M_norm(temp_vect);
      
      as_list[2*strut + end].s = 
	s_positions[st->component[end]][st->lead_vert[end]] + st->position[end] * edgelength;
      
      as_list[2*strut + end].end = end;
      as_list[2*strut + end].strut = st;
      
    }      
    
    as_list[2*strut].t = as_list[2*strut + 1].s;  /* Record t values */
    as_list[2*strut + 1].t = as_list[2*strut].s;
    
  }
  
  qsort(as_list,2*n_struts,sizeof(augmented_strut),compare_augmented_struts);
  
  if (DEBUGLEVEL > 0) { printf("Built augmented strut list...\n"); }


  /*************************************************************************/
  /*                                                                       */
  /* 2. PostScript world initialization and plot setup                     */
  /*                                                                       */
  /*************************************************************************/

  /* We are now ready to plot some boxes (!) and to identify
     potentially interesting regions of the graph
     automagically. First, we set up the graphics state for box
     plotting and have at it.*/
  
  tickplotstep = NBOXES*10.0/(double)(NTICKS);
  tick_s_step  = tot_length/(double)(NTICKS);
  
  if (pscomments) {
    fprintf(outfile,"%% Setup the graphic layout.\n");
  }
  fprintf(outfile,"%f %f translate \n",72*(3.5/2.0),72*3.0);    
  /* Center a 5"x5" box on the page */   
  
  plotscale = PLOTWIDTH*72.0/(double)(10*NBOXES);
  /* Now one (PS) unit is 1/10 box */
  
  fprintf(outfile,"%f %f scale \n",plotscale,plotscale);  

  fprintf(outfile,"/TF /Helvetica findfont %g scalefont def\n",9.0/plotscale);
  fprintf(outfile,"TF setfont\n");
  
  /* We now draw the outline box and grid. */
  
  if (pscomments) {
    fprintf(outfile,"%% Draw the outline box and grid.\n");
  }
  fprintf(outfile,"/tls { %g setgray %g setlinewidth } def\n",TEXTGRAY,1.0/plotscale);
  fprintf(outfile,"tls\n");

  fprintf(outfile,"/ct { %g %g moveto dup stringwidth pop neg 0 rmoveto "
	  " dup stringwidth %g add exch %g add exch currentpoint "
	  "4 2 roll 1.0 setgray rectfill %g setgray show"
	  " newpath 0 0 moveto %g %g lineto stroke } def \n",
	  -27.0/plotscale,((whole) ? -2.0 : 25.0)/plotscale,
          2.0/plotscale,2.0/plotscale,TEXTGRAY,
	  -25.0/plotscale,((whole) ? 0.0 : 25.0)/plotscale); 

  if (PLOT_BACKGROUND) {
    
    for(i=0;i<L->nc;i++) {
      
      cmpstart = (10.0/boxwidth)*s_positions[i][0];
      cmpwidth = (10.0/boxwidth)*s_positions[i][L->cp[i].nv] - cmpstart;

      set_ps_color(checkercolors[0],outfile); 
      if (!whole) {
        /* Draw the triangle of struts from this component to itself */    
        
        fprintf(outfile,"newpath %g dup moveto %g %g lineto %g dup lineto closepath fill\n",
	        cmpstart,
	        cmpstart+cmpwidth,
	        cmpstart,
	        cmpstart+cmpwidth);
      }
        
      /* Fill in this rectangle of the curvature plot */
      
      if (KAPPAPLOT) {
	
	if (i%2) { set_ps_color(checkercolors[1],outfile); }
	
	fprintf(outfile,"%g %g %g %g rectfill\n",
		cmpstart,-2.5*tickplotstep,cmpwidth,2*tickplotstep);
	
      }
      
      /* And finally add color-coding bars to the outside of the graph */
      
      set_ps_color(gTubeColors[i%5],outfile);
      
      if (!whole) {
        fprintf(outfile,                                   /* Main diagonal */
	        "newpath %g dup moveto %g dup lineto "
	        "%g %g lineto %g %g lineto closepath fill\n",
	        cmpstart,cmpstart+cmpwidth,
	        cmpstart+cmpwidth-4.0/(sqrt(2)*plotscale),
	        cmpstart+cmpwidth+4.0/(sqrt(2)*plotscale),
	        cmpstart-4.0/(sqrt(2)*plotscale),
	        cmpstart+4.0/(sqrt(2)*plotscale));
      }
      
      if (KAPPAPLOT) {
	
	fprintf(outfile, 		/* Bottom, under curvature plot */
		"%g %g %g %g rectfill\n",
		cmpstart,-2.5*tickplotstep-4.0/plotscale,
		cmpwidth,4.0/plotscale);
	
      } else { 

	fprintf(outfile, 		/* Bottom, under s,t plot */
		"%g %g %g %g rectfill\n",
		cmpstart,-4.0/plotscale,
		cmpwidth,4.0/plotscale);
	
      }

      if (whole) {
      cmpwidth = (10.0/boxwidth)*s_positions[i][L->cp[i].nv] - cmpstart;
        fprintf(outfile,                /* Top of the s,t plane */
		"%g %g %g %g rectfill\n",
		cmpstart,10.0*NBOXES,
		cmpwidth,4.0/plotscale);
        fprintf(outfile,		/* Left-hand side */
	        "%g %g %g %g rectfill\n",
	        -4.0/plotscale,cmpstart,4.0/plotscale,cmpwidth);
      
      }
      
      fprintf(outfile,		/* Right-hand side */
	      "%g %g %g %g rectfill\n",
	      10.0*NBOXES,cmpstart,4.0/plotscale,cmpwidth);
      
      /* Now we loop over the other components, filling in rectangles */
      
      for(j=0;j<((whole) ? L->nc : i);j++) {
	
	if ((i-j)%2) { set_ps_color(checkercolors[1],outfile); } 
	else { set_ps_color(checkercolors[0],outfile); }
	
	cmpheight = (10.0/boxwidth)*(s_positions[j][L->cp[j].nv] - 
                                     s_positions[j][0]);
	fprintf(outfile,"%g %g %g %g rectfill\n",
		cmpstart,(10.0/boxwidth)*(s_positions[j][0]),
		cmpwidth,cmpheight);
	
      }
      
      /* Last, label the component breaks with arclength values */
      
      if (i < L->nc-1) {  /* The last label would be ropelength, and that's elsewhere */
	
        if (!whole) {
	  fprintf(outfile,
            "gsave %g %g translate (%1.2f) ct grestore\n",
             cmpstart+cmpwidth,cmpstart+cmpwidth,s_positions[i][L->cp[i].nv]);
        } else {
          fprintf(outfile,
            "gsave %g %g translate (%1.2f) ct grestore\n",
             -4.0/plotscale,cmpstart+cmpwidth,s_positions[i][L->cp[i].nv]);
        }
	
      }
      
      /* We're done with the checkerboard grid on the strut plot now. */
      
    }
    
    /* The first step is to setup a procedure for drawing "label-plus-tick" */
    
    if (pscomments) {
      fprintf(outfile,"%% Now add ticks & labels.\n");
    }
    fprintf(outfile,"%g setgray\n",TEXTGRAY);
    
    /* This procedure draws labels on opaque centered white rectangles */
    
    fprintf(outfile,"/dt { %g %g moveto dup stringwidth pop neg 0 rmoveto "
	    " dup stringwidth %g add exch %g add exch "
	    " currentpoint %g sub exch %g sub exch "
	    "4 2 roll 1.0 setgray rectfill %g setgray show"
	    " newpath 0 0 moveto %g %g lineto stroke } def \n",
	    -6.0/plotscale,((whole) ? -2.0 : 4.0)/plotscale,
	    9.5/plotscale,5.0/plotscale,
	    2/plotscale,2.5/plotscale,TEXTGRAY,
	    -4.0/plotscale,((whole) ? 0.0 : 4.0)/plotscale); 
  
    /* Now drop an array of labels. */
    
    fprintf(outfile,"[ ");
    
    for(i=0,tickpos=tot_length;i<=NTICKS;i++,tickpos -= tick_s_step) {
      
      fprintf(outfile,"(%1.2f) ",(tickpos > 0 ? tickpos:0.0));
      
    }  /* Make sure that each tick position is positive. */
       /* above by printing (tickpos > 0 ? tickpos:0.0)  */
       /* Otherwise, we can print -0.00 as the start     */
       /* instead of 0.00. */
    
    fprintf(outfile," ]\n");
    
    /* Now display them. */
    
    fprintf(outfile,"aload pop \n");
    fprintf(outfile,
	    "gsave TF setfont \n"
	  "0 1 %d { pop dt %g %g translate 0 0 moveto } for grestore \n",
          NTICKS,((whole) ? 0 : tickplotstep),tickplotstep);
    
    if (KAPPAPLOT) {
      
      /* We now add tags to the curvature plot. */

      if (pscomments) {
        fprintf(outfile,"%% Add labels to the curvature plot.\n");
      }
      
      fprintf(outfile,"/lt { %g 0 moveto dup stringwidth pop neg %g rmoveto show"
	      " newpath 0 0 moveto %g 0 lineto stroke } def \n",
	      -6/plotscale,-2/plotscale,-4/plotscale); 
      
      fprintf(outfile,"gsave\n");
      fprintf(outfile,"tls TF setfont\n");
      fprintf(outfile,"0 0 moveto 0 %g translate\n",-2.5*tickplotstep);
      
      fprintf(outfile,"(2) (1) (0) 0 1 2 { pop lt 0 %g translate } for\n",
	      1.0*tickplotstep);
      
      fprintf(outfile,"grestore\n");

    }
    
    fprintf(outfile,"/rttag { " /* The "righty tag": (s) tagx tagy tickx ticky */
	    
	    "moveto "		/* (s) tagx tagy  */
	    "dup "                /* (s) tagx tagy tagy */
	    "3 2 roll dup "	/* (s) tagy tagy tagx tagx */
	    "4 1 roll exch "	/* (s) tagx tagy tagx tagy */
	    "lineto stroke "      
	    "moveto "		/* (s) */
	    "%g %g rmoveto "
	    "show } def \n",1.0/plotscale,-3/plotscale); 

    /* lshape takes 3 arguments, slo shi end and draws an "L-shape" */
    
    fprintf(outfile, "/lshp { "
	    "exch    "		/* l e h */          /* Create many copies of args */
	    "dup dup "		/* l e h h h */
	    "5 2 roll "		/* h h l e h */
	    "3 -1 roll "		/* h h e h l */
	    "dup "		/* h h e h l l */
	    "4 1 roll "           /* h h l e h l */
	    "3 -1 roll "          /* h h l h l e */
	    "dup "                /* h h l h l e e */
	    "4 1 roll "           /* h h l e h l e */
	    "exch "               /* h h l e h e l */
	    "dup "                /* h h l e h e l l */
	    "4 2 roll "           /* h h l e l l h e */
	    "3 1 roll "           /* h h l e l e l h */
	    "dup dup "            /* h h l e l e l h h h */
	    "4 3 roll "           /* h h l e l e h h h l */
	    "dup dup "            /* h h l e l e h h h l l l */
	    
	    "newpath "
	    "0 moveto "                  /* l 0 moveto */
	    "0 1 4 { pop lineto } for "
	    /* l l lineto */
	    /* h h lineto */
	    /* e h lineto */
	    /* e l lineto */
	    /* h l lineto */
	    
	    "0 lineto "          	/* 0 h lineto */ /* empty stack */
	    "closepath fill } def\n");

  }


  /************************************************************************/
  /*                                                                      */
  /*  3. Highlight (potential) straight segments                          */
  /*                                                                      */
  /************************************************************************/
   
  /* Our method is to search the list of struts for suspiciously large 
     gaps between s values. If such a gap is larger than 3*avg_edgelength,
     we will become suspicious that a straight segment is indicated. */
  

  if (SS_HIGHLIGHT) {
    
    if (pscomments) {
      fprintf(outfile,"%% Highlight straight segments.\n");
    }
    fprintf(outfile,"gsave \n");
    
    set_ps_color(gStraightSegColor,outfile);
    
    maxskip = 0; straightlength = 0; 
    
    for(i=1;i<2*n_struts;i++) {
      
      s_skip = as_list[i].s - as_list[i-1].s;
      
      maxskip = (s_skip > maxskip) ? s_skip : maxskip;
      
      if (s_skip > 3*avg_edgelength) {
	
	straightlength += s_skip;
	
	slo = (int)(10*as_list[i-1].s/boxwidth) + 5;
	shi = (int)(10*as_list[i].s/boxwidth) - 5;
	
	fprintf(outfile,"%d %d %d lshp\n",slo,shi,10*NBOXES);
	

	if (KAPPAPLOT) {
	  
	  fprintf(outfile,"%d %g %d %g rectfill\n",
		  slo,-2.5*tickplotstep,shi-slo,2*tickplotstep);

	}
	  
	/* Now add "smart" tags on the rhs. */
	
	vtags[nvtags].plotx = 10*NBOXES;
	vtags[nvtags].ploty = shi;
	
	sprintf(vtags[nvtags].text,"(%1.2f)",as_list[i].s);
	
	vtags[nvtags].tagx  = 10*NBOXES + 10.0/plotscale;
	vtags[nvtags].tagy  = shi;
	vtags[nvtags].code  = 'C';
	
	nvtags++;
	
	vtags[nvtags].plotx = 10*NBOXES;
	vtags[nvtags].ploty = slo;
	
	sprintf(vtags[nvtags].text,"(%1.2f)",as_list[i-1].s);
	
	vtags[nvtags].tagx = 10*NBOXES + 10.0/plotscale;
	vtags[nvtags].tagy = slo;
	vtags[nvtags].code = 'C';
	
	nvtags++;
	
      }
      
    }
    
    fprintf(outfile,"grestore\n");
    
    if (straightlength > 0) {
      
      printf("\tFound a total of %g units of straight segment\n",
	     straightlength);
      
    } else {
      
      printf("\tNo straight segs found. Max strut s value gap: %g\n"
	     "\t                        3*avg_edgelen        : %g\n",
	     s_skip,3*avg_edgelength);
      
    }

  }

  /*****************************************************************************/
  /*                                                                           */
  /*  4. Highlight min_rad_locs as potential kinks.                            */
  /*                                                                           */
  /*****************************************************************************/

  if (num_min_rad_locs != 0 && KINK_HIGHLIGHT) {  

    /* There may _be_ no kinks, so we might just skip */

    kinktags = calloc(num_min_rad_locs,sizeof(tag));
    nkinks = 0;

    if (pscomments) {
      fprintf(outfile,"%% Highlight kinks.\n");
    }
    fprintf(outfile,"/ht { 0 0 moveto 0 %g lineto stroke 0 %g moveto " /* Draw tick */
	    " dup stringwidth pop neg 2 div 0 rmoveto show "
	    " } def\n",-8.0/plotscale,-15.0/plotscale);

    fprintf(outfile,"/dlt { 0 0 moveto 0 %g rlineto %g 0 rlineto stroke %g %g moveto "
	    " dup stringwidth pop neg 0 rmoveto show } def \n",
	    -11.5/plotscale,-2/plotscale,-3/plotscale,-15.0/plotscale);

    fprintf(outfile,"/drt { 0 0 moveto 0 %g rlineto %g 0 rlineto stroke %g %g moveto "
	    " show } def \n",
	    -11.5/plotscale,2/plotscale,3/plotscale,-15.0/plotscale);
    
    fprintf(outfile,"gsave\n");

    set_ps_color(gKinkColor,outfile);
    /* fprintf(outfile,"1 0.2 0.2 setrgbcolor\n"); */ /* Kinks are highlighted in red. */ 
    
    qsort(min_rad_locs,num_min_rad_locs,sizeof(octrope_mrloc),compare_mrlocs);

    for(i=0;i<num_min_rad_locs;i++) {

      /* i is now at the start of a kink arc (run of consecutive mr locs) */

      for(j=i+1;
	  j<num_min_rad_locs && 
	    min_rad_locs[j].component == min_rad_locs[i].component && 
	    (min_rad_locs[j].vert == min_rad_locs[j-1].vert + 1 || min_rad_locs[j].vert == min_rad_locs[j-1].vert);
	  j++);

      j--;

      /* j is now at the end of the kink arc */

      slo=(int)(10*s_positions[min_rad_locs[i].component][min_rad_locs[i].vert]/boxwidth)-5;
      shi=(int)(10*s_positions[min_rad_locs[j].component][min_rad_locs[j].vert]/boxwidth)+5;

      /* Note that we fatten kinks a little to display one-vertex kinks with 
	 something like a responsible error bar. */

      fprintf(outfile,"%d %d %d lshp\n", /* Highlight on the main plot */
	      slo,shi,10*NBOXES);

      if (KAPPAPLOT) {

	/* Highlight kinks on the lower (curvature) plot, too. */

	fprintf(outfile,"%d %g %d %g rectfill\n", /* Highlight on the curvature plot */
		slo,-2.5*tickplotstep,shi-slo,2*tickplotstep);

      }

      /* Now display "smart" tags for the kinks on the right-hand side. */

      vtags[nvtags].plotx = NBOXES*10;
      vtags[nvtags].ploty = (shi + slo)/2.0;
      
      sprintf(vtags[nvtags].text,"((%1.3f,%1.3f))",
	      (s_positions[min_rad_locs[i].component][min_rad_locs[i].vert]), 
	       s_positions[min_rad_locs[j].component][min_rad_locs[j].vert]);

      vtags[nvtags].tagx = NBOXES*10 + 10.0/plotscale;
      vtags[nvtags].tagy = (shi + slo)/2.0;
      vtags[nvtags].code = 'C';

      nvtags++; nkinks++; i = j; /* We advance i to start looking for the next kinked region */

    }

    /* if (KAPPAPLOT) { display_htags(30/plotscale,nkinks,kinktags,outfile); } */

    printf("\tFound %d kinked regions.\n",nkinks);

  }

  /*************************************************************************/
  /*                                                                       */
  /*  5. Add a grid to the plot                                            */
  /*                                                                       */
  /*************************************************************************/

  if (DRAWGRID && PLOT_BACKGROUND) {	/* If there are no ticks, drop the grid. */
    
    if (pscomments) {
      fprintf(outfile,"%% Add a grid.\n");
    }
    fprintf(outfile,
	    "gsave "
	    "%g %g %g setrgbcolor "
	    "%g  setlinewidth \n",
	    0.8*checkercolors[0].r,0.8*checkercolors[0].g,0.8*checkercolors[0].b,
	    0.75/plotscale);
    
    for(i=1;i<NTICKS;i++) {
      
      fprintf(outfile,"newpath %g 0 moveto %g %g lineto %d %g lineto stroke\n",
	      i*tickplotstep,i*tickplotstep,i*tickplotstep,
	      NBOXES*10,i*tickplotstep);
      if (whole) {
        fprintf(outfile,
          "newpath 0 %g moveto %g %g lineto %g %d lineto stroke\n",
	  i*tickplotstep,i*tickplotstep,i*tickplotstep,
	  i*tickplotstep,NBOXES*10);
      }
    }
    
    fprintf(outfile,"grestore\n");
    
   
    if (KAPPAPLOT) {

      /* Including extending vertical tick marks onto the curvature plot below */
      
      fprintf(outfile,
	      "gsave "
	      "%g  setlinewidth \n",0.75/plotscale);
      
      for(i=1;i<NTICKS;i++) {
	
	fprintf(outfile,
		"%g %g %g  setrgbcolor newpath %g 0  moveto 0 %g rlineto stroke\n"
		"%g %g %g setrgbcolor newpath %g %g moveto 0 %g rlineto stroke\n",
		0.8*checkercolors[0].r,0.8*checkercolors[0].g,0.8*checkercolors[0].b,		
		i*tickplotstep,-0.5*tickplotstep,
		0.8*checkercolors[0].r,0.8*checkercolors[0].g,0.8*checkercolors[0].b,
		i*tickplotstep,-0.5*tickplotstep,-2*tickplotstep);
	
      }
      
      fprintf(outfile,"0.85 setgray newpath 0 %g moveto %d 0 rlineto stroke\n",
	      -1.5*tickplotstep,10*NBOXES);
      
      fprintf(outfile,"grestore\n");

    }

  }

  /******************************************************************/
  /*                                                                */
  /*  6. Plot the struts.                                           */
  /*                                                                */
  /******************************************************************/

  if (n_struts > 0 && !PLOT_COMPRESSIONS) {

    printf("\tPlotting %d struts (without compressions)...\n",n_struts);
    if (pscomments) {
      fprintf(outfile,"%% Plot the struts.\n");
    }

    set_ps_color(strutcolors[3],outfile);
    
    for(i=0;i<2*n_struts;i++) {
      
      if (as_list[i].t < as_list[i].s || /* Display lower triangle. */
          whole) {                       /* Display the whole plane. */

	if (as_list[i].s/tot_length >= sbox[0] &&
	    as_list[i].s/tot_length <= sbox[1] &&
	    as_list[i].t/tot_length >= tbox[0] &&
	    as_list[i].t/tot_length <= tbox[1]) {
	
	  fprintf(outfile, " %d %d ",
		  (int)(10*as_list[i].s/boxwidth) - (int)(5*BOXMULT),
		  (int)(10*as_list[i].t/boxwidth) - (int)(5*BOXMULT));
	  
	  struts_output++;
	  
	  if (!(struts_output%8)) { fprintf(outfile,"\n"); }

	}
	
      }
      
    }
    
    if (struts_output > 0) {

      fprintf(outfile,"0 1 %d { pop %g %g rectfill } for \n",
	      struts_output-1,10*BOXMULT,10*BOXMULT);

    } else {

      fprintf(stderr,"strutplot: Warning! No struts plotted in [%g,%g] x [%g,%g].\n",
	      sbox[0],sbox[1],tbox[0],tbox[1]);

    }
      
  } else if (n_struts > 0) {

    /* We are now plotting struts _with_ compressions. */

    mincompression = DBL_MAX; maxcompression = -1.0;

    for(i=0;i<n_struts;i++) {

      comp = strutlist[i].compression;
      mincompression = (comp < mincompression) ? comp : mincompression;
      maxcompression = (comp > maxcompression) ? comp : maxcompression;

    }

    for(i=0;i<NSTRUTCOLORS;i++) {

      comp_level[i] = log(mincompression + 1) + 
	i*(log(maxcompression + 1) - log(mincompression + 1))/(double)(NSTRUTCOLORS);
      
      comp_bins[i] = calloc(n_struts,sizeof(int));
      binsize[i] = 0;
      
    }

    /* Now sort the strut records into bins by colors. */

    for(i=0;i<2*n_struts;i++) {

      if (as_list[i].s > as_list[i].t) { /* Plot lower triangle */

	if (as_list[i].s >= sbox[0] &&   /* Clip to the sbox, tbox */
	    as_list[i].s <= sbox[1] &&
	    as_list[i].t >= tbox[0] &&
	    as_list[i].t <= tbox[1]) {

	  for(j=NSTRUTCOLORS-1;j>=0;j--) {

	    if (log(as_list[i].strut->compression + 1) > comp_level[j]) {

	      comp_bins[j][binsize[j]] = i;  /* Add to group j */
	      binsize[j]++;
	      break;

	    }
	 
	  }

	}

      }

    }

    /* We have now grouped all of the struts to be displayed into
       NSTRUTCOLORS bins based on the cutoffs in comp_level. 
       It remains to actually plot rectangles in these colors. */
    
    for(i=0;i<NSTRUTCOLORS;i++) {

      if (binsize[i] == 0) { continue; }  /* Short-circuit the loop if 
					     there's nothing to show */

      set_ps_color(strutcolors[NSTRUTCOLORS - i - 1],outfile);
    
      fprintf(outfile,"[ ");

      for(j=0;j<binsize[i];j++) {

	fprintf(outfile, " %d %d ",
		(int)(10*as_list[comp_bins[i][j]].s/boxwidth) - (int)(5*BOXMULT),
		(int)(10*as_list[comp_bins[i][j]].t/boxwidth) - (int)(5*BOXMULT));
	
	struts_output++;

	if (!(struts_output%8)) { fprintf(outfile,"\n"); }
	
      }

      fprintf(outfile,"] aload pop \n");
      fprintf(outfile,"0 1 %d { pop %g %g rectfill } for \n",
	      binsize[i]-1,(double)(10.0*BOXMULT),(double)(10.0*BOXMULT));
    
      free(comp_bins[i]);

    }

  }

  /***************************************************************/
  /*                                                             */
  /*  7. Plot curvature data.                                    */
  /*                                                             */
  /***************************************************************/
 
  if (KAPPAPLOT) {
    
    printf("\tPlotting curvature data (with subsampling)...\n");
    
    if (pscomments) {
      fprintf(outfile,"%% Plot curvature data.\n");
    }
    fprintf(outfile,"gsave\n");
    fprintf(outfile,"0 0 moveto 0 %g translate\n",
	    -2.5*tickplotstep);  

    fprintf(outfile,"%g setlinewidth\n",0.5/plotscale);
    set_ps_color(strutcolors[4],outfile);
    fprintf(outfile,"1 setlinejoin\n");
    fprintf(outfile,"1 setlinecap\n");
    fprintf(outfile,"0 0 %d %g rectclip\n",10*NBOXES,2*tickplotstep);
    
    curvpixel = 0.5/plotscale; 
    
    /* This is "one pixel" on the printed graph at 150 dpi * 1.5 inches */
    /* In order to compress the curvature data, we plot a point, 
       and then read ahead, skipping measurements as long as 
       they are within 2*curvpixel of the line joining the last 
       plotted point and the current address. */

    printed_total = 0;
    min_curv = DBL_MAX; max_curv = -1;
    
    for(cmp=0;cmp<L->nc;cmp++) {
      
      printedpt = 0;
      
      for(vert=1;vert<L->cp[cmp].nv;vert++) { 
	
	intpos = (int)(10*s_positions[cmp][vert]/boxwidth);
	yi = curv = plc_MR_curvature(L,cmp,vert)*tickplotstep;
	
	min_curv = (curv/tickplotstep < min_curv) ? curv/tickplotstep : min_curv;
	max_curv = (curv/tickplotstep > max_curv) ? curv/tickplotstep : max_curv;

	fprintf(outfile,"%d %g ",intpos,curv);      
	printedpt++;
	
	if (!(printedpt%8)) { fprintf(outfile,"\n"); }
	if (vert > L->cp[cmp].nv-3) { continue; }
	
	/* Now we attempt to read ahead */
	
	for(pt_ok=true,j=vert+2;pt_ok && j < L->cp[cmp].nv;j++) {
	  
	  yj = plc_MR_curvature(L,cmp,j)*tickplotstep;
	  mij = (yj - yi)/(s_positions[cmp][j] - s_positions[cmp][vert]);
	  
	  /* We must check every point in i+1..j-1 for closeness to the line */
	  
	  for(k=vert+1;k<j;k++) {
	    
	    predicted_curvk = mij * (s_positions[cmp][k] - s_positions[cmp][vert]);
	    
	    if (fabs(predicted_curvk - tickplotstep*plc_MR_curvature(L,cmp,k)) > 
          2*curvpixel) {
	      
	      pt_ok = false;
	      
	    }
	    
	  }
	  
	}
	
	/* We have reached the point where i and j-2 are the farthest apart we can get */
	/* without causing more than 2*curvpixel of error in the plot. We set i to j-3 */
	/* since it will be incremented when we return to the top of the loop. */
	
	vert = j-3;
	
      }
      
      fprintf(outfile,"moveto 2 1 %d { pop lineto } for stroke\n",
	      printedpt);
      
      printed_total += printedpt;
      
    }
    
    fprintf(outfile,"grestore\n");
    printf("\t Printed %d points (of %d possible).\n",printed_total,plc_num_edges(L));
    

  }

  /******************************************************************/
  /*                                                                */
  /* 9. Plot relative compressions.                                 */
  /*                                                                */
  /******************************************************************/

  if (PLOT_COMPRESSIONS) {
    
    printf("\tPlotting compression data (with subsampling)...\n");
    
    if (pscomments) {
      fprintf(outfile,"%% Plot compressions.\n");
    }
    fprintf(outfile,"gsave\n");
    fprintf(outfile,"0 0 moveto 0 %g translate\n",
	    -2.5*tickplotstep);  

    fprintf(outfile,"%g setlinewidth\n",0.75/plotscale);
    set_ps_color(compressioncolor,outfile);
    fprintf(outfile,"1 setlinejoin\n");
    fprintf(outfile,"1 setlinecap\n");
    fprintf(outfile,"0 0 %d %g rectclip\n",10*NBOXES,2*tickplotstep);
    
    curvpixel = 0.5/plotscale; 
    
    /* This is "one pixel" on the printed graph at 150 dpi * 1.5 inches */
    /* In order to compress the data, we plot a point, and then read
       ahead, skipping measurements as long as they are within
       2*curvpixel of the line joining the last plotted point and the
       current address. */

    /* We now search through our struts and add the compresions incident on 
       each edge of the graph. */

    printed_total = 0;
    
    for(cmp=0;cmp<L->nc;cmp++) {
      
      printedpt = 0;
      maxcompression = -1;

      /* Gather the aggregate compression information for each vertex */

      compressions = (double *)(calloc(L->cp[cmp].nv,sizeof(double)));
      
      for(i=0;i<n_struts;i++) {

	if (strutlist[i].component[0] == cmp) {

	  compressions[strutlist[i].lead_vert[0]] += strutlist[i].compression;

	} 
	
	if (strutlist[i].component[1] == cmp) {

	  compressions[strutlist[i].lead_vert[1]] += strutlist[i].compression;

	}

      }

      for(i=0;i<L->cp[cmp].nv;i++) { 
	
	maxcompression = (compressions[i] > maxcompression) ? 
	  compressions[i] : maxcompression; 

      }

      /* Now plot the data. */
      
      for(vert=1;vert<L->cp[cmp].nv;vert++) { 
	
	intpos = (int)(10*s_positions[cmp][vert]/boxwidth);
	yi = curv = compressions[vert]*tickplotstep/maxcompression;
	
	fprintf(outfile,"%d %g ",intpos,curv);      
	printedpt++;
	
	if (!(printedpt%8)) { fprintf(outfile,"\n"); }
	if (vert > L->cp[cmp].nv-3) { continue; }
	
	/* Now we attempt to read ahead */
	
	for(pt_ok=true,j=vert+2;pt_ok && j < L->cp[cmp].nv;j++) {
	  
	  yj = compressions[j]*tickplotstep/maxcompression;
	  mij = (yj - yi)/(s_positions[cmp][j] - s_positions[cmp][vert]);
	  
	  /* We must check every point in i+1..j-1 for closeness to the line */
	  
	  for(k=vert+1;k<j;k++) {
	    
	    predicted_curvk = mij * (s_positions[cmp][k] - s_positions[cmp][vert]);
	    
	    if (fabs(predicted_curvk - tickplotstep*compressions[k]/maxcompression) 
		> 2*curvpixel) {
	      
	      pt_ok = false;
	      
	    }
	    
	  }
	  
	}
	
	/* We have reached the point where i and j-2 are the farthest apart we can get */
	/* without causing more than 2*curvpixel of error in the plot. We set i to j-3 */
	/* since it will be incremented when we return to the top of the loop. */
	
	vert = j-3;
	
      }
      
      fprintf(outfile,"moveto 2 1 %d { pop lineto } for stroke\n",
	      printedpt);
      
      printed_total += printedpt;

      free(compressions);
      
    }
    
    fprintf(outfile,"grestore\n");
    printf("\t Printed %d points (of %d possible).\n",printed_total,plc_num_edges(L));
    

  }


  /**************************************************************/
  /*                                                            */
  /*  8. Add frames and right-hand labels and clean up memory.  */
  /*                                                            */
  /**************************************************************/

  if (PLOT_BACKGROUND) {

    if (nvtags > 0) { display_vtags(9.0/plotscale,nvtags,vtags,outfile); }
    
    if (pscomments) {
      fprintf(outfile,"%% Add frames and right-hand labels.\n");
    }
    fprintf(outfile,
	    "0 0 moveto %d 0 rlineto "
	    "0 %d rlineto -%d -%d rlineto "
	    "0.6 setgray stroke \n",
	    NBOXES*10,NBOXES*10,NBOXES*10,NBOXES*10);
    
    if (KAPPAPLOT) {
      
      fprintf(outfile,"0.6 setgray 0 %g %d %g rectstroke\n",
	      -2.5*tickplotstep,10*NBOXES,2*tickplotstep);
      
    }

  }

  /***********************************************************************/
  /*                                                                     */
  /* 2.5. Draw the key                                                   */
  /*                                                                     */
  /***********************************************************************/

  if (DATA_IN_COMMENTS) {

    fprintf(outfile,
	    "%% SPCOMMENT Filename: %s \n"
	    "%% SPCOMMENT Verts: %d \n"
	    "%% SPCOMMENT Struts: %d \n"
	    "%% SPCOMMENT Curvature Range: [%g, %g] \n"
	    "%% SPCOMMENT Kink Regions: %d \n"
	    "%% SPCOMMENT Straight Length: %g\n",
	    FILENAME,plc_num_edges(L),n_struts,
	    min_curv,max_curv,nkinks,straightlength);

  }

  if (PLOT_KEY) {

    if (pscomments) {
      fprintf(outfile,"%% Draw the key.\n");
    }
    fprintf(outfile,
	    "/textln { dup stringwidth pop exch show neg %g rmoveto } def\n",
	    -12.0/plotscale);

    fprintf(outfile,"0.1 setgray\n");

    fprintf(outfile,
	    "0 %g moveto \n"
	    "(%s / %d verts / %d struts) textln \n"
	    "(Length: %g / minstrut: %g / minrad: %g) textln \n",
	    ((whole) ? 12.0*NBOXES : 10.0*NBOXES),
	    FILENAME,plc_num_edges(L),n_struts,
	    octrope_curvelength(L),shortest,octrope_minradval(L));

    if (KAPPAPLOT) {

      fprintf(outfile,
	      "(Curvature range: [%g, %g]) textln \n",
	      min_curv,max_curv);

    }
    
    if (KINK_HIGHLIGHT) {
      
      if (nkinks > 0) {
	
	fprintf(outfile,
		"(%d kink regions) textln \n",
		nkinks);

      }
      
    }
	   
    if (SS_HIGHLIGHT) {

      if (straightlength > 0) {

	fprintf(outfile,
		"(%g units strut-free curve) textln \n",
		straightlength);

      }

    }
    
  }  
  
  /* We're done! Free memory and quit. */

  fprintf(outfile,"showpage\n"); 
  free(as_list);
  free(s_positions);  
    
}


int main(int argc,char *argv[]) {

  FILE *infile_fptr = {NULL},*outfile_fptr = {NULL},*sf_fptr = {NULL};
  plCurve *L;

  int    sl_size,n_struts;
  octrope_strut *strutlist;
  octrope_mrloc *mrloc_list;
  int num_mr_locs,mr_size;
  double eps = 1e-5;
  double tuberad = 0.5;
  double lambda = 1;
  
  int nerrors;
  char outfile_name[1000];
  
  struct arg_file *infile = arg_filen(NULL,NULL,"<VECT file>", 1, 10000, "input VECT file(s)");
  struct arg_int  *levels = arg_int0("l", "levels","<n>", "number of octree levels");
  struct arg_int  *debuglevel = arg_int0("v", "verbosity","<n>", "level of debugging information to print (0-9)");
  struct arg_dbl  *epsilon     = arg_dbl0("e","epsilon","<x>","find struts within this tolerance of minimum length");
  struct arg_dbl  *tube_radius = arg_dbl0("r","radius","<x>","plot all struts with length < x");
  struct arg_dbl  *arg_lambda = arg_dbl0("l","lambda","<x>","set stiffness lambda");
  struct arg_lit  *sciencecolors = arg_lit0(NULL,"sciencecolors","use different color scheme");

  struct arg_lit  *help        = arg_lit0("h","help","display help message");
  
  struct arg_lit  *nokinks     = arg_lit0("k","nokinks","don't highlight kinks");
  struct arg_lit  *nogrid      = arg_lit0("g","nogrid","don't draw grid");

  struct arg_lit  *noss        = arg_lit0("s","nostraight","don't highlight strut-free regions");
  struct arg_int  *maxstruts   = arg_int0("M","maxstruts","<n>","maximum number of struts to find");
  struct arg_lit  *nokappaplot = arg_lit0("c","nokappa","don't draw curvature plot");
  struct arg_dbl  *strutsize   = arg_dbl0(NULL,"strutsize","<x>","box size to plot "
					  "strut (mult. of average edgelength)");

  struct arg_lit  *epsplot     = arg_lit0(NULL,"make-eps","produce EPS output");
  struct arg_lit  *combineplots= arg_lit0(NULL,"combine","produce a single (combined) output file");
  struct arg_lit  *nokey       = arg_lit0(NULL,"nokey","don't display the key with each plot");
  struct arg_file  *strutfile   = arg_filen(NULL,"strutfile","<STRUT file>",0,10000,
					   "input STRUT file(s) to go with input VECT files");
  struct arg_int  *ticks       = arg_int0("t","ticks","<n>","number of (arclength) ticks to draw");
  struct arg_lit  *compressions = arg_lit0(NULL,"compressions","attempt to plot with compressions (experimental)");
  struct arg_dbl  *plotwidth   = arg_dbl0("w","width","<x>","plot width in inches");

  struct arg_dbl  *slo = arg_dbl0(NULL,"sl,s_clip_low","[0..1]","clip struts with (relative) s value less than this");
  struct arg_dbl  *shi = arg_dbl0(NULL,"sh,s_clip_high","[0..1]","clip struts with (relative) s value more than this");
  struct arg_dbl  *tlo = arg_dbl0(NULL,"tl,t_clip_low","[0..1]","clip struts with (relative) t value less than this");
  struct arg_dbl  *thi = arg_dbl0(NULL,"th,t_clip_high","[0..1]","clip struts with (relative) t value more than this");
  struct arg_lit  *nobg = arg_lit0(NULL,"nobackground","don't plot background");
  struct arg_lit  *nodsc = arg_lit0(NULL,"nodsc","don't use Adobe Document Structuring");
  struct arg_lit  *whole = arg_lit0(NULL,"whole","show the whole s,t plane (not just triangle)");
  struct arg_lit  *pscomments = arg_lit0(NULL,"psc","add comments to PS code");
  struct arg_lit  *dataincomments = arg_lit0(NULL,"DataComments","add data to ps file as comments");

  struct arg_file *outfile = arg_file0("o","outfile","<file>","output filename");
  struct arg_end  *end = arg_end(20);

  void *argtable[] = {help, epsilon, tube_radius, arg_lambda,levels, debuglevel, 
		      nokinks, nogrid, noss,
                      nokappaplot, ticks, plotwidth, maxstruts, outfile,
                      infile, strutfile, combineplots, nokey, compressions, sciencecolors,
                      strutsize, slo, shi, tlo, thi, nobg, nodsc, whole,
                      epsplot, pscomments, dataincomments,end}; 
  
  struct arg_end  *helpend = arg_end(20);

  void *helptable[] = {help, helpend};
  int filenum;
  char revision[20] = "$Revision: 1.20 $";
  char *dollar;

  dollar = strchr(&revision[1],'$');
  dollar[0] = '\0';
  printf("Strutplot v%s, ridgerunner v%s\n",&revision[11],PACKAGE_VERSION);
 /* The PACKAGE_VERSION preprocessor symbol is defined in "config.h" by autoconf. */

  printf("  Produce a Postscript plot of the struts for a polygonal knot.\n");
  
  /* We start by parsing the command-line arguments with argtable. */
  
  if (arg_nullcheck(argtable) != 0 || arg_nullcheck(helptable) != 0)
    printf("error: insufficient memory\n");

  nerrors = arg_parse(argc,argv,argtable);
  
  if (nerrors > 0) { /* The standard syntax didn't match-- try the help table */

    nerrors = arg_parse(argc,argv,helptable);

    if (nerrors > 0) {  /* The help table didn't match either-- the first set of */
                        /* errors was probably more helpful, so we display it. */

      arg_print_errors(stdout,helpend,"helptable");
      arg_print_errors(stdout,end,"struts");
      exit(1);

    } else {  /* The help table matched, which means we asked for help or gave nothing */
  
      printf("strutplot displays a graph of the self-contacts of a Geomview VECT file.\n"
           "usage: \n\n");
      arg_print_glossary(stdout, argtable," %-25s %s\n");
      exit(0);

    }
    
  }
  
  /* Convert command-line arguments into proper global variables */
  
  if (levels->count > 0) {
    
    octrope_set_levels(*(levels->ival));
    
  }

  if (arg_lambda->count > 0) {

    lambda = arg_lambda->dval[0];

  }

  if (tube_radius->count > 0) {

    tuberad = tube_radius->dval[0];

  }
  
  if (debuglevel->count > 0) {
    
    octrope_set_debug(*(debuglevel->ival));
    DEBUGLEVEL = *(debuglevel->ival);
    printf("Debug level set to %d.\n",*(debuglevel->ival));
    
  }

  bool use_eps = false;

  if (epsilon->count > 0) {

    eps = *(epsilon->dval);
    use_eps = true;

  }

  if (nokinks->count > 0) {

    KINK_HIGHLIGHT = false;

  }

  if (nogrid->count > 0) {
    
    DRAWGRID = false;

  }

  if (noss->count > 0) {

    SS_HIGHLIGHT = false;

  }

  if (nokappaplot->count > 0) {

    KAPPAPLOT = false;

  }

  if (nokey->count > 0) {

    PLOT_KEY = false;
    
  }

  if (ticks->count > 0) {

    NTICKS = *(ticks->ival);

    if (NTICKS < 0)    { NTICKS = 0; }
    if (NTICKS > 100)  { NTICKS = 100; }
    
  }

  if (compressions->count > 0 && strutfile->count != 0) {

    PLOT_COMPRESSIONS = true;

  } else {

    PLOT_COMPRESSIONS = false;

  }

  if (epsplot->count > 0 && combineplots->count == 0) {

    EPSPLOT = true;

  }

  if (plotwidth->count > 0) {

    PLOTWIDTH = *(plotwidth->dval);

  }

  if (strutsize->count > 0) {

    BOXMULT = *(strutsize->dval);

  } 

  if (dataincomments->count > 0) {

    DATA_IN_COMMENTS = 1;

  }

  if (sciencecolors->count > 0) {

    checkercolors[0] = sciencechecker[0];
    checkercolors[1] = sciencechecker[1];
    
    int i;
    for(i=0;i<5;i++) {strutcolors[i] = sciencestrutcolor;}

  }

  /* Set the s and t boxes. */

  if (slo->count > 0) {
    
    sbox[0] = *(slo->dval);

    sbox[0] = (sbox[0] > 0) ? sbox[0] : 0;
    sbox[0] = (sbox[0] < 1) ? sbox[0] : 1;

    if (sbox[0] == sbox[1]) {

      fprintf(stderr,"strutplot: Warning. sbox is [%g,%g].\n",sbox[0],sbox[1]);
      
    }

  }

  if (shi->count > 0) {

    sbox[1] = *(shi->dval);

    sbox[1] = (sbox[1] > 0) ? sbox[1] : 0;
    sbox[1] = (sbox[1] < 1) ? sbox[1] : 1;

    if (sbox[0] == sbox[1]) {

      fprintf(stderr,"strutplot: Warning. sbox is [%g,%g].\n",sbox[0],sbox[1]);
      
    }

  }

  if (sbox[0] > sbox[1]) {

    fprintf(stderr,"strutplot: There are no struts with s value %g <= s <= %g.\n",
	    sbox[0],sbox[1]);

    exit(1);

  }
  

  if (tlo->count > 0) {
    
    tbox[0] = *(tlo->dval);

    tbox[0] = (tbox[0] > 0) ? tbox[0] : 0;
    tbox[0] = (tbox[0] < 1) ? tbox[0] : 1;

  }

  if (thi->count > 0) {

    tbox[1] = *(thi->dval);

    tbox[1] = (tbox[1] > 0) ? tbox[1] : 0;
    tbox[1] = (tbox[1] < 1) ? tbox[1] : 1;

  }

  if (tbox[0] > tbox[1]) {

    fprintf(stderr,"strutplot: There are no struts with t value %g <= t <= %g.\n",
	    tbox[0],tbox[1]);

    exit(1);

  }

  if (nobg->count > 0) {

    PLOT_BACKGROUND = false;

  }

  if (nodsc->count > 0) {

    USE_DSC = false;

  }
  
  /* Now we begin the main loop of the program */

  for (filenum = 0;filenum < infile->count;filenum++) {

    FILENAME = infile->filename[filenum];
    
    infile_fptr = fopen(infile->filename[filenum],"r");
    
    if (infile_fptr == NULL) {
      
      fprintf(stderr,"strutplot: Couldn't open file %s.\n",infile->filename[filenum]);
      continue;
      
    }
    
    octrope_error_num = 0;
    L = plc_read(infile_fptr,&octrope_error_num,octrope_error_str,80);
    
    /* We now demonstrate the octrope library's error handling protocol: */
    
    if (octrope_error_num > 0) {   /* This is the signal for an error. */
      
      fprintf(stderr,"strutplot: link reading error\n%s\n",octrope_error_str);
      continue;
      
    }
    
    fclose(infile_fptr);

    printf("Plotting struts for %d component, %d edge link %s.\n",
	   L->nc,plc_num_edges(L),infile->filename[filenum]);
    
    /* We now need to open the output file-- there are a couple of cases here.
     * If an outfile is given, then it trumps everything (and implicitly sets
     * combineplots), so we open that output file on infile 1 (and ignore this
     * for all future infiles).  If combineplots is set, we call the output
     * "combined_plots.ps" no matter what.  In all other cases, we do
     * individual filename munging to produce a nicely named output file
     * corresponding to each input file. 
     */
    
    if (outfile->count > 0 || combineplots->count > 0) {
      
      if (filenum == 0) {
	
	if (outfile->count > 0) {
	  
	  outfile_fptr = fopen(*(outfile->filename),"w");
	  
	  if (outfile_fptr == NULL) {
	    
	    fprintf(stderr,"strutplot: Couldn't open %s for writing.\n",
		    *(outfile->filename));
	    exit(1);
	    
	  }
	  
	} else {
	  
	  outfile_fptr = fopen("combined_plots.ps","w");
	  
	  if (outfile_fptr == NULL) {
	    
	    fprintf(stderr,"strutplot: Couldn't open combined_plots.ps for writing.\n");
	    exit(1);
	    
	  }

	}

	/* This is a multi-page DSC document. */

	if (USE_DSC) {

	  fprintf(outfile_fptr,
		  "%%!PS-Adobe-3.0 \n"
		  "%%%%Pages: %d\n"
		  "%%%%Title: (%d Combined Strut Plots)\n"
		  "%%%%DocumentNeededResources: font Helvetica\n"
		  "%%%%EndComments\n",infile->count,infile->count); 

	}
      }
      
    } else {   /* Filename munging mode */

      if (outfile_fptr != NULL) { if (USE_DSC) {fprintf(outfile_fptr,"%%%%EOF\n");} fclose(outfile_fptr); }
      
      if (strlen(infile->basename[filenum]) > sizeof(outfile_name)-20) {
	
	fprintf(stderr,
          "strutplot: Ridiculously long input filename can't be parsed.\n");
	exit(1);
	
      }
      
      sprintf(outfile_name,"%s",infile->basename[filenum]);
      
      if (strstr(outfile_name,".vect") != NULL) {
	
	if (!EPSPLOT) {
	  
	  sprintf(strstr(outfile_name,".vect"),".stplot.ps");
	  
	} else {
	  
	  sprintf(strstr(outfile_name,".vect"),".stplot.eps");
	  
	}
	
      } else {
	
	if (!EPSPLOT) {
	  
	  strcat(outfile_name,".stplot.ps");
	  
	} else {
	  
	  strcat(outfile_name,".stplot.eps");
	  
	}
	
      }
      
      outfile_fptr = fopen(outfile_name,"w");
      
      if (outfile_fptr == NULL) {
	
	fprintf(stderr,"strutplot: Couldn't open %s for writing.\n",
          outfile_name);
	exit(1);
	
      }

      /* This is a single-page DSC document. */

      if (USE_DSC) {

	if (!EPSPLOT) {
	  
	  fprintf(outfile_fptr,
		  "%%!PS-Adobe-3.0 \n"
		  "%%%%Pages: %d\n"
		  "%%%%Title: (%s Strut Plot)\n"
		  "%%%%DocumentNeededResources: font Helvetica\n"
		  "%%%%EndComments\n",infile->count,infile->filename[0]); 
	  
	} else {
	  
	  /* EPS stuff goes here. */
	  
	}

      }
      
    }
    
    if (!EPSPLOT && USE_DSC) {

      fprintf(outfile_fptr,"%%%%Page: %d %d\n", filenum+1,filenum+1);

    }

    /* We have now opened an output filename one way or the other. */
    /* We must now get some struts and kinks, either from an input strutfile */
    /* or from calling struts and kinks ourselves. */
    
    if (strutfile->count > 0) {
      
      if (strutfile->count != infile->count) {
	
	fprintf(stderr,"strutplot: The number of STRUT files (%d) "
		"must match the number of VECT files (%d).\n",
		strutfile->count,infile->count);
	exit(1);
	
      }
      
      sf_fptr = fopen(strutfile->filename[filenum],"r");

      if (sf_fptr == NULL) {

	fprintf(stderr,"strutplot: Couldn't open STRUT file %s.\n",
		strutfile->filename[filenum]);
	continue;

      }
      
      octrope_strutfile_read(L,&n_struts,&strutlist,&num_mr_locs,&mrloc_list,sf_fptr);
      
      if (octrope_error_num > 0) {   /* This is the signal for an error. */
	
	fprintf(stderr,"strutplot: strutfile reading error\n%s\n",octrope_error_str);
	continue;
	
      }

      fclose(sf_fptr);

    } else {  /* We now compute the strut set. */
  
      if (maxstruts->count > 0) {
        sl_size = *(maxstruts->ival);
      } else {
        /* We allocate an extra-large buffer. */
        sl_size = 20*plc_num_edges(L);      
      }
      strutlist = (octrope_strut *)(calloc(sl_size,sizeof(octrope_strut)));

      if (use_eps) {

	n_struts = octrope_struts(L,0,eps,strutlist,sl_size,&shortest,NULL,0);
	
      } else {
	
	n_struts = octrope_struts(L,2*tuberad,0,strutlist,sl_size,&shortest,NULL,0);

      }

      if (octrope_error_num > 0) {
	
	fprintf(stderr,"%s\n",octrope_error_str);
	exit(1);
	
      }
      
      /* We now compute the set of minrad locations */
      
      mr_size = 2*plc_num_edges(L);      /* We allocate an extra-large buffer. */
      mrloc_list = (octrope_mrloc *)(calloc(mr_size,sizeof(octrope_mrloc)));
      
      if (use_eps) {

	octrope_minrad(L,(shortest+eps)/2.0,0,mrloc_list,mr_size,&num_mr_locs);

      } else {

	octrope_minrad(L,tuberad*lambda,0,mrloc_list,mr_size,&num_mr_locs);

      }
      
      if (octrope_error_num > 0) {
	
	fprintf(stderr,"%s\n",octrope_error_str);
	exit(1);
	
      }

    }

    BOXMULT *= (double)(plc_num_edges(L))/250.0;

    /* We now have a list of struts and an output file. Go ahead and plot. */
      
    create_st_plot(L,n_struts,strutlist,num_mr_locs,mrloc_list,
                   (whole->count > 0),(pscomments->count > 0),
                   outfile_fptr);       
  
    /* Now free memory. */

    free(strutlist); free(mrloc_list);
    plc_free(L);

  }
  
  /* We now close the last file */
  
  if (outfile_fptr != NULL) {

    if (USE_DSC) {fprintf(outfile_fptr,"%%%%EOF\n");}
    fclose(outfile_fptr);
    exit(0);

  } else {

    exit(1);

  }
}
