/*
 *  free_edge.c
 *  ridgerunner
 *
 *  Created by Jason Cantarella on March 18, 2009.
 *  Copyright under the usual GNU licence for rr.
 *
 */

#include "ridgerunner.h"

/* Function prototypes */

plCurve*  straighten_free_edges( plCurve* inLink, search_state* inState );

/* Subsidiary functions. */

typedef struct vertex_record_type {

  int  cp;
  int  vt; 
  bool strutfree;

} vertex_record;

int strut_free_vertices( plCurve* inLink, double tube_radius, int *cpBuf, int *vtBuf)

/* Returns two int buffers giving component and vertex numbers for all
   vertices which are known to be free of struts on adjacent
   edges. The list is sorted in dictionary order on (cp,vt) pairs. 
   The buffers cp and vt are required to be at least as large as 
   plc_num_verts(inLink). We require a run of strut-free vertices to 
   be at least 4 vertices long before we report it. */

{

  int strutStorageSize, strutCount;
  octrope_strut *strutSet;
  double shortest;
  plCurve *scratchcurve;
  octrope_strut *st;
  int vRadius = 4;
  
  if (VERBOSITY >= 10) { logprintf("\tstrut_free_vertices..."); }
  
  strutStorageSize = 6*plc_num_verts(inLink);
  strutSet = malloc_or_die(sizeof(octrope_strut)*strutStorageSize, __FILE__ , __LINE__ );
  strutCount = octrope_struts(inLink, 2*tube_radius, 0, strutSet, strutStorageSize, &shortest, NULL, 0);
  assert(strutCount < strutStorageSize);

  /* Now we can parse the results. We do this in a funny way, using a
     dummy copy of the curve as a way to index the vertices of the
     curve. */

  int i,j;

  scratchcurve = plc_copy(inLink);

  for(i=0;i<scratchcurve->nc;i++) {

    for(j=0;j<scratchcurve->cp[i].nv;j++) {

      scratchcurve->cp[i].vt[j].c[0] = 0;  /* We will store 0 for no struts, 1 for struts. */

    }
    
  }

  /* Now we read the strut list and update the dummy vertex information. */
  
  int lv,end,cp;

  for(i=0;i<strutCount;i++) {

    st = &(strutSet[i]);

    for(end=0;end<2;end++) {

      cp = st->component[end];

      if (!scratchcurve->cp[cp].open) { /* We need to wrap addressing */
	
	lv = st->lead_vert[end] + scratchcurve->cp[cp].nv; // Added to make modular arithmetic come out right
      
	for(j=0;j<vRadius;j++) {

	  scratchcurve->cp[cp].vt[(lv + j) % scratchcurve->cp[cp].nv].c[0] = 1.0;
	  scratchcurve->cp[cp].vt[(lv - j) % scratchcurve->cp[cp].nv].c[0] = 1.0;
	  
	}

      } else { /* Instead, we need to cut off addressing at ends of curve */

	lv = st->lead_vert[end];
	
	for(j=0;j<vRadius && (lv + j) <scratchcurve->cp[cp].nv;j++) {

	   scratchcurve->cp[cp].vt[(lv + j)].c[0] = 1.0;
	  
	}

	for(j=0;j<vRadius && (lv - j) >= 0;j++) {

	  scratchcurve->cp[cp].vt[(lv - j)].c[0] = 1.0;
	  
	}

      }

    }

  }
  
  /* Now we go ahead and read out the results onto the cp and vt buffers. */

  int nStrutFree = 0;

  for(i=0;i<scratchcurve->nc;i++) {
    
    for(j=0;j<scratchcurve->cp[i].nv;j++) {

      if (scratchcurve->cp[i].vt[j].c[0] < 0.5) { 

	cpBuf[nStrutFree] = i; vtBuf[nStrutFree] = j; nStrutFree++; 

      }

    }

  }

  /* Now we can go ahead and free the memory we've used. */

  plc_free(scratchcurve);
  free(strutSet);

  return nStrutFree;

}

void highlight_curve(plCurve *L, search_state *state)

/* Highlight straight segments, kinks, and other "understood" portions of a link. */

{
  int cp;

  /* First, double check that the number of colors is set large enough for every component. */

  for(cp=0;cp<L->nc;cp++) {

    if (L->cp[cp].cc != L->cp[cp].nv) { 

      FatalError("Not enough colors in a component of the final link.", __FILE__ , __LINE__ );

    }

  }

  /* Straight segments. */

  int *cpBuf,*vtBuf, nStrutFree;

  cpBuf = malloc_or_die(sizeof(int)*plc_num_verts(L), __FILE__ , __LINE__ );
  vtBuf = malloc_or_die(sizeof(int)*plc_num_verts(L), __FILE__ , __LINE__ );

  nStrutFree = strut_free_vertices(L,state->tube_radius,cpBuf,vtBuf);

  int i;

  for(i=0;i<nStrutFree;i++) {

    L->cp[cpBuf[i]].clr[vtBuf[i]] = gStraightSegColor;

  }

  free(cpBuf); free(vtBuf);

  /* Kinked regions */

  octrope_mrloc *mrBuf;
  int nKinkVerts;

  mrBuf = malloc_or_die(sizeof(octrope_mrloc)*plc_num_verts(L), __FILE__ , __LINE__ );
  octrope_minrad(L,gLambda*state->tube_radius,0,mrBuf,plc_num_verts(L),&nKinkVerts);

  for(i=0;i<nKinkVerts;i++) {

    L->cp[mrBuf[i].component].clr[mrBuf[i].vert] = gKinkColor;

  }

  free(mrBuf);

}
      
  
  
