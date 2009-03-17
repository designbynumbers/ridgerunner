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

int strut_free_vertices( plCurve* inLink, double tube_radius, int *cp, int *vt)

/* Returns two int buffers giving component and vertex numbers for all
   vertices which are known to be free of struts on adjacent
   edges. The list is sorted in dictionary order on (cp,vt) pairs. 
   The buffers cp and vt are required to be at least as large as 
   plc_num_verts(inLink). */

{

  int strutStorageSize, strutCount, nVerts;
  octrope_strut *strutSet;
  double shortest;
  plCurve *scratchcurve;
  octrope_strut *st;
  
  if (VERBOSITY >= 10) { logprintf("\tstrut_free_vertices..."); }
  
  strutStorageSize = 6*plc_num_verts(inLink);
  strutSet = malloc_or_die(sizeof(octrope_strut)*strutStorageSize);
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
  
  for(i=0;i<strutCount;i++) {

    st = &(strutSet[i]);

    scratchcurve->cp[st->component[0]].vt[st->lead_vert[0]].c[0] = 1.0;
    scratchcurve->cp[st->component[1]].vt[st->lead_vert[1]].c[0] = 1.0;
 
    /* Unfortunately, the last vertex on a closed curve is a special case. */
   
    if (!scratchcurve->cp[st->component[0]].open && st->lead_vert[0] == scratchcurve->cp[st->component[0]].nv-1) {
      
      scratchcurve->cp[st->component[0]].vt[0].c[0] = 1.0;

    } else {
      
      scratchcurve->cp[st->component[0]].vt[st->lead_vert[0]+1].c[0] = 1.0;
  
    }
	
    if (!scratchcurve->cp[st->component[1]].open && st->lead_vert[1] == scratchcurve->cp[st->component[1]].nv-1) {
       
      scratchcurve->cp[st->component[1]].vt[0].c[0] = 1.0;

    } else {
      
      scratchcurve->cp[st->component[1]].vt[st->lead_vert[1]+1].c[0] = 1.0;
  
    }
    
  }
  
  /* Now we go ahead and read out the results onto the cp and vt buffers. */

  int nStrutFree = 0;

  for(i=0;i<scratchcurve->nc;i++) {
    
    for(j=0;j<scratchcurve->cp[i].nv;j++) {

      if (scratchcurve->cp[i].vt[j].c[0] < 0.5) { 

	cp[nStrutFree] = i; vt[nStrutFree] = j; nStrutFree++; 

      }

    }

  }

  /* Now we can go ahead and free the memory we've used. */

  plc_free(scratchcurve);
  free(strutSet);

  return nFree;

}

