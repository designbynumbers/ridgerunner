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

int strut_free_vertices( plCurve* L, double tube_radius, bool *freeFlag)

/* Fills a "flat" buffer freeFlag, expected to be plc_num_verts(inLink) in size, 
   with "true" if the vertex is no closer than 4 vertices to a strut and "false" otherwise.
   The buffer refers to vertices on the plCurve inlink in dictionary order on (cp,vt). */

{

  int strutStorageSize, strutCount;
  octrope_strut *strutSet;
  double shortest;
  int vRadius = 4,tVerts,nStrutFree;
  
  tVerts = plc_num_verts(L);

  strutStorageSize = 6*tVerts;
  strutSet = malloc_or_die(sizeof(octrope_strut)*strutStorageSize, __FILE__ , __LINE__ );
  strutCount = octrope_struts(L, 2*tube_radius, 0, strutSet, strutStorageSize, &shortest, NULL, 0);
  assert(strutCount < strutStorageSize);

  int i,j;

  for(i=0;i<tVerts;i++) {
    
      freeFlag[i] = true;
      
  }
    
  /* Now we read the strut list and update the flags accordingly. */
  
  int lv,end,cp;

  for(i=0;i<strutCount;i++) {

    for(end=0;end<2;end++) {

      cp = strutSet[i].component[end];
      lv = strutSet[i].lead_vert[end];

      for(j=0;j<vRadius;j++) {
	
	freeFlag[plc_vertex_num(L,cp,lv + j)] = false;
	freeFlag[plc_vertex_num(L,cp,lv - j)] = false;

      }

    }

  }

  /* Now we can go ahead and free the memory we've used. */

  free(strutSet);

  for(i=0,nStrutFree=0;i<tVerts;i++) { if (freeFlag[i]) { nStrutFree++; } }

  return nStrutFree;

}

void highlight_curve(plCurve *L, search_state *state)

/* Highlight straight segments, kinks, and other "understood" portions of a link. */

{
  int cp,vt;

  /* First, double check that the number of colors is set large enough for every component. */

  for(cp=0;cp<L->nc;cp++) {

    if (L->cp[cp].cc != L->cp[cp].nv) { 

      FatalError("Not enough colors in a component of the final link.", __FILE__ , __LINE__ );

    }

  }

  /* Straight segments. */

  int nStrutFree;
  bool *freeFlag;

  freeFlag = malloc_or_die(sizeof(bool)*plc_num_verts(L), __FILE__ , __LINE__ );
  nStrutFree = strut_free_vertices(L,state->tube_radius,freeFlag);

  int i=0;

  for(cp=0;cp<L->nc;cp++) {

    for(vt=0;vt<L->cp[cp].nv;vt++,i++) {

      if (freeFlag[i]) { L->cp[cp].clr[vt] = gStraightSegColor; }

    }

  }

  free(freeFlag);

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
     
  
void accelerate_free_vertices( plc_vector *dLen, plCurve *L, double tube_radius)
/* Make the free vertices "faster" than the other verts. */
{

  int nStrutFree;
  bool *freeFlag;

  freeFlag = malloc_or_die(sizeof(bool)*plc_num_verts(L), __FILE__ , __LINE__ );
  nStrutFree = strut_free_vertices(L,tube_radius,freeFlag);
  
  int cp, vt, i=0;

  for(cp=0;cp<L->nc;cp++) {

    for(vt=0;vt<L->cp[cp].nv;vt++,i++) {

      if (freeFlag[i]) { dLen[i] = plc_scale_vect(5,dLen[i]); }

    }

  }

  free(freeFlag);

}

double strut_free_length( plCurve *L, search_state *state) 

/* Returns the total strut-free arclength of the curve. */

{
  int nStrutFree;
  bool *freeFlag;

  freeFlag = malloc_or_die(sizeof(bool)*plc_num_verts(L), __FILE__ , __LINE__ );
  nStrutFree = strut_free_vertices(L,state->tube_radius,freeFlag);

  int cmpItr,vItr,vt=0;
  double sfLen = 0;

  for(cmpItr=0;cmpItr < L->nc;cmpItr++) {

    for(vItr=0;vItr < L->cp[cmpItr].nv-1;vItr++,vt++) {

      if (freeFlag[vt] || freeFlag[vt+1]) {

	sfLen += plc_distance(L->cp[cmpItr].vt[vItr],L->cp[cmpItr].vt[vItr+1]);

      }

    } 

    if (L->cp[cmpItr].open == false) { /* We need to deal with the last edge */
      
      if (freeFlag[vt] || freeFlag[0]) {

	sfLen += plc_distance(L->cp[cmpItr].vt[vItr],L->cp[cmpItr].vt[0]);

      }

    } 

    vt++;

  }
      
  free(freeFlag);
   
  return sfLen;

}

