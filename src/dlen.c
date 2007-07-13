/*
 *  dlen.c
 *  ridgerunner3
 *
 *  Created by Michael Piatek on Sun Jan 18 2004.
 *  Edited/ported JHC 6/2007.

    New version doesn't attempt to use blas operations to construct
    the field, instead sticks with plCurve primitives. This should be
    slower, but more portable and maintainable.

 *  
 */

#include"ridgerunner.h"

/*
 * on input - a link for which to calculate dlen
 * on output - a vector field of the elastic force on vertices 
 *				THAT THE USER MUST DISPOSE OF

   We note that ioDL is a buffer of plc_verts(inLink) plc_vectors,
   containing velocity vectors for each of the vertices in inLink,
   in the order 

   inLink->cp[0].vt[0] 
   .
   .
   .
   inLink->cp[0].vt[inLink->cp[0].nv-1]
   inLink->cp[1].vt[0]
   .
   .
   .
   inLink->cp[1].vt[inLink->cp[1].nv-1]
   .
   .
   .
   inLink->cp[inLink->nc-1].vt[..]

   The dlenForce vector is added to ioDL. This is because the
   equilaterization force may be already added to ioDL.

 */
void
dlenForce( plc_vector* ioDL, plCurve* inLink, search_state* inState )
{
  int cItr, vItr, totalVerts=0, dlItr;
  plc_vector* diffVectors;
  char errmsg[1024],dumpname[1024];
  
  /* Allocate a "flat" buffer of vectors which will be used to construct
     the dLen vector */
  
  totalVerts = inState->totalVerts;
  diffVectors = (plc_vector*)calloc(totalVerts, sizeof(struct plc_vector_type));
  fatalifnull_(diffVectors);

  /* Now construct the vector. */
  
  int i=0,istart,iend,nv;
  bool norm1ok = true, norm2ok = true;

  for( cItr=0; cItr<inLink->nc; cItr++ ) {		   
    for( vItr=0, istart=i; vItr<inLink->cp[cItr].nv; vItr++ ) {

      /* dLen[v_i] = (v_{i+1} - v_i)/|v_{i+1} - v_i| +(v_{i-1} - v_i)/|v_{i-1} + v_i|.*/

      diffVectors[i] = plc_vect_sum(
				    
		         plc_normalize_vect(
					    plc_vect_diff(inLink->cp[cItr].vt[vItr-1],
							  inLink->cp[cItr].vt[vItr]),
					    &norm1ok),

			 plc_normalize_vect(
					    plc_vect_diff(inLink->cp[cItr].vt[vItr+1],
							  inLink->cp[cItr].vt[vItr]),
					    &norm2ok)
			 );
      i++;

      /* We actually expect some failures of the normalize calls if we are on an 
	 open component. So we can't just check that the normXok variables are both
	 true. Instead, we check something more complicated. */

      if ( (!norm1ok && (vItr != 0 || !inLink->cp[cItr].open)) || 
	   (!norm2ok && (vItr != inLink->cp[cItr].nv-1 || !inLink->cp[cItr].open)) ) {

	dumpLink(inLink,inState,dumpname);
	sprintf(errmsg,
		"ridgerunner: Zero length edge adjacent to vert %d of cmp %d \n"
		"             of link. inLink->cp[%d].open = %d. \n"
		"             Link dumped to %s.\n",
		vItr,cItr,cItr,inLink->cp[cItr].open,dumpname);

      }
							       
    }

    if (inLink->cp[cItr].open) { /* We need to make sure that end vectors are set
				    correctly for open curves. */
      iend = i-1;
      nv = inLink->cp[cItr].nv;

      diffVectors[istart] = plc_normalize_vect(plc_vect_diff(inLink->cp[cItr].vt[1],
							     inLink->cp[cItr].vt[0]),
					       &norm1ok);

      diffVectors[iend] = plc_normalize_vect(plc_vect_diff(inLink->cp[cItr].vt[nv-2],
							   inLink->cp[cItr].vt[nv-1]),
					     &norm2ok);
      assert(norm1ok && norm2ok);

    }  
  }

  /* We leave in some old code here which zeros the field if "conserveLength" is set. */
  /* This seems to be left over from some experiments with tubes which have fixed     */
  /* length and are just subject to other forces. We comment this out, since we don't */
  /* include this functionality in the current RR build. */
  
 /*  for( cItr=0; cItr<inLink->nc; cItr++ ) { */
/*     if( inState->conserveLength[cItr] != 0 ) { */
/*       for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ ) */
/* 	{ */
/* 	  diffVectors[vItr+inState->compOffsets[cItr]].c[0] = 0; */
/* 	  diffVectors[vItr+inState->compOffsets[cItr]].c[1] = 0; */
/* 	  diffVectors[vItr+inState->compOffsets[cItr]].c[2] = 0; */
/* 	} */
/*     } */
/*   } */
  
  for( dlItr=0; dlItr<totalVerts; dlItr++ ) {

    plc_M_add_vect(ioDL[dlItr],diffVectors[dlItr]); /* A += B */

  }
  
  free(diffVectors);

}

