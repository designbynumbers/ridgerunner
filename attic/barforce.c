/*

barforce.c : Code archive for ridgerunner project. 

We did some experiments with plCurves with fixed length edges acting
under imposed external forces which didn't work very well. This code
was the implementation. Apart from a mechanical port to plCurve, this
hasn't been brought into v2, and may or may not compile.

Copyright Jason Cantarella.

This file is part of ridgerunner. ridgerunner is free software: you can
redistribute it and/or modify it under the terms of the GNU General
Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version.

ridgerunner is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.  You should have received a copy of the GNU General
Public License along with ridgerunner. If not, see
<https://www.gnu.org/licenses/>.

*/

void
updateSideLengths( plCurve* inLink, search_state* inState )
{
  int cItr, vItr;
  
  inState->totalSides = 0;
  for( cItr=0; cItr<inLink->nc; cItr++ ) {

    inState->totalSides += inLink->cp[cItr].nv;
  }
  
  if( inState->sideLengths != NULL )
    free(inState->sideLengths);
  
  inState->sideLengths = (double*)malloc(inState->totalSides*sizeof(double));
  fatalifnull_(inState->sideLengths);
  
  int tot = 0;
  for( cItr=0; cItr<inLink->nc; cItr++ ) {

    for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ ) {
         
      inState->sideLengths[inState->compOffsets[cItr] + vItr] = 
	plc_M_distance(inLink->cp[cItr].vt[vItr],
		       inLink->cp[cItr].vt[vItr+1]) /*;*/
	
      inState->avgSideLength += inState->sideLengths[inState->compOffsets[cItr] + vItr];
      tot++;
    }
  }
  inState->avgSideLength /= (double)tot;
}


static void
barForce( plc_vector* dVdt, plCurve* inLink, search_state* inState )
{
  /* solve the unconstrained least squares problem Ax=b where A is the rigidity
   * matrix formed by bars on the vertices of the link, and b is the change in each
   * component's x, y, z from its original length
   */
  double*		A = NULL, *b = NULL;
  int			barVerts=0, cItr, vItr;
  int			bars=0, aIndexer;
  
  plc_fix_wrap(inLink);
  
  for( cItr=0; cItr<inLink->nc; cItr++ ) {
    
    if( inState->conserveLength[cItr] != 0 ) {
      
      barVerts += inLink->cp[cItr].nv;
      bars += plc_strand_edges(&inLink->cp[cItr]);
    }
  }
  
  A = (double*)malloc(sizeof(double)*barVerts*3*bars);
  placeVertexBars(A, inLink, 0, barVerts, bars, inState);
  b = (double*)calloc(barVerts*3, sizeof(double));
  
  /* now compute b, the change in each dimension of each vertex required to get to eq */
  plc_vector* sides;
  double*			lengths;
  int				edges;
  
  for( cItr=0; cItr<inLink->nc; cItr++ ) {
    
    if( inState->conserveLength[cItr] != 0 ) {
      
      double length;
      edges = plc_strand_edges(&inLink->cp[cItr]);
      
      for( vItr=0; vItr<edges; vItr++ ) {
	
	plc_vector s1, s2, side;
	s1.c[0] = inLink->cp[cItr].vt[vItr].c[0];
	s1.c[1] = inLink->cp[cItr].vt[vItr].c[1];
	s1.c[2] = inLink->cp[cItr].vt[vItr].c[2];
	s2.c[0] = inLink->cp[cItr].vt[vItr+1].c[0];
	s2.c[1] = inLink->cp[cItr].vt[vItr+1].c[1];
	s2.c[2] = inLink->cp[cItr].vt[vItr+1].c[2];
	side.c[0] = s2.c[0] - s1.c[0];
	side.c[1] = s2.c[1] - s1.c[1];
	side.c[2] = s2.c[2] - s1.c[2];
	
	length = plc_M_norm(side);
	
	// unit-ize sides since we will use these as the tangential motion basis below
	side.c[0] /= length;
	side.c[1] /= length;
	side.c[2] /= length;
	
	b[inState->compOffsets[cItr]+vItr + 0] = 
	  -((inState->sideLengths[inState->compOffsets[cItr]+vItr] - length) * side.c[0]);
	b[inState->compOffsets[cItr]+vItr + 1] = 
	  -((inState->sideLengths[inState->compOffsets[cItr]+vItr] - length) * side.c[1]);
	b[inState->compOffsets[cItr]+vItr + 2] = 
	  -((inState->sideLengths[inState->compOffsets[cItr]+vItr] - length) * side.c[2]);
	      
	//	b[inState->compOffsets[cItr]+vItr + 0] = 
	//         dVdt[inState->compOffsets[cItr]+vItr + 0].c[0];
	//	b[inState->compOffsets[cItr]+vItr + 1] = 
	//         dVdt[inState->compOffsets[cItr]+vItr + 0].c[1];
	//	b[inState->compOffsets[cItr]+vItr + 2] = 
	//         dVdt[inState->compOffsets[cItr]+vItr + 0].c[2];
	
	printf( "adjusting (%d) %d: %lf\n", cItr, vItr, 
		(inState->sideLengths[inState->compOffsets[cItr]+vItr] - length) );
      }
    }
  }
  
  if( barVerts != 0 ) {

    double* compressions;
    taucs_ccs_matrix* sparseA;
    
    sparseA = taucs_construct_sorted_ccs_matrix(A, bars, 3*barVerts);
    compressions = t_lsqr( sparseA, b );
    //	taucs_print_ccs_matrix(sparseA);
    taucs_ccs_free(sparseA);
    
    for( cItr=0; cItr<inLink->nc; cItr++ ) {

      if( inState->conserveLength[cItr] != 0 ) {

	double partialMult = 0;
	int totalStruts = bars;
	int dlItr, sItr;
	
	for( dlItr=inState->compOffsets[cItr]; 
	     dlItr<inState->compOffsets[cItr]+inLink->cp[cItr].nv; dlItr++ ) {

	  int cItr2;
	  /* of course, dlItr is not quite right here since we are
	   * including ALL vertices in its computation, but all
	   * vertices are not accounted for in A. We need to adjust
	   * A's indexer here to account for this, which amounts to
	   * subtracting from dlItr all the components which preceed
	   * the one we are currently processing that are NOT bar
	   * composed
					 */
	  aIndexer = dlItr;
	  for( cItr2=0; cItr2<cItr; cItr2++ ) {
	    
	    if( inState->conserveLength[cItr2] == 0 )
	      aIndexer -= inLink->cp[cItr2].nv;
	  
	  }
	  
	  for( sItr=0; sItr<totalStruts; sItr++ ) {

	    if( (totalStruts*3*aIndexer)+sItr+0 >= bars*3*barVerts )
	      printf( "*** crap!\n" );
	    
	    partialMult += compressions[sItr]*A[(totalStruts*3*aIndexer)+sItr+0];
	  
	  }
	  
	  dVdt[dlItr].c[0] = dVdt[dlItr].c[0] + partialMult;  
	  partialMult=0;
	
	  for( sItr=0; sItr<totalStruts; sItr++ ) {

	    if( (totalStruts*3*aIndexer)+sItr+0 >= bars*3*barVerts )
	      printf( "*** crap!\n" );
		
	    partialMult += compressions[sItr]*A[(totalStruts*3*aIndexer)+sItr+totalStruts];
					
	  }
	  
	  dVdt[dlItr].c[1] = dVdt[dlItr].c[1] + partialMult;
					
	  partialMult = 0;
	  for( sItr=0; sItr<totalStruts; sItr++ ) {

	    if( (totalStruts*3*aIndexer)+sItr+0 >= bars*3*barVerts )
	      printf( "*** crap!\n" );
	  
	    partialMult += 
	      compressions[sItr]*A[(totalStruts*3*aIndexer)+sItr+(2*totalStruts)];
	  }
	  dVdt[dlItr].c[2] = dVdt[dlItr].c[2] + partialMult;
	}
      }
    }
    
    free(compressions);
  }
  
  // debug code!
  char	fname[1024];
  preptmpname(fname,"dVdt.vect",inState);
  exportVect(dVdt, inLink, fname);
  
  free(A);
  free(b);
}


void
placeVertexBars( double* A, plCurve* inLink, int contactStruts, int totalBarVerts, int totalBars, search_state* inState )
{
  int totalStruts = contactStruts + totalBars;
  int cItr, vItr, sItr, next;
  
  plc_vector  strutDirection;
  plc_vector  points[2];
  double norm = 0;
  
  int thresh = ((totalBarVerts+contactStruts)*3)*(totalBars+contactStruts);
  
  octrope_strut*  vertStruts;
  
  vertStruts = (octrope_strut*)malloc(sizeof(octrope_strut)*totalBars);
  
  sItr=0;
  for( cItr=0; cItr<inLink->nc; cItr++ )
    {
      if( inState->conserveLength[cItr] != 0 )
	{
	  for( vItr=0; vItr<plc_strand_edges(&inLink->cp[cItr]); vItr++ )
	    {
	      next = ((vItr+1)%inLink->cp[cItr].nv);
	      
	      vertStruts[sItr].component[0] = cItr;
	      vertStruts[sItr].component[1] = cItr;
	      vertStruts[sItr].position[0] = 0;
	      vertStruts[sItr].position[1] = 0;
	      vertStruts[sItr].lead_vert[0] = vItr;
	      vertStruts[sItr].lead_vert[1] = next;
	      
	      sItr++;
	    }
	}
    }
  
  // and now we can place these just like minrad struts
  for( sItr=0; sItr<totalBars; sItr++ )
    {
      int entry;
      
      //	printf( "(%d) (%d) %d %d\n", vertStruts[sItr].component[0], vertStruts[sItr].component[1], vertStruts[sItr].lead_vert[0], vertStruts[sItr].lead_vert[1] );
      
      // calc the norm of the converted strut
      octrope_strut_ends( inLink, &vertStruts[sItr], points );
      
      // the normalized difference of pointOne, pointTwo is the strut force vector
      strutDirection.c[0] = points[0].c[0] - points[1].c[0];
      strutDirection.c[1] = points[0].c[1] - points[1].c[1];
      strutDirection.c[2] = points[0].c[2] - points[1].c[2];
      
      norm = plc_M_norm(strutDirection);
      strutDirection.c[0] /= norm;
      strutDirection.c[1] /= norm;
      strutDirection.c[2] /= norm;
      
      /* temporarily increment the strut's verts based on their
       * component interactions, but in the vertex bar rigidity
       * matrix, only some component's vertices are counted -- namely
       * the ones that are composed of bars. We need to increase our
       * offset, but only by those guys we have skipped who are bar
       * composed.
       */
      for( cItr=0; cItr<vertStruts[sItr].component[0]; cItr++ )
	{
	  if( inState->conserveLength[cItr] != 0 )
	    {
	      /* we update both here even though we only check first
	       * strut end because vertex bars are always on the same
	       * component */
	      vertStruts[sItr].lead_vert[0] += inLink->cp[cItr].nv;
	      vertStruts[sItr].lead_vert[1] += inLink->cp[cItr].nv;
	    }
	}
      
      // entry is the offset in A which begin this strut's influce
      // it corresponds to the x influence on the lead_vert[0]th vertex
      // after this line, entry+1 is y, +2: z.
      entry = (totalStruts*3*vertStruts[sItr].lead_vert[0])+contactStruts+sItr;
      
      if( entry+(2*totalStruts) >= thresh )
	printf( "****crap!\n" );	
      
      // the strut information includes the position from the strut.lead_vert
      // so we assign "1-position[0]" of the force to the lead vert and "position[0]"
      // of the force to lead_vert+1
      A[entry] = (1-vertStruts[sItr].position[0]) * strutDirection.c[0];
      A[entry+totalStruts] = (1-vertStruts[sItr].position[0]) * strutDirection.c[1];
      A[entry+(2*totalStruts)] = (1-vertStruts[sItr].position[0]) * strutDirection.c[2];
      
      /***************** we don't need to do this since there are no
 midpoint vertex bars *****************/
      
		// now for the next vertex, receiving "position[0]" of
		// the force, this is potential wrapping case
	/*	if(
			(vertStruts[sItr].lead_vert[0]-inState->compOffsets[vertStruts[sItr].component[0]])
			==
			(inLink->cp[vertStruts[sItr].component[0]].nv-1)
			&&
			(inLink->cp[vertStruts[sItr].component[0]].acyclic
			== 0) ) { entry = 0; } else { entry =
			(totalStruts*3*(vertStruts[sItr].lead_vert[0]+1))+contactStruts+sItr;
			}
		
		if( entry+(2*totalStruts) >= thresh ) printf(
			"****crap!\n" );
		
		A[entry] = (vertStruts[sItr].position[0]) * strutDirection.c[0];
		A[entry+totalStruts] = (vertStruts[sItr].position[0]) * strutDirection.c[1];
		A[entry+(2*totalStruts)] = (vertStruts[sItr].position[0]) * strutDirection.c[2];
		*/
		

      // we do the same thing at the opposite end of the strut, except now the 
      // force is negated
      entry = (totalStruts*3*vertStruts[sItr].lead_vert[1])+contactStruts+sItr;
      
      if( entry+(2*totalStruts) >= thresh )
	printf( "****crap!\n" );			
      
      A[entry] = (1-vertStruts[sItr].position[1]) * -strutDirection.c[0];
      A[entry+totalStruts] = (1-vertStruts[sItr].position[1]) * -strutDirection.c[1];
      A[entry+(2*totalStruts)] = (1-vertStruts[sItr].position[1]) * -strutDirection.c[2];
      
      /*if( (vertStruts[sItr].lead_vert[1]-inState->compOffsets[vertStruts[sItr].component[1]]) == (inLink->cp[vertStruts[sItr].component[1]].nv-1) &&
	(inLink->cp[vertStruts[sItr].component[1]].acyclic == 0) )
	{
	entry = 0;
	}
	else
	{
	entry = (totalStruts*3*(vertStruts[sItr].lead_vert[1]+1))+contactStruts+sItr;
	}
	
	if( entry+(2*totalStruts) >= thresh )
	printf( "****crap!\n" );	
	
	A[entry] = (vertStruts[sItr].position[1]) * -strutDirection.c[0];
	A[entry+totalStruts] = (vertStruts[sItr].position[1]) * -strutDirection.c[1];
	A[entry+(2*totalStruts)] = (vertStruts[sItr].position[1]) * -strutDirection.c[2];
      */
      /* readjusting is problematic for the same reason that adjusting
       * is (see above) and it's not necessary in this strut
       * placement, so we just don't do it
       */
      
      //vertStruts[sItr].lead_vert[0] -=
      //inState->compOffsets[vertStruts[sItr].component[0]];
      //vertStruts[sItr].lead_vert[1] -=
      //inState->compOffsets[vertStruts[sItr].component[1]];
    }
  
  free(vertStruts);
}
