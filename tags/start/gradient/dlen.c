/*
 *  dlen.c
 *  ridgerunner3
 *
 *  Created by Michael Piatek on Sun Jan 18 2004.
 *  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
 *
 *  PROBLEMS:
 *
 *  THOUGHTS:
 *      - we might be able to squeeze more speed out of a G4 version of things 
 *          by turning off Java mode
 *      - do we really even want to use double precision here?
 *
 */
#include <math.h>

#ifdef __APPLE__
#include <vecLib/vBLAS.h>
#else
#include <gsl/gsl_cblas.h>
#endif // __APPLE__

#include "dlen.h"
#include "errors.h"
#include "link.h"
#include "eqedge.h"

/*
 * on input - a link for which to calculate dlen
 * on output - a vector field of the elastic force on vertices 
 *				THAT THE USER MUST DISPOSE OF
 */
void
dlenForce( octrope_vector* ioDL, octrope_link* inLink, search_state* inState )
{
    // grab unit edges
    double* norms;
    //vector* units;
    int cItr, vItr, nextVert, totalVerts=0, dlItr;
    octrope_vector* diffVectors;
    octrope_vector* unitsCopy;
    
    totalVerts = 0;
	for( cItr=0; cItr<inLink->nc; cItr++ )
		totalVerts += inLink->cp[cItr].nv;
    
    norms = (double*)malloc(sizeof(double)*totalVerts);
    //units = (vector*)malloc(sizeof(struct vector_type)*totalVerts);
    diffVectors = (octrope_vector*)calloc(totalVerts, sizeof(struct octrope_vector_type));
    unitsCopy = (octrope_vector*)malloc(sizeof(struct octrope_vector_type)*totalVerts);
    
    fatalifnull_(norms);
    fatalifnull_(diffVectors);
    fatalifnull_(unitsCopy);
	
//	if( (rand() % 2) == 1 )
	{
		int i=0;
		inState->curvature_step = 1;
		for( cItr=0; cItr<inLink->nc; cItr++ )
		{		
			for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
			{
				nextVert = (vItr+1)%inLink->cp[cItr].nv;
				diffVectors[i].c[0] = inLink->cp[cItr].vt[nextVert].c[0] - inLink->cp[cItr].vt[vItr].c[0];
				diffVectors[i].c[1] = inLink->cp[cItr].vt[nextVert].c[1] - inLink->cp[cItr].vt[vItr].c[1];
				diffVectors[i].c[2] = inLink->cp[cItr].vt[nextVert].c[2] - inLink->cp[cItr].vt[vItr].c[2];
				
				// ddot is dot product
				norms[i] = sqrt(cblas_ddot(3, &diffVectors[i].c[0], 1, &diffVectors[i].c[0], 1));
				i++;
			}
		}
		
		// dscal is scale
		i=0;
		for( cItr=0; cItr<inLink->nc; cItr++ )
		{
			// effectively divide by norm
			for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++, i++ )
				cblas_dscal( 3, (1/norms[i]), &diffVectors[i].c[0], 1 );
		}
		
		// now diffVectors are forward units, sum with opposite of previous to get dLen field
		// duplicate the units since daxpy will operate in place
		memcpy( unitsCopy, diffVectors, sizeof(struct octrope_vector_type)*totalVerts );
		i=0;
		int prev;
		for( cItr=0; cItr<inLink->nc; cItr++ )
		{
			for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++)
			{
				// cblas_daxpy - adds a constant times a vector to another vector
				prev = (vItr == 0) ? (inLink->cp[cItr].nv - 1) : (vItr-1);
				prev += i;
				cblas_daxpy( 3, -1, &unitsCopy[prev].c[0], 1, &diffVectors[i+vItr].c[0], 1 );
			}
			i += inLink->cp[cItr].nv;
		
			// if this component is open, zero the end guys so things aren't screwy
			if( inLink->cp[cItr].acyclic != 0  )
			{
				diffVectors[i-inLink->cp[cItr].nv].c[0] = 0;
				diffVectors[i-inLink->cp[cItr].nv].c[1] = 0;
				diffVectors[i-inLink->cp[cItr].nv].c[2] = 0;
				diffVectors[i-1].c[0] = 0;
				diffVectors[i-1].c[1] = 0;
				diffVectors[i-1].c[2] = 0;
			}
		}
    }
//	else
//	{
//		printf( "*" );
//		inState->nocurvature_step = 1;
//	}
	
    free(norms);
    free(unitsCopy);
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		if( inState->conserveLength[cItr] != 0 )
		{
			for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
			{
				diffVectors[vItr+inState->compOffsets[cItr]].c[0] = 0;
				diffVectors[vItr+inState->compOffsets[cItr]].c[1] = 0;
				diffVectors[vItr+inState->compOffsets[cItr]].c[2] = 0;
			}
		}
	}

	for( dlItr=0; dlItr<inState->totalVerts; dlItr++ )
	{
		ioDL[dlItr].c[0] += diffVectors[dlItr].c[0];
		ioDL[dlItr].c[1] += diffVectors[dlItr].c[1];
		ioDL[dlItr].c[2] += diffVectors[dlItr].c[2];
	}

	free(diffVectors);
	
//    return diffVectors; // by now, the actual dlen field
}

