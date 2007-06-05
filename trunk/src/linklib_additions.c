/*
   Functions to be included in the plCurve library.
*/

#include "ridgerunner.h"

int    plc_strand_edges(plc_strand *P)

     /* Procedure computes the number of edges in a pline. */
     /* This is either equal to, or one less than, the number */
     /* of vertices, depending on whether the pline is closed. */

{
  
  return P->nv - (P->acyclic ? 1:0);
  
}

double
plCurve_torsion( plCurve* inLink, FILE* outPlot )
{
  double			totalTorsion = 0;
  plc_vector	cross1, cross2;
  double			temp;
  int				cItr, vItr;
  
  bzero(&cross1,sizeof(plc_vector));
  bzero(&cross2,sizeof(plc_vector));
  
  int				e1, e2, e3;	
  double			lengthOffset = 0;
  int				vertOffset = 0;
  
  totalTorsion = 0;
  
  for( cItr=0; cItr<inLink->nc; cItr++ )
    {
      for( vItr=0; vItr<plc_strand_edges(&inLink->cp[cItr]); vItr++)
	{
	  plc_vector s1, s2, s3;
	  
	  e1 = vItr-1;
	  e2 = vItr;
	  e3 = vItr+1;
	  
	  s1 = plc_vect_diff( inLink->cp[cItr].vt[e1+1], inLink->cp[cItr].vt[e1] );
	  s2 = plc_vect_diff( inLink->cp[cItr].vt[e2+1], inLink->cp[cItr].vt[e2] );
	  s3 = plc_vect_diff( inLink->cp[cItr].vt[(e3+1)%plc_strand_edges(&inLink->cp[cItr])], inLink->cp[cItr].vt[e3] );
	  
	  /* This is the torsion at edge vItr */
	  cross1 = plc_cross_prod(s1, s2); 
	  cross2 = plc_cross_prod(s2, s3); 
	  
	  if( (plc_M_norm(cross1)*plc_M_norm(cross2)) != 0 )
	    {
	      temp = (plc_M_dot(cross1,cross2))/(plc_M_norm(cross1)*plc_M_norm(cross2));
	      
	      if(temp >= 1)
		temp = 0;
	      else if(temp <= -1)
		temp = M_PI;
	      else
		temp = acos(temp);
	      
	      totalTorsion += temp;
	      
	      fprintf(outPlot, "%d %3.8lf %3.8lf %3.8lf\n", vertOffset, lengthOffset, temp, 2*tan(temp*.5)/plc_M_norm(s1));
	    }
	  else
	    {
	      // The norm of one of these cross products is 0.  If that's the
	      // case, the the sin of the angle between two of the edges must be 0 or
	      // Pi, which implies that they are colinear, which implies that the
	      // osculating plane isn't changing, which means that there shouldn't be
	      // any torsion at this edge.
	      
	      // then again, this probably never happens numerically.
	      totalTorsion += 0; // for emphasis!
	      
	      fprintf(outPlot, "%d %3.8lf 0.0 0.0\n", vertOffset, lengthOffset);
	    }
	  
	  lengthOffset += plc_M_norm(s2);
	  vertOffset++;
	  
	} // vert itr
    } // component itr
  return totalTorsion;
}

plCurve*
octrope_equalize_density( plCurve* inLink, search_state* inState )
{
  /*	plCurve*   doubled;
	int*			newVerts = (int*)malloc(sizeof(int)*inLink->nc);
	int*			cyclicity = (int*)malloc(sizeof(int)*inLink->nc);
	int*			color_count = (int*)malloc(sizeof(int)*inLink->nc);
	int				cItr, vItr;
	double			length_i;
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		
	//	newVerts[cItr] = inLink->cp[cItr].nv * 2;
		length_i = 0;
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
			length_i += inState->sideLengths[vItr + inState->compOffsets[cItr]];
										
		cyclicity[cItr] = inLink->cp[cItr].acyclic;
		color_count[cItr] = 1;
	}
	
	doubled = plc_new( inLink->nc,
						newVerts,
						cyclicity, 
						color_count );
	
	plc_fix_wrap(inLink);
						
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		// preserve color
		//memcpy( &doubled->cp[cItr].clr[0], &inLink->cp[cItr].clr[0], sizeof(plc_color) );
	
		doubled->cp[cItr].clr[0].r = inLink->cp[cItr].clr[0].r;
		doubled->cp[cItr].clr[0].g = inLink->cp[cItr].clr[0].g;
		doubled->cp[cItr].clr[0].b = inLink->cp[cItr].clr[0].b;
		doubled->cp[cItr].clr[0].alpha = inLink->cp[cItr].clr[0].alpha;
	
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
		{			
			if( cItr == comp )
			{
				doubled->cp[cItr].vt[vItr*2+0].c[0] = inLink->cp[cItr].vt[vItr].c[0];
				doubled->cp[cItr].vt[vItr*2+0].c[1] = inLink->cp[cItr].vt[vItr].c[1];
				doubled->cp[cItr].vt[vItr*2+0].c[2] = inLink->cp[cItr].vt[vItr].c[2];

				doubled->cp[cItr].vt[vItr*2+1].c[0] = 0.5*inLink->cp[cItr].vt[vItr].c[0] + 0.5*inLink->cp[cItr].vt[vItr+1].c[0];
				doubled->cp[cItr].vt[vItr*2+1].c[1] = 0.5*inLink->cp[cItr].vt[vItr].c[1] + 0.5*inLink->cp[cItr].vt[vItr+1].c[1];
				doubled->cp[cItr].vt[vItr*2+1].c[2] = 0.5*inLink->cp[cItr].vt[vItr].c[2] + 0.5*inLink->cp[cItr].vt[vItr+1].c[2];
			}
			else
			{
				doubled->cp[cItr].vt[vItr].c[0] = inLink->cp[cItr].vt[vItr].c[0];
				doubled->cp[cItr].vt[vItr].c[1] = inLink->cp[cItr].vt[vItr].c[1];
				doubled->cp[cItr].vt[vItr].c[2] = inLink->cp[cItr].vt[vItr].c[2];
			}
		}
	}
	
	free(newVerts);
	free(cyclicity);
	free(color_count);

	return doubled;*/
	
	return NULL;
}

plCurve*
octrope_double_component( plCurve* inLink, int comp )
{
	/* double number of sides through bisection -- keeps things eq and shouldn't munger too much with
	 * struts
	 */
	plCurve*   doubled;
	int*			newVerts = (int*)malloc(sizeof(int)*inLink->nc);
	int*			cyclicity = (int*)malloc(sizeof(int)*inLink->nc);
	int*			color_count = (int*)malloc(sizeof(int)*inLink->nc);
	int				cItr, vItr;
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		newVerts[cItr] = inLink->cp[cItr].nv * 2;
		cyclicity[cItr] = inLink->cp[cItr].acyclic;
		color_count[cItr] = 1;
	}
	
	doubled = plc_new( inLink->nc,
						newVerts,
						cyclicity, 
						color_count );
	
	plc_fix_wrap(inLink);

	cItr = comp;
//	for( cItr=inDo; cItr<inLink->nc; cItr++ )
	{
		// preserve color
		//memcpy( &doubled->cp[cItr].clr[0], &inLink->cp[cItr].clr[0], sizeof(plc_color) );
	
		doubled->cp[cItr].clr[0].r = inLink->cp[cItr].clr[0].r;
		doubled->cp[cItr].clr[0].g = inLink->cp[cItr].clr[0].g;
		doubled->cp[cItr].clr[0].b = inLink->cp[cItr].clr[0].b;
		doubled->cp[cItr].clr[0].alpha = inLink->cp[cItr].clr[0].alpha;
	
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
		{
			doubled->cp[cItr].vt[vItr*2+0].c[0] = inLink->cp[cItr].vt[vItr].c[0];
			doubled->cp[cItr].vt[vItr*2+0].c[1] = inLink->cp[cItr].vt[vItr].c[1];
			doubled->cp[cItr].vt[vItr*2+0].c[2] = inLink->cp[cItr].vt[vItr].c[2];
			
			doubled->cp[cItr].vt[vItr*2+1].c[0] = 0.5*inLink->cp[cItr].vt[vItr].c[0] + 0.5*inLink->cp[cItr].vt[vItr+1].c[0];
			doubled->cp[cItr].vt[vItr*2+1].c[1] = 0.5*inLink->cp[cItr].vt[vItr].c[1] + 0.5*inLink->cp[cItr].vt[vItr+1].c[1];
			doubled->cp[cItr].vt[vItr*2+1].c[2] = 0.5*inLink->cp[cItr].vt[vItr].c[2] + 0.5*inLink->cp[cItr].vt[vItr+1].c[2];
		}
	}
	
	free(newVerts);
	free(cyclicity);
	free(color_count);

	return doubled;
}

plCurve*
octrope_double_edges( plCurve* inLink )
{
	/* double number of sides through bisection -- keeps things eq and shouldn't munger too much with
	 * struts
	 */
	plCurve*   doubled;
	int*			newVerts = (int*)malloc(sizeof(int)*inLink->nc);
	int*			cyclicity = (int*)malloc(sizeof(int)*inLink->nc);
	int*			color_count = (int*)malloc(sizeof(int)*inLink->nc);
	int				cItr, vItr;
	
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		newVerts[cItr] = inLink->cp[cItr].nv * 2;
		cyclicity[cItr] = inLink->cp[cItr].acyclic;
		color_count[cItr] = 1;
	}
	
	doubled = plc_new( inLink->nc,
						newVerts,
						cyclicity, 
						color_count );
	
	plc_fix_wrap(inLink);
						
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		// preserve color
		//memcpy( &doubled->cp[cItr].clr[0], &inLink->cp[cItr].clr[0], sizeof(plc_color) );
	
		doubled->cp[cItr].clr[0].r = inLink->cp[cItr].clr[0].r;
		doubled->cp[cItr].clr[0].g = inLink->cp[cItr].clr[0].g;
		doubled->cp[cItr].clr[0].b = inLink->cp[cItr].clr[0].b;
		doubled->cp[cItr].clr[0].alpha = inLink->cp[cItr].clr[0].alpha;
	
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
		{
			doubled->cp[cItr].vt[vItr*2+0].c[0] = inLink->cp[cItr].vt[vItr].c[0];
			doubled->cp[cItr].vt[vItr*2+0].c[1] = inLink->cp[cItr].vt[vItr].c[1];
			doubled->cp[cItr].vt[vItr*2+0].c[2] = inLink->cp[cItr].vt[vItr].c[2];
			
			doubled->cp[cItr].vt[vItr*2+1].c[0] = 0.5*inLink->cp[cItr].vt[vItr].c[0] + 0.5*inLink->cp[cItr].vt[vItr+1].c[0];
			doubled->cp[cItr].vt[vItr*2+1].c[1] = 0.5*inLink->cp[cItr].vt[vItr].c[1] + 0.5*inLink->cp[cItr].vt[vItr+1].c[1];
			doubled->cp[cItr].vt[vItr*2+1].c[2] = 0.5*inLink->cp[cItr].vt[vItr].c[2] + 0.5*inLink->cp[cItr].vt[vItr+1].c[2];
		}
	}
	
	free(newVerts);
	free(cyclicity);
	free(color_count);

	return doubled;
}

/* 
 * The following proceedure steps around the polygon redistributing vertices as necessary
 * to make things approach equilateral. This code adapted from Eric Rawdon's TOROS 
 * polygonal runaround -- not to be confused with the spherical runaround or the tangential 
 * stepper. Contacts: rawdon@mathcs.duq.edu piatek@mathcs.duq.edu
 */
plCurve*
octrope_fixlength( plCurve* inLink )
{
	plCurve* fixed;
	int cItr, vItr;
	
	fixed = plc_copy(inLink);
	
	plc_fix_wrap(inLink);
	
	for( cItr=0; cItr<fixed->nc; cItr++ )
	{
		double goal, used, tmpgoal, left;
		int i=0, j=1, edges;
		double* lengths;
		plc_vector* sides;
		
		edges = plc_strand_edges(&inLink->cp[cItr]);
		goal = plc_strand_length(&inLink->cp[cItr]) / (double)edges;
		
		// we need the edge lengths and sides
		lengths = (double*)malloc(sizeof(double)*edges);
		sides = (plc_vector*)malloc(sizeof(plc_vector)*edges);
		for( vItr=0; vItr<edges; vItr++ )
		{
			plc_vector s1, s2;
			s1.c[0] = inLink->cp[cItr].vt[vItr].c[0];
			s1.c[1] = inLink->cp[cItr].vt[vItr].c[1];
			s1.c[2] = inLink->cp[cItr].vt[vItr].c[2];
			s2.c[0] = inLink->cp[cItr].vt[vItr+1].c[0];
			s2.c[1] = inLink->cp[cItr].vt[vItr+1].c[1];
			s2.c[2] = inLink->cp[cItr].vt[vItr+1].c[2];
			sides[vItr].c[0] = s2.c[0] - s1.c[0];
			sides[vItr].c[1] = s2.c[1] - s1.c[1];
			sides[vItr].c[2] = s2.c[2] - s1.c[2];
			lengths[vItr] = plc_M_norm(sides[vItr]);
		}
		
		left = lengths[0];
		used = 0;
		tmpgoal = goal;
		
		// commence runaround!
		while( i<edges )
		{
			while( left > tmpgoal )
			{
				fixed->cp[cItr].vt[j].c[0] = ((tmpgoal+used)/lengths[i]) * sides[i].c[0] + inLink->cp[cItr].vt[i].c[0];
				fixed->cp[cItr].vt[j].c[1] = ((tmpgoal+used)/lengths[i]) * sides[i].c[1] + inLink->cp[cItr].vt[i].c[1];
				fixed->cp[cItr].vt[j].c[2] = ((tmpgoal+used)/lengths[i]) * sides[i].c[2] + inLink->cp[cItr].vt[i].c[2];
				
				left -= tmpgoal;
				used += tmpgoal;
				tmpgoal = goal;
				j++;
			}
			
			tmpgoal -= left;
			i++;
			left = lengths[i];
			used = 0;
		}
		
		free(lengths);
		free(sides);
	}
	
	return fixed;
}

double plc_strand_length(plc_strand *P) 

     /* Procedure computes the length of P. Assumes wraparound addressing. */
     
{
  int i, edges;
  double len = 0;
  plc_vector V;

  edges = plc_strand_edges(P);

  for(i=0;i<edges;i++) {

    V = P->vt[i];
    plc_M_sub_vect(V,P->vt[i+1]);      
    len += plc_M_norm(V);

  }

  return len;
}
  


double plCurve_short_edge(plCurve *L)

     /* Procedure finds the length of the shortest edge. */

{
  int i,comp;
  double minLen = {DBL_MAX};
  plc_vector diff;

  for(comp=0;comp < L->nc;comp++) {

    for(i=0;i < L->cp[comp].nv - (L->cp[comp].acyclic ? 1:0);i++) {

      diff = L->cp[comp].vt[i];
      plc_M_sub_vect(diff,L->cp[comp].vt[i+1]);
      
      minLen =  (plc_M_norm(diff) < minLen) ? plc_M_norm(diff) : minLen;

    }

  }

  return minLen;
}

double plCurve_long_edge(plCurve *L)

     /* Procedure finds the length of the shortest edge. */

{
  int i,comp;
  double maxLen = {0};
  plc_vector diff;

  for(comp=0;comp < L->nc;comp++) {

    for(i=0;i < L->cp[comp].nv - (L->cp[comp].acyclic ? 1:0);i++) {

      diff = L->cp[comp].vt[i];
      plc_M_sub_vect(diff,L->cp[comp].vt[i+1]);
      
      maxLen =  (plc_M_norm(diff) > maxLen) ? plc_M_norm(diff) : maxLen;

    }

  }

  return maxLen;
}

void plCurve_draw(FILE *outfile, plCurve *L) 

     /* Procedure draws the link in Geomview, complete
	with red/blue alternating vertices and little spheres 
	centered on the vertices. */

{
  int i, j;
  double ballrad;

  /* We first write the polyline to the file. */

  fprintf(outfile, "LIST \n\n");
  
  fprintf(outfile, "{ = "); 
  plc_write(outfile,L); 
  fprintf(outfile, " } \n\n");

  /* Now we get ready to write the vertices. */

  ballrad = plCurve_short_edge(L)/10.0;

  /* And write them to the file, too. */

  for(i=0;i<L->nc;i++) {

    for(j=0;j<L->cp[i].nv;j++) {

      fprintf(outfile,"{ = SPHERE %g  %g %g %g }\n",ballrad, 
	      L->cp[i].vt[j].c[0],
	      L->cp[i].vt[j].c[1],
	      L->cp[i].vt[j].c[2]);

    }

  }

}

plc_vector plCurve_edge_dir(plCurve *L, int comp, int edge)

     /* Procedure returns a unit vector in the (forward) direction along
	the given edge. Assumes wraparound addressing. */

{
  plc_vector V;
  double nV;

  /* Check sanity */

  if (edge >= plc_strand_edges(&L->cp[comp]) || (edge < 0 && L->cp[comp].acyclic)) {

    fprintf(stderr,"plCurve_edge_dir: Illegal edge number %d.\n",edge);
    exit(2);

  }

  /* Now do work. */

  V = L->cp[comp].vt[edge+1];
  plc_M_sub_vect(V,L->cp[comp].vt[edge]);
  nV = plc_M_norm(V);
  plc_M_scale_vect(1/nV,V);

  return V;
}
  
plc_vector plCurve_tangent_vector(plCurve *L, int comp, int vert) 

     /* Procedure computes a (unit) tangent vector at <vert>. 
	Handles closed and open curves correctly. */

{
  plc_vector V,Tleft,Tright;
  double nV;

  /* We need some special-case code at acyclic endpoints. */

  if (L->cp[comp].acyclic) {

    if (vert == 0) return plCurve_edge_dir(L,comp,0);
    if (vert == L->cp[comp].nv-1) return plCurve_edge_dir(L,comp,vert-1);

  }

  Tleft = plCurve_edge_dir(L,comp,vert-1);
  Tright = plCurve_edge_dir(L,comp,vert);

  plc_M_vweighted(V,0.5,Tleft,Tright);
  nV = plc_M_norm(V);
  plc_M_scale_vect(1/nV,V);

  return V;
}

    

    
  
