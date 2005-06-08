/*
   Functions to be included in the octrope_link library.
*/

#include "gradient/eqedge.h"
#include "gradient/stepper.h"

int    octrope_pline_edges(octrope_pline *P)

     /* Procedure computes the number of edges in a pline. */
     /* This is either equal to, or one less than, the number */
     /* of vertices, depending on whether the pline is closed. */

{

  return P->nv - (P->acyclic ? 1:0);

}

void
link_scale( octrope_link* inLink, double factor )
{
	int cItr, vItr;
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		for( vItr=0; vItr<inLink->cp[cItr].nv; vItr++ )
		{
			inLink->cp[cItr].vt[vItr].c[0] *= factor;
			inLink->cp[cItr].vt[vItr].c[1] *= factor;
			inLink->cp[cItr].vt[vItr].c[2] *= factor;
		}
	}
} 

octrope_link*
octrope_equalize_density( octrope_link* inLink, search_state* inState )
{
/*	octrope_link*   doubled;
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
	
	doubled = octrope_link_new( inLink->nc,
						newVerts,
						cyclicity, 
						color_count );
	
	octrope_link_fix_wrap(inLink);
						
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		// preserve color
		//memcpy( &doubled->cp[cItr].clr[0], &inLink->cp[cItr].clr[0], sizeof(octrope_color) );
	
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

octrope_link*
octrope_double_component( octrope_link* inLink, int comp )
{
	/* double number of sides through bisection -- keeps things eq and shouldn't munger too much with
	 * struts
	 */
	octrope_link*   doubled;
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
	
	doubled = octrope_link_new( inLink->nc,
						newVerts,
						cyclicity, 
						color_count );
	
	octrope_link_fix_wrap(inLink);

	cItr = comp;
//	for( cItr=inDo; cItr<inLink->nc; cItr++ )
	{
		// preserve color
		//memcpy( &doubled->cp[cItr].clr[0], &inLink->cp[cItr].clr[0], sizeof(octrope_color) );
	
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

octrope_link*
octrope_double_edges( octrope_link* inLink )
{
	/* double number of sides through bisection -- keeps things eq and shouldn't munger too much with
	 * struts
	 */
	octrope_link*   doubled;
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
	
	doubled = octrope_link_new( inLink->nc,
						newVerts,
						cyclicity, 
						color_count );
	
	octrope_link_fix_wrap(inLink);
						
	for( cItr=0; cItr<inLink->nc; cItr++ )
	{
		// preserve color
		//memcpy( &doubled->cp[cItr].clr[0], &inLink->cp[cItr].clr[0], sizeof(octrope_color) );
	
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
octrope_link*
octrope_fixlength( octrope_link* inLink )
{
	octrope_link* fixed;
	int cItr, vItr;
	
	fixed = octrope_link_copy(inLink);
	
	octrope_link_fix_wrap(inLink);
	
	for( cItr=0; cItr<fixed->nc; cItr++ )
	{
		double goal, used, tmpgoal, left;
		int i=0, j=1, edges;
		double* lengths;
		octrope_vector* sides;
		
		edges = octrope_pline_edges(&inLink->cp[cItr]);
		goal = octrope_pline_length(&inLink->cp[cItr]) / (double)edges;
		
		// we need the edge lengths and sides
		lengths = (double*)malloc(sizeof(double)*edges);
		sides = (octrope_vector*)malloc(sizeof(octrope_vector)*edges);
		for( vItr=0; vItr<edges; vItr++ )
		{
			octrope_vector s1, s2;
			s1.c[0] = inLink->cp[cItr].vt[vItr].c[0];
			s1.c[1] = inLink->cp[cItr].vt[vItr].c[1];
			s1.c[2] = inLink->cp[cItr].vt[vItr].c[2];
			s2.c[0] = inLink->cp[cItr].vt[vItr+1].c[0];
			s2.c[1] = inLink->cp[cItr].vt[vItr+1].c[1];
			s2.c[2] = inLink->cp[cItr].vt[vItr+1].c[2];
			sides[vItr].c[0] = s2.c[0] - s1.c[0];
			sides[vItr].c[1] = s2.c[1] - s1.c[1];
			sides[vItr].c[2] = s2.c[2] - s1.c[2];
			lengths[vItr] = octrope_norm(sides[vItr]);
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

double octrope_pline_length(octrope_pline *P) 

     /* Procedure computes the length of P. Assumes wraparound addressing. */
     
{
  int i, edges;
  double len = 0;
  octrope_vector V;

  edges = octrope_pline_edges(P);

  for(i=0;i<edges;i++) {

    V = P->vt[i];
    octrope_vsub(V,P->vt[i+1]);      
    len += octrope_norm(V);

  }

  return len;
}
  


double octrope_link_short_edge(octrope_link *L)

     /* Procedure finds the length of the shortest edge. */

{
  int i,comp;
  double minLen = {DBL_MAX};
  octrope_vector diff;

  for(comp=0;comp < L->nc;comp++) {

    for(i=0;i < L->cp[comp].nv - (L->cp[comp].acyclic ? 1:0);i++) {

      diff = L->cp[comp].vt[i];
      octrope_vsub(diff,L->cp[comp].vt[i+1]);
      
      minLen =  (octrope_norm(diff) < minLen) ? octrope_norm(diff) : minLen;

    }

  }

  return minLen;
}

double octrope_link_long_edge(octrope_link *L)

     /* Procedure finds the length of the shortest edge. */

{
  int i,comp;
  double maxLen = {0};
  octrope_vector diff;

  for(comp=0;comp < L->nc;comp++) {

    for(i=0;i < L->cp[comp].nv - (L->cp[comp].acyclic ? 1:0);i++) {

      diff = L->cp[comp].vt[i];
      octrope_vsub(diff,L->cp[comp].vt[i+1]);
      
      maxLen =  (octrope_norm(diff) > maxLen) ? octrope_norm(diff) : maxLen;

    }

  }

  return maxLen;
}

void octrope_link_draw(FILE *outfile, octrope_link *L) 

     /* Procedure draws the link in Geomview, complete
	with red/blue alternating vertices and little spheres 
	centered on the vertices. */

{
  int i, j;
  double ballrad;

  /* We first write the polyline to the file. */

  fprintf(outfile, "LIST \n\n");
  
  fprintf(outfile, "{ = "); 
  octrope_link_write(outfile,L); 
  fprintf(outfile, " } \n\n");

  /* Now we get ready to write the vertices. */

  ballrad = octrope_link_short_edge(L)/10.0;

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

octrope_vector octrope_link_edge_dir(octrope_link *L, int comp, int edge)

     /* Procedure returns a unit vector in the (forward) direction along
	the given edge. Assumes wraparound addressing. */

{
  octrope_vector V;
  double nV;

  /* Check sanity */

  if (edge >= octrope_pline_edges(&L->cp[comp]) || (edge < 0 && L->cp[comp].acyclic)) {

    fprintf(stderr,"octrope_link_edge_dir: Illegal edge number %d.\n",edge);
    exit(2);

  }

  /* Now do work. */

  V = L->cp[comp].vt[edge+1];
  octrope_vsub(V,L->cp[comp].vt[edge]);
  nV = octrope_norm(V);
  octrope_vsmult(1/nV,V);

  return V;
}
  
octrope_vector octrope_link_tangent_vector(octrope_link *L, int comp, int vert) 

     /* Procedure computes a (unit) tangent vector at <vert>. 
	Handles closed and open curves correctly. */

{
  octrope_vector V,Tleft,Tright;
  double nV;

  /* We need some special-case code at acyclic endpoints. */

  if (L->cp[comp].acyclic) {

    if (vert == 0) return octrope_link_edge_dir(L,comp,0);
    if (vert == L->cp[comp].nv-1) return octrope_link_edge_dir(L,comp,vert-1);

  }

  Tleft = octrope_link_edge_dir(L,comp,vert-1);
  Tright = octrope_link_edge_dir(L,comp,vert);

  octrope_vweighted(V,0.5,Tleft,Tright);
  nV = octrope_norm(V);
  octrope_vsmult(1/nV,V);

  return V;
}

    

    
  
