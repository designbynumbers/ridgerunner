/*
   Functions to be included in the plCurve library.
*/

#include "ridgerunner.h"

double
plCurve_torsion( plCurve* inLink, FILE* outPlot )
     
/* Constructs a GNUplot plot of an approximation to the torsion
   of the curve. Don't take the results too seriously-- a better job
   could be done of estimating torsion using a better 5-point formula,
   and in any case this kind of computation is really sensitive to
   numerical noise in the data. */

{ 
  double totalTorsion = 0; 
  plc_vector cross1 = {{0,0,0}},cross2 = {{0,0,0}}; 
  double temp; int cItr, vItr;
			     
  int	      e1, e2, e3;	
  double      lengthOffset = 0;
  int	      vertOffset = 0;
  int         *compedge = malloc(sizeof(int)*inLink->nc);
  
  totalTorsion = 0;
  
  for( cItr=0; cItr<inLink->nc; cItr++ )    {
    for( vItr=0; vItr<compedge[cItr]; vItr++) {

      plc_vector s1, s2, s3;
      
      e1 = vItr-1;
      e2 = vItr;
      e3 = vItr+1;
      
      s1 = plc_vect_diff( inLink->cp[cItr].vt[e1+1], inLink->cp[cItr].vt[e1] );
      s2 = plc_vect_diff( inLink->cp[cItr].vt[e2+1], inLink->cp[cItr].vt[e2] );
      s3 = plc_vect_diff( inLink->cp[cItr].vt[(e3+1)%compedge[cItr]], 
			  inLink->cp[cItr].vt[e3] );
      
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
	  
	  fprintf(outPlot, "%d %3.8lf %3.8lf %3.8lf\n", 
		  vertOffset, lengthOffset, temp, 2*tan(temp*.5)/plc_M_norm(s1));
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

  free(compedge);
  return totalTorsion;
}

plCurve*
plCurve_fixresolution( plCurve* inLink, double newres)
{
  /* 
     Spline the curve so that the resolution is approximately <newres> 
     vertices per unit ropelength. Uses the global gLambda to compute the 
     thickness of the curve. 

  /* The procedure is nondestructive to the original curve. */

  plCurve*   newver;
  plc_spline *splinever;
  bool       ok;
  int        i,*nv = malloc(sizeof(int)*inLink->nc);
  double     thickness;
  double     *length = malloc(sizeof(double)*inLink->nc);

  splinever = plc_convert_to_spline(inLink,&ok);

  if (!ok) { /* Die if we can't complete the spline process */
    
    error_write(inLink);
    fatal_(kSplineProblem);

  }

  /* Now compute the number of vertices needed for the new curve. */

  thickness = octrope_thickness(inLink,NULL,0,gLambda);
  plc_arclength(inLink,length);
  for(i=0;i<inLink->nc;i++) {nv[i] = ceil(newres*length[i]/thickness);}   

  /* Generate the new curve */

  newver = plc_convert_from_spline(splinever,nv);

  /* Free memory and return. */

  free(nv);
  free(length);
  plc_spline_free(splinever);
  return newver;

}

/* 
 * The following proceedure steps around the polygon redistributing
 * vertices as necessary to make things approach equilateral. This
 * code adapted from Eric Rawdon's TOROS polygonal runaround -- not to
 * be confused with the spherical runaround or the tangential
 * stepper. Contacts: rawdon@mathcs.duq.edu piatek@mathcs.duq.edu
 */
plCurve*
octrope_fixlength( plCurve* inLink )
{
  plCurve* fixed;
  int cItr, vItr;
  double *complength = malloc(sizeof(double)*inLink->nc);
  int    *compedge = malloc(sizeof(int)*inLink->nc);
  
  fixed = plc_copy(inLink);
  plc_fix_wrap(inLink);     /* This shouldn't be needed, but... */ 

  plc_arclength(inLink,complength); /* Compute arclength for components */
  plc_edges(inLink,compedge);
  
  for( cItr=0; cItr<fixed->nc; cItr++ ) {

    double goal, used, tmpgoal, left;
    int i=0, j=1, edges;
    double* lengths;
    plc_vector* sides;
    
    edges = compedge[cItr];
    goal = complength[cItr]/(double)(edges);
    
    /* plc_strand_length(&inLink->cp[cItr]) / (double)edges; */
    
    // we need the edge lengths and sides
    lengths = (double*)malloc(sizeof(double)*edges);
    sides = (plc_vector*)malloc(sizeof(plc_vector)*edges);
    
    for( vItr=0; vItr<edges; vItr++ ) {

      /* plc_vector s1, s2; */
      
      sides[vItr] = plc_vect_diff(inLink->cp[cItr].vt[vItr+1],
				  inLink->cp[cItr].vt[vItr]);
      
      /*s1.c[0] = inLink->cp[cItr].vt[vItr].c[0];
	s1.c[1] = inLink->cp[cItr].vt[vItr].c[1];
	s1.c[2] = inLink->cp[cItr].vt[vItr].c[2];
	
	s2.c[0] = inLink->cp[cItr].vt[vItr+1].c[0];
	s2.c[1] = inLink->cp[cItr].vt[vItr+1].c[1];
	s2.c[2] = inLink->cp[cItr].vt[vItr+1].c[2];
	
	sides[vItr].c[0] = s2.c[0] - s1.c[0];
	sides[vItr].c[1] = s2.c[1] - s1.c[1];
	sides[vItr].c[2] = s2.c[2] - s1.c[2];*/
      
      lengths[vItr] = plc_M_norm(sides[vItr]);
    }
    
    left = lengths[0];
    used = 0;
    tmpgoal = goal;
    
    // commence runaround!
    while( i<edges ) {
      while( left > tmpgoal ) {
	fixed->cp[cItr].vt[j] = plc_vlincomb((tmpgoal+used)/lengths[i],sides[i],
					     1.0,inLink->cp[cItr].vt[i]);
	
	/*
	  fixed->cp[cItr].vt[j].c[0] = 
	  ((tmpgoal+used)/lengths[i]) * sides[i].c[0] + inLink->cp[cItr].vt[i].c[0];
	  
	  fixed->cp[cItr].vt[j].c[1] = 
	  ((tmpgoal+used)/lengths[i]) * sides[i].c[1] + inLink->cp[cItr].vt[i].c[1];
	  
	  fixed->cp[cItr].vt[j].c[2] = 
	  ((tmpgoal+used)/lengths[i]) * sides[i].c[2] + inLink->cp[cItr].vt[i].c[2];
	*/
	
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
  
  free(complength);
  free(compedge);
  return fixed;
}

double plCurve_short_edge(plCurve *L)

     /* Procedure finds the length of the shortest edge. */
     
{
  int i,comp;
  double minLen = {DBL_MAX};
  plc_vector diff;
  int *compedge = malloc(sizeof(int)*L->nc);

  plc_edges(L,compedge);
  
  for(comp=0;comp < L->nc;comp++) {

    for(i=0;i < compedge[comp];i++) {

      diff = L->cp[comp].vt[i];
      plc_M_sub_vect(diff,L->cp[comp].vt[i+1]);
      
      minLen =  (plc_M_norm(diff) < minLen) ? plc_M_norm(diff) : minLen;
      
    }

  }

  free(compedge);
  return minLen;
}

double plCurve_long_edge(plCurve *L)

     /* Procedure finds the length of the shortest edge. */

{
  int i,comp;
  double maxLen = {0};
  plc_vector diff;
  int *compedge = malloc(sizeof(int)*L->nc);

  plc_edges(L,compedge);

  for(comp=0;comp < L->nc;comp++) {

    for(i=0;i < compedge[comp];i++) {

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

  if (edge > L->cp[comp].nv || (edge < 0 && L->cp[comp].open)) {

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
  

    

    
  
