/*

   eqedge.c : A linear algebra solver demo showing dynamic edge 
              equilateralization using LAPACK and TAUCS. The program
	      is a bit of a monster, but is designed to check the 
	      TAUCS interface code against LAPACK.

	      We must link with octrope to get the link reading and 
	      simple vector arithmetic required.

*/

#include"eqedge.h"

#include <vecLib/vBLAS.h>
#include <vecLib/clapack.h>


void    clapack_matrix_write(double val, double *A, int LDA, int i, int j)

     /* Function writes val to position (i,j) in the LDAxN lapack matrix A. */

{
  A[j*LDA + i] = val;
}

double  clapack_matrix_read(double *A, int LDA, int i, int j)

     /* Function returns a value from position (i,j) in LDA x N lapack matrix A. */

{
  return A[j*LDA + i];
}

void    clapack_matrix_print(double *A, int M, int N) 

     /* Function prints the M x N lapack matrix A. */

{
  int i, j;

  for(i=0;i<M;i++) {
  
    for(j=0;j<N;j++) {

      printf("%3.2g ",clapack_matrix_read(A,M,i,j));

    }

    printf("\n");

  }

  printf("\n");

}

double *tangential_rigidity_matrix(octrope_link *L, int comp, int *LDA, int *N)

     /* Function returns a pointer to an LDAxN lapack-style matrix which encodes
	the change in edge lengths of component <comp> under the tangential variation 
	of its vertices. */

{

  int i;
  double *A;
  octrope_vector LeftTan,RightTan,EdgeDir;

  /* 
     
     We first allocate the matrix A. The dimensions of A are:

        <# of edges of <comp>> x <#of vertices of <comp>>  

     The number of edges is equal to the number of vertices if
     the polyline is closed (acylic == FALSE) and is one less than
     the number of vertices otherwise.

  */

  *LDA = octrope_pline_edges(&L->cp[comp]);
  *N   = L->cp[comp].nv;

  A = calloc( (*LDA) * (*N), sizeof(double) );

  /* We now fill the matrix. */

  for(i=0;i<L->cp[comp].nv-1;i++) {

    /* The i,j th entry is the effect on the length of edge i of moving vertex j. */
    /* This is nonzero only when j is an endpoint of the edge, so we first collect */
    /* the tangent vectors at the ends of the edge. */

    LeftTan = octrope_link_tangent_vector(L,comp,i);
    RightTan = octrope_link_tangent_vector(L,comp,i+1);
    EdgeDir = octrope_link_edge_dir(L,comp,i);

    clapack_matrix_write( -octrope_dot(EdgeDir,LeftTan) ,A,*LDA,i,i  );
    clapack_matrix_write( octrope_dot(EdgeDir,RightTan),A,*LDA,i,i+1);

  }

  /* If we have a wraparound edge, it has to be handled with 
     special care. */

  
  if (!L->cp[comp].acyclic) {

    LeftTan = octrope_link_tangent_vector(L,comp,i);
    RightTan = octrope_link_tangent_vector(L,comp,0);
    EdgeDir = octrope_link_edge_dir(L,comp,i);

    clapack_matrix_write( -octrope_dot(EdgeDir,LeftTan) ,A,*LDA,i,i);
    clapack_matrix_write( octrope_dot(EdgeDir,RightTan),A,*LDA,i,0);

  } 

  if (VERBOSITY > 4) {

    printf("tangential_rigidity_matrix. Constructed matrix: \n\n");
    clapack_matrix_print(A,*LDA,*N);

  }

  return A;
  
}

double *edgelength_error(octrope_link *L,int comp, double target_length)

     /* Procedure computes the difference between each edge length and
	the "target" value of octrope_pline_length/octrope_pline_edges,
	and returns it as a vector.  The argument LDB returns the number
        of rows in this vector. */

{
  double *Err;
  int    i,edges;

  octrope_vector V;

  edges = octrope_pline_edges(&(L->cp[comp]));
  Err = calloc(edges,sizeof(double));
  
  /* We now compute the raw edge lengths and total length. */

  for(i=0;i<edges;i++) {

    V = L->cp[comp].vt[i];
    octrope_vsub(V,L->cp[comp].vt[i+1]);      
    Err[i] = octrope_norm(V);
    
  }

  /* We now correct by subtracting the average (target) from each length. */

  for(i=0;i<edges;i++) {

    Err[i] -= target_length;

  }

  /* The remaining entries in "Err" are the edgelength errors. */

  return Err;

}

void lapack_eqedge(octrope_link *L, int nsteps)

     /* Procedure equilateralizes L using <nsteps> Newton steps.
	This procedure should converge exponentially in <nsteps>, with good stability
	as long as L is not far from equilateral to begin with. 

	The solution algorithm is as follows.

	For each component of L:

	1) Find the rigidity matrix A which transforms a tangential 
	   variation of the vertices of L->cp[comp] (in R^V) to a 
	   variation of the edge lengths of L->cp[comp] (in R^V or R^{V-1}). 

        2) Find the errors Err in edge lengths of L->cp[comp].

	3) Solve the (underdetermined) least-squares problem 

	       Ax = -Err,

	   using the LAPACK routine DGELSS to find the minimum norm
	   variation x which eliminates edge-length errors to first 
	   order. 

        4) Step in that direction.

	5) Repeat <nsteps> times.

     */
	
{

  __CLPK_doublereal *A, *Err, *S;
  __CLPK_integer     LDA, M,N, NRHS, LDB, RANK, LWORK, INFO;
  __CLPK_doublereal *WORK, RCOND, *B;
  int    step, comp;
  double target_length;

  octrope_vector stepV;

  int    i;
  
  // we pretend

  if (VERBOSITY > 2) {		/* Debugging code records function entry. */

    fprintf(stderr,"lapack_eqedge called. \n\n");
    fprintf(stderr,"  %d component link.\n",L->nc);

  }

  if (GRAPHICS) {

    refresh_display(L);
    i = 1;
  }

  for(comp = 0;comp < L->nc; comp++) {

    target_length = octrope_pline_length(&L->cp[comp])/(double)(octrope_pline_edges(&L->cp[comp]));

    for(step = 0; step < nsteps; step++) {

		  if (VERBOSITY > 3) {	/* Track maximum Errors before step. */

			fprintf(stderr,"    %d vertex component (# %d)\n",L->cp[comp].nv,comp);
			fprintf(stderr,"    maxErr before step %d = %g.\n",step,maxError(L,comp,target_length));
			fprintf(stderr,"    max/min before step %d = %g.\n\n",step,octrope_link_long_edge(L)/octrope_link_short_edge(L));

		  }

		  /* 1) Find the rigidity matrix. */

		  A = (__CLPK_doublereal *)(tangential_rigidity_matrix(L,comp,(int *)(&LDA),(int *)(&N)));

		  /* 2) Find the edge length errors. */

		  Err = (__CLPK_doublereal *)(edgelength_error(L,comp,target_length));

		  /* 3) Solve the least-squares problem */

		  LWORK = 20*LDA;
		  WORK = malloc(LWORK*sizeof(double));
		  M = LDA; NRHS = 1; 
		  RCOND = -1.0; /* Use machine precision */
		  S = malloc(M*sizeof(double));
		  
		  /* We now set up the RHS B. Paradoxically, though B is a vector of length M (# edges),
		 we must allocate a vector of length _N_ (3 * # verts) because this space will be 
		 overwritten by the least-squares _solution_ to the problem.

		 We similarly set LDB = N even though there are only M nonzero elements when we 
		 pass B to dgelss. */

		  B = calloc(N,sizeof(double));

		  for(i=0;i<M;i++) {

				B[i] = Err[i];

		  } 

		  LDB = N;

		  /* We now call dgelss_, the underdetermined least-squares solver. */
		  		  
		  dgelss_(&M, &N, &NRHS, A, &LDA, B, &LDB, S, &RCOND, &RANK, WORK, &LWORK, &INFO);

		  /* We now check INFO to make sure the method converged. */

		  if (INFO != 0) {

			if (INFO < 0) {

			  fprintf(stderr,"lapack_eqedge: cDGELSS reports illegal %d th argument.\n",(int)(-INFO));
			  exit(1);

			} else {

			  fprintf(stderr,"lapack_eqedge: cDGELSS couldn't find SVD for rigidity matrix.\n");
			  exit(1);

			}

		  }

		  /* 4) Take a Newton step. */

		  for(i=0;i<L->cp[comp].nv;i++) {

			stepV = octrope_link_tangent_vector(L,comp,i);
			octrope_vsmult(B[i],stepV);
			octrope_vsub((L->cp[comp].vt[i]),stepV);

		  }

      octrope_link_fix_wrap(L);

      if (VERBOSITY > 3) { /* Report on the step we just took. */

		double maxStep;

		for(maxStep = 0,i=0;i<LDB;i++) {

		  maxStep = (maxStep > fabs((double)(B[i])) ? maxStep : fabs((double)(B[i])));

		}
				 
		fprintf(stderr,"   Maximum coordinate movement %g.\n",maxStep);
		fprintf(stderr,"   Maximum error after step    %g.\n",maxError(L,comp,target_length));
		fprintf(stderr,"   max/min after step          %g.\n\n",octrope_link_long_edge(L)/octrope_link_short_edge(L));

      }

      if (GRAPHICS) {

		refresh_display(L);

      }

      /* 5) Free the memory that we allocated. */

      free(WORK); free(S); free(A); free(Err); free(B);

    }

  } /* End of main stepping loop. */

}

double maxError(octrope_link *L, int comp, double target_length) 

     /* Computes the maximum error in edge lengths. */

{
  int i;
  double maxE = 0, *Err;
  
  Err = edgelength_error(L,comp,target_length);
  
  for(i=0;i<octrope_pline_edges(&L->cp[comp]);i++) {

    maxE = (fabs(Err[i]) > maxE) ? fabs(Err[i]) : maxE;

  }

  free(Err);
  return maxE;

}
/*
int main()

{
 
  gclpipe = fopen("/tmp/OOGL","w");

  fprintf(gclpipe,"(normalization g0 none)\n");
  fprintf(gclpipe,"(bbox-draw g0 no)\n");

  VERBOSITY = 4;
  GRAPHICS = TRUE;

  segment_demo();
  ellipse_demo();

  fclose(gclpipe);


	return 0;
}
  */
  
  

  
