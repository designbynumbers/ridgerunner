
/*
 * This program is free software distributed under the GPL. A copy of the license should have been included with this 
 * archive in a file named 'LICENSE'. You can read the license there or on the web at: http://www.gnu.org/licenses/gpl.txt
 */

#include "lsqr.h"
#include "tsnnls.h"

#include <string.h>
#include <float.h>

#ifdef __APPLE__
#include <vecLib/vBLAS.h>
#include <vecLib/clapack.h>
#else
#include "gsl_cblas.h"
#endif // __APPLE__


/*

    File contains an implementation of the block-pivoting non-negative least squares
    algorithm of Portugal, Judice, and Vicente. The reference for the algorithm is their
    paper:

    @article {MR95a:90059,
    AUTHOR = {Portugal, Lu{\'{\i}}s F. and J{\'u}dice, Joaqu{\'{\i}}m J. and
              Vicente, Lu{\'{\i}}s N.},
     TITLE = {A comparison of block pivoting and interior-point algorithms
              for linear least squares problems with nonnegative variables},
   JOURNAL = {Math. Comp.},
  FJOURNAL = {Mathematics of Computation},
    VOLUME = {63},
      YEAR = {1994},
    NUMBER = {208},
     PAGES = {625--643},
      ISSN = {0025-5718},
     CODEN = {MCMPAF},
   MRCLASS = {90C20 (65K05)},
  MRNUMBER = {95a:90059},
}

   The paper is stored on JSTOR.


   Cantarella/Piatek. 5/04.

*/

void sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod );


/* The entries of setA are indices into the array varsA. Procedure finds the 
 * negative entries of varsA which correspond to indices in A, and returns a 
 * list of them in <infeas> with size <sizeInfeas>. 
 * 
 * We assume that infeas is allocated, and is large enough to hold the result. 
 */
void infeasible( const int* setA, const double* varsA, const int sizeA, 
		 int* infeas, int* sizeInfeas )
{
	int i, pItr=0;
	for( i=0; i< sizeA; i++ ) 
	{
		if( varsA[setA[i]] < 0 ) 
		{
			infeas[pItr++] = setA[i];
		}
	}
	*sizeInfeas = pItr;
}


/* Removes the elements of setB from setA, assuming 
 * that both sets are sorted ascending. 
 */
void 
int_difference( int* setA, int aSize, const int* setB, const int bSize, int* outSize )
{
	int aRead = 0, aWrite = 0, bPtr = 0;

	/* Some special cases up front */
	if( aSize == 0 && bSize == 0 )
	{
		*outSize = 0;
		return;
	}

	if( bSize == 0 )
	{
		*outSize = aSize;
		return;
	}

	if( aSize == 0 && bSize != 0 )
	{
		fprintf( stderr, "int_differe(): aSize is 0 but bSize is not!\n" );
		exit(-1);
	}

	*outSize = 0;

	for( aRead=aWrite=0; aRead<aSize; aRead++ ) 
	{
		/* We should update bPtr, if we can. */
		if (setA[aRead] > setB[bPtr] && bPtr < bSize-1) 
		{ 
			bPtr++;
		}

		/* We don't have a match. */
		if (setA[aRead] != setB[bPtr]) 
		{ 
			setA[aWrite++] = setA[aRead];
			(*outSize)++;
		}
	}
}


/*
 * performs union, assumes both setA, setB are !sorted in ascending order!.
 * writes the result to set A.
 */
void 
int_union( int* setA, int aSize, const int* setB, const int bSize, int* outSize )
{
	int aItr, bItr, uItr;
	aItr = bItr = uItr = 0;
	int i;
	int*	unionSet = NULL;

	*outSize = 0;

	if( aSize == 0 && bSize == 0 )
		return;

	unionSet = (int*)malloc(sizeof(int)*(aSize+bSize));
	
	if( aSize == 0 )
	{
		memcpy( unionSet, setB, sizeof(int)*bSize);
		*outSize = bSize;
	}
	else if( bSize == 0 )
	{
		memcpy( unionSet, setA, sizeof(int)*aSize);
		*outSize = aSize;
	}
	else
	{
		while( 1==1 )
		{
			if( setA[aItr] == setB[bItr] )
			{
				unionSet[uItr++] = setA[aItr];
				aItr++;
				bItr++;
			}
			else
			{
				if(setA[aItr] < setB[bItr] )
					unionSet[uItr++] = setA[aItr++];
				else
				{
					unionSet[uItr++] = setB[bItr++];
				}
			}

			if( aItr == aSize )
			{
				/* copy remaining b elements */
				for( ; bItr<bSize; bItr++ )
					unionSet[uItr++] = setB[bItr];
				break;
			}
			if( bItr == bSize )
			{
				/* copy remaining a elements */
				for( ; aItr<aSize; aItr++ )
					unionSet[uItr++] = setA[aItr];
				break;
			}
		} /* while */

		*outSize = uItr;
	} // else 
	
	/* Now we recopy the results to A. */
	for(i=0;i<*outSize;i++) 
	{
		setA[i] = unionSet[i];
	}

	free(unionSet);
}

#ifndef __DBL_EPSILON__
#define __DBL_EPSILON__ 2.2204460492503131e-16
#endif

static void
fix_zeros( double* values, int size, double rcond, int inPrintErrorWarnings )
{
	int i;
	double eps = __DBL_EPSILON__; // as defined in float.h
	for( i=0; i<size; i++ )
	{
		if( inPrintErrorWarnings != 0 )
		{
			double cond = (double)1/rcond;
						
			if( fabs(values[i]) < (cond*cond*eps) )
			{
				fprintf( stderr, "Variable is within (condition number)^2*eps of 0, accuracy results may be unexpected!\n" );
				inPrintErrorWarnings = 0; // It should suffice to notify the user once. Otherwise this will happen MANY times
			}
		}
		
		if( fabs(values[i]) < kZeroThreshold )
		{
			values[i] = 0.0;
		}
	}
}

taucs_ccs_matrix*
selectAprimeDotA( double* apda, int cols, int* F, int sizeF )
{
	taucs_ccs_matrix* result = NULL;
	int maxSize, i, j, valItr;
	
	result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));

	result->n = sizeF;
	result->flags = TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER;
	
	maxSize = (result->n*(result->n+1))/2;

	result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
	result->rowind = (int*)malloc(sizeof(int)*maxSize);
	result->values.d = (double*)malloc(sizeof(taucs_double)*maxSize);
	
	valItr = 0;
	for( i=0; i<sizeF; i++ )
	{
		result->colptr[i] = valItr;
		for( j=i; j<sizeF; j++ )
		{
			double val = apda[cols*F[j] + F[i]];
			if( val != 0 )
			{
				result->values.d[valItr] = val;
				result->rowind[valItr] = j;
				valItr++;
			}
		}
	}
	result->colptr[i] = valItr;

	return result;
}

void
selectAprimeDotAsparse( const taucs_ccs_matrix* apda, int* F, int sizeF, taucs_ccs_matrix* inOutApda )
{
	/* This routine presume that inOutApda has already been allocated (this save allocation and free time) 
	 * the size should be the same as apda, since this will be a submatrix, it cannot possibly be larger */
	
	int fItr, valItr, cItr, rItr;
	
	if( sizeF == 0 )
	{
		inOutApda->m = inOutApda->n = 0;
		return;
	}
	
	inOutApda->n = sizeF;
	inOutApda->flags = TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER;
	
	/* start in column F[0] and select entries that are in row F[0], F[1], F[2], ... etc */
	valItr = cItr = 0;
	for( cItr=0; cItr<sizeF; cItr++ )
	{
		inOutApda->colptr[cItr] = valItr;
		rItr=apda->colptr[F[cItr]];
		/* scan down rows until we can compare to F value */
		fItr=0;
		
		/* Since this matrix is symmetric storing the lower triangle, we need to adjust fItr to start at the main diagonal
		 */
	//	while( F[fItr] < F[cItr] )
	//		fItr++;
		fItr=cItr;
			
		while( 1==1 )
		{
			/* have we exhausted F? */
			if( fItr == sizeF )
				break;
		
			/* Continues until we can compare to F's selection or we finish this column */
			while( apda->rowind[rItr] < F[fItr] && rItr<apda->colptr[F[cItr]+1] )
				rItr++;
				
			/* Are we dont with this column? */
			if( rItr == apda->colptr[F[cItr]+1] )
				break;
			
			/* We have been been skipping along the sparse entries of apda, which may or may not correspond to 
			 * row requests from F. Once we're out of the above loop, we are either AT F[fItr] in which case 
			 * we want this entry, or we have overshot it. If we've overshot it, we need to increase fItr.
			*/
			if( apda->rowind[rItr] == F[fItr] )
			{
				inOutApda->values.d[valItr] = apda->values.d[rItr];
				inOutApda->rowind[valItr] = fItr; /* fItr here is the row in our submatrix (since it is sizeF by sizeF) */
				valItr++;
				rItr++;
				fItr++; /* we got this one */
			}
			else
			{
				fItr++;
			}
		} // while row scanning
	} // for F column selections
	
	/* fix up the last column start, valItr is one more than the index of the last guy
	 * entered, which is what we want */
	inOutApda->colptr[sizeF] = valItr;
}

taucs_ccs_matrix*
selectAprime( double* Ap, int cols, int rows, int* F, int sizeF )
{
	taucs_ccs_matrix* result = NULL;
	int maxSize, i, j, valItr;
	
	result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));

	result->m = sizeF;
	result->n = cols;
	result->flags = TAUCS_DOUBLE;
	
	maxSize = result->m*result->n;

	result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
	result->rowind = (int*)malloc(sizeof(int)*maxSize);
	result->values.d = (double*)malloc(sizeof(taucs_double)*maxSize);
	
	valItr = 0;
	for( i=0; i<result->n; i++ )
	{
		result->colptr[i] = valItr;
		for( j=0; j<result->m; j++ )
		{
			double val = Ap[cols*F[j] + i];
			if( val != 0 )
			{
				result->values.d[valItr] = val;
				result->rowind[valItr] = j;
				valItr++;
			}
		}
	}
	result->colptr[i] = valItr;

	return result;
}

/* the nature of this routine is described in lsqr.h ~line 316, primarily:
 *
 *									If MODE = 0, compute  y = y + A*x,
 *									If MODE = 1, compute  x = x + A^T*y.
 */
void 
sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod )
{
	taucs_ccs_matrix* A = (taucs_ccs_matrix*)prod;
	int i;
	int cItr, rItr;
	
	if( mode == 0 )
	{
		int i,ip,j,n, rows;
		taucs_datatype Aij;

		n = A->n;
		rows = A->m;
	
	/*	double* Ax = (double*)malloc(sizeof(double)*A->m);
		ourtaucs_ccs_times_vec( A, x->elements, Ax );
		for(i=0; i<A->m; i++ )
			y->elements[i] += Ax[i];
		free(Ax);
	*/
	//	ourtaucs_ccs_times_vec( A, x->elements, y->elements );
	
		for (j=0; j<n; j++) {
		  for (ip = (A->colptr)[j]; ip < (A->colptr[j+1]); ip++) {
			i   = (A->rowind)[ip];
			Aij = (A->values.d)[ip];
			
			y->elements[i] = y->elements[i]+(x->elements[j]*Aij);
		  }
		}
	
	}
	else if( mode == 1 )
	{
		// since A^T*y = y^T*A
	/*	double* yA = (double*)malloc(sizeof(double)*A->n);
		taucs_transpose_vec_times_matrix_nosub( y->elements, A, yA );
		for(i=0; i<A->n; i++ )
			x->elements[i] += yA[i];
		free(yA);
	*/
	//	taucs_transpose_vec_times_matrix_nosub( y->elements, A, x->elements );
		for( cItr=0; cItr<A->n; cItr++ )
		{
			//result[cItr] = 0;
			for( rItr=0; rItr<A->colptr[cItr+1]-A->colptr[cItr]; rItr++ )
			{
				int tRow = A->rowind[A->colptr[cItr] + rItr];
				x->elements[cItr] += y->elements[tRow] * A->values.d[A->colptr[cItr] + rItr];
			}
		}

	}
	else
		fprintf(stderr, "Unknown mode: %ld\n", mode );
}

#ifndef __DBL_EPSILON__
#define __DBL_EPSILON__ 2.2204460492503131e-16
#endif

taucs_double*
t_snnls( taucs_ccs_matrix *A_original_ordering, taucs_double *b, 
		 double* outResidualNorm, double inRelErrTolerance, int inPrintErrorWarnings )
{
	taucs_ccs_matrix  *Af;
	int               p, ninf, pbar = {3};
	int               m,n,i, maxSize;
	
	int				A_rows, A_cols;

	int              *F, *G, *H1, *H2, SCR[1];
	int              sizeF, sizeG, sizeH1, sizeH2, sizeSCR = {1};
	int				lsqrStep=0;
	double			rcond=1;

	/* These variables are subsets of the column indices of the matrix A, 
	 * always stored in sorted order, and sized by the corresponding
	 * "size" integers. 
	 * 
	 * Like the row indices themselves, they are 0-based. 
	 */
	taucs_double     *x,*y, *xf_raw = NULL, *yg_raw, *residual;
  
	int AprimeDotA_cols;

	taucs_ccs_matrix* AprimeDotA = taucs_ccs_aprime_times_a(A_original_ordering);
	taucs_ccs_matrix*   lsqrApA;
	
	/* create a copy of AprimeDotA memory wise to store the tlsqr submatrices */
	{
		lsqrApA = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
		lsqrApA->n = AprimeDotA->n;
		lsqrApA->flags = TAUCS_DOUBLE;
		lsqrApA->flags = lsqrApA->flags | TAUCS_SYMMETRIC;
		lsqrApA->flags = lsqrApA->flags | TAUCS_LOWER; // rep the lower half
		lsqrApA->colptr = (int*)malloc(sizeof(int)*(lsqrApA->n+1));
		/* This is the number of nonzeros in A'*A, which we cannot overflow with a submatrix */
		maxSize = AprimeDotA->colptr[AprimeDotA->n]; 
		lsqrApA->values.d = (double*)malloc(sizeof(taucs_double)*maxSize);
		lsqrApA->rowind = (int*)malloc(sizeof(int)*maxSize);
	}
	
	if( inRelErrTolerance <= 0 )
		lsqrStep = 1;
	
	A_rows = A_original_ordering->m;
	A_cols = A_original_ordering->n;
		
	AprimeDotA_cols = A_cols;

	m = A_original_ordering->m;
	n = A_original_ordering->n;
	
	// This initial values is suggested by PJV
	ninf = n+1;
  
	/* We first allocate space. */
	F   = calloc(n,sizeof(int));
	G   = calloc(n,sizeof(int));
	H1  = calloc(n,sizeof(int));
	H2  = calloc(n,sizeof(int));

	x   = calloc(n,sizeof(taucs_double));
	y   = calloc(m,sizeof(taucs_double));

	/* submatrix allocation actually takes bit of time during profiling, so we reuse
	 * an allocation that cannot be overflowed by submatrices. Note that
	 * A_original_ordering->colptr[A_original_ordering->n] is the number of
	 * nonzero entries in A
	 */
	Af = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
	Af->colptr = (int*)malloc(sizeof(int)*(A_cols+1));
	Af->rowind = (int*)malloc(sizeof(int)*(A_original_ordering->colptr[A_original_ordering->n]));
	Af->values.d = (double*)malloc(sizeof(double)*A_original_ordering->colptr[A_original_ordering->n]);

	/* Next we initialize variables. */
	for(i=0;i<n;i++) 
	{
		G[i] = i;
	}
	sizeF = 0; sizeG = n;  
	p = pbar;

	/* Set y = A'b, which is the same as y=b'A. We perform that computation as it is faster */
//	transpose_vec_times_matrix(b, A, G, A_original_ordering->n, A_original_ordering->m, A_original_ordering->n, y);
	taucs_transpose_vec_times_matrix(b, A_original_ordering, G, A_cols, y);
  
	cblas_dscal(n,-1.0,y,1);      /* Scalar multiply y *= -1. */

  /* Now we enter the main loop. */
	infeasible(F,x,sizeF,H1,&sizeH1);  
	infeasible(G,y,sizeG,H2,&sizeH2);

	for( ; sizeH1 > 0 || sizeH2 > 0 ;
		infeasible(F,x,sizeF,H1,&sizeH1),  
		infeasible(G,y,sizeG,H2,&sizeH2)  ) 
	{
		

		/* We found infeasible variables. We're going to swap them 
		   between F and G according to one of the schemes below. */

		if (sizeH1 + sizeH2 < ninf) 
		{  
			/* The number of infeasibles _decreased_ since the last try. */
			/* This is good. We reset the counters and swap all infeasibles. */
			ninf = sizeH1 + sizeH2;
			p = pbar;

			int_union(F,sizeF,H2,sizeH2,&sizeF);
			int_difference(F,sizeF,H1,sizeH1,&sizeF);

			int_union(G,sizeG,H1,sizeH1,&sizeG);
			int_difference(G,sizeG,H2,sizeH2,&sizeG);
		} 
		else 
		{  
			if( p > 0 ) 
			{  
				/* Things didn't go well last time-- we _increased_ the number */
				/* of infeasibles. But we haven't run out of chances yet, so   */
				/* we'll decrement the counter, and try the large swap again.  */
				p--;

				int_union(F,sizeF,H2,sizeH2,&sizeF);
				int_difference(F,sizeF,H1,sizeH1,&sizeF);

				int_union(G,sizeG,H1,sizeH1,&sizeG);
				int_difference(G,sizeG,H2,sizeH2,&sizeG);
			} /* (p>0) */
			else 
			{
				/* Real trouble. Large swaps aren't reducing the number of     */
				/* infeasibles, even after pbar tries. So this time, we'll     */
				/* revert to the slower "Murty's Method", and move only one    */
				/* guy-- the guy with the highest column index. */
				
				if( sizeH1 > 0 && sizeH2 > 0 )
				{
					if ( H1[sizeH1 - 1] > H2[sizeH2 - 1] ) 
					{  
						/* H1 contains the last column index. */

						SCR[0] = H1[sizeH1-1];
						
						int_difference(F,sizeF,SCR,sizeSCR,&sizeF);
						int_union(G,sizeG,SCR,sizeSCR,&sizeG);
					} 
					else 
					{
						/* H2 contains the last column index. */

						SCR[0] = H2[sizeH2-1];
						
						int_union(F,sizeF,SCR,sizeSCR,&sizeF);
						int_difference(G,sizeG,SCR,sizeSCR,&sizeG);
					}
				}
				else if( sizeH1 == 0 )
				{
					SCR[0] = H2[sizeH2-1];
					
					int_union(F,sizeF,SCR,sizeSCR,&sizeF);
					int_difference(G,sizeG,SCR,sizeSCR,&sizeG);
				}
				else
				{
					SCR[0] = H1[sizeH1-1];

					int_difference(F,sizeF,SCR,sizeSCR,&sizeF);
					int_union(G,sizeG,SCR,sizeSCR,&sizeG);
				}		
				/* p is still 0, and will remain so until we start making */
				/* progress at reducing the number of infeasibilities. */
			} /* else (p>0) */
		} /* else (sizeH1 + sizeH2 < ninf) */

		/* We have now altered F and G, and are ready to recompute everything. */
		/* We first clear x and y, since we won't need them anymore. */
		for(i=0;i<n;i++) 
		{ 
			x[i] = y[i] = 0.0; 
		}
		
		/* 
		 *
		 * By (7) in PJV, x_F should be the solution of the unconstrained
		 * least squares problem:
		 * 
		 * min || A_F x_F - b ||,
		 * 
		 * where A_F is the submatrix of A containing the columns of A 
		 * with column indices in F. 
		 * 
		 * And by (8) in PJV, the y_F should be (A_G)'*(A_F x_F - b). 
		 * 
		 * We first solve the lsqr problem.
		 * 
		 */
		taucs_ccs_submatrix(A_original_ordering, F, sizeF, Af);
		
		if( sizeF != 0 )
		{
			/* we compute these values based on selections based on F since it's faster
			 * than recalculating them in lsqr. This requires the use of a custom lsqr 
			 * solver that expects extra parameters from snnls, however.
			 */
				
			selectAprimeDotAsparse(AprimeDotA, F, sizeF, lsqrApA); 	
						
			xf_raw = NULL;
			if( inRelErrTolerance > 1 || lsqrStep != 0 && inPrintErrorWarnings == 0 )
				xf_raw = t_snnlslsqr(Af, b, lsqrApA, F, NULL);		
			else
			{
				xf_raw = t_snnlslsqr(Af, b, lsqrApA, F, &rcond );
				if( (1/rcond)*(1/rcond)*__DBL_EPSILON__ < inRelErrTolerance )
					lsqrStep = 1;
			}
			if( xf_raw == NULL )
				return NULL; // matrix probably not positive definite
																		
			/* taucs_snnls requires us to handle in some way the case where values
			 * returned from lsqr that are within error tolerance of zero get continually
			 * swapped between x and y because our comparison to zero is numerically more 
			 * precise in the infeasibles computation. In order to terminate, we must 
			 * handle this case. For our usage of snnls, it is most efficient to simply
			 * zero values within some epsilon of zero as returned from lsqr, but this 
			 * solution may not be the best for all possible applications.
			 */
			fix_zeros(xf_raw, Af->n, rcond, inPrintErrorWarnings);
			
			/* Now compute the residual A_F x_F - b. This is an m-vector. */
			residual = (taucs_double *)calloc(m,sizeof(taucs_double));
			ourtaucs_ccs_times_vec(Af,xf_raw,residual);
		}
		else
		{	  
			/* 
			 * if sizeF is 0, the meaning of residual changes (since there really is _no_ matrix),
			 * so we'll just set the residual to -b, but we still need a zeroed residual to do
			 * the below computation to make that happen, which calloc does here
			 */
			residual = (taucs_double *)calloc(m,sizeof(taucs_double));
		}
		
		cblas_daxpy(m,-1.0,b, 1, residual, 1);
		
		/* We now compute (A_G)'. */
		/* And finally take (A_G)'*residual. This is a sizeG-vector. */
		yg_raw = (taucs_double *)calloc(sizeG,sizeof(taucs_double));     

		/* 
		 * We now should compute (A_G)', and take (A_G)'*residual, but 
		 * A_G'*residual = residual'*A_G, and that's faster. 
		 * taucs_transpose_vec_times_matrix also incorporates the 
		 * selection of columns of G from which to form A_G so that 
		 * we do not have to incur the computational expense of creating
		 * a submatrix.
		 */
		taucs_transpose_vec_times_matrix(residual, A_original_ordering, G, sizeG, yg_raw);
		
		fix_zeros(yg_raw, sizeG, rcond, inPrintErrorWarnings);

		/* We're now done. It remains to scatter the entries in xf_raw
		 * and yg_raw over the x and y vectors, and free everything that
		 * we can manage to free. 
		 */
		bzero(x, sizeof(double)*n);
		bzero(y, sizeof(double)*n);

		for(i=0;i<sizeF;i++) 
		{ 
			x[F[i]] = xf_raw[i]; 
		}
		for(i=0;i<sizeG;i++) 
		{ 
			y[G[i]] = yg_raw[i]; 
		}

		free(yg_raw);
		if( sizeF != 0 )
		{
			free(xf_raw);
		}
		free(residual);

		sizeH1 = sizeH2 = 0;
	} // for
	
	if( lsqrStep != 0 )
	{
		lsqr_input   *lsqr_in;
		lsqr_output  *lsqr_out;
		lsqr_work    *lsqr_work;
		lsqr_func    *lsqr_func;
		int bItr;
						
		alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, Af->m, Af->n );
		
		/* we let lsqr() itself handle the 0 values in this structure */
		lsqr_in->num_rows = Af->m;
		lsqr_in->num_cols = Af->n;
		lsqr_in->damp_val = 0;
		lsqr_in->rel_mat_err = 0;
		lsqr_in->rel_rhs_err = 0;
		lsqr_in->cond_lim = 1e16;
		lsqr_in->max_iter = lsqr_in->num_rows + lsqr_in->num_cols + 1000;
		lsqr_in->lsqr_fp_out = NULL;	
		for( bItr=0; bItr<Af->m; bItr++ )
		{
			lsqr_in->rhs_vec->elements[bItr] = b[bItr];
		}
		/* Here we set the initial solution vector guess, which is 
		 * a simple 1-vector. You might want to adjust this value for fine-tuning
		 * t_snnls() for your application
		 */
		for( bItr=0; bItr<Af->n; bItr++ )
		{
			lsqr_in->sol_vec->elements[bItr] = 1; 
		}
		
		/* This is a function pointer to the matrix-vector multiplier */
		lsqr_func->mat_vec_prod = sparse_lsqr_mult;
		
		lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, Af );
		
		for( bItr=0; bItr<Af->n; bItr++ ) // not really bItr here, but hey
			x[F[bItr]] = lsqr_out->sol_vec->elements[bItr];
		
		free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );
	}
	
	if( outResidualNorm != NULL )
	{
		double* finalresidual = (taucs_double *)calloc(m,sizeof(taucs_double));
		ourtaucs_ccs_times_vec(A_original_ordering,x,finalresidual);
		cblas_daxpy(m,-1.0,b, 1, finalresidual, 1);
		*outResidualNorm = cblas_dnrm2(m, finalresidual, 1);
		free(finalresidual);
	}

	/* We conclude with a little more memory freeing. */
	taucs_ccs_free(Af);
	free(F); 
	free(G);     
	free(H1); 
	free(H2);
	taucs_ccs_free(AprimeDotA);
	taucs_ccs_free(lsqrApA);
  
	free(y);

	return x;
}

#pragma mark -

#ifndef taucs_add
#define taucs_add(x,y) ((x)+(y))
#endif
#ifndef taucs_mul
#define taucs_mul(x,y) ((x)*(y))
#endif

double
taucs_rcond( taucs_ccs_matrix* A )
{
	char	NORM = '1';
	int		N = 0;
	int		LDA = 0;
	double  ANORM = 0;
	double  RCOND = 0;
	double* WORK = NULL;
	int*	IWORK = NULL;
	int		INFO;
	int*	IPIV = NULL;
	double* lapackA = NULL;
	
	/* Construct LAPACK representation of A and compute the 1 norm of A */
	int vSize;
	int cItr, rItr;
	int rowCount = A->m;
	double localMax = 0;
	
	if( (A->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
	{
		vSize = A->n*A->n;
		rowCount = A->n;
	}
	else
		vSize = A->m*A->n;
	
	lapackA = (double*)malloc(sizeof(double)*vSize);
	bzero(lapackA, sizeof(double)*vSize);
		
	for( cItr=0; cItr<A->n; cItr++ )
	{
		localMax = 0;
		for( rItr=A->colptr[cItr]; rItr<A->colptr[cItr+1]; rItr++ )
		{
			int index = -1;
			index = A->rowind[rItr] + cItr*rowCount;
			if( index > vSize )
			{
				fprintf( stderr, "Rcond memory error!\n" );
				exit(-1);
			}
			lapackA[index] = A->values.d[rItr];
			localMax += fabs(A->values.d[rItr]);
		}
		if( localMax > ANORM )
			ANORM = localMax;
	}

	NORM = '1';
	N = A->n;
	LDA = A->m;
	RCOND = 0;
	WORK = (double*)malloc(sizeof(double)*4*N);
	IWORK = (int*)malloc(sizeof(int)*N);
	INFO = 0;

	IPIV = (int*)malloc(sizeof(int)*min(rowCount, A->n));
	
	dgetrf_( &rowCount, &A->n, lapackA, &rowCount, IPIV, &INFO );
	dgecon_( &NORM, &N, lapackA, &LDA, &ANORM, &RCOND, WORK, IWORK, &INFO );

	free(IPIV);
	free(IWORK);
	free(WORK);
	free(lapackA);

	return RCOND;
}

void
transpose_vec_times_matrix(double* b, double* A, int* F, int A_cols, int rows, int cols, double* result)
{
	// result = b'*A
	int cItr;
	for( cItr=0; cItr<cols; cItr++ )
		result[cItr] = cblas_ddot( rows, b, 1, &A[F[cItr]], A_cols );
}

// here cols doubles for sizeF, since F is a selection of columns from A
void
taucs_transpose_vec_times_matrix(double* b, taucs_ccs_matrix* A, int* F, int cols, double* result)
{
	int cItr, rItr;
	for( cItr=0; cItr<cols; cItr++ )
	{
		result[cItr] = 0;
		for( rItr=0; rItr<A->colptr[F[cItr]+1]-A->colptr[F[cItr]]; rItr++ )
		{
			int tRow = A->rowind[A->colptr[F[cItr]] + rItr];
				
			result[cItr] += b[tRow] * A->values.d[A->colptr[F[cItr]] + rItr];
		}
	}
}

void
transpose_vec_times_matrix_nosub(double* b, double* A, int A_cols, int rows, double* result)
{
	// result = b'*A
	int cItr;
	for( cItr=0; cItr<A_cols; cItr++ )
		result[cItr] = cblas_ddot( rows, b, 1, &A[cItr], A_cols );
}

void
taucs_transpose_vec_times_matrix_nosub(double* b, taucs_ccs_matrix* A, double* result)
{
	int cItr, rItr;
	for( cItr=0; cItr<A->n; cItr++ )
	{
		result[cItr] = 0;
		for( rItr=0; rItr<A->colptr[cItr+1]-A->colptr[cItr]; rItr++ )
		{
			int tRow = A->rowind[A->colptr[cItr] + rItr];
			result[cItr] += b[tRow] * A->values.d[A->colptr[cItr] + rItr];
		}
	}
}

void
ourtaucs_ccs_times_vec( taucs_ccs_matrix* m, 
			 taucs_datatype* X,
			 taucs_datatype* B )
{
	int i,ip,j,n, rows;
	taucs_datatype Aij;

	n = m->n;
	rows = m->m;
	if( (m->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
	{
		// this is a total hack, but taucs's thingy works when 
		// things are symmetric. keep in mind, however, that otherwise, 
		// it does not
		taucs_ccs_times_vec( m, X, B );
		return;
	}
		
	for(i=0; i < rows; i++) 
		B[i] = 0;

	for (j=0; j<n; j++) {
      for (ip = (m->colptr)[j]; ip < (m->colptr[j+1]); ip++) {
		i   = (m->rowind)[ip];
		Aij = (m->values.d)[ip];
		
		B[i] = taucs_add(B[i],taucs_mul(X[j],Aij));
      }
    }
}

void
taucs_ccs_submatrix( const taucs_ccs_matrix* A, const int* keptCols, const int inColCount, taucs_ccs_matrix* result)
{
	int cItr, colOffset, c2;
	
	result->m = A->m;
	result->n = inColCount;
	//result->flags = A->flags;
	result->flags = TAUCS_DOUBLE;
			
	colOffset = 0;
	for( cItr=0; cItr<inColCount; cItr++ )
	{
		result->colptr[cItr] = colOffset;
		for( c2=A->colptr[keptCols[cItr]]; c2<A->colptr[keptCols[cItr]+1]; c2++ )
		{
			result->rowind[colOffset] = A->rowind[c2];
			result->values.d[colOffset] = A->values.d[c2];
			colOffset++;
		}
	}
	result->colptr[cItr] = colOffset;
}



taucs_ccs_matrix*
taucs_ccs_transpose( const taucs_ccs_matrix* A )
{
  taucs_ccs_matrix* result = NULL;
  double* values = NULL;
  int cItr, rItr, colOffset;
	
  result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
	
  result->m = A->n;
  result->n = A->m;
	
  result->flags = A->flags;
	
  // we sacrifice some memory here so we don't have to do the number of 
  // nonzero entries computation
  result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
  result->rowind = (int*)malloc(sizeof(int)*result->m*result->n);
  result->values.d = (double*)malloc(sizeof(double)*result->m*result->n);
	
  values = taucs_convert_ccs_to_doubles(A);
	
  colOffset = 0;
  for( rItr=0; rItr<A->m; rItr++ )
    {
      result->colptr[rItr] = colOffset;
      for( cItr=0; cItr<A->n; cItr++ )
	{
	  double v = values[(rItr*A->n)+cItr];
	  if( v != 0 )
	    {
	      result->rowind[colOffset] = cItr;
	      result->values.d[colOffset] = v;
	      colOffset++;
	    }
	}
    }
  result->colptr[rItr] = colOffset;
	
  free(values);
	
  return result;
}

taucs_ccs_matrix*
taucs_construct_sorted_ccs_matrix( double* values, int rowsize, int rows )
{
  taucs_ccs_matrix* result = NULL;
  int nnz = 0;
  int i, rItr, cItr, colOffset;
  double v;
	
  for( i=0; i<rowsize*rows; i++ )
    {
      if( values[i] != 0 )
		nnz++;
    }
	
  result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  result->n = rowsize;
  result->m = rows;
  result->flags = TAUCS_DOUBLE;
	
  result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
  result->rowind = (int*)malloc(sizeof(int)*nnz);
  result->values.d = (double*)malloc(sizeof(taucs_double)*nnz);
	
  colOffset = 0;
	
  for( cItr=0; cItr<rowsize; cItr++ )
    {
      result->colptr[cItr] = colOffset;
      for( rItr=0; rItr<rows; rItr++ )
	{
	  v = values[rItr*rowsize+cItr];
	  if( v != 0 )
	    {
	      result->rowind[colOffset] = rItr;
	      result->values.d[colOffset++] = v;
	    }
	}
    }
  result->colptr[cItr] = colOffset;
	
  return result;
}

double*
taucs_convert_ccs_to_doubles( const taucs_ccs_matrix* A )
{
	int vSize;
	int cItr, rItr;
	double* values;
	
	if( (A->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
		vSize = A->n*A->n;
	else
		vSize = A->m*A->n;
	
	values = (double*)malloc(sizeof(double)*vSize);
	
	bzero(values, sizeof(double)*vSize);
	
	int rowSize = A->n;
	
	for( cItr=0; cItr<A->n; cItr++ )
	{
		for( rItr=A->colptr[cItr]; rItr<A->colptr[cItr+1]; rItr++ )
		{
			int index = -1;
			index = rowSize*A->rowind[rItr] + cItr;
			values[index] = A->values.d[rItr];
		}
	}
	
	return values;
}

void
taucs_print_ccs_matrix( const taucs_ccs_matrix* A )
{
  double* v = taucs_convert_ccs_to_doubles(A);
  int rowCount;
  int rItr, cItr;
	
  rowCount = A->m;
	
  if( (A->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
    {
      printf( "Matrix flagged symmetric\n" );
      rowCount = A->n;
    }
	
  for( rItr=0; rItr<rowCount; rItr++ )
    {
      for( cItr=0; cItr<A->n; cItr++ )
	printf( "%5.4g ", v[rItr*A->n + cItr] );
      printf( "\n" );
    }
	
  free(v);
}

