
/*
 * This program is free software distributed under the GPL. A copy of the license should have been included with this 
 * archive in a file named 'LICENSE'. You can read the license there or on the web at: http://www.gnu.org/licenses/gpl.txt
 */

#include "tsnnls.h"

#ifdef __APPLE__
#include <vecLib/vBLAS.h>
#include <vecLib/clapack.h>
#else
#include "gsl_cblas.h"
#endif // __APPLE__

/* This file contains code for the TAUCS version of lsqr, including
   two internal procedures called by the main routine below. 

   Cantarella/Piatek, 5/04. 
*/

/* Procedure is internal to taucs_ccs_aprime_times_a, which is 
 * presented below. 
 */
double 
taucs_dotcols( const taucs_ccs_matrix* A, int col1, int col2 )
{
	double val = 0;
	int i, j;

	i=A->colptr[col1];
	j=A->colptr[col2];

	// eliminate some loop invariant pointer derefs
	double* Avals = A->values.d;
	int* rowinds = A->rowind;
	int* colptrs = A->colptr;
	
	int stopI, stopJ;
	
	stopI = colptrs[col1+1];
	stopJ = colptrs[col2+1];

	while( 1==1 )
	{
		if( rowinds[i] == rowinds[j] )
			val += Avals[i++]*Avals[j++];
		else if( rowinds[i] < rowinds[j] )
			i++;
		else
			j++;
		if( i>=stopI )
			break;
		if( j>=stopJ )
			break;
	}
	return val;
}

/* 
 * This routine computes A'*A for full A. This is useful for comparing sparse/nonsparse
 * performance to determine if your problem can benefit from a different representation.
 * It is left in the archive as a convenience
 */
double* 
full_aprime_times_a( double* A, int rows, int cols )
{
	int rItr, cItr, colOffset;
	double* result = (double*)calloc(cols*cols, sizeof(double));
  		
	colOffset = 0;
			
	for( cItr=0; cItr<cols; cItr++ )
	{
		for( rItr=cItr; rItr<cols; rItr++ )
		{
			result[rItr*cols + cItr] = cblas_ddot( rows, &A[cItr], cols, &A[rItr], cols );
		}
	}
  	
	return result;
}

/* 
 * Procedure computes the matrix product A'*A for a sparse matrix 
 * A in compressed column storage. The result is a symmetric TAUCS
 * matrix, which is also in CCS format. 
 */
taucs_ccs_matrix* 
taucs_ccs_aprime_times_a( taucs_ccs_matrix* A )
{
	int rItr, cItr, colOffset, maxSize, scItr, srItr;
	taucs_ccs_matrix* result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
	register double v;
	int newSize, currentSize;

	result->n = A->n;

	/* we'll be symmetric */
	result->flags = TAUCS_DOUBLE;
	result->flags = result->flags | TAUCS_SYMMETRIC;
	result->flags = result->flags | TAUCS_LOWER; // rep the lower half

	result->colptr = (int*)malloc(sizeof(int)*(result->n+1));

	currentSize = A->colptr[A->n]*2; // start with the number of entries in A*2 and see if we need more
	result->values.d = (double*)malloc(sizeof(taucs_double)*currentSize);
	result->rowind = (int*)malloc(sizeof(int)*currentSize);
	
	colOffset = 0;
	double* valsPtr = result->values.d;
	int* rowptrs = result->rowind;
	int* colptrs = result->colptr;
	int Acols = A->n;
	int* Acolptrs = A->colptr;
	double* AvalsPtr = A->values.d;
	int* Arowptrs = A->rowind;
	int stopPoint, interiorStopPoint;
	
	for( cItr=0; cItr<Acols; cItr++ )
	{
		colptrs[cItr] = colOffset;
		
		for( rItr=cItr; rItr<Acols; rItr++ )
		{
			v = taucs_dotcols(A,cItr,rItr);
			
			if( v == 0.0 )
			{
				continue;
			}
			else
			{
				valsPtr[colOffset] = v;
				rowptrs[colOffset] = rItr;
				colOffset++;
				
				if( colOffset < currentSize )
					continue;
				else
				{
					/* we need to increase our allocation size.  */
					newSize = 2*currentSize;
					int* newRows = (int*)realloc(rowptrs, sizeof(int)*newSize);
					double* newVals = (double*)realloc(valsPtr, sizeof(double)*newSize);
					
					currentSize = newSize;
					
					if( newRows == NULL || newVals == NULL )
					{
						fprintf( stderr, "tsnnls: Out of memory!\n" );
					}
					
					result->values.d = newVals;
					valsPtr = newVals;
					result->rowind = newRows;
					rowptrs = newRows;
				}
			}
		}
	}
	colptrs[cItr] = colOffset;
	

	return result;
}

static void
ccs_to_lapack( taucs_ccs_matrix* L, double** lapackL, int* N, int* LDA, double* ANORM )
{	
	/* Construct LAPACK representation of A and compute the 1 norm of A */
	int vSize;
	int cItr, rItr;
	int rowCount = L->m;
	double localMax = 0;
	
	*ANORM = 0;
	
	if( (L->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
	{
		vSize = L->n*L->n;
		rowCount = L->n;
	}
	else
		vSize = L->m*L->n;
	
	*lapackL = (double*)malloc(sizeof(double)*vSize);
	bzero(*lapackL, sizeof(double)*vSize);
		
	for( cItr=0; cItr<L->n; cItr++ )
	{
		localMax = 0;
		for( rItr=L->colptr[cItr]; rItr<L->colptr[cItr+1]; rItr++ )
		{
			int index = -1;
			index = L->rowind[rItr] + cItr*rowCount;
			(*lapackL)[index] = L->values.d[rItr];
			localMax += fabs(L->values.d[rItr]);
		}
		if( localMax > *ANORM )
			*ANORM = localMax;
	}
	
	*N = L->n;
	*LDA = rowCount;
}

static double
t_condest( void* mfR )
{
	taucs_ccs_matrix* L;
	double* lapackL;
	int N, LDA, INFO;
	char	UPLO;
	double  ANORM = 0;
	double* WORK;
	int*	IWORK;
	double  RCOND;
	
	L = taucs_supernodal_factor_to_ccs(mfR);
	
	/* Construct LAPACK representation of A and compute the 1 norm of A */
	int vSize;
	int cItr, rItr;
	int rowCount = L->m;
	double localMax = 0;
	
	if( (L->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
	{
		vSize = L->n*L->n;
		rowCount = L->n;
	}
	else
		vSize = L->m*L->n;
	
	lapackL = (double*)malloc(sizeof(double)*vSize);
	bzero(lapackL, sizeof(double)*vSize);
		
	for( cItr=0; cItr<L->n; cItr++ )
	{
		localMax = 0;
		for( rItr=L->colptr[cItr]; rItr<L->colptr[cItr+1]; rItr++ )
		{
			int index = -1;
			index = L->rowind[rItr] + cItr*rowCount;
			lapackL[index] = L->values.d[rItr];
			localMax += fabs(L->values.d[rItr]);
		}
		if( localMax > ANORM )
			ANORM = localMax;
	}

	N = L->n;
	LDA = L->m;
	UPLO = 'L';
	WORK = (double*)malloc(sizeof(double)*3*N);
	IWORK = (int*)malloc(sizeof(int)*N);
	
	dpocon_( &UPLO, &N, lapackL, &LDA, &ANORM, &RCOND, WORK, IWORK, &INFO );
		
	free(WORK);
	free(IWORK);
	taucs_ccs_free(L);
	free(lapackL);
	
	return RCOND;
}

/*
 * This is a custom version of the TAUCS based lsqr solver for use with the t_snnls() block principle pivoting
 * function, which precomputes several values of interest to lsqr so as to not repeat their calcuation during
 * the multiple pivots which can occur.
 */
taucs_double*
t_snnlslsqr(taucs_ccs_matrix *A,
			taucs_double *b, 
			taucs_ccs_matrix* ApA, 
			int* F,
			double* outRcond )
{
	taucs_ccs_matrix	*ApAperm;
	void				*Apb /*n x 1*/, *ApAx /*n x 1*/, *x /* n x 1 */, *Itstep /*n x 1*/; 
	double				*x_unscrambled;
	void				*mfR; /* The Cholesky factor R, stored by TAUCS. */
	int					*perm, *invperm;
	char				*ordering;
	
	ordering = getenv("COL_ORDERING");
	if( ordering == NULL )
	{
		/* use amd ordering if the user hasn't specified anything else */
		setenv("COL_ORDERING", "amd", 0);
		ordering = getenv("COL_ORDERING");
	}
	taucs_ccs_order(ApA, &perm, &invperm, ordering);
	ApAperm = taucs_ccs_permute_symmetrically(ApA, perm, invperm);

	Apb = calloc(A->m, sizeof(taucs_double));

	mfR = taucs_ccs_factor_llt_mf(ApAperm);
	if( mfR == NULL )
	{
		// free the memory we've allocated so far and return
		taucs_ccs_free(ApAperm);
		free(Apb);
		free(perm);
		free(invperm);
		return NULL;
	}
	
    if( outRcond != NULL )
		*outRcond = t_condest(mfR);
	
	/* We now solve the first equation: x = R\(A'*A*b). */
	x = calloc(A->n,sizeof(taucs_double));
	
	taucs_transpose_vec_times_matrix_nosub(b, A, Apb);
	// we have to permute A'b to be meaningful.
	{
		double* apbperm = Apb;
		Apb = malloc(sizeof(double)*A->n);
		taucs_vec_permute(A->n, TAUCS_DOUBLE, apbperm, Apb, perm);
		free(apbperm);
	}
	
	taucs_supernodal_solve_llt(mfR,x,Apb);  /* n x n * n x 1 = n x 1 */

	/* Given the base solution, we now update it by refinement. */
	ApAx = (taucs_double*)malloc(sizeof(double)*A->n);
	Itstep = calloc(A->n,sizeof(taucs_double)); 
	
	double* scratch = (double*)malloc(sizeof(double)*ApAperm->n);
		
	memcpy(scratch, Apb, sizeof(double)*ApAperm->n);
	
	ourtaucs_ccs_times_vec(ApAperm,x,ApAx);  /* n x n * n x 1 = n * 1. */
	cblas_daxpy(A->n,-1.0,(double *)(ApAx),1,(double *)(scratch),1); /* Apb = Apb - ApAx */
	//refinementEps = cblas_dnrm2(A->n, scratch, 1);
	taucs_supernodal_solve_llt(mfR,Itstep,scratch);                 /* Itstep = R\Apb. */
	cblas_daxpy(A->n,1.0,(double *)(Itstep),1,(double *)(x),1);  /* x = x + Itstep */
			
	free(scratch);
	free(Itstep);
	free(ApAx);

	/* We now have a solution, but must clean up the memory we allocated
	 * and unscramble x 
	 */
	x_unscrambled = malloc(sizeof(double)*ApA->n);
	taucs_vec_permute(ApA->n, TAUCS_DOUBLE, x, x_unscrambled, invperm);

	taucs_ccs_free(ApAperm);
	free(Apb); 
	free(perm);
	free(invperm);
	free(x);
	taucs_supernodal_factor_free(mfR);

	return x_unscrambled;
}

taucs_double*
t_lsqr(taucs_ccs_matrix *A, taucs_double *b)
{
	int *F;
	int i;
	double* x;
	taucs_ccs_matrix* apda;
	F = malloc(sizeof(int)*A->n);
	for( i=0; i<A->n; i++ )
		F[i] = i;
	apda = taucs_ccs_aprime_times_a(A);
	x = t_snnlslsqr(A, b, apda, F, NULL);
	taucs_ccs_free(apda);
	free(F);
	
	return x;
}
