
/*
 * This program is free software distributed under the GPL. A copy of the license should have been included with this 
 * archive in a file named 'LICENSE'. You can read the license there or on the web at: http://www.gnu.org/licenses/gpl.txt
 */

/*

    tsnnls.h : Header for unconstrained and constrained sparse 
                   least-squares solvers built on the TAUCS library.

		   The unconstrained solver uses the method of normal
		   equations based on a Cholesky decomposition, as in
		   Bjork (p. 264).

		   The constrained solver is an implementation of 
		   Portugal, Judice, and Vicente's block-principal
		   pivoting algorithm.

*/

#if (__cplusplus || c_plusplus)
extern "C" {
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "taucs_basic/taucs.h"

/* 
 *Some constants that affect the operation of t_snnls() 
 */

/* t_snnls requires us to handle in some way the case where values
 * returned from lsqr that are within error tolerance of zero get continually
 * swapped between x and y because our comparison to zero is numerically more 
 * precise in the infeasibles computation. In order to terminate, we must 
 * handle this case. For our usage of snnls, it is most efficient to simply
 * zero values within some epsilon of zero as returned from lsqr, but this 
 * solution may not be the best for all possible applications. This value in the zero
 * error tolerance.
 */
#define kZeroThreshold  1e-12

/* Utility Functions for CCS matrices */
void				taucs_ccs_submatrix( const taucs_ccs_matrix* A, const int* keptCols, 
						const int inColCount, taucs_ccs_matrix* result);
 taucs_ccs_matrix*   taucs_ccs_transpose( const taucs_ccs_matrix* A );

/* taucs versions of these are sparse, otherwise full. */
void				transpose_vec_times_matrix(double* b, double* A, int* F, int A_cols, int rows, int cols, double* result);
void				taucs_transpose_vec_times_matrix(double* b, taucs_ccs_matrix* A, int* F, int cols, double* result);
void				transpose_vec_times_matrix_nosub(double* b, double* A, int A_cols, int rows, double* result);
void				taucs_transpose_vec_times_matrix_nosub(double* b, taucs_ccs_matrix* A, double* result);
void				ourtaucs_ccs_times_vec( taucs_ccs_matrix* m, taucs_datatype* X, taucs_datatype* B );

void				taucs_print_ccs_matrix( const taucs_ccs_matrix* A );
taucs_ccs_matrix*   taucs_construct_sorted_ccs_matrix( double* values, int rowsize, int rows );
double*				taucs_convert_ccs_to_doubles( const taucs_ccs_matrix* A );

/* TAUCS version of lsqr  */
double				taucs_dotcols( const taucs_ccs_matrix* A, int col1, int col2 );
double*				full_aprime_times_a(double* A, int rows, int cols);
taucs_ccs_matrix*   taucs_ccs_aprime_times_a( taucs_ccs_matrix* A );
taucs_double*		t_lsqr(taucs_ccs_matrix *A, taucs_double *b);
taucs_double*		t_snnlslsqr(taucs_ccs_matrix *A,taucs_double *b, taucs_ccs_matrix* ApA, 
							int* F, double* outRcond );
							
/* TAUCS version of snnls  */
double				taucs_rcond( taucs_ccs_matrix* A );
taucs_ccs_matrix*   selectAprimeDotA( double* apda, int cols, int* F, int sizeF );
void				selectAprimeDotAsparse( const taucs_ccs_matrix* apda, int* F, int sizeF, taucs_ccs_matrix* inOutApda );

/**
 * NOTE - general taucs_ccs_matrix types may NOT be valid in t_snnls, for instance:
 *				- A->rowind must be in SORTED ASCENDING order for each column
 *
 *  taucs_convert_ccs_to_doubles() is robust to these issues, and you can generally
 *  construct a valid t_snnls matrix by performing:
 *
 *		double* Avals = taucs_convert_ccs_to_doubles(A);
 *		taucs_ccs_matrix* fixed = taucs_construct_sorted_ccs_matrix(Avals, A->n, A->m);
 *		free(Avals);
 *		taucs_ccs_free(A);
 *		A = fixed;
 */
taucs_double*		t_snnls( taucs_ccs_matrix *A_original_ordering, taucs_double *b, 
							 double* outResidualNorm, double inRelErrTolerance, int inPrintErrorWarnings );

#if (__cplusplus || c_plusplus)
};
#endif 

