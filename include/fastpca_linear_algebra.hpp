#ifndef FASTPCA_LINEAR_ALGEBRA_H
#define FASTPCA_LINEAR_ALGEBRA_H 
#include "macros.h"
#include <time.h>

#ifdef USE_MKL
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#else
#include <lapacke.h>
#include <cblas.h>
#define USEMKL_TEST 0
#endif

#ifdef __APPLE__
//#include <Accelerate/Accelerate.h>
//#define lapack_int              long long
#include <stdlib.h>
#else
#include <malloc.h>
#endif

#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <fstream>



int fastpca_gemm(CBLAS_ORDER layout, CBLAS_TRANSPOSE transposeA,
		CBLAS_TRANSPOSE transposeB, const long long m, const long long n ,
		const long long k, const double alpha, const double * A,
		const long long lda, const double * B, const long long ldb,
		const double beta, double * C, const long long ldc);
/*********************************************************************
 *  Populate matrix with random numbers drawn from a uniform distribution
 * 
 *  USAGE 
 ** int ier = fastpca_populate_matrix_random(m,n,A) will populate matrix A
 ** of size (m,n) with values between -1 and 1 drawn from the pseudo-random number 
 ** generator rand(), with random seed set to current time
 **
 *
 *  INPUTS
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The matrix being populated
 *
 *  OUTPUTS
 ** double * A -- The matrix, populated with the random values
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully

 */
int fastpca_populate_matrix_random (long long m, long long n, double * A);

/*********************************************************************
 *  Divide matrix by a scalar.
 * 
 *  USAGE 
 ** int ier = fastpca_divide_matrix_scalar(m,n,A,k) divides the matrix A
 ** of size m,n by the scalar k
 **
 *
 *  INPUTS
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The matrix being divided by a scalar k
 ** double k --   The divisor, k
 *
 *  OUTPUTS
 ** double * A -- The matrix, divided by the scalar k.
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully

 */
int  fastpca_divide_matrix_scalar(long long m, long long n, double * A, double k);


/*********************************************************************
 *  Mean centers the rows of matrix
 * 
 *  USAGE 
 ** int ier = fastpca_mean_center_rows(m,n,A) centers the rows
 ** of A
 *
 *  INPUTS
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The matrix whose rows will be centered
 *
 *  OUTPUTS
 ** double * A -- The matrix with mean centered rows
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully
 **               ier=-30 means that memory was not allocated successfully

 */
int fastpca_mean_center_rows(long long m, long long n,double * A);
/*********************************************************************
 *  Mean centers the columns of matrix
 * 
 *  USAGE 
 ** int ier = fastpca_mean_center_columns(m,n,A) centers the columns
 ** of A
 *
 *  INPUTS
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The matrix whose columns will be centered
 *
 *  OUTPUTS
 ** double * A -- The matrix with mean centered columns
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully
 **               ier=-30 means that memory was not allocated successfully

 */
int fastpca_mean_center_columns(long long m, long long n,double * A);

/*********************************************************************
 *  Performs the LU decomposition and extracts the L matrix
 * 
 *  USAGE 
 ** int ier = fastpca_lu_l(m,n,A) factors a matrix LU as the product of a 
 ** lower triangular matrix L, and an upper triangular matrix U.  It 
 ** then replaces the values of A with L
 *
 *  INPUTS
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The matrix to be factored
 *
 *  OUTPUTS
 ** double * A -- The matrix L, where A=LU
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully
 **               ier=-30 means that memory was not allocated successfully

 */
int fastpca_lu_l(long long m, long long n, double * A);
/*********************************************************************
 *  Performs the QR decomposition and extracts the Q matrix
 * 
 *  USAGE 
 ** int ier = fastpca_qr_q(m,n,A) factors a matrix QR as the product of an
 ** orthogonal matrix Q and an upper triangular matrix R. 
 *
 *  INPUTS
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The matrix to be factored
 *
 *  OUTPUTS
 ** double * A -- The matrix Q, where A=LU
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully
 **               ier=-1 means that the dgeqrf failed
 **               ier=-2 means that the dorgqr failed
 **               ier=-30 means that memory was not allocated successfully

 */
int fastpca_qr_q(long long m, long long n, double * A);

void fastpca_square_transpose ( double * A, long long m);

void * fastpca_aligned_alloc(size_t alignment, size_t size);
void  fastpca_aligned_free(void *ptr);
#endif /* FASTPCA_LINEAR_ALGEBRA_H  */
