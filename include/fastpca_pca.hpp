#ifndef FASTPCA_PCA_H
#define FASTPCA_PCA_H 
#include "macros.h"
#include <iostream>
#include "fastpca_linear_algebra.hpp"
#include "fastpca_io.hpp"
#include "InputMatrix.h"
#include <math.h>
#include <chrono>
/*********************************************************************
 *  PCA Principal component analysis
 * 
 *  USAGE 
 ** int ier = pca(m,n,k,l,A,U,S,V,its) constructs a nearly optimal rank-k approximation
 ** USV' to to A, using its power iterations, with block size l, started with
 ** a min(m,n)*l random matrix, when A is m xn; the reference below explains
 ** "nearly optimal."  The smaller dimension of A must be >= k.  The result 
 ** of the decomposition are stored in U,S,V.
 **
 *  THEORY
 ** The low-rank approximation USV' comes in the form of a singular
 ** value decomposition (SVD) -- the columns of U are orthonormal, as
 ** are the columns of V, the entries of S are all nonnegative, and all
 ** nonzero entries of S appear in non-increasing order on its diagonal.
 ** U is m x k, V is n x k, and S is k x k, when A is m x n.
 **
 ** Increasing its or l improves the accuracy of the approximation USV';
 ** the reference below describes how the accuracy depends on its and l.
 ** Please note that even its=1 guarantees superb accuracy, whether or
 ** not there is any gap in the singular values of the matrix A being
 ** approximated, at least when measuring accuracy as the spectral norm
 ** ||A-USV'|| of A-USV' (relative to the spectral norm ||A|| of A).
 *
 *  INPUTS
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** long long  k -- rank of the approximation being constructed;
                  k must be a positive integer <= the smaller
		  dimension of A.
 ** long long  l -- block size of the normalized power iterations;
                  l must be a positive integer >=k
 ** double * A -- The matrix being approximated, preferably allocated
                  using mkl_malloc for consistent aligment.
 *
 *  OUTPUTS
 ** double * U -- (m,k) matrix in the rank-k approximation USV' to A
                  where A is (m,n); the columns of U are orthonomal
 ** double * S -- (k,k) matrix in the rank-k approximation USV' to A
                  where A is (m,n); the entires of S are all nonnegative,
		  and all nonzero entries appear in nonincreasing order
		  on the diagonal
 ** double * V -- (k,n) matrix in the rank-k approximation USV' to A
                  where A is (m,n); the columns of U are orthonomal
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the PCA function executed successfully
 **               ier=-1 means that k < 0
 **               ier=-2 means that k > min(m,n) 
 **               ier=-3 means that its < 0 
 **               ier=-4 means that l < k 
 **               ier=-5 means that A is invalid 
 **               ier=-6 means that U is invalid 
 **               ier=-7 means that S is invalid 
 **               ier=-8 means that V is invalid 
 **               ier=-20 the qr function in the power iterations failed 
 **               ier=-23 the lu function in the power iterations failed 
 **               ier=-22 the lu function before the power its 
 **               ier=-21 the qr function in the its==0 case failed 
 **               ier=-30 means an allocation failed 
 **               ier=-100 means that the input parameters are valid, but 
                  we have not yet implemented that case

 */
int fastpca_pca( long long k, long long l, 
		InputMatrix * A, double* U, double * S, double *V, long long its);


#endif /* FASTPCA_PCA_H */
