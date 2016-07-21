#ifndef DIFFSNORM_H
#define DIFFSNORM_H 
#include "macros.h"
#include <iostream>
#include "fastpca_io.hpp"
#include "fastpca_linear_algebra.hpp"
#include "InputMatrix.h"
#include <math.h>

/*********************************************************************
 *  DIFFSNORM 2-norm accuracy of an approx. to a matrix
 * 
 *  USAGE 
 ** int ier = diffsnorm(A,U,S,V,k,its,snorm) computes an estimate 
 ** snorm of the spectral norm (the operator induced by the Euclidan
 ** vector norm) of A-USV', using its iterations of the power method
 ** started with a random vector.
 **
 *
 *  INPUTS
 ** InputMatrix * A -- The matrix being approximated. First matrix in A-USV'
 ** long long  k -- rank of the approximation to A.  In other words, 
 the number of eigen values/vectors calculated in the 
 approximation.
 *
 ** double * U -- Second matrix in A-USV'
 ** double * S -- Third matrix in A-USV'
 ** double * V -- Fourth matrix in V'
 *
 *  OUTPUTS
 ** double snorm -- spectral norm of A-USV'
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the PCA function executed successfully
 **               ier=-1 means that m,n,k,or its are negative
 **               ier=-2 means that k > min(m,n) 
 **               ier=-30 means an allocation failed 
 **               ier=-100 means that the input parameters may be valid, but 
 we have not yet implemented that case

 */
int fastpca_diffsnorm(  InputMatrix* A, double* U, double* S , double* V,long long k, int its, double &snorm );

#endif /* DIFFSNORM_H */
