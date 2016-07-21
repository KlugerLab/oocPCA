#ifndef FASTPCA_IO_H
#define FASTPCA_IO_H 
#include "macros.h"
#include <time.h>
#include "fastpca_linear_algebra.hpp"
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include "InputMatrix.h"
#include <cstring>
/*********************************************************************
 *  For diagnostic purposes. Prints to standard output the contents of a matrix. 
 *  Uncomment the content of the function for diagnostic output.
 * 
 *  USAGE 
 ** int ier = fastpca_print_matrix(desc, m,n a, lda) outputs the contents
 ** of a matrix labelled desc, of size m,n, and stored in a, with a
 ** leading edge lda. On row-major machines, lda==n
 *
 *  INPUTS
 ** const char *  desc  -- description or name of matrix
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The (m,n) matrix A, to print to standard output
 ** long long  lda -- the leading edge of matrix A
 *
 *  OUTPUTS
 ** None.
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully

 */
int fastpca_print_matrix( const char* desc, long long m, long long n, double* A );
/*********************************************************************
 *  Reads the principal components of the matrix in plink BED format
 * 
 *  USAGE 
 ** int ier = fastpca_read_bed_format(filename, m,n,A) reads input 
 *  filename into a matrix size (m,n) pointed to by A.  Note that 
 *  m,n,and A are all outputs, as they are determined from the file.
 *
 *  INPUTS
 ** const char *  filename  -- filename of output file
 *
 *  OUTPUTS
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The (m,n) input matrix A
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully
 **               ier=-1 means could not open the .bim file
 **               ier=-2 means could not open the .bed file
 **               ier=-3 means could not open the .fam file
 **               ier=-20 means that the binary file was not in proper BED format
 **               ier=-21 means that the BED file was not in row-major format
 **               ier=-30 means that memory was not allocated successfully
 **               ier=-400 means that one of the values was missing

 */
int fastpca_read_bed_format(const char * filename, long long &m, long long &n, double *& A);
/*********************************************************************
 *  Impute the missing data by assigning it the average value for the SNP
 * 
 *  USAGE 
 ** int ier = fastpca_impute_by_average(A,m,n, missingValue) 
 *
 *  INPUTS
 ** double * A -- The (m,n) input matrix A
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double  missingValue -- The missing value to replace with the average of each SNP
 *
 *  OUTPUTS
 ** double * A -- The (m,n) input matrix A
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully
 **               ier=-400 there are no non-missing individuals for at least one SNP

 */
int fastpca_impute_missing_average(double *& A, long long m, long long n, double missingValue);
/********************************************************************* 
 *  Outputs matrix in binary format 
 *
 *  USAGE 
 ** int ier = fastpca_write_binary_format(filename, m,n,A) outputs
 *
 *  INPUTS
 ** const char *  filename  -- filename of output file
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The (m,n) matrix A to be output
 *
 *  OUTPUTS
 ** None.
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully
 **               ier=-1 means that the function could not open the output file
 **               ier=-30 means that memory was not allocated successfully

 */
int fastpca_write_binary_format(const char * filename, long long m, long long n,double * A);
int fastpca_write_input_matrix_binary_format(const char * filename, InputMatrix * A) ;
/********************************************************************* 
 *  Outputs the principal components of the matrix in eigenstrat format * 
 *  USAGE 
 ** int ier = fastpca_write_eigenstrat_format(filename, m,n,k, S,V) outputs
 * S and V in the eigenstrat format
 *
 *  INPUTS
 ** const char *  filename  -- filename of output file
 ** long long  m -- number of rows in the input matrix
 ** long long  n -- number of columns  in the input matrix
 ** long long  k -- rank of the decomposition 
 ** double * S -- The (k,k) matrix S, of the SVD A=USV'
 ** double * V -- The (n,k) matrix V, of the SVD A=USV
 *
 *  OUTPUTS
 ** None.
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully
 **               ier=-1 means that the function could not open the output file
 **               ier=-30 means that memory was not allocated successfully

 */
int fastpca_write_eigenstrat_format(const char * filename, long long m, long long n,long long k,double *S, double *V);
int fastpca_save_bin( const char * basename, InputMatrix* A, double* U, double* S , double* V,long long k );
#endif /* FASTPCA_IO_H */
