#ifndef FASTPCA_INTERFACE_H
#define FASTPCA_INTERFACE_H
/*
 * These functions provide an interface for the R C++ code to call.  The reason it is separated is so that the user does not have to have Intel MKL installed to compile the R package.  
 *
 * This file is actually copied then into the fastRPCA folder as it needs to be included by the R package.
 */
extern "C" {
int fastPCAMemory( double * A, double **U, double **S, double **V, long long m, long long n, long long k, long long l, long long its,  int centering_row, int centering_column, int dits, int diffsnorm,double &snorm );
int fastPCAFile ( int ,  const char * inputFileName, double **U, double **S, double **V, long long &m, long long &n, long long k, long long l, long long its, long long memory, int centering_row, int centering_column, int dits, int diffsnorm, double& snorm);
}



#endif /* FASTPCA_IO_H */
