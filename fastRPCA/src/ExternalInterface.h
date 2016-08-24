#ifndef FASTPCA_INTERFACE_H
#define FASTPCA_INTERFACE_H
extern "C" {
int fastPCAMemory( double * A, double **U, double **S, double **V, long long m, long long n, long long k, long long l, long long its, long long memory, int centering_row, int centering_column, int dits, int diffsnorm,double &snorm );
int fofo();
}



#endif /* FASTPCA_IO_H */
