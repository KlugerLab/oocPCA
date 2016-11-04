#include "fastpca_linear_algebra.hpp"
#include <math.h>
#include <dlfcn.h>

#ifdef __APPLE__
#include <stdlib.h>
#else
#include <malloc.h>
#endif



int fastpca_gemm(CBLAS_ORDER layout, CBLAS_TRANSPOSE transposeA,
		CBLAS_TRANSPOSE transposeB, const long long m, const long long n ,
		const long long k, const double alpha, const double * A,
		const long long lda, const double * B, const long long ldb,
		const double beta, double * C, const long long ldc){

	void (*dgemm_fxn)(char, char,
			char, const long long, const long long,
			const long long, const double, const double *,
			const long long, const double *, const long long,
			const double, double *, const long long);
	void 	*hndl;
	//hndl = dlopen("libcustom_mkl.so", RTLD_LOCAL | RTLD_LAZY);
	//hndl = dlopen("/lib/libcustom_mkldist.so", RTLD_LOCAL | RTLD_LAZY);
	#ifdef __APPLE__
	hndl = dlopen("libfastpca_custommkl.dylib", RTLD_LOCAL | RTLD_LAZY);
	#else
	hndl = dlopen("libfastpca_custommkl.so", RTLD_LOCAL | RTLD_LAZY);
	#endif
	//hndl = dlopen(MKL_CUSTOM_LIBRARY, RTLD_LOCAL | RTLD_LAZY);
	if (!hndl) {
		fprintf(stderr, "%s\n", dlerror());
		exit(EXIT_FAILURE);
		return -1;
	}
	*(void **) (&dgemm_fxn) = dlsym(hndl,"cblas_dgemm");
	dgemm_fxn(layout, transposeA, transposeB, 
			m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
	return 1;
}
void * fastpca_aligned_alloc(size_t alignment, size_t size) {
	void * ptr;
	posix_memalign(&ptr, alignment, size);
	return ptr;
}
void  fastpca_aligned_free(void* ptr) {
	free(ptr);
}
int fastpca_populate_matrix_random (long long m, long long n, double * A) {
	srand(time(NULL));
	//srand(3);
	for (int i=0; i<m*n; i++) {
		A[i] = ((double) 2.0*rand())/((double) RAND_MAX) - 1;

	}
	return 0;
}

int  fastpca_divide_matrix_scalar(long long m, long long n, double * A, double k) {
	for (int i=0; i<m*n; i++ ){
		A[i] = A[i] / k;
	}
	return 0;
}

int fastpca_mean_center_rows(long long m, long long n,double * A) {

	/*********************************************************************
	 *  In order to take advantage of the highly paralellized matrix multiplication
	 *  the row sums are found by multiplying A by a vector of size (1,m) 
	 *  populated with the value 1
	 */
	double * ones = (double *)fastpca_aligned_alloc(64, n*sizeof( double ));
	if (NULL == ones) {
		return -30;
	}
	for (int i=0; i<n; i++){
		ones[i] = 1;
	}
	double * rowsums = (double *)fastpca_aligned_alloc(64, m*sizeof( double ));
	fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m,1,n, 1.0, A,n, ones,1, 0.0, rowsums,1);

	/*********************************************************************
	 *  Determine the mean of each row, and then subtract the mean from each row 
	 */
	for (int i=0; i<m; i++ ){
		rowsums[i] /= n;
	}

	for (long long i=0; i<m; i++){
		for (long  long j=0; j<n; j++ ){
			A[i*n+j] -= rowsums[i];
		}
	}
	return 1;
}
int fastpca_mean_center_columns(long long m, long long n,double * A) {

	/*********************************************************************
	 *  In order to take advantage of the highly paralellized matrix multiplication
	 *  the column sums are found by multiplying A by a vector of size (1,m) 
	 *  populated with the value 1
	 */
	double * ones = (double *)fastpca_aligned_alloc(64, m*sizeof( double ));
	if (NULL == ones) {
		return -30;
	}
	for (int i=0; i<m; i++){
		ones[i] = 1;
	}
	double * colsums = (double *)fastpca_aligned_alloc(64, n*sizeof( double ));
	fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 1,n,m, 1.0, ones,m, A,n, 0.0, colsums,n);

	/*********************************************************************
	 *  Determine the mean of each column, and then subtract the mean from each column 
	 */
	for (int i=0; i<n; i++ ){
		colsums[i] /= m;
	}

	for (long long i=0; i<m; i++){
		for (long  long j=0; j<n; j++ ){
			A[i*n+j] -= colsums[j];
		}
	}
	return 1;
}

int fastpca_lu_l(long long m, long long n, double * A) {

	/*********************************************************************
	 * Perform the LU decomposition using the ?GETRF function.  P is a 
	 * permutation matrix that is necessary to pass, but will not be used
	 */
	 //TODO:MACROFY? Change the lapack_int, me thinks
	//long long * P = (long long *)fastpca_aligned_alloc(64, m*sizeof( long long ));
	//lapack_int * P = (lapack_int*)fastpca_aligned_alloc(64, m*sizeof( lapack_int));
	lapack_int * P = (lapack_int*)fastpca_aligned_alloc(64, m*sizeof( lapack_int));
	//int * P = (int  *)fastpca_aligned_alloc(64, m*sizeof( long long ));
	//int test =2;
	//TODO: undo this test change
	//long long info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, test, test, A, test,P);
	long long info = LAPACKE_dgetrf( LAPACK_ROW_MAJOR, m, n, 
		A, n, 
		P);
	if (info < 0) {

		return info;
	}

	/*********************************************************************
	 * Note that the GETRF function overwrites A with L and U.  We only want
	 * L, and in its proper form.  The following loop sets the upper triang-
	 * ular portion to 0, and the diagonal to 1
	 */
	for (long long i=0; i<m; i++) {
		for(long long j=i; j<n; j++){
			double newVal = 0.0;
			if (j==i){
				newVal=1;
			}
			A[i*n+j] =newVal;
		}
	}
	fastpca_aligned_free(P);
	return info;
}
int fastpca_qr_q(long long m, long long n, double * A) {

	/*********************************************************************
	 * the ?GEQRF family of functions produces a QR factorization,
	 * but does not explicitly form the matrix Q, which is represented as a 
	 * product of elementary reflectors 
	 */

	double * tau = (double  *)fastpca_aligned_alloc(64, std::min(m,n)*sizeof( double ));
	

	//TODO:MACROFY it
	int info = LAPACKE_dgeqrf( LAPACK_ROW_MAJOR, m,n, A, n, tau);
	//int info = _dgeqrf( LAPACK_ROW_MAJOR, m,n, A, n, tau);
	if ( info <0 ) {
		fastpca_debug_print ("LAPACKE_dgeqrf: %d\n ",info);
		return -1;
	}
	//lwork is recommended to be set to n*blocksize
	long long lwork = n*64;
	info = LAPACKE_dorgqr(LAPACK_ROW_MAJOR, m,n, n, A, n, tau);
	if ( info <0 ) {
		fastpca_debug_print("LAPACKE_dorgqr: %d ",info);
		return -2;
	}

	return 0;
}

void fastpca_square_transpose(double *A, long long m){

//#ifdef USEMKL
		//mkl_dimatcopy('R', 'T', m, m, 1.0, A,m, m);
//#else
//#endif
	for (int i = 0; i < m; i++){
		for (int j =i+1; j<m; j++){
			double temp = A[j*m +i];
			A[j*m+i] = A[i*m+j];
			A[i*m+j] = temp;
		}
	}

}
