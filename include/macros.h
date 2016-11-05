#ifndef MACROS_H 
#define MACROS_H
/*
#ifdef PCA_DATATYPE_SINGLE
#define PCA_DATATYPE float
#define PCA_GEMM cblas_sgemm
#define PCA_GEQRF LAPACKE_sgeqrf
#define PCA_ORGQR LAPACKE_sorgqr
#define PCA_GESDD LAPACKE_sgesdd
#define PCA_GETRF LAPACKE_sgetrf
#define PCA_NRM2  cblas_snrm2
#define PCA_IMATCOPY mkl_simatcopy
#else
#define PCA_DATATYPE double
#define PCA_GEMM cblas_dgemm
#define PCA_GEQRF LAPACKE_dgeqrf
#define PCA_ORGQR LAPACKE_dorgqr
#define PCA_GESDD LAPACKE_dgesdd
#define PCA_GETRF LAPACKE_dgetrf
#define PCA_NRM2  cblas_dnrm2
#define PCA_IMATCOPY mkl_dimatcopy
#endif
*/

/*
#ifdef USEMKL
#ifdef OSX
#define MKL_CUSTOM_LIBRARY "libcustom_mkl.dylib"
#else
#define MKL_CUSTOM_LIBRARY "libcustom_mkl.so"
#endif
#else
#define MKL_CUSTOM_LIBRARY "/usr/lib/libblas.so"
#endif
#pragma message "Library is " MKL_CUSTOM_LIBRARY
*/

#ifdef DEBUG
#define DEBUG_TEST 1
#else
#define DEBUG_TEST 0
#endif
#define fastpca_debug_print(fmt, ...) \
	do { if (DEBUG_TEST) fprintf(stdout, "%s:%d: " fmt, __FILE__, \
                                __LINE__, __VA_ARGS__);} while (0)

#endif
