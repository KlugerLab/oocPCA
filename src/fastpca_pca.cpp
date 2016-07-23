#include <chrono>
#include "fastpca_pca.hpp"
//#include "mkl_cblas_gcl.h"
int fastpca_pca( long long k, long long l, 
		InputMatrix * A, double* U, double * S, double *V, long long its){
using namespace std::chrono;
 
	int m = A->m;
	int n = A->n;

	/*********************************************************************
	 *  Check input.  Please see function comment for more information.
	 */
	if (k< 0 ){
		return -1;
	}
	if ( (k> m) || (k>n) ){
		return -2;
	}
	if ( its <0){
		return -3;
	}
	if (l <k){
		return -4;
	}
	if (l >n ){
		return -99;
		//make new one
	}
	if (l >m ){
		return -99;
		//make new one
	}
	if (A==NULL ) {
		return -5;
	}
	if (U==NULL ) {
		return -6;
	}
	if (S==NULL ) {
		return -7;
	}
	if (V==NULL ) {
		return -8;
	}



	/*********************************************************************
	 *  Apply A to a random (n,l) matrix, Ran, obtaining Q (m,l)
	 */

	double * Ran = (double *)fastpca_aligned_alloc(64, m*l*sizeof( double ));
	//if (NULL == Ran) {
		//return -30;
	//}
	fastpca_populate_matrix_random (m,l,Ran);
	//Matrix product: Q=A'*Ran
	//Dimensions: (n,l)=(m,n)'*(m,l)
	double * Q = (double *)fastpca_aligned_alloc(64, n*l*sizeof( double ));
	for (int i = 0; i < n*l; i++ ){
		Q[i] = 0;
	}
	if (NULL == Q) {
		return -30;
	}
	while (A->hasNext()){
		fastpca_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, A->n,l,A->blockSize, 1.0, A->block,A->n,&(Ran[A->blockStart*l]),l, 1.0,Q ,l);

	}


	/*********************************************************************
	 *  Form a matrix Q whose columns constitute a well-conditioned basis 
	 *  for the columns of the earlier Q
	 */
	if (its>0) {
		int info = fastpca_lu_l(n,l,Q);
		if (info < 0 ){
			fastpca_debug_print("We got an error of %d\n", info);
			return -22;
		}
	}else {
		int info = fastpca_qr_q(n,l,Q);
		if (info < 0 ){
			fastpca_debug_print("We got an error of %d\n", info);
			return -21;
		}
	}



	/*********************************************************************
	 *  Conduct its normalized power iteration.  Each LU is a normalization operation.
	 */
	for (int i=0; i<its; i++) {
		//Matrix product: Q=A'*AQ
		//Dimensions: (m,l)=(m,n)*(n,l)
		double * Q_temp = (double *)fastpca_aligned_alloc(64, m*l*sizeof( double ));
		for (int i = 0; i < m*l; i++ ){
			Q_temp[i] = 0;
		}
		if (NULL == Q_temp) {
			return -30;
		}
		double * Q_temp2 = (double *)fastpca_aligned_alloc(64, n*l*sizeof( double ));
		for (int i = 0; i < n*l; i++ ){
			Q_temp2[i] = 0;
		}
		if (NULL == Q_temp2) {
			return -30;
		}
		while (A->hasNext()){
			//A(blockSize,n)*Q (n,l)
			//printf("%d, %d,%d\n", A->blockStart, A->blockSize, l); 
			fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->blockSize,l,A->n, 1.0, A->block,A->n,Q,l, 0.0,&(Q_temp[A->blockStart*l]) ,l);
			//fastpca_print_matrix("Error Debug2", A->blockSize,l, Q_temp+A->blockStart*l);
			//fastpca_print_matrix("Error Debug3", m,l, Q_temp);
		}
		int info = fastpca_lu_l(m,l,Q_temp);
		if (info < 0 ){
			//printf("m:%d, n:%d, l:%d, blockSize %d\n", m,n,l,A->blockSize);
			fastpca_debug_print("We got an error of %d\n", info);
			//fastpca_print_matrix("Error Debug2", A->n,l, Q);
			//fastpca_print_matrix("Error Debug3", m,l, Q_temp);
			return -23;
		}
		while (A->hasNext()){
			//A(blockSize,n)*Q (n,l)
			fastpca_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, A->n,l,A->blockSize, 1.0, A->block,A->n,&(Q_temp[A->blockStart*l]),l, 1.0,Q_temp2 ,l);
		}
		/*
		while (A->hasNext()){
			//A(blockSize,n)*Q (n,l)
			fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->blockSize,l,A->n, 1.0, A->block,A->n,Q,l, 0.0,Q_temp ,l);
			int info = fastpca_lu_l(A->blockSize,l,Q_temp);
			if (info < 0 ){
				//printf("m:%d, n:%d, l:%d, blockSize %d\n", m,n,l,A->blockSize);
				fastpca_debug_print("We got an error of %d\n", info);
				//fastpca_print_matrix("Error Debug1", A->blockSize,A->n,A->block);
				//fastpca_print_matrix("Error Debug2", A->n,l, Q);
				//fastpca_print_matrix("Error Debug3", m,l, Q_temp);
				return -23;
			}
			
			fastpca_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, A->n,l,A->blockSize, 1.0, A->block,A->n,Q_temp,l, 1.0,Q_temp2 ,l);
		}
		*/
		fastpca_aligned_free(Q);
		Q=Q_temp2;

		
		if (i +1 == its){
			int info = fastpca_qr_q(n,l,Q);
			if (info < 0 ){
				return -20;
			}
		}else{
			int info = fastpca_lu_l(n,l,Q);
			if (info < 0 ){
				fastpca_debug_print("We got an error of %d\n", info);
				return -23;
			}
		}
	}

	/*********************************************************************
	 *  SVD A*Q to obtain approximations to the singular values 
	 *  and left signular values of A;  The right singular vectors
	 *  will be adjusted in the next step
	 */

	//Matrix product: AQ=A*Q
	//Dimensions: (m,l)=(m,n)*(n,l)
	double * AQ = (double *)fastpca_aligned_alloc(64, m*l*sizeof( double ));
	if (NULL == AQ) {
		return -30;
	}

	while (A->hasNext()){
		fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A->blockSize,l,A->n, 1.0, A->block,A->n,Q,l, 0.0,&(AQ[A->blockStart * l]) ,l);
	}

	//U*s*Rt = AQ
	//U is m,l, Rt is l,l
	double * s = (double *)fastpca_aligned_alloc(64, l*sizeof( double ));
	if (NULL == s) {
		return -30;
	}
	double * Rt = (double *)fastpca_aligned_alloc(64, l*l*sizeof( double ));
	if (NULL == Rt) {
		return -30;
	}
	double * U_l = (double *)fastpca_aligned_alloc(64, m*l*sizeof( double ));
	if (NULL == U_l) {
		return -30;
	}
	int info = LAPACKE_dgesdd( LAPACK_ROW_MAJOR, 'S', m, l, AQ, l, s, U_l, l, Rt, l );

	if( info > 0 ) {
		fastpca_debug_print( "%s", "The algorithm computing SVD failed to converge.\n" );
		exit( 1 );
	}

	//Transpose Vt, as the output of blas' SVD is V', not V
	fastpca_square_transpose(Rt, l);

	/*********************************************************************
	 *  The right singular vectors of AQ, R, must be adjusted in order to be equal to
	 *  V, the right singular vectors of A.  V= Q*R.
	 */

	//Matrix product: V_l=Q*R
	//Dimensions: (n,l)=(n,l)*(l,l)
	double * V_l = (double *)fastpca_aligned_alloc(64, n*l*sizeof( double ));

	fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n,l,l, 1.0, Q,l, Rt,l, 0.0, V_l,l);

	/*********************************************************************
	 * Retain only the leftmost k columns of U,
	 * the leftmost k columns of V, and the
	 * uppermost leftmost k x k block of S
	 */
	for (int i=0; i<m; i++){
		for (int j=0; j<k; j++) {
			U[i*k+j] = U_l[i*l+j];
		}
	}

	for (int i=0; i<k; i++){
		for (int j=0; j<k; j++) {
			if (i==j){
				S[i*k+j] = s[i];
			}else{
				S[i*k+j] = 0;
			}
		}
	}
	for (int i=0; i<n; i++){
		for (int j=0; j<k; j++) {
			V[i*k+j] = V_l[i*l+j];
		}
	}

	return 1;



}


