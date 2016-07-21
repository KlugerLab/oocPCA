#include "fastpca_diffsnorm.hpp"

int fastpca_diffsnorm(  InputMatrix* A, double* U, double* S , double* V,long long k, int its, double &snorm ) {

	//fastpca_print_matrix("Test", A->blockSize,A->n,A->block);
	//For convenience
	int m = A->m;
	int n = A->n;
	/*********************************************************************
	 *  Check input.  Please see function comment for more information.
	if ((m <0) || ( n<0) || (k<0) || (its<0)) {
		return -1;
	}
	 */

	if ( k > m || k > n) {
		return -2;
	}

//	if (m>=n) {

		/*********************************************************************
		 *  Generate a random vector x of size n
		 */
		double * x = (double *)fastpca_aligned_alloc(64, n*1*sizeof( double ));
		if (NULL == x) {
			return -30;
		}
		fastpca_populate_matrix_random(n,1,x);

		double norm = cblas_dnrm2(n, x, 1);

		fastpca_divide_matrix_scalar(n,1,x,norm);

		/*********************************************************************
		 *  Run iterations of the power method, starting with the randm x
		 */
		fastpca_debug_print("%s", "Begin calculating the L2 norm of the difference between A and USV'\n");
		for (int i=0; i<its; i++) {
			//printf("Power iteration i=%d\n", i);

			/*********************************************************************
			 *  Set y = (A-USV')x
			 */

			//Matrix product: y=V'*x
			//Dimensions: (k,1)=(n,k)'*(n,1)
			double * y = (double *)fastpca_aligned_alloc(64, 1*n*sizeof( double ));
			if (NULL == y) {
				return -30;
			}
			fastpca_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, k,1,n, 1.0, V,k, x,1, 0.0, y,1);

			//Matrix product: y=S*y
			//Dimensions: (k,1)=(k,k)*(k,1)
			double * y_temp = (double *)fastpca_aligned_alloc(64, m*1*sizeof( double ));
			if (NULL == y_temp) {
				return -30;
			}
			fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, k,1,k, 1.0, S,k, y,1, 0.0, y_temp,1);
			fastpca_aligned_free(y);
			y = y_temp;

			//Matrix product: y=U*y
			//Dimensions: (m,1)=(m,k)*(k,1)
			y_temp = (double *)fastpca_aligned_alloc(64, m*1*sizeof( double ));
			if (NULL == y_temp) {
				return -30;
			}
			fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m,1,k, 1.0, U,k, y,1, 0.0, y_temp,1);
			fastpca_aligned_free(y);
			y = y_temp;

			//Matrix product: y=A*x -y
			//Dimensions: (m,1)=(m,n)*(n,1) - (m,1)
			while (A->hasNext()){
				fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasTrans, A->blockSize,1,n, 1.0, A->block,n, x,n, -1.0, &(y[A->blockStart]),1);
			}

			/*********************************************************************
			 * Set x = (A'-VS'U')y
			 */

			//Matrix product: x=U'*y 
			//Dimensions: (k,1)=(m,k)'*(m,1)
			x = (double *)fastpca_aligned_alloc(64, 1*m*sizeof( double ));
			fastpca_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, k,1,m, 1.0, U,k, y,1, 0.0, x,1);

			//TODO: Make sure this is okay.
			double* x_temp = (double *)fastpca_aligned_alloc(64, 1*k*sizeof( double ));
			if (NULL == x_temp) {
				return -30;
			}

			//Matrix product: x=S'*x 
			//Dimensions: (k,1)=(k,k)'*(k,1)
			fastpca_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, k,1,k, 1.0, S,k, x,1, 0.0, x_temp,1);
			fastpca_aligned_free(x);
			x = x_temp;

			//Matrix product: x=V*x 
			//Dimensions: (n,1)=(n,k)*(k,1)
			x_temp = (double *)fastpca_aligned_alloc(64, n*1*sizeof( double ));
			if (NULL == x_temp) {
				return -30;
			}
			fastpca_gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n,1,k, 1.0, V,k, x,1, 0.0, x_temp,1);
			fastpca_aligned_free(x);
			x = x_temp;


			//Matrix product: x = A'*y -x
			//Dimensions: (n,1)=(m,n)'*(m,1) - (n,1)

			//Use a temporary matrix, because we need to add the
			//results of each block, so we cannot directly be
			//subtracting from x
			double* Aty_temp = (double *)fastpca_aligned_alloc(64, 1*n*sizeof( double ));
			for (int i=0; i<n; i++ ){
				Aty_temp[i] = 0;
			}
			if (NULL == Aty_temp) {
				return -30;
			}
			while (A->hasNext()){
				fastpca_gemm(CblasRowMajor, CblasTrans, CblasNoTrans, A->n,1,A->blockSize, 1.0, A->block,A->n,&(y[A->blockStart]),1, 1.0,Aty_temp ,1);
			}
			for (int i=0; i<n; i++ ){
				x[i] = Aty_temp[i] - x[i];
			}

			norm = cblas_dnrm2(n, x, 1);
			fastpca_divide_matrix_scalar(n,1,x,norm);

		}
		snorm = sqrt(norm);
	//}
	return 0;
}

