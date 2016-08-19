#include "omp.h"
#include "ExternalInterface.h"
#include "fastpca_io.hpp"
#include "fastpca_pca.hpp"
#include "fastpca_diffsnorm.hpp"
#include "InputMatrix.h"
#include "InputMatrixMemory.h"
#include "InputMatrixEigenstrat.h"
#include "InputMatrixGeneralBinary.h"
#include "InputMatrixBedInCore.h"
#include "InputMatrixEigenInCore.h"
#include "InputMatrixCSV.h"

//int fastPCAFile ( std::string inputFormat, std::string inputFileName, long long m, long long n, long long k, long long l, long long its, long long memory, int centering_row, int centering_column, int its, int diffsnorm ) {
/*(
		fastpca_debug_print("%s", "Input type is a file path\n");
		std::string inputFile = Rcpp::as<std::string>(inputMatrix_temp);
		std::cout <<"Reading file "<< inputFile << std::endl;
		std::cout <<k<<"," << l << ","<<memory << "," << its<< inputFile << std::endl;
		if ("bed" == inputFormat) {
		
			inputMatrix = new InputMatrixBedInCore(inputFile, memory);
			fastpca_debug_print("%s","BED input\n");
			inputMatrix->imputeMissingOn();	


		}else if("eigen" == inputFormat) {
			fastpca_debug_print("%s","Eigen input\n");
			inputMatrix = new InputMatrixEigenstrat(inputFile, memory);
			inputMatrix->imputeMissingOn();	

		}else if("csv" == inputFormat) {
			fastpca_debug_print("%s","CSV input\n");
			inputMatrix = new InputMatrixCSV(inputFile, memory);
			
		}else {
			char error [100];
			sprintf(error,"The following input format is not recognized, please use either 'eigen' or 'bed' or 'csv': %s", inputFormat.c_str()); 
			::Rf_error(error);
		}
		*/

 int fastPCAMemory ( double * A, double **U, double **S, double **V, long long m, long long n, long long k, long long l, long long its, long long memory, int centering_row, int centering_column, int dits, int diffsnorm, double& snorm) {

	fastpca_debug_print("%s","We are in\n");
	InputMatrix * inputMatrix;

	//double * inputMatrixRaw = A;
	double * inputMatrixRaw = (double *)fastpca_aligned_alloc( 64, (long long) m*n*sizeof( double ) );

	for (int i =0; i<m; i++){
		for (int j=0;j<n;j++){
			inputMatrixRaw[i*n+j] = A[i*n+j];
		}
	}
	std::string p ("na");
	fastpca_debug_print("%s","realloced\n");

	fastpca_debug_print("Maximum memory set to be %e\n", (double)memory);
	//Make sure memory is checked
	inputMatrix = new InputMatrixMemory(inputMatrixRaw, memory,m,n,p );
	if (centering_row == 1) {
		inputMatrix->centerRowsOn();	
	}


	int info = inputMatrix->init();
	fastpca_debug_print("%s"," Inited it");
	if (centering_column == 1) {
		inputMatrix->centerColumns();	
	}
	if (info < 0 ) {
	//TODO: Add error
		char error [100];
		sprintf(error,"An unknown error has occured in inputMatrix->init(): %d", info); 
		return -100;
//		::Rf_error(error);
	}

	for (int test = 0; test <0; test ++ ){
		fastpca_debug_print("%d iteration of Test!\n", test);
		while (inputMatrix->hasNext()){
			fastpca_print_matrix("Test", inputMatrix->blockSize,inputMatrix->n,inputMatrix->block);
		}
	}
	fastpca_debug_print("%s", "\n\nBegin PCA...\n");
	fastpca_debug_print("Running PCA with k=%lld, m=%e, n=%e, l=%lld, its=%d\n", k,(double)inputMatrix->m,(double)inputMatrix->n,l,its);
	//high_resolution_clock::time_point t_pca = high_resolution_clock::now();
	*U = (double *)fastpca_aligned_alloc( 64, (long long) inputMatrix->m*k*sizeof( double ) );
	*S = (double *)fastpca_aligned_alloc( 64, (long long) k*k*sizeof( double ));
	*V = (double *)fastpca_aligned_alloc( 64, (long long) k*inputMatrix->n*sizeof( double ));

	info = fastpca_pca(k,l, inputMatrix, *U,*S,*V,its);
	if (info <0 ) {
	//TODO: error here
		char error [100];
		sprintf(error,"An unknown error has occured in fastpca_pca: %d", info); 
		return -101;
		//#::Rf_error(error);
	}
	//duration<double> time_span_pca = duration_cast<duration<double>>(high_resolution_clock::now() - t_pca);
	//fastpca_debug_print("The PCA function took %lf seconds\n", time_span_pca.count());
	if (diffsnorm == 1) {
		fastpca_debug_print("%s", "\n\nBegin diffsnorm...\n");
		info = fastpca_diffsnorm(inputMatrix, *U, *S,*V,k,dits, snorm);
		if (info <0 ) {
			char error [100];
			sprintf(error,"An unknown error has occured in diffsnorm: %d", info); 
			return -102;
			//::Rf_error(error);
		}
		fastpca_debug_print("The norm is %le\n", snorm);
		fastpca_debug_print("%s", "End diffsnorm\n\n");
	}

	//fastpca_print_matrix("U", inputMatrix->m,k,U);
	//fastpca_print_matrix("S", k,k,S);
	//fastpca_print_matrix("V", inputMatrix->n,k,V);


return 1;

}
