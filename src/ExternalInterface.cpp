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

int fastPCAFile (int inputFormat, const char *inputFileName, double **U, double **S, double **V, long long &m, long long &n, long long k, long long l, long long its,long long memory, int centering_row, int centering_column, int dits, int diffsnorm, double& snorm) {
		fastpca_debug_print("Input type is file path: %s", inputFileName);
//fastpca_debug_print("%s", "Input type is file path: ");
		//std::cout <<"Reading file "<< inputFileName << std::endl;
		//std::cout <<k<<"," << l << ","<<memory << "," << its<< inputFileName << std::endl;
		InputMatrix * inputMatrix;
		if (2 == inputFormat) {
		
			inputMatrix = new InputMatrixBedInCore(inputFileName, memory);
			fastpca_debug_print("%s","BED input\n");
			inputMatrix->imputeMissingOn();	


		}else if(3 == inputFormat) {
			fastpca_debug_print("%s","Eigen input\n");
			inputMatrix = new InputMatrixEigenstrat(inputFileName, memory);
			inputMatrix->imputeMissingOn();	

		}else if(1 == inputFormat) {
			fastpca_debug_print("%s","CSV input\n");
			inputMatrix = new InputMatrixCSV(inputFileName, memory);
			
		}else {
			char error [100];
			sprintf(error,"The following input format is not recognized, please use either 'eigen' or 'bed' or 'csv': %d", inputFormat); 
			//fastpca_debug_print("%s", "error");
			//fastpca_debug_print("%s", inputFormat.c_str());
			//::Rf_error(error);
			return(-105);
		}



		fastpca_debug_print("Maximum memory set to be %e\n", (double)memory);
	//Make sure memory is checked
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
			//exit(-1);
		}
	}
	fastpca_debug_print("%s", "\n\nBegin PCA...\n");
	m = (long long) inputMatrix->m;
	n = (long long) inputMatrix->n;
	fastpca_debug_print("Running PCA with k=%lld, m=%e, n=%e, l=%lld, its=%lld\n", k,(double)inputMatrix->m,(double)inputMatrix->n,l,its);
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

 int fastPCAMemory ( double * A, double **U, double **S, double **V, long long m, long long n, long long k, long long l, long long its,  int centering_row, int centering_column, int dits, int diffsnorm, double& snorm) {

	 
	InputMatrix * inputMatrix;

	double * inputMatrixRaw = (double *)fastpca_aligned_alloc( 64, (long long) m*n*sizeof( double ) );

	for (int i =0; i<m; i++){
		for (int j=0;j<n;j++){
			inputMatrixRaw[i*n+j] = A[i*n+j];
		}
	}
	std::string p ("na");

	//You have to pass this so that all the InputMatrix objects are
	//consistent. Obviously, we aren't doing any blocking on the matrices
	//that are passed by memory
	long long memory = m*n*8;
	inputMatrix = new InputMatrixMemory(inputMatrixRaw, memory,m,n,p );
	if (centering_row == 1) {
		inputMatrix->centerRowsOn();	
	}


	int info = inputMatrix->init();
	if (centering_column == 1) {
		inputMatrix->centerColumns();	
	}
	fastpca_debug_print("Centering columns: %d\n\n",centering_column);
	if (info < 0 ) {
		char error [100];
		sprintf(error,"An unknown error has occured in inputMatrix->init(): %d", info); 
		std::cerr <<  error << std::endl;
		return -100;
	}

	for (int test = 0; test <0; test ++ ){
		fastpca_debug_print("%d iteration of Test!\n", test);
		while (inputMatrix->hasNext()){
			fastpca_print_matrix("Test", inputMatrix->blockSize,inputMatrix->n,inputMatrix->block);
		}
	}
	fastpca_debug_print("%s", "Begin PCA...\n");
	fastpca_debug_print("Running PCA with k=%lld, m=%e, n=%e, l=%lld, its=%lld\n", k,(double)inputMatrix->m,(double)inputMatrix->n,l,its);
	//high_resolution_clock::time_point t_pca = high_resolution_clock::now();
	*U = (double *)fastpca_aligned_alloc( 64, (long long) inputMatrix->m*k*sizeof( double ) );
	*S = (double *)fastpca_aligned_alloc( 64, (long long) k*k*sizeof( double ));
	*V = (double *)fastpca_aligned_alloc( 64, (long long) k*inputMatrix->n*sizeof( double ));

	info = fastpca_pca(k,l, inputMatrix, *U,*S,*V,its);
	if (info <0 ) {
		char error [100];
		sprintf(error,"An unknown error has occured in fastpca_pca: %d", info); 
		std::cerr <<  error << std::endl;
		return -101;
	}
	//duration<double> time_span_pca = duration_cast<duration<double>>(high_resolution_clock::now() - t_pca);
	//fastpca_debug_print("The PCA function took %lf seconds\n", time_span_pca.count());
	if (diffsnorm == 1) {
		fastpca_debug_print("%s", "Begin diffsnorm...\n\n");
		info = fastpca_diffsnorm(inputMatrix, *U, *S,*V,k,dits, snorm);
		if (info <0 ) {
			char error [100];
			sprintf(error,"An unknown error has occured in diffsnorm: %d", info); 
			std::cerr <<  error << std::endl;
			return -102;
		}
		fastpca_debug_print("The norm is %le\n", snorm);
		fastpca_debug_print("%s", "End diffsnorm\n\n");
	}

	//fastpca_print_matrix("U", inputMatrix->m,k,U);
	//fastpca_print_matrix("S", k,k,S);
	//fastpca_print_matrix("V", inputMatrix->n,k,V);


return 1;

}
