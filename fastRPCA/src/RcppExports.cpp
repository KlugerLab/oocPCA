

#include <string>
#include <Rcpp.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include "ExternalInterface.h"
//void dgemm_test();
/*********************************************************************
 *  Function to provide an interface for performing PCA and calculating the accuracy
 *  of said PCA, on GWAS data sets.
 *  The logic for the mathematical operations are, of course, not
 *  in this file, but rather called from the fastpca_ooc_pca() and fastpca_ooc_diffsnorm()
 *  functions, in their respective files.
 * 
 *  In this file:
 *  1) Parse the R  arguments using the ezOptionParser
 *  3) Calculate the SVD/PCA of the matrix A, such that A= USV', using the fastpca_ooc_pca() function
 *  4) Calculate the accuracy of the decomposition as the L2 norm of A-USV'.
 *  5) Output the resulting matrices as eigenstrat files
 *  
 *  USAGE: result <- fastRPCA("/data/Linderman/bu/FastPCA/LindermanCImplementation_ooc/run/experiments/012115_1/example", k=5, mem=1e+10,l=5)
 
 *  
RcppExport SEXP fastRPCA(SEXP inputFile_temp,SEXP m, SEXP n, SEXP k_temp, SEXP l_temp, SEXP its_temp, SEXP memory_temp){
 *  INPUTS
 ** SEXP  inputFile_temp 
 ** long long  n -- number of columns  in the input matrix
 ** double * A -- The matrix being populated
 *
 *  OUTPUTS
 ** double * A -- The matrix, populated with the random values
 *
 *  RETURN VALUE
 ** int    ier --   error return code
 **               ier=0 means that the function executed successfully


 */  

using namespace Rcpp;
RcppExport SEXP fastRPCA(SEXP inputFormat_temp, SEXP inputMatrix_temp,SEXP m_temp, SEXP n_temp, SEXP k_temp, SEXP l_temp, SEXP its_temp, SEXP memory_temp,SEXP centering_row_temp, SEXP centering_column_temp, SEXP diffsnorm_temp){

	std::string inputFormat = Rcpp::as<std::string>(inputFormat_temp);

	//Rank of the decomposition
	long int k = Rcpp::as<long int>(k_temp);

	//Block size of the normalized power iterations
	long int l = Rcpp::as<long int>(l_temp);

	long int m = Rcpp::as<long int>(m_temp);
	long int n = Rcpp::as<long int>(n_temp);

	//Number of normalized power iterations to conduct
	int its = Rcpp::as<int>(its_temp);

	//Maximum amount of memory
	long int  memory = Rcpp::as<long int >(memory_temp);
	//fastpca_debug_print("Maximum memories set to be %e\n", (double)memory);
	//To calculate the diffsnorm or not
	int diffsnorm = Rcpp::as<int>(diffsnorm_temp);

	//To center or not
	int centering_row = Rcpp::as<int>(centering_row_temp);
	int centering_column = Rcpp::as<int>(centering_column_temp);

	//fastpca_debug_print("Input format: %s", inputFormat.c_str());

	//fastpca_debug_print("%s", "Checking input formats\n");
	//Number of power iterations to conduct for the diffsnorm
	int dits;
	dits = 20;
	std::vector<double> U_out;
	std::vector<double> S_out;
	std::vector<double> V_out;
	std::vector<double> dimensions;
	double *U;
	double * V;
	double *S;
	double snorm;
	long long m_,n_;

	using namespace Rcpp;
	int returnCode;
	if (inputFormat == "memory"){
		Rcpp::NumericVector X(inputMatrix_temp);
		double * inputMatrixRaw = X.begin();
		returnCode = fastPCAMemory(inputMatrixRaw, &U, &S, &V, m,n,k, l, its, centering_row, centering_column, dits, diffsnorm, snorm);
		m_ = m;
		n_ = n;

	}else if (inputFormat=="csv"){
		std::string inputMatrixPath= Rcpp::as<std::string>(inputMatrix_temp);
		returnCode=fastPCAFile(1, inputMatrixPath.c_str(), &U, &S, &V, m_,n_,k, l, its,memory, centering_row, centering_column, dits, diffsnorm, snorm);

	}else {
		std::string inputMatrixPath= Rcpp::as<std::string>(inputMatrix_temp);
		returnCode=fastPCAFile(2, inputMatrixPath.c_str(), &U, &S, &V, m_,n_,k, l, its,memory, centering_row, centering_column, dits, diffsnorm, snorm);
	}

	switch (returnCode){
		case -100:
			::Rf_error("An unknown error has occurred in inputMatrix->init()");
		break;
		case -101:
			::Rf_error("An unknown error has occurred in fastpca_pca()");
		break;
		case -102:
			::Rf_error("An unknown error has occurred in diffsnorm()");
		break;
				
	}


	std::map<std::string,std::vector<double> > outs;

	dimensions.push_back((double)m_);
	dimensions.push_back((double)n_);
	for (int p=0; p<k*m_; p++){
		U_out.push_back(U[p]);
	}
	for (int p=0; p<k*k; p++){
		S_out.push_back(S[p]);
	}
	for (int p=0; p<k*n_; p++){
		V_out.push_back(V[p]);
	}
	
	outs["U"] = U_out;
	outs["S"] = S_out;
	outs["V"] = V_out;
	outs["dims"] = dimensions;
	if (diffsnorm == 1) {
		std::vector<double> diffsnorm_vec;
		diffsnorm_vec.push_back(snorm);
		outs["diffsnorm"] = diffsnorm_vec;
	}

	return Rcpp::wrap(outs);

/*
	InputMatrix * inputMatrix;

	if ("memory" == inputFormat) {
		Rcpp::NumericVector X(inputMatrix_temp);
		double * inputMatrixRaw = X.begin();
		std::string p ("%s", "na");
		//Make sure memory is checked
		inputMatrix = new InputMatrixMemory(inputMatrixRaw, memory,m,n,p );

	} else{
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
	}

	if (centering_row == 1) {
		inputMatrix->centerRowsOn();	
	}


	int info = inputMatrix->init();
	if (centering_column == 1) {
		inputMatrix->centerColumns();	
	}
	if (info < 0 ) {
		char error [100];
		sprintf(error,"An unknown error has occured in inputMatrix->init(): %d", info); 
		::Rf_error(error);
	}

	for (int test = 0; test <0; test ++ ){
		fastpca_debug_print("%d iteration of Test!\n", test);
		while (inputMatrix->hasNext()){
			fastpca_print_matrix("Test", inputMatrix->blockSize,inputMatrix->n,inputMatrix->block);
		}
	}
	fastpca_debug_print("%s", "\n\nBegin PCA...\n");
	fastpca_debug_print("Running PCA with k=%lld, m=%e, n=%e, l=%lld, its=%d\n", k,(double)inputMatrix->m,(double)inputMatrix->n,l,its);
	high_resolution_clock::time_point t_pca = high_resolution_clock::now();
	double * U = (double *)fastpca_aligned_alloc( 64, (long long) inputMatrix->m*k*sizeof( double ) );
	double * S = (double *)fastpca_aligned_alloc( 64, (long long) k*k*sizeof( double ));
	double * V = (double *)fastpca_aligned_alloc( 64, (long long) k*inputMatrix->n*sizeof( double ));

	info = fastpca_pca(k,l, inputMatrix, U,S,V,its);
	if (info <0 ) {
		char error [100];
		sprintf(error,"An unknown error has occured in fastpca_pca: %d", info); 
		::Rf_error(error);
	}
	duration<double> time_span_pca = duration_cast<duration<double>>(high_resolution_clock::now() - t_pca);
	fastpca_debug_print("The PCA function took %lf seconds\n", time_span_pca.count());
	double snorm;
	if (diffsnorm == 1) {
		fastpca_debug_print("%s", "\n\nBegin diffsnorm...\n");
		info = fastpca_diffsnorm(inputMatrix, U, S,V,k,dits, snorm);
		if (info <0 ) {
			char error [100];
			sprintf(error,"An unknown error has occured in diffsnorm: %d", info); 
			::Rf_error(error);
		}
		fastpca_debug_print("The norm is %le\n", snorm);
		fastpca_debug_print("%s", "End diffsnorm\n\n");
	}

	//fastpca_print_matrix("U", inputMatrix->m,k,U);
	//fastpca_print_matrix("S", k,k,S);
	//fastpca_print_matrix("V", inputMatrix->n,k,V);
	*/
} 


