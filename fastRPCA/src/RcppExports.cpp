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


//fastRPCA is the function that is called by using the R .Call interface in the RcppExports.R file.  This function then proceeds to call the respective function in the ExternalInterace file
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

} 


