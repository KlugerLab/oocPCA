#include "macros.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fastpca_io.hpp"
#include "fastpca_pca.hpp"
#include "fastpca_diffsnorm.hpp"
#include <chrono>
//#include "omp.h"
#include "ezOptionParser.hpp"
#include <iostream>
#include <fstream>
#include <limits.h>
#include "InputMatrix.h"
#include "InputMatrixEigenstrat.h"
#include "InputMatrixGeneralBinary.h"
#include "InputMatrixBedInCore.h"
#include "InputMatrixEigenInCore.h"
#include "InputMatrixCSV.h"
#include <iostream>
#include <fstream>
//export DYLD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib/:/opt/intel/mkl/lib
/*********************************************************************
 *  main() function to provide an interface for performing PCA and calculating the accuracy
 *  of said PCA, on GWAS data sets.
 *  The logic for the mathematical operations are, of course, not
 *  in this file, but rather called from the fastpca_ooc_pca() and fastpca_ooc_diffsnorm()
 *  functions, in their respective files.
 * 
 *  In this file:
 *  1) Parse the command line arguments using the ezOptionParser
 *  3) Calculate the SVD/PCA of the matrix A, such that A= USV', using the fastpca_ooc_pca() function
 *  4) Calculate the accuracy of the decomposition as the L2 norm of A-USV'.
 *  5) Output the resulting matrices as eigenstrat files
 *  
 *  USAGE: fastpca_ooc.x [OPTIONS]
 *  
 *  PARAMETERS:
 *  -block, -l, --block ARG       Block size of the normalized power iterations; l
 *                                must be positive integer >=k.  Defaults to k+2
 *  
 *  -c, -center, --center         Center the rows of the matrix
 *
 *  -d, -diffnorm, --diffnorm     Calculate the L2 norm of A-USV as a measure of the
 *                                accuracy of the decomposition
 *  
 * 
 *  -k, -rank, --rank ARG         The rank of the decomposition
 *  
 *  -p, -pits, --pits ARG         Number of normalized power iterations iterations
 *                                for the PCA
 *  
 *  -h, -help, --help, --usage    Display usage instructions.
 * 
 *  OUTPUT FORMATS:
 *  -binaryO, -binaryOut, --binaryOut ARG       Output results as ARG.U, ARG.S, ARG.V
 *  
 *  -eigenO, -eigenOutput, --eigenOutput ARG     Output file ARG, in eigenstrat format
 *
 **  INPUT FORMATS:
 *  
 *  -eigenI, -eigenInput, --eigenInput ARG       Input file, in eigenstrat format
 *  
 *  -bedI, -bedInput, --bedInput ARG       Binary Input files, as ARG.bed,ARG.bim,ARG.fam mode
 *  
 *  -binaryI, -binaryInput, --binaryInput ARG    General binary input ARG.  Row-major output of doubles.  Most provide -m and -n parameters
 *  
 *  -m, -rows, --rows ARG         Number of rows (only applicable for -binaryInput run)
 *  
 *  -cols, -n, --cols ARG         Number of columns (only applicable for -binaryInput run)
 *
 * 
 *   DIAGNOSTIC: Not required for typical usage
 * 
 *  -bina, -binaryA, --binaryA ARG       Output the preprocessed matrix as ARG
 *  -time, --time ARG       Output the time for PCA to file ARG
 *  -t, -threads, --threads ARG   Number of threads to be used by Intel MKL. 
 *                                Recommended to not set, and Intel MKL will
 *                                dynamically determine optimal thread number
 *
 *  EXAMPLES:
 *  
 *  ../bin/fastpca_ooc.x -k 5 -l 5 -eigenInput /home/george/software/EIG5.0.2/EIGENSTRAT/example.geno -eigenOutput test.pca -d -mem 140 
 * 
 *  ../bin/fastpca_ooc.x -k 20 -bedI /data/GERA_DECRYPTED/LindermanAnalysis/EUR/Benchmark/eur_1000_62318/eur_1000_62318 -mem 140000000 -c -binaryOut eur_1000_62318_k20.9.bim.binmatrix -d
 * 
 *  ../bin/generate_test_matrix.x -k 20 -m 1000 -n 3000 -o /data/BigDataSandbox/1000x3000.binmatrix
 *  ../bin/fastpca_ooc.x -m 1000 -n 3000 -k 20 -binaryI /data/BigDataSandbox/1000x3000.binmatrix -mem 140000000 -binaryOutput 1000x3000.binmatrix -d
 * 
 *  
 */
using namespace std::chrono;
   int main(int argc, const char * argv[]) {
   
   	/***********************************************************
   	 *  Declare variables for the input parameters. (Not all must be
   	 *  passed.  For example, m and n are calculated from the input
   	 *  matrix.
   	 */
   
	//Maximum amount of memory
	long long memory;

	//Rank of the decomposition
	long long k;

	InputMatrix * inputMatrix;

	//Block size of the normalized power iterations
	long long l; 


	//Number of normalized power iterations to conduct
	int its;

	//Number of power iterations to conduct for the diffsnorm
	int dits;
	

	/***********************************************************
	 *  Define each input parameters with ezOptionParser
	 */
	using namespace ez;
	ezOptionParser opt;

	opt.overview = "Fast PCA implementation";
	opt.syntax = "fastpca.x [OPTIONS]";
	opt.example = "./fastpca.x -k 5 -eigenInput /home/george/software/EIG5.0.2/EIGENSTRAT/example.geno  -o test.pca -d\n\n./fastpca.x -e -k 10 -m 400000 -n 10000 -d\n\n";
	opt.footer = "This program is free and without warranty.\n";

	opt.add(
			"", // Default.
			0, // Required?
			0, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Display usage instructions.", // Help description.
			"-h",     // Flag token. 
			"-help",  // Flag token.
			"--help", // Flag token.
			"--usage" // Flag token.
	       );
	opt.add(
			"10", // Default.
			1, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"The rank of the decomposition", // Help description.
			"-k", // Flag token.
			"-rank", // Flag token.
			"--rank" // Flag token.
	       );
	opt.add(
			"1000", // Default.
			1, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"The max amount of memory the software can use", // Help description.
			"-mem", // Flag token.
			"--memory" // Flag token.
	       );
	opt.add(
			"12", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Number of rows (only applicable for example run)", // Help description.
			"-m", // Flag token.
			"-rows", // Flag token.
			"--rows" // Flag token.
	       );
	opt.add(
			"10", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Number of columns (only applicable for example run)", // Help description.
			"-n", // Flag token.
			"-cols", // Flag token.
			"--cols" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Block size of the normalized power iterations; l must be positive integer >=k.  Defaults to k+2", // Help description.
			"-l", // Flag token.
			"-block", // Flag token.
			"--block" // Flag token.
	       );
	opt.add(
			"2", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Number of normalized power iterations iterations for the PCA", // Help description.
			"-p", // Flag token.
			"-pits", // Flag token.
			"--pits" // Flag token.
	       );
	opt.add(
			"20", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Number of power iterations iterations for the diffsnorm", // Help description.
			"-dits", // Flag token.
			"--dits" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Number of threads to be used by Intel MKL.  Recommended to not set, and Intel MKL will dynamically determine optimal thread number", // Help description.
			"-t", // Flag token.
			"-threads", // Flag token.
			"--threads" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Output as matlab file", // Help description.
			"-binaryO",     // Flag token. 
			"-binaryOutput",  // Flag token.
			"--binaryOutput" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Output the preprocessed matrix as a binary file", // Help description.
			"-binA",     // Flag token. 
			"-binaryA",  // Flag token.
			"--binaryA" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			0, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Center the columns of the matrix", // Help description.
			"-cc",     // Flag token. 
			"-ccenter",  // Flag token.
			"--ccenter" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			0, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Center the rows of the matrix", // Help description.
			"-c",     // Flag token. 
			"-center",  // Flag token.
			"--center" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			0, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Calculate the L2 norm of A-USV as a measure of the accuracy of the decomposition", // Help description.
			"-d",     // Flag token. 
			"-diffnorm",  // Flag token.
			"--diffnorm" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Input file, in eigenstrat format", // Help description.
			"-eigenI",     // Flag token. 
			"-eigenInput",  // Flag token.
			"--eigenInput" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Input file, in eigenstrat plaintext format, calculations in core", // Help description.
			"-eigenInCoreI",     // Flag token. 
			"-eigenInCoreInput",     // Flag token. 
			"--eigenInCoreInput"     // Flag token. 
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Input file, in eigenstrat binary format", // Help description.
			"-bedI",     // Flag token. 
			"-bedInput",     // Flag token. 
			"--bedInput"     // Flag token. 
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Input file, in general binary format", // Help description.
			"-binaryI",     // Flag token. 
			"-binaryInput",  // Flag token.
			"--binaryInput" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Input file, in CSV format", // Help description.
			"-csvI",     // Flag token. 
			"-csvInput",  // Flag token.
			"--csvInput" // Flag token.
	       );
	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Output file, in eigenstrat format", // Help description.
			"-eigenO",     // Flag token. 
			"-eigenOutput",  // Flag token.
			"--eigenOutput" // Flag token.
	       );

	opt.add(
			"", // Default.
			0, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Output the elapsed time for PCA to this file", // Help description.
			"-time",     // Flag token. 
			"--time" // Flag token.
	       );

	/***********************************************************
	 *  Parse the command line arguments.  
	 */
	opt.parse(argc, argv);

	// -h flag is for help
	if (opt.isSet("-h")) {
		std::string usage;
		opt.getUsage(usage);
		std::cout << usage;
		return -1;
	}

	std::vector<std::string> badOptions;
	int i;
	if(!opt.gotRequired(badOptions)) {
		for(i=0; i < badOptions.size(); ++i)
			std::cerr << "ERROR: Missing required option " << badOptions[i] << ".\n\n";
		std::string usage;
		opt.getUsage(usage);
		std::cout << usage;
		return -1;
	}

	if(!opt.gotExpected(badOptions)) {
		for(i=0; i < badOptions.size(); ++i)
			std::cerr << "ERROR: Got unexpected number of arguments for option " << badOptions[i] << ".\n\n";
		std::string usage;
		opt.getUsage(usage);
		std::cout << usage;
		return -1;
	}

	if (opt.isSet("-mem")) {
		double memory_temp;
		opt.get("-mem")->getDouble(memory_temp);
		memory=(long long) memory_temp;
		if ( memory < 0 ) {
			std::cerr << "ERROR: memory cannot be negative ";
			return -1;
		}
		fastpca_debug_print("Maximum memory set to be %e\n", (double)memory);
	}
	if (opt.isSet("-k")) {
		int k_temp;
		opt.get("-k")->getInt(k_temp);
		k=k_temp;
		if ( k < 0 ) {
			std::cerr << "ERROR: k cannot be negative ";
			return -1;
		}
	}
	if (opt.isSet("-l")) {
		int l_temp;
		opt.get("-l")->getInt(l_temp);
		l = l_temp;
		if ( l < 0 ) {
			std::cerr << "ERROR: l cannot be negative ";
			return -1;
		}
	}else {
		l = k+2;
	}

	//For general binary input, m and n must be passed.  Otherwise, the m and n are
	//determined by the size of the input matrix
	if (opt.isSet("-binaryI")) {
		if (!opt.isSet("-m") || !opt.isSet("-n")){
			std::cerr << "ERROR: -m and -n must be set when passing a general binary file! " << std::endl;
			return -1;
		}
	}
	if (opt.isSet("-p")) {
		opt.get("-p")->getInt(its);
	}else {
		its = 2;
	}
	if (opt.isSet("-dits")) {
		opt.get("-dits")->getInt(dits);
	}else {
		dits = 20;
	}

	if (opt.isSet ("-e") && opt.isSet("-i")){
		std::cerr << "ERROR: cannot run an example and provide input! ";
		return -1;
	}

	if (opt.isSet ("-e") && opt.isSet("-binaryI")){
		std::cerr << "ERROR: cannot run an example and provide input! ";
		return -1;
	}

	if (opt.isSet ("-i") && opt.isSet("-binaryI")){
		std::cerr << "ERROR: cannot provide both an EIGENSTRAT GENO input and a binary input! ";
		return -1;
	}

	//Manually set the number of threads used by MKL
	if (opt.isSet("-t")) {
		/*
		int t;
		opt.get("-t")->getInt(t);
		mkl_set_num_threads(t);
		mkl_set_dynamic(0);
		fastpca_debug_print("The number of threads has been switched from dynamic to %d.  Please use with caution.\n", t);
		*/
		fastpca_debug_print("%s", "Thread setting is not yet available.");
	}



	/*********************************************************************
	 *  Read Eigenstrat Input.  
	 */
	if (opt.isSet("-eigenI")) {

		//Read the file
		std::string inputFile;
		opt.get("-eigenI")->getString(inputFile);
		fastpca_debug_print("Reading file: %s\n", inputFile.c_str());

		inputMatrix = new InputMatrixEigenstrat(inputFile, memory);
		inputMatrix->imputeMissingOn();	
	/*********************************************************************
	 *  Read Eigenstrat Input.  
	 */
	}else if (opt.isSet("-csvI")) {

		//Read the file
		std::string inputFile;
		opt.get("-csvI")->getString(inputFile);
		fastpca_debug_print("Reading file: %s\n", inputFile.c_str());

		inputMatrix = new InputMatrixCSV(inputFile, memory);

	/*********************************************************************
	 *  Read text input, in core 
	 */
	}else if (opt.isSet("-eigenInCoreI")) {
		std::string inputFile;
		opt.get("-eigenInCoreI")->getString(inputFile);
		fastpca_debug_print("Reading file: %s\n", inputFile.c_str());

		inputMatrix = new InputMatrixEigenInCore(inputFile, memory);
		inputMatrix->imputeMissingOn();	

	/*********************************************************************
	 *  Read Binary Input.  
	 */
	}else if (opt.isSet("-bedI")) {
		std::string inputFile;
		opt.get("-bedI")->getString(inputFile);
		fastpca_debug_print("Reading file: %s\n", inputFile.c_str());

		inputMatrix = new InputMatrixBedInCore(inputFile, memory);
		inputMatrix->imputeMissingOn();	


	/*********************************************************************
	 *  Read Binary Input.  
	 */
	}else if (opt.isSet("-binaryI")) {
		int m_arg,n_arg;
		opt.get("-m")->getInt(m_arg);
		opt.get("-n")->getInt(n_arg);
		//Read the file
		std::string inputFile;
		opt.get("-binaryI")->getString(inputFile);
		fastpca_debug_print("Reading file: %s\n", inputFile.c_str());

		inputMatrix = new InputMatrixGeneralBinary(inputFile, memory, m_arg,n_arg);
		

	}else{
		fprintf(stderr,"No input matrix.  Please set -eigenI, -binaryI, eigenInCoreI, or -bedI.\n");
		return -1;
	}
	/*********************************************************************
	 *  Initialize the InputMatrix, and output for diagnositcs (-binA)  
	 */

	if (opt.isSet("-c")) {
		inputMatrix->centerRowsOn();	
		fastpca_debug_print("%s", "Centering rows: True\n");
	}
	int info = inputMatrix->init();
	if (opt.isSet("-cc")){
		inputMatrix->centerColumns();
		fastpca_debug_print("%s", "Centering Columns: True\n");
	}
	if (info < 0 ) {
		std::cerr << "An unknown error has occured in inputMatrix->init(): " << info << std::endl;
		return -1;
	}
	if (opt.isSet("-binA")) {
		std::string outputFile;
		opt.get("-binA")->getString(outputFile);
		std::cout <<"Outputting preprocessed input matrix to binary file "<< outputFile << std::endl;
		fastpca_debug_print("Load as\n Ag= loadFastPcaBinary('%s', %lld,%lld);\n", outputFile.c_str(), inputMatrix->m, inputMatrix->n);
		fastpca_write_input_matrix_binary_format(outputFile.c_str(), inputMatrix);
	}


	/*********************************************************************
	 *  Calculate the rank k SVD of A as defined as USV'=A.  Sizes of U,S,V
	 *  are (m,k), (k,k), and (k,n) respectively.
	 */
	for (int test = 0; test <0; test ++ ){
		printf("%d iteration of Test!\n", test);
		while (inputMatrix->hasNext()){
			printf("Begin: %d, End: %d\n", inputMatrix->blockStart, inputMatrix->blockEnd);
			fastpca_print_matrix("TestAfter", inputMatrix->blockSize,inputMatrix->n,inputMatrix->block);
		}
	}
	//fastpca_debug_print("%s", "Should end here...\n");
	//exit(-5);
	fastpca_debug_print("%s", "Begin PCA...\n");
	fastpca_debug_print("Running PCA with k=%lld, m=%e, n=%e, l=%lld, its=%d\n", k,(double)inputMatrix->m,(double)inputMatrix->n,l,its);
	high_resolution_clock::time_point t_pca = high_resolution_clock::now();
	double * U = (double *)fastpca_aligned_alloc( 64, (long long) inputMatrix->m*k*sizeof( double ) );
	double * S = (double *)fastpca_aligned_alloc(64, (long long) k*k*sizeof( double ));
	double * V = (double *)fastpca_aligned_alloc( 64, (long long) k*inputMatrix->n*sizeof( double ));

	info = fastpca_pca(k,l, inputMatrix, U,S,V,its);
	if (info <0 ) {
		std::cerr << "pca() returned an error code of "<< info << std::endl;
		return -1;
	}

	duration<double> time_span_pca = duration_cast<duration<double>>(high_resolution_clock::now() - t_pca);
	if (opt.isSet("-time")) {
		std::string timeFile;
		opt.get("-time")->getString(timeFile);
		std::ofstream outputFile( timeFile, std::ios::out | std::ios::app	);
		if (!outputFile.is_open()){
			return -1;
		}
		outputFile << inputMatrix->m << "," << inputMatrix->n << "," << k << "," << time_span_pca.count() <<std::endl;
	}
	fastpca_debug_print("The PCA function took %lf seconds.\n", time_span_pca.count());

	//fastpca_print_matrix("S", k,k,S);
//	fastpca_print_matrix("U", inputMatrix->m,k,U);
//fastpca_print_matrix("V", inputMatrix->n,k,V);
	/*********************************************************************
	 *  Print the results of the SVD for diagnostic purposes
	 */
	
	fastpca_debug_print("%s", "End PCA\n\n");
	if (opt.isSet("-eigenOutput")) {
		std::string outputFile;
		opt.get("-eigenOutput")->getString(outputFile);

		fastpca_debug_print ("Outputting to file %s\n ",outputFile.c_str());

		fastpca_write_eigenstrat_format(outputFile.c_str(), inputMatrix->m, inputMatrix->n, k,S,V);
	}
	if (opt.isSet("-binaryOutput")) {
		std::string outputFile;
		opt.get("-binaryOutput")->getString(outputFile);

		fastpca_debug_print ("Outputting to binary file %s\n ",outputFile.c_str());
		fastpca_save_bin(outputFile.c_str(),inputMatrix,U,S,V,k);
	}
	if (! opt.isSet("-eigenOutput") && ! opt.isSet("-binaryOutput")){
		fastpca_debug_print("%s",   "-eigenOutput and -binaryOutput not set, no file outputting");
	}


	/*********************************************************************
	 *  Calculate the L2 norm of A-USV' as a measure of accuracy of the 
	 *  decomposition
	 */
	if (opt.isSet("-d")) {
		fastpca_debug_print("%s", "Begin diffsnorm...\n");
		double snorm;
		int info = fastpca_diffsnorm(inputMatrix, U, S,V,k,dits, snorm);
		if (info <0 ) {
			std::cerr << "diffsorm() returned an error code of "<< info << std::endl;
			return -1;
		}
		fastpca_debug_print("The norm is %le\n", snorm);
		fastpca_debug_print("%s", "End diffsnorm\n\n");
	}
	return 1;
} 

