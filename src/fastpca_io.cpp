#include "fastpca_io.hpp"
int fastpca_print_matrix( const char* desc, long long m, long long n, double* A ) {
	   //printf("%s\n", desc);
	   fastpca_debug_print("m:%lld, n:%lld\n", m,n);
	   long long i, j;
	   long long n_long = n;
	   printf( "\n \n%s = [", desc );
	   long long int totalCount = 0;
	   for( i =0; i < m; i++ ) {
		   for( j = 0; j < n; j++ ) {
		   
//		   printf( " %6.12f", A[i*n+j] );
		   printf( " %6.2f", A[i*n+j] );
		   totalCount++;
		if (isnan(A[i*n_long+j] )) {
			fastpca_debug_print ("A DOES HAVE A NAN AT %lld %lld", i,j);
			exit(-1);
		}

		   }
		   printf( ";\n" );
	   }
	   printf("]\n");
	return 0;
}
int fastpca_read_bed_format(const char * filename, long long &m, long long &n, double *& A) {
	/*********************************************************************
	 *  First determine the number of rows and columns from the .bim and .fam
	 *  files
	 */
	std::string baseFileName(filename);
	std::string bimExtension(".bim");
	std::string bimFileName = baseFileName + bimExtension;
	std::ifstream bimFile;
	bimFile.open(bimFileName, std::ios::in);
	if (!bimFile.is_open()){
		return -1;

	}

	std::string bedExtension(".bed");
	std::string bedFileName = baseFileName + bedExtension;
	std::ifstream bedFile;
	bedFile.open(bedFileName, std::ios::in | std::ios::binary);
	if (!bedFile.is_open()){
		return -2;

	}


	std::string famExtension(".fam");
	std::string famFileName = baseFileName + famExtension;
	std::ifstream famFile;
	famFile.open(famFileName, std::ios::in);
	if (!famFile.is_open()){
		return -3;

	}


	//Count the number of rows
	m =0;
	char ch;
	while(bimFile >> std::noskipws >> ch ) {
		if (ch == '\n'){
			m++;
		}
	}
	fastpca_debug_print("Identified m=%lld snps\n", m);

	//Count the number of rows
	n=0;
	while(famFile >> std::noskipws >> ch ) {
		if (ch == '\n'){
			n++;
		}
	}
	fastpca_debug_print("Identified n=%lld patients\n", n);

	fastpca_debug_print("Allocating a matrix of size %le\n",(double) m*n*sizeof(double));
	A = (double *)fastpca_aligned_alloc(64, (long long) m*n*sizeof( double ));
	if (NULL == A) {
		std::cerr << "ERROR: A was not successfully allocated.  Is "
			<< (long long) m*n*sizeof(double) << 
			" bytes larger than the virtual memory?" << std::endl ;
		return -1;
	}
	fastpca_debug_print("%s", "Matrix successfully allocated\n");

	/*********************************************************************
	 *  Read Binary Input from the BED file, first checking its format
	 * 
	 */
	char magicNumberBuffer[3];

	bedFile.read(magicNumberBuffer,3);


	if ( magicNumberBuffer[0] != 108 ||  magicNumberBuffer[1] != 27 ) {
		return -20;
	}

	if (magicNumberBuffer[2] != 1) {
		return -21;
	}

	fastpca_debug_print("Read:%d_%d_%d\n", magicNumberBuffer[0], magicNumberBuffer[1], magicNumberBuffer[2]);

	//Read in the data to the passed array
	//char * genotypeBuffer = malloc(sizeof(char) * 
	bedFile.seekg(3);
	long long int i=0;
	long long int j=0;
	char genotypeByte;
	while (bedFile.get(genotypeByte))
	{
		for (int jump=0; jump<8; jump += 2 ){
			if (j>=n){
				//printf("i:%lld\n", i);
				i++;
				j=0;
				break;
			}else{
				char ones = ((genotypeByte>>jump) & 1);
				char twos = ((genotypeByte>>(jump+1)) & 1);
				//Homozygote genotype, encoded as 0
				if ( ones == 0 && twos == 0 ) {
					(A)[(long long)i*n+j] = 0;
				
				//Other homozygote, encoded as 2
				}else if(ones == 1 && twos == 1 ){
					(A)[(long long)i*n+j] = 2;

				//Heterozygote
				}else if ( ones == 0 && twos == 1 ){
					(A)[(long long)i*n+j] = 1;
				}else if (ones== 1 && twos == 0) {
					(A)[(long long)i*n+j] = 9;
				}
				j++;
			}

		}

	}


	return 1;
}
int fastpca_impute_missing_average(double *& A, long long m, long long n, double missingValue) {
	for (int i=0; i<m; i++ ){
		//Find the average (not including the missing values)
		int sum = 0;
		int countOfNonMissing = 0;
		for (int j=0; j<n; j++ ){
			double val = A[i*n+j];
			if (val != missingValue) {
				sum += val;
				countOfNonMissing++;
			}
		}
		if (countOfNonMissing ==0 ){
			return -400;
		}
		double average = sum/(double)countOfNonMissing;

		for (int j=0; j<n; j++ ){
			double val = A[i*n+j];
			if (val == missingValue) {
				A[i*n+j] = average;
			}
		}

	}
	 return 1;
}
int fastpca_write_binary_format(const char * filename, long long m, long long n,double * A) {
	std::ofstream outputFile( filename, std::ios::out | std::ios::binary | std::ios::trunc	);
	//outputFile.open(filename);
	if (!outputFile.is_open()){
		return -1;
	}
	//outputFile.precision(4);
	outputFile.write((char*) A,sizeof (double)*m*n );
	outputFile.close();


	return 1;
}
int fastpca_write_input_matrix_binary_format(const char * filename, InputMatrix * A) {

	std::ofstream outputFile( filename, std::ios::out | std::ios::binary | std::ios::trunc	);
	//outputFile.open(filename);
	if (!outputFile.is_open()){
		return -1;
	}
	while (A->hasNext()){
		//fastpca_print_matrix("Test", A->m, A->n, A->block);
		outputFile.write((char*) A->block,sizeof (double) * A->blockSize*A->n );
	}
	outputFile.close();


	 return 1;
}
int fastpca_write_eigenstrat_format(const char * filename, long long m, long long n,long long k,double *S, double *V) {
	std::ofstream outputFile;	
	outputFile.open(filename);
	if (!outputFile.is_open()){
		return -1;
	}
	outputFile.precision(4);
	outputFile << k << "\n";
	for(int i=0; i <k; i++ ){
		outputFile << S[i*k+i] << std::endl;
	}
	outputFile << "  ";
	for (int i=0; i<n; i++ ) {
		for (int j=0; j<k; j++){
			outputFile << std::fixed << V[i*k+j] << "  ";
			//printf("%6.5f  ",V[i*n+j]);

		}
		//printf("  \n");
		outputFile << std::endl << "  ";
	}


	return 1;
}
int fastpca_save_bin( const char * basename, InputMatrix* A, double* U, double* S , double* V,long long k ) {
		std::string baseFileName(basename);
		std::string U_suffix(".U");
		std::string S_suffix(".S");
		std::string V_suffix(".V");
		std::string UFileName = baseFileName +  U_suffix;  
		std::string SFileName = baseFileName +  S_suffix;  
		std::string VFileName = baseFileName +  V_suffix;  
		fastpca_write_binary_format(UFileName.c_str(), A->m,k, U);
		fastpca_write_binary_format(SFileName.c_str(), k,k, S);
		fastpca_write_binary_format(VFileName.c_str(), k,A->n, V);
		return 1;
}

