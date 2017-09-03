#include <iostream>
#include <fstream>
#include "ezOptionParser.hpp"
using namespace std::chrono;
   int main(int argc, const char * argv[]) {
		using namespace ez;
		ezOptionParser opt;
	opt.add(
			"", // Default.
			1, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Input file, in CSV format", // Help description.
			"-csvI",     // Flag token. 
			"-csvInput",  // Flag token.
			"--csvInput" // Flag token.
	       );
	opt.add(
			"", // Default.
			1, // Required?
			1, // Number of args expected.
			0, // Delimiter if expecting multiple args.
			"Output file, in binary format", // Help description.
			"-binO",     // Flag token. 
			"-binOutput",  // Flag token.
			"--binOutput" // Flag token.
	       );


	opt.parse(argc, argv);
	std::string inputFile;
	std::string outputFile;
	if (opt.isSet("-csvI")) {
		//Read the file
		opt.get("-csvI")->getString(inputFile);
	}
	if (opt.isSet("-binO")) {
		//Read the file
		opt.get("-binO")->getString(outputFile);
	}
	//First traversal to determine m and n.  eachLineN is a variable updated after every line
	//to contain the number of columns.  If the columns on a given line do not equal those
	// on the line before it, return an error.
	//printf("Converting CSV to binary");
	std::ifstream inCSV;
	   inCSV.open(inputFile.c_str(), std::ios::in);
	int eachLineN=0;
	int m =0;
	int n=0;
	char ch;
	while(inCSV >> std::noskipws >> ch ) {
	//	std::cout << ch << std::endl;
		if (ch == '\n'){
			if (m==0) {
				n = eachLineN+1;
			}
			if (eachLineN+1 != n) {
				printf("ERROR: The number of columns on row %d does not match  \
						the number of columns on the previous row.\n", m);
				return -2;
			}
			n = eachLineN+1;
			eachLineN = 0;
			m++;
		}else if (ch == ','){
			eachLineN++;
		}
	}
	inCSV.clear();
	inCSV.seekg(0, std::ios::beg);
	//printf("Converting CSV file with %d rows and %d columns into binary format.\n", m,n);

	double * block = (double*) malloc(m*n*sizeof(double));
	int i = 0;
	int j = 0;
	double num;
	while(inCSV>> std::noskipws >> num ) {
		block[(long long)i*n+j] = num;
		//printf("%d,%d = %.12lf\n", i, j, block[(long long)i*n+j]);
		//If reach end of file, start from the beginning
		inCSV >> std::noskipws >> ch;
		if (ch == '\n'){
			i++;
			j=0;
		} else if (ch == ','){
			j++;
		}
	}

	std::ofstream outBin( outputFile.c_str(), std::ios::out | std::ios::binary | std::ios::trunc	);
	//outputFile.open(filename);
	if (!outBin.is_open()){
		return -1;
	}
	//outputFile.precision(4);
	outBin.write((char*) block,sizeof (double)*m*n );
	outBin.close();
	/*
	   double num;
		inCSV >> std::noskipws >> ch;
		while(inCSV >> std::noskipws >> num ) {

		}
		*/

} 

