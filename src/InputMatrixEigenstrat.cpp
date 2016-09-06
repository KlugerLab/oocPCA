#include "InputMatrix.h"
#include "InputMatrixEigenstrat.h"
#include <stdlib.h>
#include <chrono>
#include <math.h>
#include "fastpca_io.hpp"
/*********************************************************************
 *  Implementation file for InputMatrixEigenstrat.  
 *  Please see the header file for extensive class documentation
 */ 
InputMatrixEigenstrat::InputMatrixEigenstrat( std::string p, long long maxMemoryBytes )
	: InputMatrix(p,maxMemoryBytes){

}

int InputMatrixEigenstrat::init() {
	eigenstratFile.open(this->path.c_str(), std::ios::in);
	if (!eigenstratFile.is_open()){
		return -11;

	}

	//First traversal to determine m and n.  eachLineN is a variable updated after every line
	//to contain the number of columns.  If the columns on a given line do not equal those
	// on the line before it, return an error.
	int eachLineN=0;
	this->m =0;
	this->n=0;
	char ch;
	while(eigenstratFile >> std::noskipws >> ch ) {
		//printf("%d,%d=%c\n",m,eachLineN,ch);
		if (ch == '\n'){
			if (m==0) {
				this->n = eachLineN;
			}
			if (eachLineN != this->n) {
				fastpca_debug_print("ERROR: The number of columns on row %lld does not match  \
						the number of columns on the previous row.\n", m,eachLineN,ch);
				return -2;
				
			}
			this->n = eachLineN;
			eachLineN = 0;
			this->m++;
		}else{
			eachLineN++;
		}
	}
	eigenstratFile.clear();
	eigenstratFile.seekg(0, std::ios::beg);
	fastpca_debug_print("The number of rows: %lld and columns: %lld \n", m,n);

	return InputMatrix::init();

}

int InputMatrixEigenstrat::loadNextBlock() {
	long long i= 0;
	long long j = 0;
	char ch;

	//load this->blockSize rows
	while(i < ( this->blockSize) && eigenstratFile >> std::noskipws >> ch ) {

		//If reach end of file, start from the beginning
		if (eigenstratFile.eof() || eigenstratFile.peek() == EOF) {
			eigenstratFile.clear();
			eigenstratFile.seekg(0, std::ios::beg);
		}

		if (ch == '\n'){
			i++;
			j=0;
		}else {
			double temp;
			temp  = ch - '0';

			//Eigenstrat files' genotypes can only take on values 0,1,2, or 9 (for missing data)
			if (temp != 0 && temp != 1 && temp != 2 && temp != 9) {
				fastpca_debug_print("ERROR: The %lld,%lld element is invalid: %c\n", i,j,ch);
				return -1;
			}
			(this->block)[(long long)i*n+j] = temp;
//			printf("i,j: %i,%i=%lf (%lf)\n", i,j, temp,(this->block)[(long long)i*n+j] );
			j++;
		}
	}
	//printf("----------Loaded a block--------------\n");
	return InputMatrix::loadNextBlock();

}
