#include "InputMatrix.h"
#include "InputMatrixCSV.h"
#include <stdlib.h>
#include <chrono>
#include <math.h>
/*********************************************************************
 *  Implementation file for InputMatrixCSV.  
 *  Please see the header file for extensive class documentation
 */ 
InputMatrixCSV::InputMatrixCSV( std::string p, long long maxMemoryBytes )
	: InputMatrix(p,maxMemoryBytes){

}

int InputMatrixCSV::init() {
	eigenstratFile.open(this->path.c_str(), std::ios::in);
	if (!eigenstratFile.is_open()){
	fastpca_debug_print("The file %s could not open\n", this->path.c_str());
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
		if (ch == '\n'){
			if (m==0) {
				this->n = eachLineN+1;
			}
			if (eachLineN+1 != this->n) {
				fastpca_debug_print("ERROR: The number of columns on row %lld does not match  \
						the number of columns on the previous row.\n", m,eachLineN,ch);
				return -2;
				
			}
			this->n = eachLineN+1;
			eachLineN = 0;
			this->m++;
		}else if (ch == ','){
			eachLineN++;
		}
	}
	eigenstratFile.clear();
	eigenstratFile.seekg(0, std::ios::beg);
	fastpca_debug_print("The number of rows: %lld and columns: %lld \n", m,n);

	return InputMatrix::init();

}

int InputMatrixCSV::loadNextBlock() {
	long long i= 0;
	long long j = 0;
	double num;
	char ch;

	//fastpca_debug_print("%s\n", "Load Next Block");
	int loc = eigenstratFile.tellg();
	//fastpca_debug_print("Current location %d\n", loc);
	//load this->blockSize rows
	while(i < ( this->blockSize) && eigenstratFile >> std::noskipws >> num ) {


		(this->block)[(long long)i*n+j] = num;
		//fastpca_debug_print("%d,%d = %lf\n", i*n, j, (this->block)[(long long)i*n+j]);
		//If reach end of file, start from the beginning
		eigenstratFile >> std::noskipws >> ch;
		if (ch == '\n'){
			i++;
			j=0;
		} else if (ch == ','){
			j++;
		}
		int loc = eigenstratFile.tellg();
		//fastpca_debug_print("In loop%d,%d ,%d\n",i,j, loc);
		if (eigenstratFile.eof() || eigenstratFile.peek() == EOF) {
			eigenstratFile.clear();
			eigenstratFile.seekg(0, std::ios::beg);
			i++;
			j=0;
			continue;
		}
		//Check it's a number.
		/*
			double temp;
			temp  = ch - '0';
			//printf("i,j: %i,%i=%lf\n", i,j, temp);

			//Eigenstrat files' genotypes can only take on values 0,1,2, or 9 (for missing data)
			if (temp != 0 && temp != 1 && temp != 2 && temp != 9) {
				fastpca_debug_print("ERROR: The %lld,%lld element is invalid: %c\n", i,j,ch);
				return -1;
			}
			*/
		}
	//printf("----------Loaded a block--------------\n");
	return InputMatrix::loadNextBlock();

}
