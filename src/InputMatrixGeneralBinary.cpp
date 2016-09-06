#include "InputMatrix.h"
#include "InputMatrixGeneralBinary.h"
#include <stdlib.h>
#include <chrono>
#include <math.h>
/*********************************************************************
 *  Implementation file for InputMatrixGeneralBinary.  
 *  Please see the header file for extensive class documentation
 */ 
InputMatrixGeneralBinary::InputMatrixGeneralBinary( std::string p, long long maxMemoryBytes, int m, int n)
	: InputMatrix(p,maxMemoryBytes) {
		this->m = m;
		this->n = n;

	}

int InputMatrixGeneralBinary::init() {
	inputFile = new std::ifstream(this->path.c_str(), std::ios::in | std::ios::binary);
	if (NULL == inputFile || !inputFile->is_open()) {
		return -11;
	}
	return InputMatrix::init();
	return 1;
}

int InputMatrixGeneralBinary::loadNextBlock() {
	if (inputFile->eof() || inputFile->peek() == EOF) {
		inputFile->clear();
		inputFile->seekg(0, std::ios::beg);
	}
	inputFile->read((char *) this->block, sizeof(double)*this->blockSize*this->n);
	InputMatrix::loadNextBlock();
	return 1;
}
