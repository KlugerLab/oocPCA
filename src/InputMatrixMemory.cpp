#include "InputMatrix.h"
#include "InputMatrixMemory.h"
#include <stdlib.h>
#include <chrono>
#include <math.h>
/*
void printBits(unsigned int num)
{
   for(int bit=0;bit<(8); bit++)
      {
            printf("%i ", num & 0x01);
	          num = num >> 1;
      }
      printf("\n");
}
*/
/*********************************************************************
 *  Implementation file for InputMatrixMemory.  
 *  Please see the header file for extensive class documentation
 */ 
InputMatrixMemory::InputMatrixMemory( double * A, long long maxMemoryBytes, int m, int n, std::string p)
	: InputMatrix(p,maxMemoryBytes){
	this->block = A;
	this->m = m;
	this->n = n;
	//fastpca_debug_print ("m %lld n %lld", this->m,this->n);

}

int InputMatrixMemory::init() {
	this->calculatedBlockSize   = this->m;
	this->calculatedRemainderBlockSize	 = 0;
	this->blockSize   = this->m;
	this->blockStart = 0;
	this->blockEnd = this->blockStart + this->blockSize;
	this->traversalStart = 0;
	this->firstIt = true;
	this->remainderIt = false;
	fastpca_debug_print ("%s", "Init complete");
	return 1;
}
int InputMatrixMemory::hasNext() {

	if (this->firstIt){
		this->firstIt = false;
		return true;
	}else{
		this->firstIt = true;
		return false;
	}

}

int InputMatrixMemory::loadNextBlock() {
	return 1;
}
