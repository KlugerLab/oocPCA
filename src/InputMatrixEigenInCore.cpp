#include "InputMatrix.h"
#include "InputMatrixEigenInCore.h"
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
 *  Implementation file for InputMatrixEigenInCore.  
 *  Please see the header file for extensive class documentation
 */ 
InputMatrixEigenInCore::InputMatrixEigenInCore( std::string p, long long maxMemoryBytes)
	: InputMatrix(p,maxMemoryBytes){

}

int InputMatrixEigenInCore::init() {
	std::ifstream eigenstratFile;
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

	long long int bytesPerSnp = ceil(n / 4.0);

	this->matrixSize = bytesPerSnp * this->m;

	fastpca_debug_print("Bytes per snp: %le, Allocating a matrix of size %le\n",(double) bytesPerSnp, (double) matrixSize);
	//TODO: make sure this can exceed integer (long long)
	this->source = (char *)fastpca_aligned_alloc(64,  matrixSize);
	if (NULL == this->source) {
		std::cerr << "ERROR: source was not successfully allocated.  Is "
			<< matrixSize << 
			" bytes larger than the virtual memory?" << std::endl ;
		return -1;
	}
	fastpca_debug_print("%s", "Matrix successfully allocated\n");


	long long int byteCount = 0;
	for (long long int i=0; i<this->m; i++ ) {
		long long int j = 0;
		//this->source[i*bytesPerSnp + j] = 0;
		this->source[byteCount] = 0;
		//int k =0;
		//for (int j = 0;  j<this->n; j++) {
		//for (int j = 0;  j<this->n; j++) {
		while(j<this->n) {
			for (int offset=0; offset< 8; offset +=2 ){
				if (j>= this->n){
					break;
				}
				char ch;
				eigenstratFile >> std::skipws >> ch;

				int toAssign;
				if (ch == '0'){
					toAssign = 0;
				}else if (ch=='1'){
					toAssign = 2;
				}else if (ch=='2') {
					toAssign = 3;
				}else if (ch=='9'){
					toAssign = 1;
				}else{
					return -21;
				}

				toAssign = toAssign << offset;
					
				this->source[byteCount] = this->source[byteCount] | toAssign;
				//printf("Assiging %c-> %d,%d,%d, at byte: %d\n", ch, i,j,toAssign,byteCount);
				j++;
			}
			//printBits(this->source[byteCount]);
			byteCount++;
		}
	}

	//Since we are allocating all this memory, we need to subtract it from the total
	//memory allowed
	this->maxMemoryBytes -= matrixSize;
	byteCount = 0;
	return InputMatrix::init();

}

int InputMatrixEigenInCore::loadNextBlock() {

	printf("loading next block at %d with %d blockSize",byteCount, this->blockSize );
	long long int i=0;
	long long int j=0;
	char genotypeByte;
	long long int index = 0;
	while ( i < this->blockSize ) {
		

		//if (fmod((double)index, (double) 1E8)==0) {
			//printf("%lld...\n", index);

		//}	
		if ( byteCount == this->matrixSize  ) {
			byteCount = 0;
		}
		genotypeByte = this->source [byteCount];
		for (int jump=0; jump<8; jump += 2 ){
			index = (long long)i*(long long)this->n+j;
			//printf("reading %d,%d, at byte: %d\n",  i,j,byteCount);
			char ones = ((genotypeByte>>jump) & 1);
			char twos = ((genotypeByte>>(jump+1)) & 1);
			//Homozygote genotype, encoded as 0
			if ( ones == 0 && twos == 0 ) {
				(this->block)[index] = 0;
			
			//Other homozygote, encoded as 2
			}else if(ones == 1 && twos == 1 ){
				(this->block)[index] = 2;

			//Heterozygote
			}else if ( ones == 0 && twos == 1 ){
				(this->block)[index] = 1;
			}else if (ones== 1 && twos == 0) {
				(this->block)[index] = 9;
			}
			j++;
			if (j>=n){
				//printf("i:%lld\n", i);
				i++;
				j=0;
				break;
			}
			//printf("%lld:%lf\n", index, this->block[index]);
		}

	byteCount++;
	}
	printf("Successfully loaded a block\n");
	InputMatrix::loadNextBlock();
	printf("Finished loaded a block!!\n");
	return 1;
}
