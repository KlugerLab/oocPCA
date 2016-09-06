#include "InputMatrix.h"
#include "InputMatrixBedInCore.h"
#include <stdlib.h>
#include <chrono>
#include <math.h>
/*********************************************************************
 *  Implementation file for InputMatrixBedInCore.  
 *  Please see the header file for extensive class documentation
 */ 
InputMatrixBedInCore::InputMatrixBedInCore( std::string p, long long maxMemoryBytes)
	: InputMatrix(p,maxMemoryBytes){

}

int InputMatrixBedInCore::init() {

	std::string baseFileName(this->path.c_str());
	std::string bimExtension(".bim");
	std::string bimFileName = baseFileName + bimExtension;
	std::ifstream bimFile;
	bimFile.open(bimFileName, std::ios::in);
	if (!bimFile.is_open()){
		return -101;

	}


	std::string famExtension(".fam");
	std::string famFileName = baseFileName + famExtension;
	std::ifstream famFile;
	famFile.open(famFileName, std::ios::in);
	if (!famFile.is_open()){
		return -103;

	}

	std::string bedExtension(".bed");
	std::string bedFileName = baseFileName + bedExtension;
	std::ifstream bedFile;
	bedFile.open(bedFileName, std::ios::in | std::ios::binary);
	if (!bedFile.is_open()){
		return -102;

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

	fastpca_debug_print("Read magic numbers:%d_%d_%d\n", magicNumberBuffer[0], magicNumberBuffer[1], magicNumberBuffer[2]);

	//Read data into the memory
	std::streampos fsize = bedFile.tellg();
	bedFile.seekg(0, bedFile.end);
	fsize = bedFile.tellg();
	this->matrixSize = (long long)fsize - 3;

	fastpca_debug_print("Allocating a matrix of size %le\n",(double) matrixSize);
	//TODO: make sure this can exceed integer (long long)
	this->source = (char *)fastpca_aligned_alloc(64,  matrixSize);
	if (NULL == this->source) {
		std::cerr << "ERROR: source was not successfully allocated.  Is "
			<< matrixSize << 
			" bytes larger than the virtual memory?" << std::endl ;
		return -1;
	}
	bedFile.seekg(3);
	fastpca_debug_print("%s","Matrix successfully allocated\n");

	//Since we are allocating all this memory, we need to subtract it from the total
	//memory allowed
	this->maxMemoryBytes -= matrixSize;
	bedFile.read(this->source, matrixSize);
	bedFile.close();
	byteCount = 0;
	return InputMatrix::init();

}

int InputMatrixBedInCore::loadNextBlock() {

	//printf("loading next block at %d with %d blockSize",byteCount, this->blockSize );
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
			}
			//printf("%lld:%lf\n", index, this->block[index]);
		}

	byteCount++;
	}
	//printf("Successfully loaded a block\n");
	return InputMatrix::loadNextBlock();

}
