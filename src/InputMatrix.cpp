#include <climits>
#include "InputMatrix.h"
#include "fastpca_io.hpp"
/*********************************************************************
 *  Implementation file for InputMatrix.  
 *  Please see the header file for extensive class documentation
 */ 

InputMatrix::InputMatrix(std::string p, long long maxMemoryBytes)
	: path(p), maxMemoryBytes(maxMemoryBytes)
{
	centerRows = false;
	imputeMissing = false;
}
int InputMatrix::loadNextBlock () {
	//fastpca_debug_print("%s", "Preprocess new block\n");
	preprocessBlock();
	return 1;
}
int InputMatrix::init () {

	//Calculate the appropriate block size
	long long m_long, n_long, double_size;
	m_long =m;
	n_long = n;
	double_size = sizeof(double);
	long long totalMem = (long long)m*(long long)n*(long long)double_size;
	//long long totalMem = this->m*this->n*(long long)sizeof( double );
	int blockNum = ceil(totalMem/(double)this->maxMemoryBytes);
	double blockNum_test = totalMem/(double)this->maxMemoryBytes;
	fastpca_debug_print("totalMem: %lld, this->maxMemoryBytes: %lld, blockNum: %d, blockNum_test: %e\n",  totalMem, this->maxMemoryBytes, blockNum, blockNum_test);
	if (this->maxMemoryBytes < 1 || blockNum > this->m) {
		std::cerr << "ERROR: must at least allocate enough memory for a single row (and in-core storage, if applicable): "
			<< (long long) this->n*sizeof(double) <<  std::endl;
		return -2;

	}

	this->calculatedBlockSize   = floor(m/(double)blockNum);
	this->calculatedRemainderBlockSize	 = m % this->calculatedBlockSize;
	this->blockSize   = this->calculatedBlockSize;

	fastpca_debug_print("m:%lld, n:%lld, blockSize:%d,  totalMem: %lld, blockNum:%d, calculated block size %lld,  remainder: %d\n", m,n, blockSize, totalMem, blockNum, this->calculatedBlockSize, this->calculatedRemainderBlockSize);
	long long blockSizeBytes = (long long) this->calculatedBlockSize*(long long) this->n* (long long) sizeof( double );
	//this->block    = (double *)mkl_malloc( (long long) this->calculatedBlockSize*this->n*sizeof( double ), 64 );
	//fastpca_debug_print("Would you like to proceed in allocating %le bytes? \n<enter> to continue, ctrl+c to cancel.\n", (double)blockSizeBytes);
	this->block    = (double *)fastpca_aligned_alloc(64, blockSizeBytes);
	if (NULL == this->block) {
		std::cerr << "ERROR: Block was not successfully allocated.  Is "
			<< (long long) this->calculatedBlockSize*this->n*sizeof(double) << 
			" bytes larger than the virtual memory?" << std::endl ;
		return -1;
	}
	if (this->maxMemoryBytes < this->n*sizeof(double)) {
		std::cerr << "ERROR: must at least allocate enough memory for a single row: "
			<< (long long) this->n*sizeof(double) <<  std::endl;
		return -2;

	}
	this->blockStart = 0;
	this->blockEnd = this->blockStart + this->blockSize;
	this->traversalStart = 0;
	this->firstIt = true;
	this->remainderIt = false;


	//Setting centers to zero
	this->rowMeans = (double *) malloc(this->m*sizeof(double));
	this->rowMeansFlags = (bool *)malloc(this->m*sizeof(bool));
	for (int i=0; i<this->m; i++) {
		this->rowMeansFlags[i] = false;
	}
	this->colMeans = (double *) malloc(this->n*sizeof(double));
	for (int i=0; i<this->n; i++) {
		this->colMeans[i] = 0;
	}
	this->centerColumnsFlag = false;

	fastpca_debug_print("%s", "The block was allocated, proceeding to load the first block...\n");
	loadNextBlock();
	fastpca_debug_print("%s", "Loaded the first block successfully\n");
	
	return 0;

}
//TODO: It would be nice to refactor this into a child class. Because you do NOT want to do column centering if you have to impute missing data on a GWAS dataset.
bool InputMatrix::centerColumns() {
	while (this->hasNext()){
		for (int j=0; j<this->n; j++){
			for (int i=0; i<this->blockSize; i++){
				this->colMeans[j] += this->block[i*this->n+j];	
				//fastpca_debug_print("%d,%d = %lf, then added to %lf\n",i, j,this->block[i*this->n+j], this->colMeans[j]);
				//this->colMeans[j] = 1;	
			}
			
		}
	}
	for (int j=0; j<this->n; j++){
		this->colMeans[j] /= this->m;	
		//fastpca_debug_print("Mean for %d, is %lf\n", j, this->colMeans[j]);
	}

	this->centerColumnsFlag = true;
	this->preprocessBlock();
	return true;

}

int InputMatrix::preprocessBlock() {
	
	if (!this->imputeMissing && !this->centerRows && !this->centerColumnsFlag){
		return 1;
	}
	//fastpca_debug_print("%s", "Begin preprocessing\n");

	//double * rowsums = (double *)mkl_malloc( this->blockSize*sizeof( double ), 64 );
	int countOfNonMissing = 0;
	double missingValue = 9;
	for (long long i=0; i<this->blockSize; i++ ){
		if (this->rowMeansFlags[this->blockStart + i] == false) {
			double sum =0;
			int countOfNonMissing =0;
			for ( long long j=0; j<n; j++ ){
				double val = this->block[i*n+j];
				/*
				if (this->blockStart + i < 697 + 5 && j < 100) {
					fastpca_debug_print("%lf,",val);
				}
				*/
				if (!imputeMissing || val != missingValue) {
					sum += val;
					countOfNonMissing++;
				}
			}
			if (countOfNonMissing == 0 ) {
				this->rowMeans[this->blockStart + i] =0;
				fastpca_debug_print("Line %d is all missing! Will fill it with zeros.\n", i);
			}else {
				this->rowMeans[this->blockStart + i] = sum/(double)countOfNonMissing;
			}
			this->rowMeansFlags[this->blockStart + i] = true;
		}
		double average = this->rowMeans[this->blockStart + i];	

		for (long long j=0; j<n; j++ ){
			double val = this->block[i*n+j];
			if (imputeMissing && val == missingValue) {
				this->block[i*n+j] = average;
			}
			if (centerRows) {
				this->block[i*n+j] -= average;
			}
			if (centerColumnsFlag) {
				this->block[i*n+j] -= this->colMeans[j];
			//fastpca_debug_print("Changing %lf by %lf \n", this->block[i*n+j], this->colMeans[j]);
			}
		}

	}
	//fastpca_debug_print("%s", "End preprocessing\n");

	 return 1;
}

bool InputMatrix::hasNext() {

	//fastpca_debug_print("blockStart %d,blockStop: %d,  blockSize %d, this->m %d\n", this->blockStart, this->blockEnd, this->blockSize, this->m);
	//Case 1: This is the first iteration of a traversal, use the already
	//loaded chunk
	if (firstIt) {
		//fastpca_debug_print("%s", "case 1");
		firstIt = false;
		return true;
	}

	//Case 2: This is the end of the traversal.  Prepare for the next
	//traversal, and end this one by returning false
	if((this->blockEnd % this->m) == this->traversalStart) {
		//fastpca_debug_print("%s", "case 2");
		firstIt = true;
		this->traversalStart = this->blockStart;
		return false;
	}

	//Case 3: Continue from the beginning.  Since we are reusing blocks,
	//blockEnd == this->m does not mean the end of the matrix, since a
	//traversal does not always start from zero. 
	if (this->blockEnd == this->m) {
		//fastpca_debug_print("%s", "case 3");
		this->blockSize = this->calculatedBlockSize;
		this->blockStart = 0;

	//Case 4: Use the "Remainder Block"
	}else if ( (this->m - this->blockEnd  ) == this->calculatedRemainderBlockSize ) {
		//fastpca_debug_print("%s", "case 4");
		this->blockSize = this->calculatedRemainderBlockSize;
		this->blockStart = this->blockEnd;

	//Case 5: Every other iteration, just add the normal block size
	}else {
		//fastpca_debug_print("%s", "case 5");
		this->blockSize = this->calculatedBlockSize;
		this->blockStart = this->blockEnd;
	}
	this->blockEnd = this->blockStart + this->blockSize;
	loadNextBlock();
	//fastpca_print_matrix("HasNextEnd", this->blockSize,this->n,this->block);
	return true;
}

void InputMatrix::imputeMissingOn() {
	this->imputeMissing = true;
}
void InputMatrix::centerRowsOn() {
	this->centerRows = true;
}
