#ifndef INPUTMATRIXBEDINCORE_H 
#define INPUTMATRIXBEDINCORE_H
#include "InputMatrix.h"
/*********************************************************************
 *  Input Matrix Class for Eigenstrat .bed files
 * 
 *  
 **  Use an InputMatrixEigensrat class to represent a matrix stored in the
 **  eigenstrat .bed format.
 **  
 **  Unlike the general case, SNPs can be stored very effectively as 2 bits/snp.
 **  However, they cannot be manipulated in this format, they must be converted to double/single precision.
 **  As such, we read the entire matrix from the hard drive and store it in the memory in the 2 bit/snp
 **  format, and then we pull blocks from this source.
 **
 *  USAGE 
** 	//Instantiate the matrix
**	inputMatrix = new InputMatrixEigenstrat(inputFile, memory);
**	//Initialize
**	int info = inputMatrix->init();
**	if (info < 0 ) {
**	//Error handling
**	}
**	//Loop through blocks
**	while (inputMatrix->hasNext()){
**		//Perform operation on each block as inputMatrix->block
**	}
*
 */
class InputMatrixBedInCore : public InputMatrix {
	public:
		//Constructor.  Initializes the class variables
		InputMatrixBedInCore( std::string p, long long memory);


		/*
		**Initializes the object with the first block.  Only run once
		**per object.  
		*
		* RETURN
		 ** int   ier --   error return code
		 *   	  ier=0 The matrix initialized without error
		 *    	  ier=-1 to -10 check InputMatrix::init() error codes
		 *       ier=-101 means could not open the .bim file
		 *       ier=-102 means could not open the .bed file
		 *       ier=-103 means could not open the .fam file
		 *       ier=-20 means that the binary file was not in proper BED format
		 *       ier=-21 means that the BED file was not in row-major format
		 *       ier=-30 means that memory was not allocated successfully
		 *       ier=-400 means that one of the values was missing

		 */
		int  init();

	protected:
		/*
		** Load the next block. 
		*
		* RETURN
		** int ier -- error return code
		** 	ier=0 The matrix initialized without error
		** 	ier=-1 Invalid element encountered
		*/
		int loadNextBlock();

	private:
		char * source;
		long long matrixSize;
		long long byteCount;
};
#endif
