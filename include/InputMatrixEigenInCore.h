#ifndef INPUTMATRIXEIGENINCORE_H 
#define INPUTMATRIXEIGENINCORE_H
#include "InputMatrix.h"
/*********************************************************************
 *  Input Matrix Class for Eigenstrat .eigenstratgeno files
 * 
 *  
 **  Use an InputMatrixEigensrat class to represent a matrix stored in the
 **  eigenstrat .eigenstratgeno format.
 **  
 **  Unlike the general case, SNPs can be stored very effectively as 2 bits/snp.
 **  However, they cannot be manipulated in this format, they must be converted to double/single precision.
 **  As such, we read the entire matrix from the hard drive and store it in the memory in the 2 bit/snp
 **  format, and then we pull blocks from this source.
 **
*
 */
class InputMatrixEigenInCore : public InputMatrix {
	public:
		//Constructor.  Initializes the class variables
		InputMatrixEigenInCore( std::string p, long long memory);


		/*
		**Initializes the object with the first block.  Only run once
		**per object.  
		*
		* RETURN
		 ** int   ier --   error return code
		 *   	  ier=0 The matrix initialized without error
		 *    	  ier=-1 to -10 check InputMatrix::init() error codes
		 *       ier=-101 means could not open the .eigenstrageno file
		 *       ier=-21 means that an invalid character was encountered
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
