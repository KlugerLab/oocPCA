#ifndef INPUTMATRIXGENERALBINARY_H 
#define INPUTMATRIXGENERALBINARY_H
#include "InputMatrix.h"
/*********************************************************************
 *  Input Matrix Class for GeneralBinary .geno files
 * 
 *  
 **  Use an InputMatrixEigensrat class to represent a matrix stored in the
 **  eigenstrat .geno format.
 **
 *  USAGE 
** 	//Instantiate the matrix
**	inputMatrix = new InputMatrixGeneralBinary(inputFile, memory);
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
class InputMatrixGeneralBinary : public InputMatrix {
	public:
		//Constructor.  Initializes the class variables
		InputMatrixGeneralBinary( std::string p, long long memory, int, int);


		/*
		**Initializes the object with the first block.  Only run once
		**per object.  It is separate from the constructor so that
		**error codes can be returned
		*
		* RETURN
		** int ier -- error return code
		** 	ier=0 The matrix initialized without error
		 *    	ier=-1 to -10 check InputMatrix::init() error codes
		* 	ier=-11 File cannot be opened
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
		std::ifstream * inputFile;
};
#endif
