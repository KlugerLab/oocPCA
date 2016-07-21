#ifndef INPUTMATRIXEIGENSTRAT_H 
#define INPUTMATRIXEIGENSTRAT_H
#include "InputMatrix.h"
/*********************************************************************
 *  Input Matrix Class for Eigenstrat .geno files
 * 
 *  
 **  Use an InputMatrixEigensrat class to represent a matrix stored in the
 **  eigenstrat .geno format.  The difference between the "InCore" is that it does not use
 **  the 2bit per snp representation. 
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
class InputMatrixEigenstrat : public InputMatrix {
	public:
		//Constructor.  Initializes the class variables
		InputMatrixEigenstrat( std::string p, long long memory);


		/*
		**
		*  Initializes the object with the first block.  Only run once
		*  per object. 
		*  The matrix is first traversed in order to determine the number
		*  of rows (m) and number of columns (n) and also check the consistency
		*  of the input.  Then it is traversed a second time, to read the values
		*  into A.  Please see README for optimization ideas.  
		* It is separate from the constructor so that error codes can be returned
		*
		* RETURN
		** int ier -- error return code
		** 	ier=0 The matrix initialized without error
		 *    	ier=-1 to -10 check InputMatrix::init() error codes
		* 	ier=-11 File cannot be opened
		* 	ier=-12 File integrity error
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
		//A stream of the actual input file
		std::ifstream eigenstratFile;
};
#endif
