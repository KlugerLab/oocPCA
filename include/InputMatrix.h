#ifndef INPUTMATRIX_H 
#define INPUTMATRIX_H
#include "macros.h"
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "fastpca_linear_algebra.hpp"
/*********************************************************************
 *  Input Matrix General Class
 * 
 *  
 **  InputMatrix is an abstract class based on the Iterator pattern,
 ** representing the matrix, A, which is to be analyzed.  It is written as an
 ** iterator so that one can easily iterate through blocks of the matrix.  All
 ** logic for iterating through the blocks of the matrix is restricted to the
 ** InputMatrix class, while its subclasses contain implementations specific to
 ** each file format.
 **
 *  THEORY
 *
 ** A major goal of this class is to reduce the number of times a matrix must
 ** be read from the source (e.g. disk).  In order to do so, it reuses the last
 ** block loaded from the previous traversal.  For example, consider a row-major
 ** matrix with 8 rows, and a block size of 2 rows.  After the first traversal,
 ** rows 7:8 are already loaded, so it makes sense that the next traversal begins
 ** with this block.  The alternative, of starting from the first row, will
 ** unnecessarily cause an extra block to be loaded for every traversal.
 ** InputMatrix begins with the last block of the previous traversal (if it
 ** exists), and iterates through the rest of the matrix.
 **
 ** A nice result of our "reusing" strategy is that when the memory is large
 ** enough to load the entire matrix, then the block is simply the entire matrix.
 ** As such, every traversal will "reuse" the entire matrix, and it will only be
 ** loaded once.  In other words, an InputMatrix object is the most optimal
 ** representation of a matrix, regardless of whether it is "out-of-core", or if
 ** the entire matrix can fit in the memory.  The user can be totally oblivious. 
 *  USAGE 
 *
 **  InputMatrix is an abstract class that must be implemented for various
 **  input formats. InputMatrix itself cannot be instantiated.
 *

 */
class InputMatrix {
	public:
		//The number of rows in the whole matrix
		long long m;	

		//The number of columns in the whole matrix
		long long n;	

		//The size of the block.  (e.g. If iterating through rows, this is
		//the number of rows in a block)
		int blockSize;	

		//The start of the block. (e.g. if iterating through rows, this
		//is the starting row of the block
		int blockStart;	

		//The (inclusive) end of the block. (e.g. if iterating through rows, this
		//is the last row of the block)
		int blockEnd;	

		//Pointer to the currently loaded block
		double * block;

		/*
		** Initializes the object with the first block.  Only run once
		** per object. Intended to be overriden in derived classes.
		*
		* RETURN
		 ** int   ier --   error return code
		 *   	  ier=0 The matrix initialized without error
		 *    	  ier=-1 block allocation was unsuccessful
		 *        ier=-2 Insufficient memory set
		 *        ier>10 derived init() functions for  error codes

		 */
		virtual int init() = 0;

		//hasNext() checks if there is another block to be loaded, and
		//if so, calls next() to load it (unless it is the first
		//iteration, in which case the previous block is reused)
		bool hasNext();

		//loadNextBlock() loads the block.  Only called by hasNext().
		//Must be implemented in the derived class, as loading data is
		//specific to the format being loaded.  
		virtual int loadNextBlock() = 0;

		void imputeMissingOn();
		void centerRowsOn();
		bool centerColumns();


	protected:
		
		//Constructor, only used to initialize the path and memory variables
		InputMatrix(std::string p, long long maxMemoryBytes);
		double * rowMeans;
		double * colMeans;
		bool * rowMeansFlags;

		int preprocessBlock();

		//Path to the file containing the full matrix
		std::string path;	
		

		//Indicates the beginning of a matrix traversal.  When true,
		//hasNext() will not load another block, as it will reuse the
		//block from the previous iteration
		bool firstIt;

		//Number of rows in a normal (not a remainder) block. For
		//example, in a row-major implementation, if there are 9 rows,
		//and block size is 4, then the calculatedBlockSize is 2  
		int calculatedBlockSize;

		//Size of the "remainder" block. For example, in a row-major
		//implementation, if there are 9 rows and block size is 4, then
		//calculatedRemainderBlockSize is 1
		int calculatedRemainderBlockSize;


		//The row/column at which the traversal started.  Most
		//importantly, it is the start of the block at which the last
		//iteration finished.  In other words, if there is a row-major
		//matrix of size 8, and the blockSize is 2, then for the first
		//traversal traversalStart == 0, but for the second traversal,
		//it will be 6. 
		int traversalStart;


		//True if this iteration will be the "remainder iteration"
		bool remainderIt;

		//Max memory to allocate for a given block
		long long maxMemoryBytes;	

	private:
		bool centerRows;
		bool centerColumnsFlag;
		bool imputeMissing;
};
#endif
