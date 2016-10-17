#' Prepare the fastPCA result object 
#'
#' This is an internal function that formats the resulting fastPCA object so that it can be handled appropriately by R
#'
#' @useDynLib fastRPCA
#' @importFrom Rcpp evalCpp
#' @param result Value of one of the other functions.
#' @param k  Rank of the transformation.
#' @return A FastPCA object, containing the decomposition.
#' @examples
fastPCA_base <- function(result, k) { result$U <- t(matrix(result$U, ncol = result$dim[1], nrow = k))
	result$V <- t(matrix(result$V, ncol = result$dim[2], nrow = k))
	result$S <- matrix(result$S, ncol = k, nrow = k)
	return (result);

}
#' Perform fast SVD for a matrix in the memory
#'
#' This function performs a nearly optimal rank-k approximation to the singular value decomposition inputMatrix = USV' on a matrix that is passed via the memory.  Please see references for explanation of 'nearly optimal.'
#'
#' @param inputMatrix Matrix to decompose.
#' @param k Rank of decomposition. Default: 5
#' @param l Block size. Default k+2
#' @param its Number of normalized power iterations. Default: 2
#' @param diffsnorm Calculate 2-norm accuracy, i.e. ||A-USV||_2. 
#' @param centeringRow Center the rows prior to decomposition.
#' @param centeringRow Center the columns prior to decomposition.
#' @return A FastPCA object, containing the decomposition.
#' @examples
#' 
#' k_ <- 20;
#' m = 9E2;
#' n = 10E3;
#' B <- matrix(rexp(m*k_), m)
#' C <- matrix(rexp(k_*n), k_)
#' D <- B %*%C;
#' dim(D)
#' fastDecomp <- fastPCA(D, k=k_, diffsnorm=TRUE)
#' str(fastDecomp)
#' norm(D - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')
#' @export
fastPCA<- function (inputMatrix,k=5, l, its=2,diffsnorm=0,centeringRow=0, centeringColumn = 0) {
	if (missing(l)){
		l = k+2;
	}
	#Swapping n and m becaue it is column major.  So, the U will be the V and vice versa.
	n = nrow(inputMatrix);
	m = ncol(inputMatrix);
	
	mem= m*n*8;

	print("About to call");
	result = .Call( 'fastRPCA', 'memory', as.matrix(inputMatrix), m, n, k,l,its, mem,centeringRow, centeringColumn,diffsnorm, PACKAGE = 'fastRPCA');
	V<- t(matrix(result$U, ncol = m, nrow = k))
	result$U <- t(matrix(result$V, ncol = n, nrow = k))
	result$V <-  V;
	result$S <- matrix(result$S, ncol = k, nrow = k)
	return (result)

}
#
fastPCA_CSV <- function (inputFile,k=5, mem=144, l=5, its=2,diffsnorm=0,centeringRow=0, centeringColumn = 0) {
	result = .Call( 'fastRPCA','csv', inputFile, -1, -1, k,l,its, mem,centeringRow, centeringColumn,diffsnorm);
	return (fastPCA_base(result, k));

}
#fastPCA_Phenos <- function (inputFile) {
#
#
#}
#fastPCA_BED <- function (inputFile,k=5, mem=144, l=5, its=2,phenoFile=-1, diffsnorm=0,centeringRow=0, centeringColumn = 0) {
#
#	#result = .Call( 'fastRPCA','bed', inputFile, -1, -1, k,l,its, mem,centeringRow, centeringColumn,diffsnorm);
#	result$U <- t(matrix(result$U, ncol = result$dim[1], nrow = k))
#	result$V <- t(matrix(result$V, ncol = result$dim[2], nrow = k))
#	result$S <- matrix(result$S, ncol = k, nrow = k)
#	fam <- read.table(paste(inputFile,'.fam', sep=""), sep=" ", header=FALSE);
#	row.names(result$V) <- fam[,2];
#
#	#TODO: make this an optional argument to this function
#	if (phenoFile != -1) {
#		phenos <- read.table(phenoFile, header=TRUE);
#		row.names(phenos) <- phenos$SUBJID;
#		result$phenos = phenos[row.names(result$V),];
#	}
#	return (result);
#
#}
#
#fastPCA_EIGEN <- function (inputFile,k=5, mem=144, l=5, its=2,diffsnorm=0,centering=0) {
#
#	#result = .Call( 'fastRPCA','eigen', inputFile, -1, -1, k,l,its, mem,centering,diffsnorm);
#	result$U <- t(matrix(result$U, ncol = result$dim[1], nrow = k))
#	result$V <- t(matrix(result$V, ncol = result$dim[2], nrow = k))
#	result$S <- matrix(result$S, ncol = k, nrow = k)
#	return (result);
#
#}
