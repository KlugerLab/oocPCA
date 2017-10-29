#' Prepare the oocPCA result object 
#'
#' This is an internal function that formats the resulting oocPCA object so that it can be handled appropriately by R
#'
#' @param base_outname The basename for the matrices U, V, and S from the fastpca.xx call
#' @param k  Rank of the transformation.
#' @return A list containing the decomposition.
oocPCA_base <- function( base_outname, k) {
	m = file.size(sprintf("%s.U", base_outname))/(8*k);
	n = file.size(sprintf("%s.V", base_outname))/(8*k);

	Ufile = file(sprintf("%s.U", base_outname), "rb")
	U = t(matrix(readBin(Ufile, numeric(),size=8, n=(m*k*8) ),k));
	close(Ufile)
	Sfile = file(sprintf("%s.S", base_outname), "rb")
	S = matrix(readBin(Sfile, numeric(),size=8, n=(k*k*8) ),k);
	close(Sfile)
	Vfile = file(sprintf("%s.V", base_outname), "rb")
	V = t(matrix(readBin(Vfile, numeric(),size=8, n=(n*k*8) ),k));
	close(Vfile)
	unlink(sprintf("%s.U", base_outname));
	unlink(sprintf("%s.S", base_outname));
	unlink(sprintf("%s.V", base_outname));

	result <- list()
	result$U <- U	
	result$V <- V
	result$S <- S
	result$dims <- c(m,n);
	return (result);

}

#' Perform out-of-core SVD for a matrix in CSV format
#'
#' This function performs a nearly optimal rank-k approximation to the singular value decomposition inputMatrix = USV' on a matrix that is passed via CSV format.  Please see references for explanation of 'nearly optimal.'
#'
#' @param inputFile csv containing matrix to decompose
#' @param k Rank of decomposition. Default: 5
#' @param l Block size. Default k+2
#' @param mem Amount of memory the algorithm is allowed to use in bytes.  Must be greater than the number of columns * 8 (i.e. the memory to store one row).
#' @param its Number of normalized power iterations. Default: 2
#' @param diffsnorm Calculate 2-norm accuracy, i.e. ||A-USV||_2. 
#' @param centeringRow Center the rows prior to decomposition.
#' @param centeringColumn Center the columns prior to decomposition.
#' @param logTransform add 1 to every value and then take the log
#' @return A list containing the decomposition.
#' @examples
#' 
#' k_ <- 10;
#' m = 30;
#' n = 40;
#' B <- matrix(rexp(m*k_), m)
#' C <- matrix(rexp(k_*n), k_)
#' D <- B %*%C;
#' dim(D)
#' fn = "test_csv.csv"
#' write.table(D,file=fn,sep=',',col.names=FALSE, row.names=FALSE)
#' library(oocRPCA);
#' fastDecomp <- oocPCA_CSV(fn, k=k_, mem=n*8*100, diffsnorm=TRUE)
#' norm( D - fastDecomp$U %*% fastDecomp$S %*% t(fastDecomp$V))
#' unlink(fn)
#' @export
oocPCA_CSV<- function (inputFile,k=5, l, mem=2e9, its=2,diffsnorm=FALSE,centeringRow=FALSE, centeringColumn = FALSE, logTransform = FALSE) {
	if (!file.exists(inputFile)) {
		stop("File does not exist");
	}
	if (missing(l)){
		l = k+2;
	}
	
	if (k <=0 || l<=0  || mem <=0  ||  its <=0 ){
		stop("k,l,mem, its all must be positive numbers");
	}		

	if (!identical(TRUE,diffsnorm)  && !identical(FALSE, diffsnorm)) {
		    stop("diffsnorm must be TRUE or FALSE");
	}
	if (!identical(TRUE,centeringRow)  && !identical(FALSE, centeringRow)) {
		    stop("centeringRow must be TRUE or FALSE");
	}
	if (!identical(TRUE,centeringColumn)  && !identical(FALSE, centeringColumn)) {
		    stop("centeringColumn must be TRUE or FALSE");
	}
	if (identical(TRUE,centeringRow)  && identical(TRUE, centeringColumn)) {
		    stop("Both centeringRow and centeringColumn cannot be TRUE");
	}

	if (Sys.info()["sysname"] == 'Darwin') {
		executable <- 'fastpca.osx'
	}else{
		executable <- 'fastpca.linux'
		system(sprintf("chmod u+x %s/%s",system.file("build", package="oocRPCA"), executable));
	}

	base_outname = "oocRPCA.binmatrix"
	pcacall <- sprintf('%s/%s -k %d -csvI %s -binaryOutput %s -mem %0.f -l %d -its %d', system.file("build", package="oocRPCA"), executable, k, inputFile, base_outname, mem , l, its )
	if (identical(TRUE, centeringRow)) {
		pcacall <- paste(pcacall, " -center",  sep="")
	}else if (identical(TRUE, centeringColumn)) {
		pcacall <- paste(pcacall, " -ccenter",  sep="")
	}
	if (identical(TRUE, logTransform)) {
		pcacall <- paste(pcacall, " -log",  sep="")
	}
	
	pcareturn <- system(pcacall);
	if (pcareturn != 1) {
		stop(sprintf("An unkonwn error has occurred during call: %s", pcacall));
	}
	return (oocPCA_base(base_outname, k));

}
#' Convert CSV to binary format
#'
#' The the computation time of an out-of-core algorithm is typically determined by the cost of loading 
#' the matrix from the slow memory (i.e. the hard disk).  It is much faster to load a matrix from a binary format
#' than CSV, because the former represents it in much more compressed format.  This function converts 
#' a CSV to the compressed binary format.
#'
#' @param inputFile file containing matrix to convert
#' @param outputFile name of output file
#' @export
oocPCA_csv2binary<- function (inputFile,outputFile) {
	if (!file.exists(inputFile)) {
		stop("File does not exist");
	}
	if (Sys.info()["sysname"] == 'Darwin') {
		executable <- 'csv2binary.osx'
	}else{
		executable <- 'csv2binary.linux'
		system(sprintf("chmod u+x %s/%s",system.file("build", package="oocRPCA"), executable));
	}
	system(sprintf('%s/%s -csvI %s -binO %s',system.file("build", package="oocRPCA"), executable, inputFile, outputFile ))
}


#' Perform out-of-core SVD for a matrix in binary format
#'
#' This function performs a nearly optimal rank-k approximation to the singular value decomposition inputMatrix = USV' on a matrix that is passed via binary format.  Please see references for explanation of 'nearly optimal.'
#'
#' @param inputFile file containing matrix to decompose
#' @param m Number of rows in matrix
#' @param n Number of columns in matrix
#' @param k Rank of decomposition. Default: 5
#' @param l Block size. Default k+2
#' @param mem Amount of memory the algorithm is allowed to use in bytes.  Must be greater than the number of columns * 8 (i.e. the memory to store one row).
#' @param its Number of normalized power iterations. Default: 2
#' @param diffsnorm Calculate 2-norm accuracy, i.e. ||A-USV||_2. 
#' @param centeringRow Center the rows prior to decomposition.
#' @param centeringColumn Center the columns prior to decomposition.
#' @param logTransform add 1 to every value and then take the log
#' @return A list containing the decomposition.
#' @examples
#' 
#'k_ <- 10;
#'m = 1e3;
#'n = 1e2;
#'B <- matrix(rexp(m*k_), m)
#'C <- matrix(rexp(k_*n), k_)
#'D <- B %*%C;
#'dim(D)
#'fn = "test_csv.csv"
#'fnb = "test_csv.bin"
#'write.table(D,file=fn,sep=',',col.names=FALSE, row.names=FALSE)
#'oocPCA_csv2binary(fn, fnb)
#'library(oocRPCA);
#'fastDecomp <- oocPCA_BIN(fnb, m,n, k=k_, mem=floor(n*8*m/3), diffsnorm=TRUE)
#'norm( D - fastDecomp$U %*% fastDecomp$S %*% t(fastDecomp$V))
#'unlink(fn)
#'unlink(fnb)
#' @export
oocPCA_BIN<- function (inputFile,m, n,k=5, l, mem=2e9, its=2,diffsnorm=FALSE,centeringRow=FALSE, centeringColumn = FALSE, logTransform = FALSE) {
	if (!file.exists(inputFile)) {
		stop("File does not exist");
	}
	if (missing(l)){
		l = k+2;
	}
	
	if (k <=0 || l<=0  || mem <=0  ||  its <=0 ){
		stop("k,l,mem, its all must be positive numbers");
	}		

	if (!identical(TRUE,diffsnorm)  && !identical(FALSE, diffsnorm)) {
		    stop("diffsnorm must be TRUE or FALSE");
	}
	if (!identical(TRUE,centeringRow)  && !identical(FALSE, centeringRow)) {
		    stop("centeringRow must be TRUE or FALSE");
	}
	if (!identical(TRUE,centeringColumn)  && !identical(FALSE, centeringColumn)) {
		    stop("centeringColumn must be TRUE or FALSE");
	}
	if (identical(TRUE,centeringRow)  && identical(TRUE, centeringColumn)) {
		    stop("Both centeringRow and centeringColumn cannot be TRUE");
	}
	#TODO: check if m*n*8 = filesize

	if (Sys.info()["sysname"] == 'Darwin') {
		executable <- 'fastpca.osx'
	}else{
		executable <- 'fastpca.linux'
		system(sprintf("chmod u+x %s/%s",system.file("build", package="oocRPCA"), executable));
	}

	base_outname = "oocRPCA.binmatrix"
	pcacall <- sprintf('%s/%s -m %d -n %d -k %d -binaryI %s -binaryOutput %s -mem %0.f -l %d -its %d', system.file("build", package="oocRPCA"), executable, m, n, k, inputFile, base_outname, mem , l, its )
	if (identical(TRUE, centeringRow)) {
		pcacall <- paste(pcacall, " -center",  sep="")
	}else if (identical(TRUE, centeringColumn)) {
		pcacall <- paste(pcacall, " -ccenter",  sep="")
	}
	if (identical(TRUE, logTransform)) {
		pcacall <- paste(pcacall, " -log",  sep="")
	}
	
	pcareturn <- system(pcacall);
	if (pcareturn != 1) {
		stop(sprintf("An unkonwn error has occurred during call: %s", pcacall));
	}
	return (oocPCA_base(base_outname, k));

}
