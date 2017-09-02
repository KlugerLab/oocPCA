#' Prepare the fastPCA result object 
#'
#' This is an internal function that formats the resulting fastPCA object so that it can be handled appropriately by R
#'
#' @useDynLib fastRPCA
#' @importFrom Rcpp evalCpp
#' @param result Value of one of the other functions.
#' @param k  Rank of the transformation.
#' @return A FastPCA object, containing the decomposition.
fastPCA_base <- function( base_outname, k) {
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
	reformed <- U %*% S %*% t(V);
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

#' Perform fast SVD for a matrix in CSV format
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
#' @return A FastPCA object, containing the decomposition.
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
#' fastDecomp <- fastPCA_CSV(fn, k=k_, mem=n*8*100, diffsnorm=TRUE)
#' norm( D - fastDecomp$U %*% fastDecomp$S %*% t(fastDecomp$V))
#' @export
fastPCA_CSV<- function (inputFile,k=5, l, mem, its=2,diffsnorm=FALSE,centeringRow=FALSE, centeringColumn = FALSE) {
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


	base_outname = "oocRPCA.binmatrix"
	pcacall <- sprintf('%s/fastpca.xx -k %d -csvI %s -binaryOutput %s -mem %d -l %d -its %d', system.file("build", package="oocRPCA"), k, inputFile, base_outname, mem , l, its )
	if (identical(TRUE, centeringRow)) {
		pcacall <- paste(pcacall, " -center",  sep="")
	}else if (identical(TRUE, centeringColumn)) {
		pcacall <- paste(pcacall, " -ccenter",  sep="")
	}
	#print(pcacall);
	system(pcacall);
	return (fastPCA_base(base_outname, k));

}
