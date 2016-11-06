#install_github("KlugerLab/FastPCA",subdir="fastRPCA", host="git.yale.edu/api/v3", auth_token="3135ea8797c3ae6471fc96fc92cd382d0ab330f5")
#testthat::test_dir(sprintf("%s/testthat", system.file("tests", package="fastRPCA"))
#' Prepare the fastPCA result object 
#'
#' This is an internal function that formats the resulting fastPCA object so that it can be handled appropriately by R
#'
#' @useDynLib fastRPCA
#' @importFrom Rcpp evalCpp
#' @param result Value of one of the other functions.
#' @param k  Rank of the transformation.
#' @return A FastPCA object, containing the decomposition.
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
#' @param centeringColumn Center the columns prior to decomposition.
#' @return A FastPCA object, containing the decomposition.
#' @examples
#' 
#' k_ <- 20;
#' m = 9E2;
#' n = 10E3;
#' B <- matrix(rexp(m*k_), m)
#' C <- matrix(rexp(k_*n), k_)
#' D <- B %*% C;
#' dim(D)
#' fastDecomp <- fastPCA(D, k=k_, diffsnorm=TRUE)
#' str(fastDecomp)
#' norm(D - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')
#' @export
fastPCA<- function (inputMatrix,k=5, l, its=2,diffsnorm=FALSE,centeringRow=FALSE, centeringColumn = FALSE) {

	if (missing(l)){
		l = k+2;
	}
	if (!"matrix" == class(inputMatrix)) {
		stop("inputMatrix needs to be a matrix");
	}
	m <- nrow(inputMatrix); n <- ncol(inputMatrix);
	if (k <=0 || l<=0  ||  its <=0 ){
		stop("k,l,mem, its all must be positive numbers");
	}		

	if (k>l ){
		stop("k must be less than or equal to l");
	}

	if (l> min(m,n) ){
		stop("l should not exceed the smallest dimension of the matrix");
	}

	if (k > min(m,n) ){
		stop("k must be less than or equal to the smallest dimension of the matrix");
	}		
	if (!identical(TRUE,diffsnorm)  && !identical(FALSE, diffsnorm)) {
		    stop("diffsnorm must be TRUE or FALSE");
	}
	if (!identical(TRUE,centeringColumn)  && !identical(FALSE, centeringColumn)) {
		    stop("centeringColumn must be TRUE or FALSE");
	}
	if (identical(TRUE,centeringRow)  && identical(TRUE, centeringColumn)) {
		    stop("Both centeringRow and centeringColumn cannot be TRUE");
	}
	if (identical(TRUE,centeringRow)  ||  identical(TRUE, centeringColumn)) {
		    stop("For in-memory matrices, centering matrices is not supported ");
	}
	#Swapping n and m becaue it is column major.  So, the U will be the V and vice versa.
	n = nrow(inputMatrix);
	m = ncol(inputMatrix);
	
	result = .Call( 'fastRPCA', 'memory', as.matrix(inputMatrix), m, n, k,l,its, -1, FALSE, FALSE,diffsnorm, PACKAGE = 'fastRPCA');
	V<- t(matrix(result$U, ncol = m, nrow = k))
	result$U <- t(matrix(result$V, ncol = n, nrow = k))
	result$V <-  V;
	result$S <- matrix(result$S, ncol = k, nrow = k)
	return (result)

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
#' k_ <- 20;
#' m = 200;
#' n = 100;
#' B <- matrix(rexp(m*k_), m)
#' C <- matrix(rexp(k_*n), k_)
#' D <- B %*%C;
#' dim(D)
#' fn = "test_csv.csv"
#' write.table(D,file=fn,sep=',',col.names=FALSE, row.names=FALSE)
#' fastDecomp <- fastPCA_CSV(fn, k=k_, mem=n*8*100, diffsnorm=TRUE)
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

	result = .Call( 'fastRPCA', 'csv', inputFile, -1, -1, k,l,its, mem,centeringRow, centeringColumn,diffsnorm, PACKAGE = 'fastRPCA');
	m = result$dims[1]
	n = result$dims[2]

	result$U<- matrix(result$U, nrow = m, ncol = k)
	result$V <- matrix(result$V, nrow = n, ncol = k)
	result$S <- matrix(result$S, ncol = k, nrow = k)
	return (fastPCA_base(result, k));

}
#' Perform fast SVD for a matrix in BED format
#'
#' This function performs a nearly optimal rank-k approximation to the singular value decomposition inputMatrix = USV' on a GWAS matrix that is passed via BED format.  Please see references for explanation of 'nearly optimal.'
#' 
#' Missing values will be imputed by replacing the missing value with the row average.  If centeringRow is set to TRUE, then this just means replacing the missing values with 0.
#' Note that unlike in the other functions, centeringRow defaults to true here
#'
#' @param inputFile Base file name for the .bed,.fam, .bim files.  All three of these files must be present in the same folder: <inputFile>.bed <inputFile>.fam <inputFile>.bim.
#' @param k Rank of decomposition. Default: 5
#' @param l Block size. Default k+2
#' @param mem Amount of memory the algorithm is allowed to use in bytes.  Must be greater than the number of columns * 8 (i.e. the memory to store one row).
#' @param phenoFile Path to a phenotype file, the attributes of which will be returned in the FastPCA object, in the same order as the V matrix
#' @param its Number of normalized power iterations. Default: 2
#' @param diffsnorm Calculate 2-norm accuracy, i.e. ||A-USV||_2. 
#' @param centeringRow Center the rows prior to decomposition. Default: TRUE
#' @param centeringColumn Center the columns prior to decomposition. Default: FALSE
#' @return A FastPCA object, containing the decomposition.  Also, has the ordered phenotype information
#' @examples
#' 
#' fn <-sprintf("%s/example_with_impute", system.file("tests", package="fastRPCA")) ;
#' fastDecomp <- fastPCA_BED(fn, k=5,l=5, mem=1E3, diffsnorm=TRUE)
#' str(fastDecomp)
#' 
#' #Now, to check, let's take the same matrix and do the imputation and row centering in R
#' Dm <- matrix(c( 1,1,1,NA,2,2, 1, 0,NA,0, 0,1,1,2,1,2,2,1,0,0,0,1,1,2,2,1,1,NA,1,2,2,2,NA,1,1),nrow=7,ncol=5, byrow=TRUE);
#' D <- Dm;
#' k <- which(is.na(D), arr.ind=TRUE)
#' D[k] <- rowMeans(D, na.rm=TRUE)[k[,1]]
#' Dt <- t(scale(t(D), center=TRUE, scale=FALSE));
#' norm(Dt - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')
#' 
#' 
#' #Another example
#' #Before publication, this example will be replaced by test data
#' #fastDecompBed10E3 <- fastPCA_BED("/data/GERA_DECRYPTED/LindermanAnalysis/EUR/Benchmark/eur_100000_62318/eur_100000_62318",k = 20, mem=8e+10, l=22,centeringRow=TRUE, phenoFile="/data/GERA_DECRYPTED/38852/PhenoGenotypeFiles/RootStudyConsentSet_phs000674.GERA.v1.p1.c1.HROO/PhenotypeFiles/phs000674.v1.pht003641.v1.p1.c1.Survey.HROO.txt")
#' #str(fastDecompBed10E3)

#' #plot(fastDecompBed10E3$V[,1], fastDecompBed10E3$V[,2], pch='.', col=factor(fastDecompBed10E3$phenos$RACE))
#'
#' @export
fastPCA_BED<- function (inputFile,k=5, l, mem, its=2,phenoFile=-1,diffsnorm=FALSE,centeringRow=TRUE, centeringColumn = FALSE) {

	bedFile <- sprintf("%s.bed", inputFile);
	famFile <- sprintf("%s.fam", inputFile);
	bimFile <- sprintf("%s.bim", inputFile);

	if (!file.exists(bedFile)) {
		stop("Bed file does not exist");
	}

	if (!file.exists(famFile)) {
		stop("Fam file does not exist");
	}


	if (!file.exists(bimFile)) {
		stop("Bim file does not exist");
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
	result = .Call( 'fastRPCA', 'bed', inputFile, -1, -1, k,l,its, mem,centeringRow, centeringColumn,diffsnorm, PACKAGE = 'fastRPCA');
	m = result$dims[1]
	n = result$dims[2]

	result$U <- t(matrix(result$U, ncol = result$dim[1], nrow = k))
	result$V <- t(matrix(result$V, ncol = result$dim[2], nrow = k))
	result$S <- matrix(result$S, ncol = k, nrow = k)
	fam <- read.table(paste(inputFile,'.fam', sep=""), sep="", header=FALSE);
	row.names(result$V) <- fam[,2];

	if (phenoFile != -1) {
		phenos <- read.table(phenoFile, header=TRUE);
		row.names(phenos) <- phenos$SUBJID;
		result$phenos = phenos[row.names(result$V),];
	}
	return (result);

}
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
