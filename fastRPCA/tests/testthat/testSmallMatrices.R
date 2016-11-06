require(fastRPCA)

context("Small matrices")
testDim <- c(10E1, 10E2)
for (m in testDim){
	for (n in testDim){
		k_ <- 20;
		B <- matrix(rexp(m*k_), m)
		C <- matrix(rexp(k_*n), k_)
		D <- B %*%C;
		dim(D)
		fastDecomp <- fastPCA(D, k=k_, diffsnorm=TRUE)
		diffNorm <- norm(D - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')
		test_that( sprintf("Memory - m: %d, n: %d, k: %d", fastDecomp$dim[0], fastDecomp[1], k_), {
				  expect_that(diffNorm,equals(0));
})

	}
}


testDim <- c(50, 100)
for (m in testDim){
	for (n in testDim){
		for (centering in c(0,1,2)){



			k_ <- 20;
			B <- matrix(rexp(m*k_), m)
			C <- matrix(rexp(k_*n), k_)
			D <- B %*%C;
			fn = "test_csv.csv"
			write.table(D,file=fn,sep=',',col.names=FALSE, row.names=FALSE)

			memories <- c(5*n*8 +2342 , 100*n*8, 1E9)


			if (centering == 0 ){
				rowCentering_ <- FALSE;
				colCentering_ <- FALSE;
				Dt <- D;
			}else if (centering ==1){
				rowCentering_ <- TRUE;
				colCentering_ <- FALSE;
				Dt <- t(scale(t(D), center=TRUE, scale=FALSE));
			}else {
				rowCentering_ <- FALSE;
				colCentering_ <- TRUE;
				Dt <- scale(D, center=TRUE, scale=FALSE);
			}

			for (mem in memories) {
				fastDecomp <- fastPCA_CSV(fn, k=k_, mem=mem, diffsnorm=TRUE, centeringColumn = colCentering_, centeringRow = rowCentering_)
				diffNorm <- norm(Dt - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')
				test_that( sprintf("CSV - m: %d, n: %d, k: %d, mem %d, centeringRow :%d, centeringCol: %d", fastDecomp$dims[1], fastDecomp$dims[2], k_, mem, rowCentering_, colCentering_), {
						  expect_that(diffNorm,equals(0));
})
			}
		}
	}
}





context("Small BED files")
D <- matrix(c(1,1,1,2,2,2,1,0,1,0,0,1,1,2,1,2,2,1,0,0,0,1,1,2,2,2,2,1,1,1,0,0,1,1,2),nrow=7,ncol=5, byrow=TRUE);
fn <- "/Users/george/Research_Local/FastPCA4/FastPCA/fastRPCA/inst/tests/example";
for (centering in c(0,1,2)){
	k_ = 5;
	l_ =5;
	n = 5;
	memories <- c(2*n*8,3*n*8, 4*n*8, 1E9 )
	if (centering == 0 ){
		rowCentering_ <- FALSE;
		colCentering_ <- FALSE;
		Dt <- D;
	}else if (centering ==1){
		rowCentering_ <- TRUE;
		colCentering_ <- FALSE;
		Dt <- t(scale(t(D), center=TRUE, scale=FALSE));
	}else {
		rowCentering_ <- FALSE;
		colCentering_ <- TRUE;
		Dt <- scale(D, center=TRUE, scale=FALSE);
	}

	for (mem in memories) {
		fastDecomp <- fastPCA_BED(fn, k=k_,l=l_, mem=mem, diffsnorm=TRUE, centeringColumn = colCentering_, centeringRow = rowCentering_)
		diffNorm <- norm(Dt - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')
		test_that( sprintf("CSV - m: %d, n: %d, k: %d, mem %d, centeringRow :%d, centeringCol: %d", fastDecomp$dims[1], fastDecomp$dims[2], k_, mem, rowCentering_, colCentering_), {
				  expect_that(diffNorm,equals(0));
})
	}
}



