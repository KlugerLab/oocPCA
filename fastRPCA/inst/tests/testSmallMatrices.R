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
		k_ <- 20;
		B <- matrix(rexp(m*k_), m)
		C <- matrix(rexp(k_*n), k_)
		D <- B %*%C;
		fn = "test_csv.csv"
		write.table(D,file=fn,sep=',',col.names=FALSE, row.names=FALSE)

		memories <- c(5*n*8 +2342 , 100*n*8, 1E9)
		for (mem in memories) {
			fastDecomp <- fastPCA_CSV(fn, k=k_, mem=mem, diffsnorm=TRUE)
		}
	diffNorm <- norm(D - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')
	test_that( sprintf("CSV - m: %d, n: %d, k: %d", fastDecomp$dim[0], fastDecomp[1], k_), {
			  expect_that(diffNorm,equals(0));
		})
	}
}

