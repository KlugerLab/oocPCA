\name{fastRPCA-package}
\alias{fastRPCA-package}
\alias{fastRPCA}
\docType{package}
\title{
\packageTitle{fastRPCA}
}
\description{
\packageDescription{fastRPCA}
}
\details{

\packageDESCRIPTION{fastRPCA}
\packageIndices{fastRPCA}
}
\author{
\packageAuthor{fastRPCA}

Maintainer: \packageMaintainer{fastRPCA}
}
\references{
}
\keyword{ package }
\seealso{
}
\examples{
 k_ <- 20;
 m = 9E2;
 n = 10E3;
 B <- matrix(rexp(m*k_), m)
 C <- matrix(rexp(k_*n), k_)
 D <- B \%*\%C;
 dim(D)
 fastDecomp <- fastPCA(D, k=k_, diffsnorm=TRUE)
 str(fastDecomp)
}
