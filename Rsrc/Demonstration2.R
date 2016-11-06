#setwd('/data/Linderman/FastPCA4/Rsrc/')
#install_github("KlugerLab/FastPCA",subdir="fastRPCA", host="git.yale.edu/api/v3", auth_token="3135ea8797c3ae6471fc96fc92cd382d0ab330f5") 
#testthat::test_dir(sprintf("%s/testthat", system.file("tests", package="fastRPCA"))

source('fastRPCA.R') 

require(fastRPCA)
#####Small example
k_ <- 20;
m = 9E2;
n = 10E3;
B <- matrix(rexp(m*k_), m)
C <- matrix(rexp(k_*n), k_)
D <- B %*%C;
dim(D)


#require(fastRPCA)
#m=12;
#n =10;
#k_ = 5
#B <- matrix(rexp(m*k_), m)
#C <- matrix(rexp(k_*n), k_)
#D <- B %*%C;

require(fastRPCA)
#####Small example
k_ <- 20;
m = 9E2;
n = 10E3;
B <- matrix(rexp(m*k_), m)
C <- matrix(rexp(k_*n), k_)
D <- B %*%C;
dim(D)

Dcc <- scale(D, center=TRUE, scale=FALSE);
Drc <- t(scale(t(D), center=TRUE, scale=FALSE));
fastDecomp <- fastPCA(D, k=k_, diffsnorm=TRUE, centeringRow=TRUE)
#fastDecomp <- fastPCA(D, k=k_, diffsnorm=TRUE, centeringRow = TRUE)

str(fastDecomp)
norm(Dcc - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')

slowDecomp <- svd(D)
norm(D - slowDecomp$u %*% diag(slowDecomp$d) %*%t(slowDecomp$v), type='2')


#####Small example csv
require(fastRPCA)
k_ <- 20;
m = 97E1;
n = 93E0;
#m=12;
#n =10;
#k_ =5
B <- matrix(rexp(m*k_), m)
C <- matrix(rexp(k_*n), k_)
D <- B %*%C;
Dcc <- scale(D, center=TRUE, scale=FALSE);
Drc <- t(scale(t(D), center=TRUE, scale=FALSE));
fn = "/Users/george/Research_Local/FastPCA4/FastPCA/fastRPCA/test_csv.csv"
write.table(D,file=fn,sep=',',col.names=FALSE, row.names=FALSE)
#fastDecomp <- fastPCA_CSV(fn, k=k_, mem=(10E6 + m* n*8), diffsnorm=TRUE)
fastDecomp <- fastPCA_CSV(fn, k=k_, mem=( 12312 + 8*m*8), diffsnorm=TRUE,centeringColumn=TRUE)
norm(Dcc - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')
require(fastRPCA)
k_ <- 20;
m = 100;
n = 50;
#m=12;
#n =10;
#k_ =5
B <- matrix(rexp(m*k_), m)
C <- matrix(rexp(k_*n), k_)
D <- B %*%C;
fn = "/Users/george/Research_Local/FastPCA4/FastPCA/fastRPCA/test_csv.csv"
write.table(D,file=fn,sep=',',col.names=FALSE, row.names=FALSE)
Dcc <- scale(D, center=TRUE, scale=FALSE);
Drc <- t(scale(t(D), center=TRUE, scale=FALSE));

#fastDecomp <- fastPCA_CSV(fn, k=k_, mem=(4342), diffsnorm=TRUE,centeringCol=TRUE)
fastDecomp <- fastPCA_CSV(fn, k=k_, mem=(4342), diffsnorm=TRUE,centeringCol=TRUE)
norm(Dcc - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')

####Big example: in core R
B <- matrix(rexp(1E4*20), 1E4)
C <- matrix(rexp(20*2E5), 20)
D <- B %*%C;
dim(D)
fastDecomp <- fastPCA(D, k=20, diffsnorm=TRUE)
str(fastDecomp)

##################
#A <- matrix(c(1,1,1,0,0,0,1,2,1,2,2,1,1,0,1,0,0,1,2,2,2,1,1,0,0,0,0,1,1,1,2,2,1,1,0),nrow=7,ncol=5, byrow=TRUE);
A <- matrix(c(1,1,1,2,2,2,1,0,1,0,0,1,1,2,1,2,2,1,0,0,0,1,1,2,2,2,2,1,1,1,0,0,1,1,2),nrow=7,ncol=5, byrow=TRUE);
fastDecomp <- fastPCA_BED("/Users/george/Research_Local/FastPCA4/FastPCA/fastRPCA/inst/tests/example", k=5, l=5, mem=2*5*8,diffsnorm=TRUE, centeringRow=TRUE);

Acc <- scale(A,center=TRUE, scale=FALSE)
Arc <- t(scale(t(A),center=TRUE, scale=FALSE))

norm(Arc - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')
reconA <- fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V)


fn <-  "/Users/george/Research_Local/FastPCA4/FastPCA/fastRPCA/inst/tests/example.fam";
fam <- read.table(fn, sep=" ", header=FALSE);





###################
fastDecompBed10E3 <- fastPCA_BED("/data/GERA_DECRYPTED/LindermanAnalysis/EUR/Benchmark/eur_100000_62318/eur_100000_62318",k = 20, mem=8e+10, l=22,centeringRow=TRUE, phenoFile="/data/GERA_DECRYPTED/38852/PhenoGenotypeFiles/RootStudyConsentSet_phs000674.GERA.v1.p1.c1.HROO/PhenotypeFiles/phs000674.v1.pht003641.v1.p1.c1.Survey.HROO.txt")
str(fastDecompBed10E3)

#plot(fastDecompBed10E3$V[,1], fastDecompBed10E3$V[,2], pch='.', col=factor(fastDecompBed10E3$phenos$RACE))

fastDecompBed <- fastPCA_BED("/data/GERA_DECRYPTED/LindermanAnalysis/EUR/Benchmark/eur_20000_62318/eur_20000_62318",k = 20, mem=1e+9, l=22, phenoFile="/data/GERA_DECRYPTED/38852/PhenoGenotypeFiles/RootStudyConsentSet_phs000674.GERA.v1.p1.c1.HROO/PhenotypeFiles/phs000674.v1.pht003641.v1.p1.c1.Survey.HROO.txt",centeringRow = TRUE)

#plot(fastDecompBed$V[,1], fastDecompBed$V[,2], pch='.', col=factor(fastDecompBed$phenos$GENERALHEALTHCAT))


######CSV

fastDecompDropSeqCSV <- fastPCA_CSV('/data/Linderman/DropSeq/GSE63472_P14Retina_merged_digital_expression_headerless.csv',k=5, l=7, mem=50E9 )
