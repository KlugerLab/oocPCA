#setwd('/data/Linderman/FastPCA4/Rsrc/')
#install_github("KlugerLab/FastPCA",subdir="fastRPCA", host="git.yale.edu/api/v3", auth_token="3135ea8797c3ae6471fc96fc92cd382d0ab330f5") 

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
fastDecomp <- fastPCA(D, k=k_, diffsnorm=TRUE)
str(fastDecomp)
norm(D - fastDecomp$U %*% fastDecomp$S %*%t(fastDecomp$V), type='2')

slowDecomp <- svd(D)
norm(D - slowDecomp$u %*% diag(slowDecomp$d) %*%t(slowDecomp$v), type='2')

####Big example: in core R
B <- matrix(rexp(1E4*20), 1E4)
C <- matrix(rexp(20*2E5), 20)
D <- B %*%C;
dim(D)
fastDecomp <- fastPCA(D, k=20, diffsnorm=TRUE)
str(fastDecomp)


###################
fastDecompBed10E3 <- fastPCA_BED("/data/GERA_DECRYPTED/LindermanAnalysis/EUR/Benchmark/eur_100000_62318/eur_100000_62318",k = 20, mem=8e+10, l=22,centeringRow=TRUE, phenoFile="/data/GERA_DECRYPTED/38852/PhenoGenotypeFiles/RootStudyConsentSet_phs000674.GERA.v1.p1.c1.HROO/PhenotypeFiles/phs000674.v1.pht003641.v1.p1.c1.Survey.HROO.txt")
str(fastDecompBed10E3)

#plot(fastDecompBed10E3$V[,1], fastDecompBed10E3$V[,2], pch='.', col=factor(fastDecompBed10E3$phenos$RACE))

fastDecompBed <- fastPCA_BED("/data/GERA_DECRYPTED/LindermanAnalysis/EUR/Benchmark/eur_20000_62318/eur_20000_62318",k = 20, mem=1e+9, l=22, phenoFile="/data/GERA_DECRYPTED/38852/PhenoGenotypeFiles/RootStudyConsentSet_phs000674.GERA.v1.p1.c1.HROO/PhenotypeFiles/phs000674.v1.pht003641.v1.p1.c1.Survey.HROO.txt",centeringRow = TRUE)

#plot(fastDecompBed$V[,1], fastDecompBed$V[,2], pch='.', col=factor(fastDecompBed$phenos$GENERALHEALTHCAT))


######CSV

fastDecompDropSeqCSV <- fastPCA_CSV('/data/Linderman/DropSeq/GSE63472_P14Retina_merged_digital_expression_headerless.csv',k=5, l=7, mem=50E9 )
