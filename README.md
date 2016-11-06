# FastPCA
FastPCA is a Intel MKL-based, out-of-core C++ implementation of randomized
SVD for rank-k approximation of matrices that are 1) too large to fit into
memory and/or 2) would take too long to compute using traditional methods.

Two interfaces are available, via an R wrapper called fastRPCA and by
commandline.   OS X and Linux users can install the pre-compiled binaries,
following the processes outlined below:

##R Package Installation
1. Install devtools: `install.packages('devtools')`
2. Install fastRPCA: `devtools::install_github("KlugerLab/FastPCA",subdir="fastRPCA",
   host="git.yale.edu/api/v3", auth_token="<>")`

OR:

1.  Clone this git repository
2. `cd fastRPCA`
3. `R CMD INSTALL .`

Please see the documentation for usage: `?fastPCA ?fastPCA_CSV  ?fastPCA_BED `

###R Testing Suite
Test cases for this software use the popular testing package `testthat`:

1. `install.packages('testthat')`
2. `testthat::test_dir(sprintf("%s/testthat", system.file("tests", package="fastRPCA")))`




##Installing Command-line Implementation
1. Clone this git repository
2. Export the `LD_LIBRARY_PATH` or (on OS X)  `DYLD_LIBRARY_PATH` so that the new executable can find the
   necessary dynamic libraries in the `lib` folder.  ` export
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/user/Downloads/fastPCA/lib"`. 

##Features
* Variety of input formats and use cases
  * From memory in R
  * CSV format
  * plink's binary [BED format] (http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml) for GWAS data. 
* All matrix algebra is done with Intel MKL (pre-compiled version already linked) making it extremely fast
* The calculations are 'blocked' allowing it to be 'out-of-core' when necessary, so that the user to specify the maximum amount of memory to be used.
  * CSV files: when too large for the memory, read block by block from the hard drive
  * BED files: when too large for the memory, stored in a compressed 2 bit-per-SNP format, and then decompressed block by block for calculations
* Row-centering and column-centering
* Imputation by averaging of missing data for GWAS


##Development
###Compiling from source
####OSX
#####Prerequisities:
* Intel Math Kernel Library: Highly optimized implementations of BLAS and LAPACK (Free download [here](https://software.intel.com/sites/campaigns/nest/) ).  The `lib` folder contains a [custom built shared library](https://software.intel.com/en-us/node/528690), but the headers cannot be distributed.  As such, to compile from source, Intel MKL must be installed on your machine.

#####How to compile from source
1. Use a terminal to cd into the src folder, and run the command `make` which will build the software in that folder.
2. Export the `LD_LIBRARY_PATH` or (on OS X)  `DYLD_LIBRARY_PATH` so that the new executable can find the
   necessary dynamic libraries in the `lib` folder. For example:
  ` export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/user/Downloads/fastPCA/lib"`. 

####Note
The code base uses the LAPACKE interface for LAPACK as opposed to the f2c
generated LAPACK functions.  The OS X Accelerate Framework does not include the
former, and hence, we use Intel MKL in lieu of OS X Accelerate Framework. 

##TODO
* Transition from LAPACKE to the CLAPACK f2c'd functions, allowing linkage with OS X Accelerate Framework
