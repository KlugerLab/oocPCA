# FastPCA
oocPCA is a Intel MKL-based, out-of-core C++ implementation of randomized
SVD for rank-k approximation of matrices that are too large to fit into
memory.

Two interfaces are available, via an R wrapper called oocRPCA and by
commandline.   OS X and Linux users can install the pre-compiled binaries,
following the processes outlined below:

## R Package Installation
```R
if(!require(devtools)) install.packages("devtools") # If not already installed
devtools::install_github("linqiaozhi/fmmRtsne", auth_token = "df6fb169997c1c59d34c6dc7254657cdf54ae8f1")
```

OR:

1.  Clone this git repository
2. `cd fastRPCA`
3. `R CMD INSTALL .`

Please see the documentation for usage: `?oocPCA_CSV`

## R Testing Suite
Test cases for this software use the popular testing package `testthat`:

```R
if(!require(testthat)) install.packages('testthat')
testthat::test_dir(sprintf("%s/testthat", system.file("tests", package="oocRPCA")))
```

## Features
* Variety of input formats and use cases
  * CSV format
  * plink's binary [BED format] (http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml) for GWAS data. 
* All matrix algebra is done with Intel MKL (pre-compiled version already linked) making it extremely fast
* The calculations are 'blocked' allowing it to be 'out-of-core' when necessary, so that the user to specify the maximum amount of memory to be used.
  * CSV files: when too large for the memory, read block by block from the hard drive
  * BED files: when too large for the memory, stored in a compressed 2 bit-per-SNP format, and then decompressed block by block for calculations
* Row-centering and column-centering
* Imputation by averaging of missing data for GWAS


## Development
### Compiling from source
This implementation relies heavily on a highly optimized implementations of BLAS and LAPACK called Intel Math Kernel Library (MKL) (Free download [here](https://software.intel.com/sites/campaigns/nest/) ).  The `lib` folder contains a [custom built shared library](https://software.intel.com/en-us/node/528690), but the headers cannot be distributed.  As such, to compile from source, Intel MKL must be installed on your machine.  To recompile the custom built shared library, follow the instructions in lib/generate_custom_mkl.sh

#### Note
The code base uses the LAPACKE interface for LAPACK as opposed to the f2c generated LAPACK functions.  The OS X Accelerate Framework does not include the former, and hence, we use Intel MKL in lieu of OS X Accelerate Framework. 

## TODO
* Windows Support

