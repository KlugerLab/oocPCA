# FastPCA
FastPCA is a C++ implementation of randomized SVD. Two interfaces are
available, via an R wrapper called fastRPCA and by commandline.   OS X and
Linux users can install the pre-compiled binaries, following the processes
outlined below:

##R Package Installation
1. Install devtools: `install.packages('devtools')`
2. Install fastRPCA: `install_github("KlugerLab/FastPCA",subdir="fastRPCA",
   host="git.yale.edu/api/v3", auth_token="<>")`

OR:

1.  Clone this git repository
2. `cd fastRPCA`
3. `R CMD INSTALL .`



###R Testing
Run the testing suite using the following command: `devtools::test('fastRPCA')`

##Installing Command-line Implementation
1. Clone this git repository
2. Export the `LD_LIBRARY_PATH` or (on OS X)  `DYLD_LIBRARY_PATH` so that the new executable can find the
   necessary dynamic libraries in the `lib` folder.  ` export
LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/user/Downloads/fastPCA/lib"`. 


##Compiling from source
###OSX
####Prerequisities:
* Intel Math Kernel Library: Highly optimized implementations of BLAS and LAPACK (Free download [here](https://software.intel.com/sites/campaigns/nest/) ).  The `lib` folder contains a [custom built shared library](https://software.intel.com/en-us/node/528690), but the headers cannot be distributed.  As such, to compile from source, Intel MKL must be installed on your machine.

####How to compile from source
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
