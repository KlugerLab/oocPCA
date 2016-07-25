# FastPCA

##Installation
###OSX
####Prerequisities:
* Cmake: Open source tool for building software cross-platforms (Download here: https://cmake.org/download/)
* Intel Math Kernel Library: Highly optimized implementations of BLAS and LAPACK (Free download here: https://software.intel.com/sites/campaigns/nest/)

####Installation process
1. Open Cmake, and set the "Source Directory" to be the root folder of the project, and the the build directory to be a directory of your choosing (you can just let it be a folder called "build" in the FastPCA folder).
2. Click "Configure," and Cmake will find the location of your Intel MKL libraries and include files.  If it does not find them, you can specify them yourself
3. Click "Generate," which will generate a Makefile
4. Use a terminal to cd into the build folder, and run the command "make" which will build the software in that folder.
5. Export the DYLD_LIBRARY_PATH so that the new executable can find the necessary dynamic libraries.  The folder containing libmkl_rt.dylib and the folder containing libiomp5.dylib must be in the path.  For example, : export DYLD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib/:/opt/intel/mkl/lib"

####Note
The code base uses the LAPACKE interface for LAPACK as opposed to the f2c
generated LAPACK functions.  The OS X Accelerate Framework does not include the
former, and hence, we require users to install Intel MKL to use FastPCA on
their Macs.  

##TODO
* Transition from LAPACKE to the CLAPACK f2c'd functions, allowing linkage with OS X Accelerate Framework
