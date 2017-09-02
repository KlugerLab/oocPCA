make -C /opt/intel/mkl/tools/builder/ libuni export=`pwd`/mkl_fxn_list name=`pwd`/libfastpca_custommkl
install_name_tool -id "@rpath/libfastpca_custommkl.dylib" libfastpca_custommkl.dylib
#install_name_tool -id "@rpath/libiomp5.dylib" libiomp5.dylib

#install_name_tool -change libiomp5.dylib @rpath/libiomp5.dylib libfastpca_custommkl.dylib
#Now, need to change the rpath of libfast_custommkl.dylib to also be @executablelocal../
