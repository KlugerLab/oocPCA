cmake_minimum_required(VERSION 2.8.9)
#SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#SET(BUILD_SHARED_LIBRARIES OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")

set(Intel_Threading_Library           "/opt/intel/mkl/lib/libmkl_intel_thread.dylib")
set(IOMP5           "/opt/intel/compilers_and_libraries_2016.1.111/mac/compiler/lib/libiomp5.dylib")
set(ARMA_USE_LAPACK           false)
set(ARMA_USE_BLAS             false)
set(ARMA_USE_ATLAS            false)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake_aux/Modules/")
project(FastPCA)
add_definitions(-DDEBUG)



#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp")


#Can manually add the sources using the set command as follows:
#TODO: If BUIDLR
#set(SOURCES ${SOURCES} Rsrc/fastRPCA.cpp)


include(CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
	CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

##
## Find LAPACK and BLAS libraries, or their optimised versions
##

if(APPLE)
  set(ARMA_OS macos)
  include(ARMA_FindMKL)
  message(STATUS "     MKL_FOUND = ${MKL_FOUND}"     )
    add_definitions(-D__APPLE__)
  if(MKL_FOUND )
    
    set(ARMA_LAPACKE_MKL_INCLUDE_DIR ${LAPACKE_MKL_INCLUDE_DIR})
    message(STATUS "")
    message(STATUS "ARMA_LAPACKE_MKL_INCLUDE_DIR = ${ARMA_LAPACKE_MKL_INCLUDE_DIR}")
    
    set(ARMA_LIBS ${ARMA_LIBS} ${MKL_LIBRARIES})
    add_definitions(-DMKL_CUSTOM_LIBRARY="${MKL_LIBRARIES}")
    add_definitions(-DUSE_MKL)


    #add_definitions(-DMKL_CUSTOM_LIBRARY=fofo)
  else()
  
	  set(ARMA_LIBS ${ARMA_LIBS} "-framework Accelerate")  # or "-framework accelerate" ?
	  add_definitions(-D__APPLE__)
	  message(STATUS "MacOS X detected with no MKL. Added '-framework accelerate' to compiler flags")
  
	  add_definitions(-DMKL_CUSTOM_LIBRARY="fofo")
  endif()
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    message(STATUS "Clang compiler on MacOS X detected. Added '-stdlib=libc++' to compiler flags")
  endif()
  
else()
  
  set(ARMA_OS unix)
  
  include(ARMA_FindMKL)
  include(ARMA_FindOpenBLAS)
  include(ARMA_FindATLAS)
  include(ARMA_FindBLAS)
  include(ARMA_FindLAPACK)
  
  #set(MKL_FOUND           false)
  message(STATUS "     MKL_FOUND = ${MKL_FOUND}"     )
  message(STATUS "OpenBLAS_FOUND = ${OpenBLAS_FOUND}")
  message(STATUS "   ATLAS_FOUND = ${ATLAS_FOUND}"   )
  message(STATUS "    BLAS_FOUND = ${BLAS_FOUND}"    )
  message(STATUS "  LAPACK_FOUND = ${LAPACK_FOUND}"  )
  
  if(MKL_FOUND )
    
    set(ARMA_USE_LAPACK true)
    set(ARMA_USE_BLAS   true)
    
    message(STATUS "")
    message(STATUS "*** If the MKL or ACML libraries are installed in non-standard locations such as")
    message(STATUS "*** /opt/intel/mkl, /opt/intel/composerxe/, /usr/local/intel/mkl")
    message(STATUS "*** make sure the run-time linker can find them.")
    message(STATUS "*** On Linux systems this can be done by editing /etc/ld.so.conf")
    message(STATUS "*** or modifying the LD_LIBRARY_PATH environment variable.")
    message(STATUS "")
    message(STATUS "*** On systems with SELinux enabled (eg. Fedora, RHEL),")
    message(STATUS "*** you may need to change the SELinux type of all MKL/ACML libraries")
    message(STATUS "*** to fix permission problems that may occur during run-time.")
    message(STATUS "*** See README.txt for more information")
    message(STATUS "")
    
#    set(ARMA_LIBS ${ARMA_LIBS} ${MKL_LIBRARIES})
      message(STATUS "NOT LINKING!!")

      #list(APPEND CMAKE_PREFIX_PATH "./")

      #make libintel64 export=mkl_fxn_list name=/data/Linderman/FastPCA4/lib/libfastpca_custommkl
    set(ARMA_LIBS ${ARMA_LIBS} /data/Linderman/FastPCA4/lib/libfastpca_custommkl.so)
    #set(ARMA_LIBS ${ARMA_LIBS} libcustom_mkldist.so)
    #set(ARMA_LIBS ${ARMA_LIBS})
    add_definitions(-DMKL_CUSTOM_LIBRARY="libfastpca_custommkl.so")
    #add_definitions(-DMKL_CUSTOM_LIBRARY=fofo)
    
  else()
    
    if(OpenBLAS_FOUND AND ATLAS_FOUND)
      message(STATUS "")
      message(STATUS "*** WARNING: found both OpenBLAS and ATLAS; ATLAS will not be used")
    endif()
    
    if(OpenBLAS_FOUND AND BLAS_FOUND)
      message(STATUS "")
      message(STATUS "*** WARNING: found both OpenBLAS and BLAS; BLAS will not be used")
    endif()
    
    if(OpenBLAS_FOUND)
      
      set(ARMA_USE_BLAS true)
      set(ARMA_LIBS ${ARMA_LIBS} ${OpenBLAS_LIBRARIES})
      add_definitions(-DMKL_CUSTOM_LIBRARY="${OpenBLAS_LIBRARIES}")
      
      message(STATUS "")
      message(STATUS "*** If the OpenBLAS library is installed in")
      message(STATUS "*** /usr/local/lib or /usr/local/lib64")
      message(STATUS "*** make sure the run-time linker can find it.")
      message(STATUS "*** On Linux systems this can be done by editing /etc/ld.so.conf")
      message(STATUS "*** or modifying the LD_LIBRARY_PATH environment variable.")
      message(STATUS "")
      
    else()
      
      if(ATLAS_FOUND)
        set(ARMA_USE_ATLAS true)
        set(ARMA_ATLAS_INCLUDE_DIR ${ATLAS_INCLUDE_DIR})
        set(ARMA_LIBS ${ARMA_LIBS} ${ATLAS_LIBRARIES})
        add_definitions(-DMKL_CUSTOM_LIBRARY="${ATLAS_LIBRARIES}")
        
        message(STATUS "ATLAS_INCLUDE_DIR = ${ATLAS_INCLUDE_DIR}")
      endif()
      
      if(BLAS_FOUND)
        set(ARMA_USE_BLAS true)
        set(ARMA_LIBS ${ARMA_LIBS} ${BLAS_LIBRARIES})
        add_definitions(-DMKL_CUSTOM_LIBRARY="${BLAS_LIBRARIES}")
      endif()
      
    endif()
    
    if(LAPACK_FOUND)
      set(ARMA_USE_LAPACK true)
      set(ARMA_LIBS ${ARMA_LIBS} ${LAPACK_LIBRARIES})
      message(STATUS "Foudn LAPACK, ARMA_LIBS = ${ARMA_LIBS}")
    endif()
      
  endif()
  
endif()



#Bring the headers, such as Student.h into the project
include_directories(include ${ARMA_LAPACKE_MKL_INCLUDE_DIR})

add_executable(fastpca.xx ${SOURCES})
add_library(fastpca SHARED ${SOURCES})
#set(ARMA_LIBS ${ARMA_LIBS} ${IOMP5})
#set(ARMA_LIBS ${ARMA_LIBS} ${Intel_Threading_Library})

#find_package(Threads REQUIRED)
#if(THREADS_HAVE_PTHREAD_ARG)
#set_property(TARGET fastpca.xx PROPERTY COMPILE_OPTIONS "-pthread")
#set_property(TARGET fastpca.xx PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
#endif()
#if(CMAKE_THREAD_LIBS_INIT)
#target_link_libraries(fastpca.xx "${CMAKE_THREAD_LIBS_INIT}")
#endif()

#set_property(TARGET fastpca.xx PROPERTY COMPILE_OPTIONS "-Wl,-rpath,/opt/intel/mkl/lib")
#set_property(TARGET fastpca.xx PROPERTY INTERFACE_COMPILE_OPTIONS "-Wl,-rpath,/opt/intel/mkl/lib")

message(STATUS "ARMA_LIBS = ${ARMA_LIBS}")
target_link_libraries( fastpca.xx ${ARMA_LIBS}  ${CMAKE_DL_LIBS})
target_link_libraries( fastpca ${ARMA_LIBS}  ${CMAKE_DL_LIBS})
message(STATUS "Test  ${ARMA_LIBS}")
message(STATUS "Test  ${CMAKE_DL_LIBS}")


#install(FILES lib/libcustommkldist_local.so lib/libiomp5.so DESTINATION /lib/libfastpca.so)

