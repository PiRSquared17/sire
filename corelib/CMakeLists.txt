################################
#
#  CMake Project file for Sire
#
#   (C) Christopher Woods
#
################################

############
############ SETTING INITIAL SIRE VARIABLES
############

# require cmake >= 2.8.12
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
cmake_policy(VERSION 2.8)

# use loose loop syntax in CMakeLists files
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

# get the number of cores on this machine - useful for
# compiling dependencies
include(ProcessorCount)
ProcessorCount(NCORES)
if (NCORES EQUAL 0)
    set( NCORES 1 )
endif()
message(STATUS "Configuring on a machine with number of cores = ${NCORES}")

# Option to set the install directory of Sire
set( SIRE_INSTALL_PREFIX "$ENV{HOME}/sire.app" CACHE PATH
        "Install directory for Sire. This directory can be moved after compilation" )

set(CMAKE_INSTALL_PREFIX "${SIRE_INSTALL_PREFIX}" CACHE INTERNAL
         "Prefix prepended to install directories" FORCE ) 

message( STATUS "Sire will be compiled and installed to directory ${CMAKE_INSTALL_PREFIX}" )

# Absolutely can't run cmake in the source directory!
if ( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
  message( FATAL_ERROR "You must run CMake in a different directory to the source!" ) 
endif()

# name the project
project (SIRE CXX C)

# Whether or not to build the parts of Sire that depend on
# a graphical environment (e.g. inquire)
option ( SIRE_BUILD_GUI "Build graphical applications" OFF )

# Create a file in which we can save the values of all useful variables.
# This will mean that projects using Sire won't have to set these variables
set (SIRE_VARIABLES_FILE "${CMAKE_BINARY_DIR}/SireCompileVariables.cmake")
file ( REMOVE ${SIRE_VARIABLES_FILE} )

# This macro is used to save the names of variable to SireCompileVariables.cmake
macro( SAVE_SIRE_VARIABLE _var _value )
  file(APPEND ${SIRE_VARIABLES_FILE}
        "set( ${_var} \"${_value}\" )\n")
endmacro (SAVE_SIRE_VARIABLE)

# project version
set (SIRE_VERSION_MAJOR "0")
set (SIRE_VERSION_MINOR "0")
set (SIRE_VERSION_PATCH "1")

set (SIRE_VERSION "${SIRE_VERSION_MAJOR}.${SIRE_VERSION_MINOR}.${SIRE_VERSION_PATCH}")

math (EXPR SIRE_VERSION_NUMBER 
      "${SIRE_VERSION_MAJOR}*100000 + ${SIRE_VERSION_MINOR}*100 + ${SIRE_VERSION_PATCH}")

set (SIRE_VERSION_STRING "${SIRE_VERSION_MAJOR}_${SIRE_VERSION_MINOR}_${SIRE_VERSION_PATCH}")

save_sire_variable( "SIRE_VERSION_MAJOR" "${SIRE_VERSION_MAJOR}" )
save_sire_variable( "SIRE_VERSION_MINOR" "${SIRE_VERSION_MINOR}" )
save_sire_variable( "SIRE_VERSION_PATCH" "${SIRE_VERSION_PATCH}" )
save_sire_variable( "SIRE_VERSION" "${SIRE_VERSION}" )
save_sire_variable( "SIRE_VERSION_STRING" "${SIRE_VERSION_STRING}" )
save_sire_variable( "SIRE_APP_DIR" "${CMAKE_INSTALL_PREFIX}" )

# Ensure that the Sire's directories are searched for header files
# before the system directories - this prevents problems in case
# an older version of Sire is already installed on the system
INCLUDE_DIRECTORIES( BEFORE ${SIRE_SOURCE_DIR} )

# Sire must be built as a set of shared libraries
# (otherwise the modules cannot be loaded in python)
set (BUILD_SHARED_LIBS ON)

# This is where you set your own compiler flags
set (SIRE_USER_CFLAGS "" CACHE STRING "User supplied C compiler flags" )
set (SIRE_USER_CPPFLAGS "" CACHE STRING "User supplied C++ compiler flags" )

save_sire_variable( "SIRE_USER_CFLAGS" "${SIRE_USER_CFLAGS}" )
save_sire_variable( "SIRE_USER_CPPFLAGS" "${SIRE_USER_CPPFLAGS}" )

# Name of executable install dir
set (SIRE_BIN "bin")
# Shared library install dir
set (SIRE_LIBS "lib")
# Where to install all bundled components
set (SIRE_BUNDLED_DIR "bundled" )
set (SIRE_BUNDLED_LIBS "${SIRE_BUNDLED_DIR}/lib")
# Static library install dir
set (SIRE_ARCHIVES "lib/static")
# Include (header) files install dir
set (SIRE_INCLUDES "include/Sire")
# CMake install files install dir
set (SIRE_CMAKEFILES "include/Sire/cmake")
# Sire share directory (for parameters etc.)
set (SIRE_SHARE "share/Sire")
# Where will the test files be saved?
set (SIRE_TEST "test")

# Add these directories to the RPATH variable for the libraries / executables,
# so that they can find each other when loading
set( CMAKE_SKIP_RPATH FALSE )
set (CMAKE_SKIP_BUILD_RPATH FALSE)

if (APPLE)
  # make sure that all libraries are named "@rpath/library.dylib". This will allow
  # the library to be found as long as the executable can substitute @rpath with the 
  # path to the library
  set (CMAKE_INSTALL_NAME_DIR "@rpath")

  # Add the two passed directories to the @rpaths to search for libraries. @executable_path
  # is the directory containing the executable, so the libraries are included in 
  # @exeutable_path/../${SIRE_LIBS} and @exeutable_path/../${SIRE_BUNDLED_LIBS}
  set (SIRE_INSTALL_RPATH "@executable_path/.:@executable_path/../${SIRE_LIBS}:@executable_path/../${SIRE_BUNDLED_LIBS}")
else()
  set (SIRE_INSTALL_RPATH "\$ORIGIN/.:\$ORIGIN/../${SIRE_LIBS}:\$ORIGIN/../${SIRE_BUNDLED_LIBS}")
endif()

set (CMAKE_INSTALL_RPATH "${SIRE_INSTALL_RPATH}")
set (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

save_sire_variable( "SIRE_INSTALL_PREFIX" "${CMAKE_INSTALL_PREFIX}" )
save_sire_variable( "SIRE_BIN" "${SIRE_BIN}" )
save_sire_variable( "SIRE_LIBS" "${SIRE_LIBS}" )
save_sire_variable( "SIRE_BUNDLED_DIR" "${SIRE_BUNDLED_DIR}" )
save_sire_variable( "SIRE_BUNDLED_LIBS" "${SIRE_BUNDLED_LIBS}" )
save_sire_variable( "SIRE_ARCHIVES" "${SIRE_ARCHIVES}" )
save_sire_variable( "SIRE_INCLUDES" "${SIRE_INCLUDES}" )
save_sire_variable( "SIRE_SHARE" "${SIRE_SHARE}" )
save_sire_variable( "SIRE_CMAKEFILES" "${SIRE_CMAKEFILES}" )
save_sire_variable( "SIRE_INSTALL_RPATH" "${SIRE_INSTALL_RPATH}" )
save_sire_variable( "SIRE_TEST" "${SIRE_TEST}" )

# Uncomment for more verbose compiling/linking
# (or run 'make VERBOSE=1')
# set (CMAKE_VERBOSE_MAKEFILE ON)

# Path to the project's extra cmake files
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/build/cmake ${CMAKE_SOURCE_DIR}/build/cmake/cmake_cxx11)

# Options to pass to 'strip' when stripping the libraries
# (on linux, aix and mac need to use '-x' to ensure we leave
#  the global symbols)
set ( SIRE_STRIP_OPTIONS "-x" CACHE STRING "Options to pass to 'strip'" )
save_sire_variable( "SIRE_STRIP_OPTIONS" "${SIRE_STRIP_OPTIONS}" )

set ( SIRE_STRIP_COMMAND "${CMAKE_STRIP} ${SIRE_STRIP_OPTIONS}" )

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_strip.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip"
  IMMEDIATE @ONLY)

set ( CMAKE_STRIP "${CMAKE_CURRENT_BINARY_DIR}/cmake_strip" )

# Add option to turn on or off symbol hiding
option ( SIRE_SYMBOL_HIDING "Turn on library symbol hiding" ON )
save_sire_variable( "SIRE_SYMBOL_HIDING" "${SIRE_SYMBOL_HIDING}" )

# Add option to turn on or off vectorisation of the code
option ( SIRE_VECTORISE "Turn on code vectorisation" ON )
save_sire_variable( "SIRE_VECTORISE" "${SIRE_VECTORISE}" )

# Add option to disable use of SSE
option ( SIRE_DISABLE_SSE "Prevent use of SSE when compiling" OFF )

# Add option to disable use of AVX
option ( SIRE_DISABLE_AVX "Prevent use of AVX when compiling" ON )

# Add option to enable use of SSE4 (only works if SSE is not disabled)
# This is on by default as most processors now support SSE4
option ( SIRE_ENABLE_SSE4 "Enable use of SSE4 when compiling" OFF )

# Add option to enable use of AVX2 (only works if AVX is not disabled)
# This is off by default as AVX2 is not widespread
option ( SIRE_ENABLE_AVX2 "Enable use of AVX2 when compiling" OFF )

# Add option to enable use of AVX512 (only works if AVX is not disabled)
# This is off by default as AVX512 is not widespread
option ( SIRE_ENABLE_AVX512 "Enabled use of AVX512 when compiling" OFF )

# Add an option to force off Fortran compilation and linking
option ( SIRE_DISABLE_FORTRAN "Turn off Fortran compilation and linking" ON )
save_sire_variable( "SIRE_DISABLE_FORTRAN" "${SIRE_DISABLE_FORTRAN}" )

if (SIRE_DISABLE_FORTRAN)
  message( STATUS "Disabling any code that requires compilation or linking to Fortran" )
endif()

# Add an option to force the compilation to ignore any MPI libraries
# (this is useful on clusters where the MPI-enabled build cannot run
#  on the head node, and yet you still want to run tests!)
option( SIRE_DISABLE_MPI "Turn off MPI compilation and linking" ON )

#############
############# PLATFORM TESTS
#############

# Now run some platform tests - the results will be placed in config.h
include (CheckIncludeFiles)
include (CheckSymbolExists)
include (CheckFunctionExists)

check_include_files( unistd.h HAVE_UNISTD_H )
check_function_exists( lseek64 HAVE_LSEEK64 )

check_include_files( malloc/malloc.h HAVE_MALLOC_MALLOC_H )
check_include_files( emmintrin.h HAVE_EMMINTRIN_H )
check_include_files( smmintrin.h HAVE_SMMINTRIN_H )
check_include_files( immintrin.h HAVE_IMMINTRIN_H )
check_include_files( execinfo.h HAVE_EXECINFO_H )
check_include_files( stdint.h HAVE_STDINT_H )
check_function_exists( mstats HAVE_MSTATS )
check_function_exists( mallinfo HAVE_MALLINFO )
check_function_exists( sysctl HAVE_SYSCTL )
check_function_exists( cbrt HAVE_CUBEROOT )
check_symbol_exists( backtrace "execinfo.h" HAVE_BACKTRACE )

#create config.h
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/sire_config.h.in 
               ${CMAKE_CURRENT_BINARY_DIR}/sire_config.h)

#  We need to get the compiler flags for different options
set( SIRE_SMALL_FLAGS "" CACHE INTERNAL "Compiler flags to produce small code" )
set( SIRE_WARNALL_FLAGS "" CACHE INTERNAL "Compiler flags to turn on all warnings" )
set( SIRE_DEBUG_FLAGS "" CACHE INTERNAL "Compiler flags to produce debugging code" )
set( SIRE_RELEASE_FLAGS "" CACHE INTERNAL "Compiler flags to produce optimised release code" )
set( SIRE_OPENMP_FLAGS "" CACHE INTERNAL "Compiler flags to activate OpenMP" )
set( SIRE_VECTOR_FLAGS "" CACHE INTERNAL "Compiler flags to vectorise code (e.g. SSE)" )
set( SIRE_VISIBILITY_FLAGS "" CACHE INTERNAL "Compiler flags to activate symbol visibility" )
set( SIRE_PLATFORM_FLAGS "" CACHE INTERNAL "Additional flags needed by the compiler" )
set( SIRE_SHARE_LINK_FLAGS "" CACHE INTERNAL "Flags needed to link shared libraries" )
set( SIRE_STATIC_LINK_FLAGS "" CACHE INTERNAL "Flags needed to link static libraries" )
set( SIRE_EXE_LINK_FLAGS "" CACHE INTERNAL "Flags needed to link executables" )

#############
############# DETECTING AND CONFIGURING THE COMPILER
#############

include (CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG( "-stdlib=libc++" SIRE_HAS_CPP_LIB )
CHECK_CXX_COMPILER_FLAG( "-std=c++11" SIRE_HAS_CPP_11 )

save_sire_variable( "SIRE_HAS_CPP_LIB" "${SIRE_HAS_CPP_LIB}" )
save_sire_variable( "SIRE_HAS_CPP_11" "${SIRE_HAS_CPP_11}" )

if (SIRE_HAS_CPP_LIB)
  message(STATUS "Using libc++ in place of libstdc++")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  set(SIRE_CXX11_FLAGS "${SIRE_CXX11_FLAGS} -stdlib=libc++")
endif()

set(NEED_UNDEF_STRICT_ANSI FALSE)

if (SIRE_HAS_CPP_11)
  message(STATUS "Compiler supports C++ 2011")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(SIRE_CXX11_FLAGS "${SIRE_CXX11_FLAGS} -std=c++11")

  if (NEED_UNDEF_STRICT_ANSI)
    message( STATUS "Undefining __STRICT_ANSI__ because of compile problems")
    set(SIRE_CXX11_FLAGS "${SIRE_CXX11_FLAGS} -U__STRICT_ANSI__")
  endif()
endif()

save_sire_variable( "SIRE_CXX11_FLAGS" "${SIRE_CXX11_FLAGS}" )

# Find which C++ 2011 features are supported
include (CheckCXX11Features)

message( STATUS "C++11 supported features: ${CXX11_FEATURE_LIST}" )

set( SIRE_COMPILER "UNKNOWN" )

if (NOT ${CMAKE_C_COMPILER_ID} STREQUAL ${CMAKE_CXX_COMPILER_ID})
  message( WARNING "It is not sensible to compile Sire using different C and C++ compilers. "
                   "You are using ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID}) and "
                   "${CMAKE_CXX_COMPILER} (${CMAKE_CXX_COMPILER_ID})")
endif()

string( TOUPPER ${CMAKE_CXX_COMPILER_ID} SIRE_CXX_COMPILER_ID )

# Add compiler specific flags
if ( ${SIRE_CXX_COMPILER_ID} MATCHES "GNU" )

    message( STATUS "Detected and using GCC C and C++ compilers" )
    set( SIRE_COMPILER "GCC" )

    # get the compiler version (via gcc -dumpversion)
    # (I'll assume that we are using the same version
    # of compiler for C and C++...)
    EXEC_PROGRAM( ${CMAKE_CXX_COMPILER}
                  ARGS -dumpversion
                  OUTPUT_VARIABLE exe_output )

    # These two regexps are copied from the Chicken CMake
    # files (which provide a very useful cmake tutorial - thanks :-)                  
    SET( VERSION_BUILD_REGEX "^([0-9]+)[.]([0-9]+).*$" )
    STRING( REGEX REPLACE ${VERSION_BUILD_REGEX} "\\1"
            GCC_MAJOR_VERSION ${exe_output} )
    STRING( REGEX REPLACE ${VERSION_BUILD_REGEX} "\\2"
            GCC_MINOR_VERSION ${exe_output} )
        
    message( STATUS "Using GCC version "
                    "${GCC_MAJOR_VERSION}.${GCC_MINOR_VERSION}" )

    # First, we need at least GCC 3.0, as I don't think that
    # GCC 2.x can handle the template code that I use
    if ( GCC_MAJOR_VERSION LESS 3 )
       message( FATAL_ERROR "Sire requires GCC >= 3.0. Please upgrade "
                      "your version of GCC." )
    endif()

    set ( SIRE_SMALL_FLAGS  "-Os" )
    set ( SIRE_WARNALL_FLAGS  "-Wall" )
    set ( SIRE_DEBUG_FLAGS   "-g" )

    set ( SIRE_RELEASE_FLAGS "-O3 -ffast-math" )

    # -fomit-frame-pointer breaks backtrace on an apple
    # and can cause segfaults when exceptions are not caught
    # set ( SIRE_RELEASE_FLAGS "${SIRE_RELEASE_FLAGS} -fomit-frame-pointer" )

    set ( SIRE_OPENMP_FLAGS "-fopenmp" )
    set ( SIRE_PLATFORM_FLAGS "-pipe" )
    set ( SIRE_STATIC_LINK_FLAGS "-static" )

    CHECK_CXX_COMPILER_FLAG( "-msse2" HAVE_SSE_FLAG )
    CHECK_CXX_COMPILER_FLAG( "-mavx" HAVE_AVX_FLAG )

    if ( NOT SIRE_DISABLE_AVX)

      if ( HAVE_IMMINTRIN_H AND HAVE_AVX_FLAG )

        # we have immintrin.h so we have AVX
        message( STATUS "Compiling with AVX enabled. Resulting binary will only work on newer Intel CPUs" )
        
        CHECK_CXX_COMPILER_FLAG( "-mavx2" HAVE_AVX2_FLAG )

        if ( SIRE_ENABLE_AVX2 AND HAVE_AVX2_FLAG )
          message( STATUS "Compiling with AVX2 enabled. Resulting binary will only work on really new Intel CPUs" )
          set( SIRE_VECTOR_FLAGS "-mavx2 -DSIRE_USE_AVX -DSIRE_USE_AVX2" )
        else()
          set ( SIRE_VECTOR_FLAGS "-mavx -DSIRE_USE_AVX" )
        endif()
      elseif ( HAVE_EMMINTRIN_H AND HAVE_SSE_FLAG )
        # we have emmintrin.h so we have SSE
        message( STATUS "Compiling with SSE2 enabled. Resulting binary will not work on really old Intel CPUs" )

        CHECK_CXX_COMPILER_FLAG( "-msse4.1" HAVE_SSE4_FLAG )

        if ( SIRE_ENABLE_SSE4 AND HAVE_SMMINTRIN_H AND HAVE_SSE4_FLAG )
          message( STATUS "Compiling with SSE4 enabled. Resulting binary will only work on newer Intel CPUs" )
          set( SIRE_VECTOR_FLAGS "-msse4.1 -DSIRE_USE_SSE -DSIRE_USE_SSE4" )
        else()
          set ( SIRE_VECTOR_FLAGS "-msse2 -DSIRE_USE_SSE" )
        endif()
      endif()
    elseif (NOT SIRE_DISABLE_SSE)   
        # we have emmintrin.h so we have SSE
        message( STATUS "Compiling with SSE2 enabled. Resulting binary will not work on really old Intel CPUs" )

        CHECK_CXX_COMPILER_FLAG( "-msse4.1" HAVE_SSE4_FLAG )

        if ( SIRE_ENABLE_SSE4 AND HAVE_SMMINTRIN_H AND HAVE_SSE4_FLAG )
          message( STATUS "Compiling with SSE4 enabled. Resulting binary will only work on newer Intel CPUs" )
          set( SIRE_VECTOR_FLAGS "-msse4.1 -DSIRE_USE_SSE -DSIRE_USE_SSE4" )
        else()
          set ( SIRE_VECTOR_FLAGS "-msse2 -DSIRE_USE_SSE" )
        endif()  
    endif()

    # Now gcc 4 specific options
    if ( GCC_MAJOR_VERSION GREATER 3 )
      # remove warnings about strict aliasing
      set( SIRE_PLATFORM_FLAGS "${SIRE_PLATFORM_FLAGS} -Wno-strict-aliasing" )

      # gcc-4.2 and above has openMP - turn it on (except for Mac as libgomp is broken!)
      if (GCC_MINOR_VERSION GREATER 1)
          if (SIRE_OPENMP)
              if (APPLE)
                message( STATUS "OpenMP support is unavailable as OpenMP + pthreads "
                                "is broken in Apple gcc4" )
                set( SIRE_OPENMP_FLAGS "-fno-openmp" )
              endif()
          endif()
      endif()

      set ( SIRE_VISIBILITY_FLAGS "-DSIRE_VISIBILITY_AVAILABLE -fvisibility=hidden -fvisibility-inlines-hidden" )
  
    else()
      set ( SIRE_VISIBILITY_FLAGS "-DSIRE_NO_VISIBILITY_AVAILABLE" ) 
    endif()

    # Ensure that the libraries don't contain any undefined symbols
    if(APPLE)
        set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} -Wl" )
    else(UNIX)
        set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} -rdynamic -Wl,--no-undefined" )
    endif()

elseif ( ${SIRE_CXX_COMPILER_ID} MATCHES "CLANG" )

    message( STATUS "Detected and using clang C and C++ compilers" )
    set( SIRE_COMPILER "CLANG" )

    set ( SIRE_SMALL_FLAGS  "-Os" )
    set ( SIRE_WARNALL_FLAGS  "-Wall" )
    set ( SIRE_DEBUG_FLAGS   "-g" )

    set ( SIRE_RELEASE_FLAGS "-O3 -ffast-math" )

    # -fomit-frame-pointer breaks backtrace on an apple
    # set ( SIRE_RELEASE_FLAGS "${SIRE_RELEASE_FLAGS} -fomit-frame-pointer" )

    set ( SIRE_OPENMP_FLAGS "-fopenmp" )
    set ( SIRE_PLATFORM_FLAGS "-pipe" )
    set ( SIRE_STATIC_LINK_FLAGS "-static" )

    CHECK_CXX_COMPILER_FLAG( "-msse2" HAVE_SSE_FLAG )
    CHECK_CXX_COMPILER_FLAG( "-mavx" HAVE_AVX_FLAG )

    if ( NOT SIRE_DISABLE_AVX )
      if ( HAVE_IMMINTRIN_H AND HAVE_AVX_FLAG )
        # we have immintrin.h so we have AVX
        message( STATUS "Compiling with AVX enabled. Resulting binary will only work on newer Intel CPUs" )

        CHECK_CXX_COMPILER_FLAG( "-mavx2" HAVE_AVX2_FLAG )

        if ( SIRE_ENABLE_AVX2 AND HAVE_AVX2_FLAG )
          message( STATUS "Compiling with AVX2 enabled. Resulting binary will only work on really new Intel CPUs")
          set( SIRE_VECTOR_FLAGS "-mavx2 -DSIRE_USE_AVX -DSIRE_USE_AVX2" )
        else()
          set ( SIRE_VECTOR_FLAGS "-mavx -DSIRE_USE_AVX" )
        endif() 
     elseif ( HAVE_EMMINTRIN_H AND HAVE_SSE_FLAG )
        # we have emmintrin.h so we have SSE
        message( STATUS "Compiling with SSE2 enabled. Resulting binary will not work on really old Intel CPUs" )

        CHECK_CXX_COMPILER_FLAG( "-msse4.1" HAVE_SSE4_FLAG )
        
        if ( SIRE_ENABLE_SSE4 AND HAVE_SMMINTRIN_H AND HAVE_SSE4_FLAG )
          message( STATUS "Compiling with SSE4.1 enabled. Resulting binary will only work on newer Intel CPUs" )
          set( SIRE_VECTOR_FLAGS "-msse4.1 -DSIRE_USE_SSE -DSIRE_USE_SSE4" )
        else()
          set ( SIRE_VECTOR_FLAGS "-msse2 -DSIRE_USE_SSE" )
        endif() 
     endif()
    elseif( NOT SIRE_DISABLE_SSE )
        # we have emmintrin.h so we have SSE
        message( STATUS "Compiling with SSE2 enabled. Resulting binary will not work on really old Intel CPUs" )

        CHECK_CXX_COMPILER_FLAG( "-msse4.1" HAVE_SSE4_FLAG )

        if ( SIRE_ENABLE_SSE4 AND HAVE_SMMINTRIN_H AND HAVE_SSE4_FLAG )
          message( STATUS "Compiling with SSE4.1 enabled. Resulting binary will only work on newer Intel CPUs" )
          set( SIRE_VECTOR_FLAGS "-msse4.1 -DSIRE_USE_SSE -DSIRE_USE_SSE4" )
        else()
          set ( SIRE_VECTOR_FLAGS "-msse2 -DSIRE_USE_SSE" )
        endif()    
    endif()
       
    set ( SIRE_VISIBILITY_FLAGS "-DSIRE_VISIBILITY_AVAILABLE -fvisibility=hidden -fvisibility-inlines-hidden" )

    if (APPLE)
    elseif(UNIX)
        set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} -rdynamic -Wl,--no-undefined" )
    endif()

elseif ( ${CMAKE_CXX_COMPILER} MATCHES "xlC" )
   
    message( STATUS "Detected the xlc and xlC C and C++ compilers" ) 
    set( SIRE_COMPILER "XLC" )

    set ( SIRE_SMALL_FLAGS  "-qcompact" )
    set ( SIRE_WARNALL_FLAGS  "-qwarn64 -qflag=w:w" )
    set ( SIRE_RELEASE_FLAGS "-O3 -qstrict" )
    set ( SIRE_DEBUG_FLAGS   "-g" )
    set ( SIRE_VISIBILITY_FLAGS "-DSIRE_NO_VISIBILITY_AVAILABLE" )
    set ( SIRE_PLATFORM_FLAGS "-qrtti=all -qthreaded" )

    set ( SIRE_SHARE_LINK_FLAGS "-qmkshrobj" )

    # Turn off Fortran code because this produces weird link errors!
    message( STATUS "Disabling Fortran code as this leads to link errors!" )
    message( STATUS "This will disable sire_blas, sire_lapack and sire_linpack" )
    set ( SIRE_DISABLE_FORTRAN ON )

    # This is a hack to automatically turn on MPI support on HPCx
    if ( NOT SIRE_DISABLE_MPI )
      if ( NOT SIRE_USE_MPI )
        set( SIRE_USE_MPI ON )
        set( MPI_INCLUDE_PATH "/usr/lpp/ppe.poe/include;/opt/rsct/lapi/include" )
        set( MPI_LIBRARY "-L/usr/lpp/ppe.poe/lib -L/usr/lpp/ppe.poe/lib/ip -lmpi_r -lvtd_r -llapi_r -binitfini:poe_remote_main" )
        set( MPI_EXTRA_LIBRARY "-L/usr/lpp/ppe.poe/lib/threads -lpthreads" )
        set( MPI_DEFINITIONS "-D_MPI_CPP_BINDINGS" )
      endif()
    endif()

elseif( ${CMAKE_CXX_COMPILER} MATCHES "icpc" )

    message( STATUS "Detected the Intel C++ and C compilers" )
    set( SIRE_COMPILER ICPC )

    set ( SIRE_SMALL_FLAGS  "-Os" )
    set ( SIRE_WARNALL_FLAGS  "-w1" )
    set ( SIRE_RELEASE_FLAGS  "-O3 -inline-level=2" )
    set ( SIRE_VECTOR_FLAGS  "-DSIRE_USE_SSE -msse3" )
    set ( SIRE_DEBUG_FLAGS   "-g" )
    set ( SIRE_VISIBILITY_FLAGS "-DSIRE_NO_VISIBILITY_AVAILABLE" )    
    set ( SIRE_OPENMP_FLAGS "-openmp" )
    set ( SIRE_SHARE_LINK_FLAGS "-shared" )
    set ( SIRE_STATIC_LINK_FLAGS "-static" )

else()
    message( STATUS "CMAKE_SYSTEM_NAME == ${CMAKE_SYSTEM_NAME}" )
    message( STATUS "CMAKE_C_COMPILER == ${CMAKE_C_COMPILER}" )
    message( STATUS "CMAKE_CXX_COMPILER == ${CMAKE_CXX_COMPILER}" )
    message( STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}" )
    message( FATAL_ERROR "Could not detect compiler type - cannot set compiler options." )
endif()

save_sire_variable( "SIRE_COMPILER" "${SIRE_COMPILER}" )

################
################ SORTING OUT ALL OF SIRE'S DEPENDENCIES
################

# Tell CMake to run moc when necessary
set(CMAKE_AUTOMOC OFF)

# As moc files are generated in the binary dir, tell CMake
# to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Option to disable all bundles
option( SIRE_ENABLE_BUNDLING "Allow Sire to use bundled libraries" ON )

# Fine-grain control of bundling each dependency
option( SIRE_BUNDLE_TBB "Whether Sire should use its bundled copy of TBB" ON )
option( SIRE_BUNDLE_GSL "Whether Sire should use its bundled copy of GSL" ON )
option( SIRE_BUNDLE_CPUID "Whether Sire should use its bundled copy of libcpuid" ON )
option( SIRE_BUNDLE_BOOST "Whether Sire should use its bundled copy of Boost" ON )
option( SIRE_BUNDLE_BLASLAPACK "Whether Sire should use its bundled copy of BLAS and LAPACK" ON )
option( SIRE_BUNDLE_QT "Whether Sire should use its bundled copy of Qt" ON )

save_sire_variable( "SIRE_ENABLE_BUNDLING" "${SIRE_ENABLE_BUNDLING}" )
save_sire_variable( "SIRE_BUNDLE_TBB" "${SIRE_BUNDLE_TBB}" )
save_sire_variable( "SIRE_BUNDLE_CPUID" "${SIRE_BUNDLE_CPUID}" )
save_sire_variable( "SIRE_BUNDLE_GSL" "${SIRE_BUNDLE_GSL}" )
save_sire_variable( "SIRE_BUNDLE_BOOST" "${SIRE_BUNDLE_BOOST}" )
save_sire_variable( "SIRE_BUNDLE_BLASLAPACK" "${SIRE_BUNDLE_BLASLAPACK}" )
save_sire_variable( "SIRE_BUNDLE_QT" "${SIRE_BUNDLE_QT}" )

set (SIRE_FOUND_TBB FALSE)
set (SIRE_FOUND_GSL FALSE)
set (SIRE_FOUND_BOOST FALSE)
set (SIRE_FOUND_BLASLAPACK FALSE)
set (SIRE_FOUND_QT FALSE)

if ( ${SIRE_ENABLE_BUNDLING} )
    message( STATUS "Sire will be compiled to use its own copy of bundled libraries." )

    set( BUNDLE_BUILDDIR "${CMAKE_BINARY_DIR}/build_bundled" )
    set( BUNDLE_STAGEDIR "${CMAKE_INSTALL_PREFIX}/${SIRE_BUNDLED_DIR}" )

    if (NOT EXISTS "${BUNDLE_BUILDDIR}")
      file(MAKE_DIRECTORY ${BUNDLE_BUILDDIR})
    endif()

    if (NOT EXISTS "${BUNDLE_STAGEDIR}")
      file(MAKE_DIRECTORY ${BUNDLE_STAGEDIR})
    endif()
    
    if (NOT EXISTS "${BUNDLE_STAGEDIR}/lib")
      file(MAKE_DIRECTORY ${BUNDLE_STAGEDIR}/lib)
    endif()
          
    if (NOT EXISTS "${BUNDLE_STAGEDIR}/include")
      file(MAKE_DIRECTORY ${BUNDLE_STAGEDIR}/include)
    endif()

    if ( ${SIRE_BUNDLE_TBB} )
      include( "src/bundled/install_tbb.cmake" )
    endif()

    if ( ${SIRE_BUNDLE_CPUID} )
      include( "src/bundled/install_cpuid.cmake" )
    endif()

    if ( ${SIRE_BUNDLE_GSL} )
      include( "src/bundled/install_gsl.cmake" )

      if (${SIRE_FOUND_GSL})
        list(APPEND SIREMATHS_EXTRA_LIBRARIES ${GSL_LIBRARY})
        list(APPEND SIREMATHS_EXTRA_LIBRARIES ${GSL_CBLAS_LIBRARY})
      endif()

    endif()

    if ( ${SIRE_BUNDLE_BOOST} )
      include( "src/bundled/install_boost.cmake" )
    endif()

    if ( ${SIRE_BUNDLE_QT} )
      include( "src/bundled/install_qt5.cmake" )
    endif()

    if ( NOT ${BUNDLE_STAGEDIR} STREQUAL "${CMAKE_INSTALL_PREFIX}/${SIRE_BUNDLED_DIR}" )
      # All of the bundled objects are built into the directory ${CMAKE_BINARY_DIR}/bundled
      # This directory must be installed into the sire.app/ directory
      install( DIRECTORY ${BUNDLE_STAGEDIR} 
               DESTINATION ${CMAKE_INSTALL_PREFIX}
               USE_SOURCE_PERMISSIONS )
    endif()
endif()

if ( NOT ${SIRE_FOUND_QT} )
    # look for Qt5 in the system
    find_package( Qt5Core REQUIRED )
endif()

get_target_property(QT5_LOCATION Qt5::Core LOCATION)
get_filename_component(SIRE_QT5CORE_DIR ${QT5_LOCATION} PATH)
message( STATUS "Qt5Core location ${SIRE_QT5CORE_DIR}" )
save_sire_variable( "SIRE_QT5CORE_DIR" "${SIRE_QT5CORE_DIR}" )

if (NOT ${SIRE_FOUND_BOOST})
  # Sire depends on >= boost 1.31 - all libraries must be dynamically linked
  set (BOOST_ALL_DYN_LINK "YES")
  set (Boost_ADDITIONAL_VERSIONS "1.39" "1.39.0")
  FIND_PACKAGE ( Boost 1.31 REQUIRED )

  if (Boost_FOUND)
    message(STATUS "Boost paths ${Boost_LIBRARY_DIRS} | ${Boost_INCLUDE_DIR}" )
    set ( BOOST_INCLUDE_DIRS "${Boost_INCLUDE_DIR}" )
    include_directories( ${Boost_INCLUDE_DIR} )

    #save the path to this include directory so that it can be
    #used by anything compiling against Sire
    save_sire_variable( "SIRE_BOOST_INCLUDE_DIR" "${Boost_INCLUDE_DIR}" )
  else()
    message(FATAL_ERROR "Cannot find the boost libraries.")
  endif()  
endif()

if (NOT ${SIRE_FOUND_GSL})
  # Sire depends on >= GSL 1.6
  find_package(GSL 1.6 REQUIRED)

  message( STATUS "GSL paths ${GSL_LINK_DIRECTORIES} "
                  "${GSL_LIBRARIES} | ${GSL_INCLUDE_DIR}" )
         
  list(APPEND SIREMATHS_EXTRA_LIBRARIES ${GSL_LIBRARIES})
                                  
  include_directories( BEFORE ${GSL_INCLUDE_DIR})
                                  
  save_sire_variable( "SIRE_GSL_INCLUDE_PATH" "${GSL_INCLUDE_DIR}" )
endif()

if (NOT ${SIRE_FOUND_BLASLAPACK})
  # Now look for BLAS and LAPACK - if these are not found then
  # a local version will have to be used instead
  set( BLA_STATIC OFF )
  find_package( BLAS )
  find_package( LAPACK )

  if (BLAS_FOUND)
    message( STATUS "Found BLAS library ${BLAS_LIBRARIES}" )
    set( SIRE_USE_BLAS ON )
    save_sire_variable( "SIRE_BLAS_LIBRARIES" "${BLAS_LIBRARIES}" )
  else()
    message( STATUS "Using internal (non-optimal) BLAS functions." )
    save_sire_variable( SIRE_USE_BLAS OFF )
  endif()

  save_sire_variable( "SIRE_USE_BLAS" "${SIRE_USE_BLAS}" )

  if (LAPACK_FOUND)
    message( STATUS "Found LAPACK library ${LAPACK_LIBRARIES}" )
    set( SIRE_USE_LAPACK ON )
    save_sire_variable( "SIRE_LAPACK_LIBRARIES" "${LAPACK_LIBRARIES}" )
  else()
    message( STATUS "Using internal (non-optimal) LAPACK functions." )
    set( SIRE_USE_LAPACK OFF )
  endif()

  save_sire_variable( "SIRE_USE_LAPACK" "${SIRE_USE_LAPACK}" )
endif()

if (NOT ${SIRE_FOUND_CPUID})
  message( STATUS "Cannot find libcpuid. Will disable CPU detection code." )
endif()

###########
########### PRINTING OUT SETTINGS AND FINALISING COMPILE FLAGS
###########

set( SIRE_DISABLE_WARNINGS ON )

if (SIRE_DISABLE_WARNINGS)
  # Disable some of the unnecessary compiler warnings...
  CHECK_CXX_COMPILER_FLAG( "-Wno-attributes" SIRE_HAS_WNO_ATTRIBUTES )

  if (SIRE_HAS_WNO_ATTRIBUTES)
    list( APPEND SIRE_DISABLE_WARNING_FLAGS "-Wno-attributes" )
  endif()

  set( SIRE_WARNALL_FLAGS "${SIRE_WARNALL_FLAGS} ${SIRE_DISABLE_WARNING_FLAGS}" )

endif()

message( STATUS  "CMAKE_SYSTEM_NAME      == ${CMAKE_SYSTEM_NAME}" )
message( STATUS  "CMAKE_C_COMPILER       == ${CMAKE_C_COMPILER}" )
message( STATUS  "CMAKE_CXX_COMPILER     == ${CMAKE_CXX_COMPILER}" )
message( STATUS  "CMAKE_Fortran_COMPILER == ${CMAKE_Fortran_COMPILER}" )

message ( STATUS "SIRE_CXX11_FLAGS       == ${SIRE_CXX11_FLAGS}" )
message ( STATUS "SIRE_SMALL_FLAGS       == ${SIRE_SMALL_FLAGS}" )
message ( STATUS "SIRE_WARNALL_FLAGS     == ${SIRE_WARNALL_FLAGS}" )
message ( STATUS "SIRE_DEBUG_FLAGS       == ${SIRE_DEBUG_FLAGS}" )
message ( STATUS "SIRE_RELEASE_FLAGS     == ${SIRE_RELEASE_FLAGS}" )
message ( STATUS "SIRE_OPENMP_FLAGS      == ${SIRE_OPENMP_FLAGS}" )
message ( STATUS "SIRE_VECTOR_FLAGS      == ${SIRE_VECTOR_FLAGS}" )
message ( STATUS "SIRE_VISIBILITY_FLAGS  == ${SIRE_VISIBILITY_FLAGS}" )
message ( STATUS "SIRE_PLATFORM_FLAGS    == ${SIRE_PLATFORM_FLAGS}" )
message ( STATUS "SIRE_SHARE_LINK_FLAGS  == ${SIRE_SHARE_LINK_FLAGS}" )
message ( STATUS "SIRE_STATIC_LINK_FLAGS == ${SIRE_STATIC_LINK_FLAGS}" )
message ( STATUS "SIRE_EXE_LINK_FLAGS    == ${SIRE_EXE_LINK_FLAGS}" )

set( SIRE_COMPILE_FLAGS "${SIRE_WARNALL_FLAGS} ${SIRE_PLATFORM_FLAGS}" )

if ( SIRE_SYMBOL_HIDING )
  set ( SIRE_COMPILE_FLAGS "${SIRE_COMPILE_FLAGS} ${SIRE_VISIBILITY_FLAGS}" )
  set ( SIRE_SHARE_LINK_FLAGS "${SIRE_SHARE_LINK_FLAGS} ${SIRE_VISIBILITY_FLAGS}" )
endif()

if ( SIRE_VECTORISE )
  set ( SIRE_COMPILE_FLAGS "${SIRE_COMPILE_FLAGS} ${SIRE_VECTOR_FLAGS}" )
endif()

set( SIRE_C_FLAGS_RELEASE "${SIRE_COMPILE_FLAGS} ${SIRE_RELEASE_FLAGS} ${SIRE_USER_CFLAGS}" )
set( SIRE_CXX_FLAGS_RELEASE "${SIRE_COMPILE_FLAGS} ${SIRE_CXX11_FLAGS} ${SIRE_RELEASE_FLAGS}  ${SIRE_USER_CPPFLAGS}" )

set( SIRE_C_FLAGS_DEBUG "${SIRE_COMPILE_FLAGS} ${SIRE_DEBUG_FLAGS} ${SIRE_USER_CFLAGS}" )
set( SIRE_CXX_FLAGS_DEBUG "${SIRE_COMPILE_FLAGS} ${SIRE_CXX11_FLAGS} ${SIRE_DEBUG_FLAGS} ${SIRE_USER_CPPFLAGS}" )

set( SIRE_C_FLAGS_SMALL "${SIRE_COMPILE_FLAGS} ${SIRE_SMALL_FLAGS} ${SIRE_USER_CFLAGS}" )
set( SIRE_CXX_FLAGS_SMALL "${SIRE_COMPILE_FLAGS} ${SIRE_CXX11_FLAGS} ${SIRE_SMALL_FLAGS} ${SIRE_USER_CPPFLAGS}" )

# Copy the Sire compile flags to CMake
set( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${SIRE_C_FLAGS_RELEASE}" )
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${SIRE_CXX_FLAGS_RELEASE}" )
set( CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${SIRE_C_FLAGS_DEBUG}" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${SIRE_CXX_FLAGS_DEBUG}" )
set( CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL} ${SIRE_C_FLAGS_SMALL}" )
set( CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL} ${SIRE_CXX_FLAGS_SMALL}" )

if ( CMAKE_BUILD_TYPE )
    message( STATUS "Compiling Sire using ${CMAKE_BUILD_TYPE}" )
else()
    message( STATUS "Compiling Sire using a default (RELEASE) build" )
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_RELEASE} ${SIRE_C_FLAGS_RELEASE}" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${SIRE_CXX_FLAGS_RELEASE}" )
endif()

message ( STATUS "C compiler flags       == ${CMAKE_C_FLAGS}" )
message ( STATUS "C++ compiler flags     == ${CMAKE_CXX_FLAGS}" )

# Add this to all shared libraries
set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${SIRE_SHARE_LINK_FLAGS}" )
set( CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${SIRE_STATIC_LINK_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${SIRE_EXE_LINK_FLAGS}" )

message( STATUS "Shared library flags    == ${CMAKE_SHARED_LINKER_FLAGS}" )
message( STATUS "Static library flags    == ${CMAKE_STATIC_LINKER_FLAGS}" )
message( STATUS "Executable link flags   == ${CMAKE_EXE_LINKER_FLAGS}" )

# Now save the compile variables so that they can be re-used by libraries that use Sire
save_sire_variable( "SIRE_C_COMPILER" "${CMAKE_C_COMPILER}" )
save_sire_variable( "SIRE_CXX_COMPILER" "${CMAKE_CXX_COMPILER}" )
save_sire_variable( "SIRE_C_FLAGS_RELEASE" "${CMAKE_C_FLAGS_RELEASE}" )
save_sire_variable( "SIRE_CXX_FLAGS_RELEASE" "${CMAKE_CXX_FLAGS_RELEASE}" )
save_sire_variable( "SIRE_C_FLAGS_SMALL" "${CMAKE_C_FLAGS_MINSIZEREL}" )
save_sire_variable( "SIRE_CXX_FLAGS_SMALL" "${CMAKE_CXX_FLAGS_MINSIZEREL}" )
save_sire_variable( "SIRE_C_FLAGS_DEBUG" "${CMAKE_C_FLAGS_DEBUG}" )
save_sire_variable( "SIRE_CXX_FLAGS_DEBUG" "${CMAKE_CXX_FLAGS_DEBUG}" )
save_sire_variable( "SIRE_SHARE_LINK_FLAGS" "${SIRE_SHARE_LINK_FLAGS}" )
save_sire_variable( "SIRE_STATIC_LINK_FLAGS" "${SIRE_STATIC_LINK_FLAGS}" )
save_sire_variable( "SIRE_EXE_LINK_FLAGS" "${SIRE_EXE_LINK_FLAGS}" )
save_sire_variable( "SIRE_SMALL_FLAGS" "${SIRE_SMALL_FLAGS}" )
save_sire_variable( "SIRE_WARNALL_FLAGS" "${SIRE_WARNALL_FLAGS}" )
save_sire_variable( "SIRE_DEBUG_FLAGS" "${SIRE_DEBUG_FLAGS}" )
save_sire_variable( "SIRE_RELEASE_FLAGS" "${SIRE_RELEASE_FLAGS}" )
save_sire_variable( "SIRE_OPENMP_FLAGS" "${SIRE_OPENMP_FLAGS}" )
save_sire_variable( "SIRE_VECTOR_FLAGS" "${SIRE_VECTOR_FLAGS}" )
save_sire_variable( "SIRE_VISIBILITY_FLAGS" "${SIRE_VISIBILITY_FLAGS}" )
save_sire_variable( "SIRE_PLATFORM_FLAGS" "${SIRE_PLATFORM_FLAGS}" )

#include the build directory so that we can get sire_config.h and sire_version.h
include_directories(${CMAKE_BINARY_DIR})

################
################ TELLING CMAKE WHERE THE SOURCE IS LOCATED
################

# Add the compiler test directory - this allows us to run some compiler
# tests to make sure that Sire will compile properly
add_subdirectory (build/test_compiler)

# Recurse into the src/libs and src/apps subdirectories
add_subdirectory (src/libs)
add_subdirectory (src/apps)

################
################ CONFIG FILES AND INSTALLATION INSTRUCTIONS
################

# Finally, write the version header file and finish off the component dependencies file
include( GetSvnVersion )

# install the extra includes file in the Sire/cmake directory
install (FILES ${SIRE_VARIABLES_FILE} 
         DESTINATION ${SIRE_CMAKEFILES}
        )

install (FILES ${CMAKE_CURRENT_BINARY_DIR}/sire_config.h
         DESTINATION ${SIRE_INCLUDES})

# install the parameters file (containing forcefield parameters) into
# share/Sire
install (DIRECTORY ${CMAKE_SOURCE_DIR}/parameters
         DESTINATION ${SIRE_SHARE}
         USE_SOURCE_PERMISSIONS
         PATTERN ".svn" EXCLUDE
        )

install (DIRECTORY ${CMAKE_SOURCE_DIR}/share/build
         DESTINATION ${SIRE_SHARE}
         USE_SOURCE_PERMISSIONS
         PATTERN ".svn" EXCLUDE
        )

# These commands are used to create an 'uninstall' target
# (this is copied from the vtk.org wiki)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/build/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")
