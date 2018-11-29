#=============================================================================
# CMake configuration file for the Chrono CVXOPT module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================


option(ENABLE_MODULE_CVXOPT "Enable the Chrono CVXOPT module" OFF)

IF(NOT ENABLE_MODULE_CVXOPT)
    MARK_AS_ADVANCED(FORCE BLAS_LIBRARIES_USER)
    MARK_AS_ADVANCED(FORCE CVXOPT_DIR)
    RETURN()
ENDIF()

MARK_AS_ADVANCED(CLEAR CVXOPT_DIR)

message(STATUS "==== Chrono Engine CVXOPT module ====")
cmake_minimum_required(VERSION 3.12 FATAL_ERROR) # find_package (Python2 COMPONENTS Development REQUIRED) requires >=3.12

SET(CVXOPT_DIR "" CACHE PATH  "CVXOPT folder (must contain src/C folder")

find_path (
          CVXOPT_INCLUDE_DIRS
          "cvxopt.h"
          HINTS "${CVXOPT_DIR}/src/C"
          DOC "Path to the CVXOPT C headers folder"
          NO_DEFAULT_PATH
         )
         
### find python
find_package (Python2 COMPONENTS Development REQUIRED)
IF(Python2_FOUND)
    MESSAGE(STATUS "Python2_INCLUDE_DIRS is ${Python2_INCLUDE_DIRS}")
    MESSAGE(STATUS "Python2_LIBRARIES is ${Python2_LIBRARIES}")
ELSE()
    MESSAGE(ERROR "Python2 required.")
ENDIF()

################# BLAS library ##########################
#BLAS configuration
SET(BLAS_VENDOR "Intel10_64lp" CACHE STRING "BLAS library vendor")
SET_PROPERTY(CACHE BLAS_VENDOR PROPERTY STRINGS  All Goto OpenBLAS ATLAS PhiPACK CXML DXML SunPerf SCSL SGIMATH IBMESSL
												Intel10_32 #(intel mkl v10 32 bit)
												Intel10_64lp #(intel mkl v10 64 bit lp thread model  lp64 model)
												Intel10_64lp_seq #(intel mkl v10 64 bit sequential code  lp64 model)
												Intel #( older versions of mkl 32 and 64 bit) 
												ACML ACML_MP ACML_GPU Apple NAS)

UNSET(BLAS_FOUND)
UNSET(BLAS_LIBRARIES)
SET(CMAKE_LIBRARY_PATH_BKP ${CMAKE_LIBRARY_PATH} )

# Automatic inclusion of MKL folders
if(CMAKE_SIZEOF_VOID_P MATCHES 4)
	SET(MKL_SUBFOLDERARCH "ia32")
elseif (CMAKE_SIZEOF_VOID_P MATCHES 8)
	SET(MKL_SUBFOLDERARCH "intel64")
endif()

if (WIN32)
	SET(PROGRAM_FILE_ENVVAR "PROGRAMFILES(x86)")
	FILE(TO_CMAKE_PATH "$ENV{${PROGRAM_FILE_ENVVAR}}" PRG_FOLD)
	list(APPEND INTEL_LIBRARIES_DIRS "${PRG_FOLD}/IntelSWTools/compilers_and_libraries/windows/mkl/lib/${MKL_SUBFOLDERARCH}") # default for ParallelStudioXE2016 and later
	list(APPEND INTEL_LIBRARIES_DIRS "${PRG_FOLD}/IntelSWTools/compilers_and_libraries/windows/compiler/lib/${MKL_SUBFOLDERARCH}") # default for ParallelStudioXE2016 and later
	
	list(APPEND INTEL_LIBRARIES_DIRS "${PRG_FOLD}/Intel/Composer XE/mkl/lib/${MKL_SUBFOLDERARCH}") # default until ParallelStudioXE2015
	list(APPEND INTEL_LIBRARIES_DIRS "${PRG_FOLD}/Intel/Composer XE/compiler/lib/${MKL_SUBFOLDERARCH}") # default until ParallelStudioXE2015
elseif(UNIX AND NOT APPLE)
	list(APPEND INTEL_LIBRARIES_DIRS "/opt/intel/compilers_and_libraries/linux/mkl/lib/${MKL_SUBFOLDERARCH}") # default for ParallelStudioXE2016 and later (root permissions)
	list(APPEND INTEL_LIBRARIES_DIRS "/opt/intel/compilers_and_libraries/linux/compiler/lib/${MKL_SUBFOLDERARCH}") # default for ParallelStudioXE2016 and later (root permissions)
	
	list(APPEND INTEL_LIBRARIES_DIRS "$ENV{HOME}/intel/compilers_and_libraries/linux/mkl/lib/${MKL_SUBFOLDERARCH}") # default for ParallelStudioXE2016 and later (no root permissions)
	list(APPEND INTEL_LIBRARIES_DIRS "$ENV{HOME}/intel/compilers_and_libraries/linux/compiler/lib/${MKL_SUBFOLDERARCH}") # default for ParallelStudioXE2016 and later (no root permissions)
endif()

SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${BLAS_LIBRARIES_USER} ${INTEL_LIBRARIES_DIRS} )


find_package(BLAS)
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH_BKP}) # restore original CMAKE_LIBRARY_PATH

if (NOT BLAS_FOUND)
	SET(BLAS_LIBRARIES_USER "" CACHE PATH "Locate BLAS library directory") # let the user specifies its own BLAS library
	MARK_AS_ADVANCED(CLEAR BLAS_LIBRARIES_USER)
else()
	MESSAGE(STATUS "Linking BLAS library: ${BLAS_LIBRARIES}")
	MARK_AS_ADVANCED(FORCE BLAS_LIBRARIES_USER)
endif()


string(FIND BLAS_VENDOR "Intel" USING_INTEL)

if (NOT USING_INTEL EQUAL -1)
    SET(CH_CVXOPT_LINK_FLAGS "${IOMP5_LINKERFLAG}")
endif()

LINK_DIRECTORIES(${INTEL_LIBRARIES_DIRS})

				
############# Fortran runtime libraries


IF (CVXOPT_INCLUDE_DIRS)
    
  SET(CVXOPT_DEPENDENCY_LIBRARIES
							${BLAS_LIBRARIES} # provided by FindBlas
                            ${Python2_LIBRARIES}
							)
	
  SET(CVXOPT_DEPENDENCY_LIBRARIES ${CVXOPT_DEPENDENCY_LIBRARIES} PARENT_SCOPE)
							
  SET(CH_CVXOPT_INCLUDE_DIRS
                            ${CVXOPT_INCLUDE_DIRS}
                            ${Python2_INCLUDE_DIRS}
                            )
                            
  SET(CH_CVXOPT_INCLUDE_DIRS ${CH_CVXOPT_INCLUDE_DIRS} PARENT_SCOPE)
							
  SET(CH_CVXOPT_LINK_FLAGS ${CH_CVXOPT_LINK_FLAGS} PARENT_SCOPE)
  

  check_symbol_exists(snprintf stdio.h HAVE_SNPRINTF)
  IF(HAVE_SNPRINTF)
    SET(CH_CVXOPT_DEFS "HAVE_SNPRINTF")
  ENDIF()
  
  SET(CH_CVXOPT_DEFS ${CH_CVXOPT_DEFS} PARENT_SCOPE)
  
  
ENDIF()

# ------------------------------------------------------------------------------
# List all files in the Chrono cvxopt module
# ------------------------------------------------------------------------------

SET(ChronoEngine_CVXOPT_HEADERS
  ChApiCvxopt.h
  ChCvxoptConeQpEngine.h
  # ChCvxoptLinearEngine.h
  ChSolverCvxoptConeQp.h
  ChCSMatrixIndexType.h
)

SET(ChronoEngine_CVXOPT_SOURCES
  ChCvxoptConeQpEngine.cpp
  # ChCvxoptLinearEngine.cpp
  ChSolverCvxoptConeQp.cpp
)

SOURCE_GROUP("" FILES ${ChronoEngine_CVXOPT_HEADERS} ${ChronoEngine_CVXOPT_SOURCES})


# ------------------------------------------------------------------------------
# Add the ChronoEngine_cvxopt library
# ------------------------------------------------------------------------------



ADD_LIBRARY(ChronoEngine_cvxopt SHARED
            ${ChronoEngine_CVXOPT_SOURCES}
            ${ChronoEngine_CVXOPT_HEADERS})

SET_TARGET_PROPERTIES(ChronoEngine_cvxopt PROPERTIES
                      COMPILE_FLAGS "${CH_CXX_FLAGS}"
                      LINK_FLAGS "${CH_LINKERFLAG_SHARED}"
                      COMPILE_DEFINITIONS "CH_API_COMPILE_CVXOPT;${CH_CVXOPT_DEFS}")

TARGET_LINK_LIBRARIES(ChronoEngine_cvxopt
                      ChronoEngine
                      ${CVXOPT_DEPENDENCY_LIBRARIES}
                      )
					  
TARGET_INCLUDE_DIRECTORIES(ChronoEngine_cvxopt PUBLIC
							${CH_CVXOPT_INCLUDE_DIRS})
                            
INSTALL(TARGETS ChronoEngine_cvxopt
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib64
        ARCHIVE DESTINATION lib64)

INSTALL(FILES ${ChronoEngine_CVXOPT_HEADERS} 
        DESTINATION include/chrono_cvxopt)