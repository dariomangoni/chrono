#=============================================================================
# CMake configuration file for the Chrono ModelReduction module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================

option(ENABLE_MODULE_MODELREDUCTION "Enable the Chrono ModelReduction module" OFF)

IF(NOT ENABLE_MODULE_MODELREDUCTION)
    MARK_AS_ADVANCED(FORCE SPECTRA_ROOT)
    RETURN()
ENDIF()

MESSAGE(STATUS "==== Chrono ModelReduction module ====")

MARK_AS_ADVANCED(CLEAR SPECTRA_ROOT)

# ------------------------------------------------------------------------------
# Dependencies for Chrono ModelReduction module
# ------------------------------------------------------------------------------
# ------------------ Spectra
message(STATUS "Find Spectra library")

UNSET(SPECTRA_INCLUDE_DIR CACHE)

find_path(SPECTRA_INCLUDE_DIR
          "Spectra/GenEigsSolver.h"
          HINTS ${SPECTRA_ROOT}
          PATH_SUFFIXES "include"
         )

message(STATUS "   SPECTRA_INCLUDE_DIR dir:  ${SPECTRA_INCLUDE_DIR}")

if(NOT SPECTRA_INCLUDE_DIR)
	set(SPECTRA_ROOT "C:/workspace/libraries/spectra" CACHE PATH "Location of Spectra library")
	# RETURN()
else()
  MARK_AS_ADVANCED(FORCE SPECTRA_ROOT)
	MARK_AS_ADVANCED(FORCE SPECTRA_INCLUDE_DIR)
endif()

	
# ------------------------------------------------------------------------------
# Collect all additional include directories necessary for the ModelReduction module
# ------------------------------------------------------------------------------

SET(CH_MODELREDUCTION_INCLUDES ${SPECTRA_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${CH_MODELREDUCTION_INCLUDES})
SET(CH_MODELREDUCTION_INCLUDES "${CH_MODELREDUCTION_INCLUDES}" PARENT_SCOPE)

# ------------------------------------------------------------------------------
# Add the ChronoEngine_modelreduction library
# ------------------------------------------------------------------------------

if(ENABLE_MODULE_PARDISO_MKL)
  include_directories(${CH_MKL_INCLUDES})
  set(COMPILER_FLAGS "${COMPILER_FLAGS} ${CH_MKL_CXX_FLAGS}")
  set(LINKER_FLAGS "${LINKER_FLAGS} ${CH_MKL_LINK_FLAGS}")
  list(APPEND ADDITIONAL_LIBS "ChronoEngine_pardisomkl")
  SET(MKL_FOUND true)
  MESSAGE(STATUS "ENABLE_MODULE_PARDISO_MKL in ModalAnalysis")
else()
  message(STATUS "Find MKL libraries")
  MESSAGE(STATUS "Find MKL libraries in ModalAnalysis")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
  #set(MKL_USE_STATIC_LIBS ON)
  #set(MKL_FIND_DEBUG ON)
  find_package(MKL)

  message(STATUS "   MKL include dirs:   ${MKL_INCLUDE_DIRS}")
  message(STATUS "   MKL libraries:      ${MKL_LIBRARIES}")
  message(STATUS "   IOMP5 library:      ${IOMP5_LIBRARY}")
  message(STATUS "   MATH library:       ${MATH_LIBRARY}")
  message(STATUS "   MKL library dirs:   ${MKL_LIBRARY_DIRS}")

  if(NOT MKL_LIBRARIES)
    set(MKL_ROOT "" CACHE PATH "Location of MKL installation")
  endif()

  # Make required libraries visible from outside current directory
  # set(CH_MKL_LIBRARIES ${MKL_LIBRARIES} ${IOMP5_LIBRARY})
  # set(CH_MKL_LIBRARIES "${CH_MKL_LIBRARIES}" PARENT_SCOPE)


  # ------------------------------------------------------------------------------
  # Set up linking line for Intel OpenMP library
  # ------------------------------------------------------------------------------

  if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
    if (NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
      set(CH_MKL_LINK_FLAGS	"-Wl,--as-needed")
    endif()
  endif ()
endif()

if(MKL_FOUND)
  MESSAGE(STATUS "MKL_FOUND in ModalAnalysis")
  find_package(arpack-ng)
  if (arpack-ng_FOUND)
    SET(ARPACK::ARPACK PARENT_SCOPE)
    SET(ENABLE_ARPACK true)
    SET(ENABLE_ARPACK true PARENT_SCOPE)
  endif()
endif()

if(ENABLE_MODULE_IRRLICHT)
    include_directories(${CH_IRRLICHTINC})
    set(COMPILER_FLAGS "${COMPILER_FLAGS} ${CH_IRRLICHT_CXX_FLAGS}")
    list(APPEND ADDITIONAL_LIBS "ChronoEngine_irrlicht")
    list(APPEND ADDITIONAL_LIBS "${CH_IRRLICHTLIB}")
endif()

if(ENABLE_MODULE_MUMPS)
  include_directories(${CH_MUMPS_INCLUDES})
  list(APPEND ADDITIONAL_LIBS "ChronoEngine_mumps")
endif()

# ------------------------------------------------------------------------------
# List all files in the Chrono ModelReduction module
# ------------------------------------------------------------------------------

SET(ChronoEngine_ModelReduction_HEADERS
  ChApiModelReduction.h
  ChModelReduction.h
	ChEigenAnalysis.h
)

SET(ChronoEngine_ModelReduction_SOURCES
  ChModelReduction.cpp
)

if(ENABLE_ARPACK)
  SET(ChronoEngine_ModelReduction_HEADERS
    ${ChronoEngine_ModelReduction_HEADERS}
    ChModalAnalysisArpack.h
  )
  SET(ChronoEngine_ModelReduction_SOURCES
    ${ChronoEngine_ModelReduction_SOURCES}
    ChModalAnalysisArpack.cpp
  )
endif()

SOURCE_GROUP("" FILES ${ChronoEngine_ModelReduction_HEADERS} ${ChronoEngine_ModelReduction_SOURCES})

ADD_LIBRARY(ChronoEngine_modelreduction SHARED
            ${ChronoEngine_ModelReduction_SOURCES}
            ${ChronoEngine_ModelReduction_HEADERS})

SET_TARGET_PROPERTIES(ChronoEngine_modelreduction PROPERTIES
                      COMPILE_DEFINITIONS "CH_API_COMPILE_MODELREDUCTION"
                      COMPILE_FLAGS "${COMPILER_FLAGS}"
                      LINK_FLAGS "${LINKER_FLAGS}")


TARGET_LINK_LIBRARIES(ChronoEngine_modelreduction
                      ChronoEngine
					            ${ADDITIONAL_LIBS}
					  )

INSTALL(TARGETS ChronoEngine_modelreduction
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib64
        ARCHIVE DESTINATION lib64)

INSTALL(FILES ${ChronoEngine_MODELREDUCTION_HEADERS} 
        DESTINATION include/chrono_modelreduction)


if (ENABLE_ARPACK)
MESSAGE(STATUS "MKL_LIBRARY_DIRS is ${MKL_LIBRARY_DIRS}")
  set_target_properties(ChronoEngine_modelreduction PROPERTIES
    LINK_FLAGS "${CH_MKL_LINK_FLAGS}"
    LINK_DIRECTORIES "${MKL_LIBRARY_DIRS}")
  target_include_directories(ChronoEngine_modelreduction PUBLIC ARPACK::ARPACK "${MKL_INCLUDE_DIRS}")
  target_compile_definitions(ChronoEngine_modelreduction PRIVATE "MKL_LP64")
  target_link_libraries(ChronoEngine_modelreduction ARPACK::ARPACK ${MKL_LIBRARIES} ${IOMP5_LIBRARY})
  target_link_directories(ChronoEngine_modelreduction INTERFACE ARPACK::ARPACK "${MKL_LIBRARY_DIRS}")
endif()

