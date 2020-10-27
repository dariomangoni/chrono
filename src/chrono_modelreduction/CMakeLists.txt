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
	RETURN()
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
# List all files in the Chrono ModelReduction module
# ------------------------------------------------------------------------------

SET(ChronoEngine_ModelReduction_HEADERS
    ChApiModelReduction.h
    ChModelReduction.h
	# ChSolverMR.h
	ChEigenAnalysis.h
)

SET(ChronoEngine_ModelReduction_SOURCES
    ChModelReduction.cpp
	# ChSolverMR.cpp
)

SOURCE_GROUP("" FILES ${ChronoEngine_ModelReduction_HEADERS} ${ChronoEngine_ModelReduction_SOURCES})

# ------------------------------------------------------------------------------
# Add the ChronoEngine_modelreduction library
# ------------------------------------------------------------------------------

if(ENABLE_MODULE_PARDISO_MKL)
  include_directories(${CH_MKL_INCLUDES})
  set(COMPILER_FLAGS "${COMPILER_FLAGS} ${CH_MKL_CXX_FLAGS}")
  set(LINKER_FLAGS "${LINKER_FLAGS} ${CH_MKL_LINK_FLAGS}")
  list(APPEND ADDITIONAL_LIBS "ChronoEngine_pardisomkl")
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
