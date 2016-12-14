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
          GenEigsSolver.h
          HINTS ${SPECTRA_ROOT}
          PATH_SUFFIXES "include"
         )

message(STATUS "   Spectra include dir:  ${SPECTRA_INCLUDE_DIR}")

if(NOT SPECTRA_INCLUDE_DIR)
	set(SPECTRA_ROOT "C:/workspace/libraries/spectra" CACHE PATH "Location of Spectra library")
	RETURN()
else()
    MARK_AS_ADVANCED(FORCE SPECTRA_ROOT)
	MARK_AS_ADVANCED(FORCE SPECTRA_INCLUDE_DIR)
endif()

# ------------------ Eigen

message(STATUS "Find Eigen library")

UNSET(EIGENHEADER_DIR CACHE)

find_path(EIGENHEADER_DIR # path that contains Eigen library
          Eigen
          HINTS ${EIGEN_ROOT_CUSTOM}
          PATH_SUFFIXES Eigen
         )

if(NOT EIGENHEADER_DIR)
	SET(EIGEN_ROOT_CUSTOM "C:/workspace/libraries/eigen" CACHE PATH "Location of Eigen main folder (in which Eigen subfolder can be found)")
	MARK_AS_ADVANCED(CLEAR EIGEN_ROOT_CUSTOM)
	RETURN()
else()
    SET(EIGEN_INCLUDE_DIR "${EIGENHEADER_DIR}/../")
	MARK_AS_ADVANCED(FORCE EIGEN_ROOT_CUSTOM)
    MARK_AS_ADVANCED(FORCE EIGENHEADER_DIR)
endif()

message(STATUS "   Eigen include dir:  ${EIGEN_INCLUDE_DIR}")

	
# ------------------------------------------------------------------------------
# Collect all additional include directories necessary for the ModelReduction module
# ------------------------------------------------------------------------------

SET(CH_MODELREDUCTION_INCLUDES ${SPECTRA_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})

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

ADD_LIBRARY(ChronoEngine_modelreduction SHARED
            ${ChronoEngine_ModelReduction_SOURCES}
            ${ChronoEngine_ModelReduction_HEADERS})

SET_TARGET_PROPERTIES(ChronoEngine_modelreduction PROPERTIES
                      COMPILE_DEFINITIONS "CH_API_COMPILE_MODELREDUCTION")

TARGET_LINK_LIBRARIES(ChronoEngine_modelreduction
                      ChronoEngine
					  )

INSTALL(TARGETS ChronoEngine_modelreduction
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib64
        ARCHIVE DESTINATION lib64)

INSTALL(FILES ${ChronoEngine_MODELREDUCTION_HEADERS} 
        DESTINATION include/chrono_modelreduction)
