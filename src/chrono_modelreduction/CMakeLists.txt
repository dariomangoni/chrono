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

message(STATUS "Find Spectra library")

UNSET(SPECTRA_INCLUDE_DIRS CACHE)
find_path(SPECTRA_INCLUDE_DIRS
          GenEigsSolver.h
          HINTS ${SPECTRA_ROOT} "C:/workspace/libraries/spectra"
          PATH_SUFFIXES "include"
         )

message(STATUS "   Spectra include dirs:  ${SPECTRA_INCLUDE_DIRS}")

if(NOT SPECTRA_INCLUDE_DIRS)
  set(SPECTRA_ROOT "" CACHE PATH "Location of Spectra library")
else()
    MARK_AS_ADVANCED(FORCE SPECTRA_ROOT)
	MARK_AS_ADVANCED(FORCE SPECTRA_INCLUDE_DIRS)
endif()



message(STATUS "Find Eigen library")

UNSET(EIGEN_INCLUDE_DIRS CACHE)
find_path(EIGEN_INCLUDE_DIRS
          Eigen
          HINTS ${EIGEN_ROOT} "C:/workspace/libraries/eigen"
          PATH_SUFFIXES "Eigen"
         )

message(STATUS "   Eigen include dirs:  ${EIGEN_INCLUDE_DIRS}")

if(NOT EIGEN_INCLUDE_DIRS)
  set(EIGEN_ROOT "" CACHE PATH "Location of Eigen library (in which Eigen subfolder can be found)")
else()
    SET(EIGEN_INCLUDE_DIRS "${EIGEN_INCLUDE_DIRS}/../")
    MARK_AS_ADVANCED(FORCE EIGEN_ROOT)
	MARK_AS_ADVANCED(FORCE EIGEN_INCLUDE_DIRS)
endif()

	
# ------------------------------------------------------------------------------
# Collect all additional include directories necessary for the MODELREDUCTION module
# ------------------------------------------------------------------------------

SET(CH_MODELREDUCTION_INCLUDES ${SPECTRA_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(${CH_MODELREDUCTION_INCLUDES})
SET(CH_MODELREDUCTION_INCLUDES "${CH_MODELREDUCTION_INCLUDES}" PARENT_SCOPE)

# ------------------------------------------------------------------------------
# List all files in the Chrono ModelReduction module
# ------------------------------------------------------------------------------

SET(ChronoEngine_ModelReduction_HEADERS
    ChApiModelReduction.h
    ChModelReduction.h
)

SET(ChronoEngine_ModelReduction_SOURCES
    ChModelReduction.cpp
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
