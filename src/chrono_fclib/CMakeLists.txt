#=============================================================================
# CMake configuration file for the Chrono fclib module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================

option(ENABLE_MODULE_FCLIB "Enable the Chrono fclib module" OFF)

IF(NOT ENABLE_MODULE_FCLIB)
    RETURN()
ENDIF()

MESSAGE(STATUS "==== Chrono FCLIB module ====")
# Find HDF5 library
find_package(HDF5 REQUIRED COMPONENTS HL)
include_directories(${HDF5_INCLUDE_DIRS})

# Find fclib-cpp
find_package(FClib REQUIRED)

# Set ChronoEngine fclib variables
set(CH_FCLIB_INCLUDE_DIRS ${FCLib_INCLUDE_DIRS} ${FCLib_EXTRA_INCLUDE_DIRS} PARENT_SCOPE)
set(CH_FCLIB_LIBRARY_DIRS ${FCLib_LIBRARY_DIRS} PARENT_SCOPE)
set(CH_FCLIB_LIBS ${FCLib_LIBRARIES} PARENT_SCOPE)
set(CH_FCLIB_INCLUDES ${FClib_INCLUDES} PARENT_SCOPE)

