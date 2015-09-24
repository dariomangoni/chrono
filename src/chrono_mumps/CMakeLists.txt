#=============================================================================
# CMake configuration file for the Chrono MUMPS module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================

message(STATUS "==== Chrono Engine Mumps module ====")

option(ENABLE_MODULE_MUMPS "Enable the Chrono MUMPS module" OFF)

# Arithmetic configuration
SET(MUMPS_ARITHMETIC "double" CACHE STRING "MUMPS arithmetic")
SET_PROPERTY(CACHE MUMPS_ARITHMETIC PROPERTY STRINGS float double complex complex-double)


IF (MUMPS_ARITHMETIC MATCHES "float")
	SET(MUMPS_ARITH_PREFIX "s")
	MESSAGE("No float library implemented yet")
ELSEIF (MUMPS_ARITHMETIC MATCHES "double")
	SET(MUMPS_ARITH_PREFIX "d")
ELSEIF (MUMPS_ARITHMETIC MATCHES "complex")
	SET(MUMPS_ARITH_PREFIX "c")
	MESSAGE("No complex-float library implemented yet")
ELSEIF (MUMPS_ARITHMETIC MATCHES "complex-double")
	SET(MUMPS_ARITH_PREFIX "z")
	MESSAGE("No complex-double library implemented yet")
ENDIF ()

MESSAGE(STATUS "Mumps arithmetic set to '${MUMPS_ARITHMETIC}'")


find_package(mumps)
find_package(BLAS)
#TODO: do I need to link to BLAS???	



IF (MUMPS_FOUND)
SET (CH_MUMPS_LIBRARIES
						"${MUMPS_LIBRARIES}"
						"${BLAS_LIBRARIES}" #aka LIBBLAS
						)
						
SET (CH_MUMPS_INCLUDES
						"${MUMPS_INCLUDE_DIRS}"
						"${BLAS_INCLUDE_DIRS}" #aka LIBBLAS
						)
ENDIF()


# ------------------------------------------------------------------------------
# Collect all additional include directories necessary for the Mumps module
# ------------------------------------------------------------------------------
INCLUDE_DIRECTORIES(${CH_MUMPS_INCLUDES})

SET(CH_MUMPS_INCLUDES "${CH_MUMPS_INCLUDES}" PARENT_SCOPE)

# ------------------------------------------------------------------------------
# List all files in the Chrono mumps module
# ------------------------------------------------------------------------------

SET(ChronoEngine_MUMPS_HEADERS
  ChApiMumps.h
  ChSolverMumps.h
  ChMumpsEngine.h
  ChCOOMatrix.h
)

SET(ChronoEngine_MUMPS_SOURCES
  ChSolverMumps.cpp
  ChMumpsEngine.cpp
  ChCOOMatrix.cpp
)

SOURCE_GROUP("Headers" FILES ${ChronoEngine_MUMPS_HEADERS})
SOURCE_GROUP("Sources" FILES ${ChronoEngine_MUMPS_SOURCES})


# ------------------------------------------------------------------------------
# Add the ChronoEngine_mumps library
# ------------------------------------------------------------------------------
#TODO
LINK_DIRECTORIES("C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows/compiler/lib/intel64")

ADD_LIBRARY(ChronoEngine_mumps SHARED
            ${ChronoEngine_MUMPS_SOURCES}
            ${ChronoEngine_MUMPS_HEADERS})

SET_TARGET_PROPERTIES(ChronoEngine_mumps PROPERTIES
                      COMPILE_FLAGS "${CH_CXX_FLAGS}"
                      LINK_FLAGS "${CH_LINKERFLAG_SHARED}"
                      COMPILE_DEFINITIONS "CH_API_COMPILE_MUMPS")

TARGET_LINK_LIBRARIES(ChronoEngine_mumps
                      ChronoEngine
                      ${CH_MUMPS_LIBRARIES})
					  
TARGET_INCLUDE_DIRECTORIES(ChronoEngine_mumps PUBLIC
							${CH_MUMPS_INCLUDES})

INSTALL(TARGETS ChronoEngine_mumps
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

INSTALL(FILES ${ChronoEngine_MUMPS_HEADERS} 
        DESTINATION include/chrono_mumps)