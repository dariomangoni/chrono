#=============================================================================
# CMake configuration file for the Chrono MUMPS module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================


option(ENABLE_MODULE_MUMPS "Enable the Chrono MUMPS module" OFF)

IF(NOT ENABLE_MODULE_MUMPS)
    MARK_AS_ADVANCED(FORCE BLAS_LIBRARIES_USER)
    RETURN()
ENDIF()

message(STATUS "==== Chrono Engine Mumps module ====")

# Arithmetic configuration
SET(MUMPS_ROOT "C:/Program Files/Mumps535_libseq" CACHE PATH "MUMPS folder")

# find_package(BLAS REQUIRED)
find_package(MKL REQUIRED)
# MESSAGE(STATUS "BLAS_LIBRARIES are: ${BLAS_LIBRARIES}")
MESSAGE(STATUS "MKL_FOUND is: ${MKL_FOUND}")

			
############# Fortran runtime libraries

SET(MUMPS_LIBRARIES
      "${MUMPS_ROOT}/lib/dmumps.lib;${MUMPS_ROOT}/lib/mumps_common.lib;${MUMPS_ROOT}/lib/pord.lib;${MUMPS_ROOT}/lib/mpiseq.lib" # provided by FindMumps
      ${MKL_LIBRARIES}
      ${MATH_LIBRARY}
      ${IOMP5_LIBRARY}
            )

SET(MUMPS_LIBRARIES ${MUMPS_LIBRARIES} PARENT_SCOPE)
            
SET(CH_MUMPS_INCLUDES
          "${MUMPS_ROOT}/include" # provided by FindMumps
            )
            
SET(CH_MUMPS_INCLUDES ${CH_MUMPS_INCLUDES} PARENT_SCOPE)

SET(CH_MUMPS_LINK_FLAGS
      ${MKL_LINKER_FLAGS}
            )
            
SET(CH_MUMPS_LINK_FLAGS ${CH_MUMPS_LINK_FLAGS} PARENT_SCOPE)

SET(CH_MUMPS_LIBRARY_DIRS
      ${MKL_LIBRARY_DIRS}
            )
            
SET(CH_MUMPS_LIBRARY_DIRS ${CH_MUMPS_LIBRARY_DIRS} PARENT_SCOPE)



# ------------------------------------------------------------------------------
# Collect all additional include directories necessary for the Mumps module
# ------------------------------------------------------------------------------
INCLUDE_DIRECTORIES(${CH_MUMPS_INCLUDES})

# ------------------------------------------------------------------------------
# List all files in the Chrono mumps module
# ------------------------------------------------------------------------------

SET(ChronoEngine_MUMPS_HEADERS
  ChApiMumps.h
  ChSolverMumps.h
  ChMumpsEngine.h
)

SET(ChronoEngine_MUMPS_SOURCES
  ChSolverMumps.cpp
  ChMumpsEngine.cpp
)

SOURCE_GROUP("" FILES ${ChronoEngine_MUMPS_HEADERS} ${ChronoEngine_MUMPS_SOURCES})


# ------------------------------------------------------------------------------
# Add the ChronoEngine_mumps library
# ------------------------------------------------------------------------------
ADD_DEFINITIONS("-D_OPENMP_NOFORCE_MANIFEST")


ADD_LIBRARY(ChronoEngine_mumps SHARED
            ${ChronoEngine_MUMPS_SOURCES}
            ${ChronoEngine_MUMPS_HEADERS})

SET_TARGET_PROPERTIES(ChronoEngine_mumps PROPERTIES
                      COMPILE_FLAGS "${CH_CXX_FLAGS}"
                      LINK_FLAGS "${CH_LINKERFLAG_SHARED} ${CH_MUMPS_LINK_FLAGS}")

target_compile_definitions(ChronoEngine_mumps PRIVATE "CH_API_COMPILE_MUMPS")
target_compile_definitions(ChronoEngine_mumps PRIVATE "CH_IGNORE_DEPRECATED")

TARGET_LINK_LIBRARIES(ChronoEngine_mumps
                      ChronoEngine
                      ${MUMPS_LIBRARIES})

TARGET_LINK_DIRECTORIES(ChronoEngine_mumps PRIVATE
                      ${CH_MUMPS_LIBRARY_DIRS})
					  
TARGET_INCLUDE_DIRECTORIES(ChronoEngine_mumps PUBLIC
							${CH_MUMPS_INCLUDES})
                            
# IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#   ADD_CUSTOM_COMMAND(
#     TARGET ChronoEngine_mumps POST_BUILD
#     COMMAND ${CMAKE_COMMAND} -E make_directory
#             "${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>"
#     COMMAND ${CMAKE_COMMAND} -E copy_if_different
#             "${MUMPS_SHARED_LIBRARIES}"
#             "${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>"
#     MAIN_DEPENDENCY ChronoEngine_mumps
#     COMMENT "Copying Mumps DLL..."
#   )

  # # Make variable visible from outside this directory
  # SET(MUMPS_SHARED_LIBRARIES ${MUMPS_SHARED_LIBRARIES} PARENT_SCOPE)

# ENDIF()

INSTALL(TARGETS ChronoEngine_mumps
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

INSTALL(FILES ${ChronoEngine_MUMPS_HEADERS} 
        DESTINATION include/chrono_mumps)