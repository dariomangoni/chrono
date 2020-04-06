#=============================================================================
# CMake configuration file for the Chrono OSQP module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================

option(ENABLE_MODULE_OSQP "Enable the Chrono OSQP module" OFF)

if(NOT ENABLE_MODULE_OSQP)
    mark_as_advanced(FORCE osqp_DIR)
    return()
endif()

message(STATUS "==== Chrono OSQP module ====")

mark_as_advanced(CLEAR osqp_DIR)

# ------------------------------------------------------------------------------
# Dependencies for OSQP module
# ------------------------------------------------------------------------------

message(STATUS "Find OSQP libraries")
find_package(osqp REQUIRED)

if(NOT osqp_FOUND)
  set(osqp_DIR "" CACHE PATH "Location of OSQP installation")
endif()

# ------------------------------------------------------------------------------
# List all files in the Chrono osqp module
# ------------------------------------------------------------------------------

set(ChronoEngine_OSQP_HEADERS
  ChApiOsqp.h
  ChSolverOSQP.h
)

set(ChronoEngine_OSQP_SOURCES
  ChSolverOSQP.cpp
)

source_group("" FILES ${ChronoEngine_OSQP_HEADERS} ${ChronoEngine_OSQP_SOURCES})

# ------------------------------------------------------------------------------
# Add the ChronoEngine_osqp library
# ------------------------------------------------------------------------------

add_library(ChronoEngine_osqp SHARED
            ${ChronoEngine_OSQP_SOURCES}
            ${ChronoEngine_OSQP_HEADERS})

target_link_libraries(ChronoEngine_osqp
            ChronoEngine
            osqp::osqp
            )

set_target_properties(ChronoEngine_osqp PROPERTIES
                      COMPILE_FLAGS "${CH_CXX_FLAGS}"
                      LINK_FLAGS "${CH_LINKERFLAG_SHARED}")

target_compile_definitions(ChronoEngine_osqp PRIVATE "CH_API_COMPILE_OSQP")



install(TARGETS ChronoEngine_osqp
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(FILES ${ChronoEngine_OSQP_HEADERS}
        DESTINATION include/chrono_osqp)
