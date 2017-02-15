#=============================================================================
# CMake configuration file for the Chrono InteriorPoint module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#=============================================================================

option(ENABLE_MODULE_INTERIORPOINT "Enable the Chrono InteriorPoint module" OFF)

# Return now if this module is not enabled
IF(NOT ENABLE_MODULE_INTERIORPOINT)
  RETURN()
ENDIF()

IF (ENABLE_MODULE_MUMPS)
	SET(KERNEL_SOLVER_MODULE_LIB "ChronoEngine_mumps")
	INCLUDE_DIRECTORIES( ${CH_MUMPS_INCLUDES} )
ELSEIF (ENABLE_MODULE_MKL)
	SET(KERNEL_SOLVER_MODULE_LIB "ChronoEngine_mkl")
	INCLUDE_DIRECTORIES( ${CH_MKL_INCLUDES} )
	MESSAGE(SEND_ERROR "InteriorPoint with MKL not implemented yet.")
	RETURN()
ELSE()
	MESSAGE(SEND_ERROR "ChronoEngine InteriorPoint module needs MKL or MUMPS modules to be built.")
	RETURN()
ENDIF()

MESSAGE(STATUS "==== Chrono InteriorPoint module ====")

#-----------------------------------------------------------------------------

SET(ChronoEngine_IP_SOURCES
    ChInteriorPoint.cpp
)

SET(ChronoEngine_IP_HEADERS
    ChInteriorPoint.h
	ChApiInteriorPoint.h
)

SOURCE_GROUP("" FILES 
    ${ChronoEngine_IP_SOURCES} 
    ${ChronoEngine_IP_HEADERS})

#-----------------------------------------------------------------------------	
# Create the ChronoEngine_interiorpoint library

ADD_LIBRARY(ChronoEngine_interiorpoint SHARED 
    ${ChronoEngine_IP_SOURCES}
    ${ChronoEngine_IP_HEADERS})

SET_TARGET_PROPERTIES(ChronoEngine_interiorpoint PROPERTIES
                      COMPILE_FLAGS "${CH_CXX_FLAGS}"
                      LINK_FLAGS "${CH_LINKERFLAG_SHARED}" 
                      COMPILE_DEFINITIONS "CH_API_COMPILE_INTERIORPOINT")
                          
TARGET_LINK_LIBRARIES(ChronoEngine_interiorpoint
                      ChronoEngine
					  ${KERNEL_SOLVER_MODULE_LIB})

ADD_DEPENDENCIES(ChronoEngine_interiorpoint
	ChronoEngine
	${KERNEL_SOLVER_MODULE_LIB}
)

INSTALL(TARGETS ChronoEngine_interiorpoint
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib64
        ARCHIVE DESTINATION lib64)

INSTALL(FILES ${ChronoEngine_IP_HEADERS}
        DESTINATION include/chrono_interiorpoint)
