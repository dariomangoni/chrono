#===============================================================================
# CMake configuration file for the Chrono irrlicht module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#===============================================================================

option(ENABLE_MODULE_IRRLICHT "Enable the Chrono Irrlicht module" OFF)

# Return now if this module is not enabled
IF(NOT ENABLE_MODULE_IRRLICHT)
  RETURN()
ENDIF()

MESSAGE(STATUS "\n==== Chrono Irrlicht module ====\n")

#-------------------------------------------------------------------------------
# List all the files in the ChronoEngine_irrlicht lib
SET(ChronoEngine_IRRLICHT_SOURCES
  ChIrrCamera.cpp
  ChIrrMeshTools.cpp
  ChIrrTools.cpp
  ChIrrSkyBoxSceneNode.cpp
  ChIrrGUI.cpp
  ChVisualSystemIrrlicht.cpp
  ChIrrNodeModel.cpp
  ChIrrNodeShape.cpp
)

SET(ChronoEngine_IRRLICHT_HEADERS
  ChApiIrr.h
  ChIrrTools.h
  ChIrrCamera.h
  ChIrrMeshTools.h
  ChIrrEffects.h
  ChIrrSkyBoxSceneNode.h
  ChIrrGUI.h
  ChVisualSystemIrrlicht.h
  ChIrrNodeModel.h
  ChIrrNodeShape.h
)

SOURCE_GROUP("" FILES
  ${ChronoEngine_IRRLICHT_SOURCES}
  ${ChronoEngine_IRRLICHT_HEADERS}
)

#-------------------------------------------------------------------------------
# Find the IRRLICHT library

find_package(IRRLICHT)

if (IRRLICHT_FOUND)
    message(STATUS "IRRLICHT found")
else()
    message(ERROR "IRRLICHT not found. Please set IRRLICHT_DIR to the Irrlicht SDK installation directory.")
endif()

# Let some variables be visible also from outside this directory
set(CH_IRRLICHT_INCLUDES  "${IRRLICHT_INCLUDE_DIR}" PARENT_SCOPE)
set(CH_IRRLICHT_CXX_FLAGS "${IRRLICHT_CXX_FLAGS}"   PARENT_SCOPE)
set(CH_IRRLICHT_C_FLAGS   "${IRRLICHT_C_FLAGS}"     PARENT_SCOPE)

#-------------------------------------------------------------------------------
# Create the ChronoEngine_irrlicht library

include_directories(${IRRLICHT_INCLUDE_DIR})

# set(COMPILER_FLAGS "${CH_CXX_FLAGS} ${IRRLICHT_CXX_FLAGS}")
# set(LINKER_FLAGS "${CH_LINKERFLAG_LIB}")
set(LIBRARIES 
    ChronoEngine
    Irrlicht::Irrlicht
    )

if (ENABLE_MODULE_POSTPROCESS)
    set(LIBRARIES ${LIBRARIES} ChronoEngine_postprocess)
endif()
if (ENABLE_MODULE_MODAL)
    set(LIBRARIES ${LIBRARIES} ChronoEngine_modal)
endif()

add_library(ChronoEngine_irrlicht
  ${ChronoEngine_IRRLICHT_SOURCES}
  ${ChronoEngine_IRRLICHT_HEADERS}
)

add_library(ChronoEngine::ChronoEngine_irrlicht ALIAS ChronoEngine_irrlicht)

#-----------------------------------------------------------------------------
# Select proper MSVC Runtime Library
#-----------------------------------------------------------------------------
set_target_properties(ChronoEngine PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})

target_compile_definitions(ChronoEngine_irrlicht PRIVATE "CH_API_COMPILE_IRRLICHT")
target_compile_definitions(ChronoEngine_irrlicht PRIVATE "CH_IGNORE_DEPRECATED")
target_compile_definitions(ChronoEngine_irrlicht PRIVATE "IGNORE_DEPRECATED_WARNING")

target_link_libraries(ChronoEngine_irrlicht PUBLIC
  ${LIBRARIES}
)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_link_libraries(ChronoEngine_irrlicht
    ${MAC_LIBS}
  )
ENDIF()

ADD_DEPENDENCIES(ChronoEngine_irrlicht ChronoEngine)

if (ENABLE_MODULE_POSTPROCESS)
    ADD_DEPENDENCIES(ChronoEngine_irrlicht ChronoEngine_postprocess)
endif()

if (ENABLE_MODULE_MODAL)
    ADD_DEPENDENCIES(ChronoEngine_irrlicht ChronoEngine_modal)
endif()


#-------------------------------------------------------------------------------
# On Windows, copy the Irrlicht DLL to the bin/ directory
# We do this through a custom command so that the Irrlich DLL is copied to the
# appropriate directory (depending on the build type); however, we use
# copy_if_different.


IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

  # Make variable visible from outside this directory
  SET(CH_IRRLICHT_DLL "${IRRLICHT_DLL}")
  SET(CH_IRRLICHT_DLL "${CH_IRRLICHT_DLL}" PARENT_SCOPE)

  add_custom_command(
    TARGET ChronoEngine_irrlicht POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CH_IRRLICHT_DLL}"
            "${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>"
    MAIN_DEPENDENCY ChronoEngine_irrlicht
    COMMENT "Copying Irrlicht DLL..."
  )

  # copy the IRRLICHT dll into the installed package to be on the safe side 
  install(FILES "${CH_IRRLICHT_DLL}" DESTINATION bin)

ENDIF()

if(CH_STATIC)
  SET(CH_IRRLICHT_LIB "${IRRLICHT_LIBRARY}" PARENT_SCOPE)
endif()



#-------------------------------------------------------------------------------
# Install the ChronoEngine_irrlicht library

install(TARGETS ChronoEngine_irrlicht
  EXPORT ChronoEngineTargets
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION include/chrono_irrlicht
)

# # TODO: DARIOM check if this is already covered by the install TARGETS above
# install(FILES ${ChronoEngine_IRRLICHT_HEADERS}
#         DESTINATION include/chrono_irrlicht)

file(COPY "${CMAKE_SOURCE_DIR}/cmake/FindIrrlicht.cmake"
     DESTINATION ${CMAKE_BINARY_DIR}/cmake/)

install(FILES "${CMAKE_SOURCE_DIR}/cmake/FindIrrlicht.cmake"
        DESTINATION ${CH_CONFIG_INSTALL_PATH}
        )