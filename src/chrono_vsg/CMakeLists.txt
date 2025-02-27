#===============================================================================
# CMake configuration file for the Chrono vsg module
# 
# Cannot be used stand-alone (it's loaded by CMake config. file in parent dir.)
#===============================================================================

option(ENABLE_MODULE_VSG "Enable the Chrono VSG module" OFF)

# Return now if this module is not enabled
if(NOT ENABLE_MODULE_VSG)
   mark_as_advanced(FORCE vsg_DIR)
   mark_as_advanced(FORCE vsgXchange_DIR)
   mark_as_advanced(FORCE vsgImGui_DIR)
   return()
endif()

message(STATUS "\n==== Chrono VSG module ====\n")

# VSG requires C++17 support. Disable module if not available
if(NOT CH_CXX17)
   message(WARNING "Chrono::VSG requires C++17 support which is not available; disabling Chrono::VSG")
   set(ENABLE_MODULE_VSG OFF CACHE BOOL "Enable the Chrono VSG module" FORCE)
   mark_as_advanced(FORCE vsg_DIR)
   mark_as_advanced(FORCE vsgXchange_DIR)
   mark_as_advanced(FORCE vsgImGui_DIR)
   return()
endif()

mark_as_advanced(CLEAR VSG_LIBRARY)
mark_as_advanced(CLEAR VSG_ROOT)

#-------------------------------------------------------------------------------
# Find the VSG library

if(VULKAN_SDK)
    set(ENV{VULKAN_SDK} ${VULKAN_SDK})
endif()

find_package(vsg 1.1.0 REQUIRED)
find_package(vsgXchange 1.1.0 REQUIRED)
find_package(vsgImGui REQUIRED)

message(STATUS "vsg version:        ${vsg_VERSION}")
message(STATUS "vsgXchange version: ${vsgXchange_VERSION}")
message(STATUS "vsgImGui version:   ${vsgImGui_VERSION}")

# Check if vsgXchange is a static library. 
# If yes, we assume it includes STB and we do not add those files in the the Chrono library
# (else we'd get multiple defined symbols at link time).
get_target_property(vsg_type vsgXchange::vsgXchange TYPE)
if (vsg_type STREQUAL STATIC_LIBRARY)
    set(INCLUDE_STB OFF)
    message(STATUS "vsgXchange library: static")
else()
    set(INCLUDE_STB ON)
    message(STATUS "vsgXchange library: dynamic")
endif()

#-------------------------------------------------------------------------------
# List all the files in the ChronoEngine_vsg lib

set(CE_VSG_BASE_FILES
    ChApiVSG.h
    ChVisualSystemVSG.cpp
    ChVisualSystemVSG.h
    ChGuiComponentVSG.h
    ChEventHandlerVSG.h
)

set(CE_VSG_SHAPES_FILES
    shapes/GetSurfaceShapeData.cpp
    shapes/GetSurfaceShapeData.h
    shapes/ShaderUtils.cpp
    shapes/ShaderUtils.h
    shapes/ShapeBuilder.cpp
    shapes/ShapeBuilder.h
)

set(CE_VSG_RESOURCES_FILES
    resources/chronoLineShader_frag.h
    resources/chronoLineShader_vert.h
    resources/chronoPbrShader_frag.h
    resources/chronoPbrShader_vert.h
    resources/pcShader_frag.h
    resources/pcShader_vert.h
    resources/cubemapShaders.h
)

set(CE_VSG_UTILS_FILES
    utils/ChUtilsVSG.cpp
    utils/ChUtilsVSG.h
    utils/ChConversionsVSG.cpp
    utils/ChConversionsVSG.h)

source_group("" FILES ${CE_VSG_BASE_FILES})
source_group("shapes" FILES ${CE_VSG_SHAPES_FILES})
source_group("resources" FILES ${CE_VSG_RESOURCES_FILES})
source_group("utils" FILES ${CE_VSG_UTILS_FILES})

if(INCLUDE_STB)
    set(CE_VSG_STB_FILES
        ${CMAKE_SOURCE_DIR}/src/chrono_thirdparty/stb/stb.h
        ${CMAKE_SOURCE_DIR}/src/chrono_thirdparty/stb/stb_image.h
        ${CMAKE_SOURCE_DIR}/src/chrono_thirdparty/stb/stb_image.cpp
        ${CMAKE_SOURCE_DIR}/src/chrono_thirdparty/stb/stb_image_write.h
        ${CMAKE_SOURCE_DIR}/src/chrono_thirdparty/stb/stb_image_write.cpp
    )
    
    source_group("utils" FILES ${CE_VSG_STB_FILES})
else()
    set(CE_VSG_STB_FILES "")
endif()

#-------------------------------------------------------------------------------
# Let some variables be visible also from outside this directory

#-------------------------------------------------------------------------------
# Create the ChronoEngine_vsg library

add_library(ChronoEngine_vsg
            ${CE_VSG_BASE_FILES}
            ${CE_VSG_SHAPES_FILES}
            ${CE_VSG_RESOURCES_FILES}
            ${CE_VSG_UTILS_FILES}
            ${CE_VSG_STB_FILES})

target_compile_definitions(ChronoEngine_vsg PRIVATE $<$<COMPILE_LANGUAGE:CXX>:CH_API_COMPILE_VSG>)

target_link_libraries(ChronoEngine_vsg PRIVATE ChronoEngine)
target_link_libraries(ChronoEngine_vsg PUBLIC
                      vsg::vsg
                      vsgImGui::vsgImGui
                      vsgXchange::vsgXchange)

# TODO DARIOM: verify if needed
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_link_options(ChronoEngine_vsg PUBLIC "-framework IOKit -framework Cocoa -framework OpenGL")
endif()

#-------------------------------------------------------------------------------
# Install the ChronoEngine_vsg library

install(TARGETS ChronoEngine_vsg
        EXPORT ChronoEngineTargets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        INCLUDES DESTINATION include/chrono_vsg)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
        DESTINATION include/chrono_vsg
        FILES_MATCHING PATTERN "*.h")
