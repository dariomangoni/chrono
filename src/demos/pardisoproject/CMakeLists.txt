set(PARDISOPROJECT_DEMOS
    demo_PP_PardisoProject
    demo_PP_ChPardisoProjectEngine
)

# Additional include paths
INCLUDE_DIRECTORIES( ${CH_PARDISOPROJECT_INCLUDES} )
 
MESSAGE(STATUS "Demo programs for PardisoProject module...")

# ---------------------------------------------------------------------

FOREACH(PROGRAM ${PARDISOPROJECT_DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")
    SOURCE_GROUP("" FILES ${PROGRAM}.cpp)

    ADD_EXECUTABLE(${PROGRAM} ${PROGRAM}.cpp)
    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES 
                FOLDER demos
                COMPILE_FLAGS "${CH_CXX_FLAGS} ${CH_IRRLICHT_CXX_FLAGS}"
                LINK_FLAGS "${CH_LINKERFLAG_EXE}")
    SET_PROPERTY(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")
    TARGET_LINK_LIBRARIES(${PROGRAM}
                            ${PARDISOPROJECT_LIBRARIES}
                            ChronoEngine
                            ChronoEngine_pardisoproject
                            )

    INSTALL(TARGETS ${PROGRAM} DESTINATION ${CH_INSTALL_DEMO})
ENDFOREACH()
