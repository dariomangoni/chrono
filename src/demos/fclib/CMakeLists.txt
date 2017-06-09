set(FCLIB_DEMOS
    demo_FCLIB_collisionNSC
 )

# Additional include paths
INCLUDE_DIRECTORIES( ${CH_IRRLICHTINC} )


# ---------------------------------------------------------------------
MESSAGE(STATUS "Demo programs for FClib module...")

FOREACH(PROGRAM ${FCLIB_DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")
    SOURCE_GROUP("" FILES ${PROGRAM}.cpp)

    ADD_EXECUTABLE(${PROGRAM} ${PROGRAM}.cpp)
    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES 
                          FOLDER demos
                          LINK_FLAGS "${CH_LINKERFLAG_EXE}")
	TARGET_INCLUDE_DIRECTORIES(${PROGRAM} PUBLIC ${CH_FCLIB_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES(${PROGRAM} ChronoEngine ChronoEngine_irrlicht ${CH_FCLIB_LIBS} ${CH_IRRLICHTLIB})
    INSTALL(TARGETS ${PROGRAM} DESTINATION ${CH_INSTALL_DEMO})
ENDFOREACH()
