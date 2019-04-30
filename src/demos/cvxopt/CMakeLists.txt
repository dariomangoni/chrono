INCLUDE_DIRECTORIES( ${CH_CVXOPT_INCLUDE_DIRS} )
set(COMPILER_FLAGS "${CH_CXX_FLAGS}")
set(LINKER_FLAGS "${CH_LINKERFLAG_EXE} ${CH_CVXOPT_LINK_FLAGS}" )
list(APPEND LIBRARIES "${CVXOPT_DEPENDENCY_LIBRARIES};ChronoEngine;ChronoEngine_cvxopt;")

MESSAGE(STATUS "Demo programs for CVXOPT module...")

set(CVXOPT_DEMOS
    demo_CVXOPT_ConeQpEngine
)

FOREACH(PROGRAM ${CVXOPT_DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")
    SOURCE_GROUP("" FILES ${PROGRAM}.cpp)

    ADD_EXECUTABLE(${PROGRAM} ${PROGRAM}.cpp)
    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES 
                          FOLDER demos
                          COMPILE_DEFINITIONS "${CH_CVXOPT_DEFS}"
                          LINK_FLAGS ${LINKER_FLAGS})
    TARGET_LINK_LIBRARIES(${PROGRAM} ${LIBRARIES})
    INSTALL(TARGETS ${PROGRAM} DESTINATION ${CH_INSTALL_DEMO})
    
ENDFOREACH()

IF (ENABLE_MODULE_IRRLICHT)

    set(CVXOPT_DEMOS_IRRLICHT
        demo_CVXOPT_SolverConeQp
        demo_CVXOPT_contacts
        demo_CVXOPT_friction
        demo_CVXOPT_beamsIGA
        demo_CVXOPT_torqueIGA
        demo_CVXOPT_ThreePointBending
    )
    
    FOREACH(PROGRAM ${CVXOPT_DEMOS_IRRLICHT})
        MESSAGE(STATUS "...add ${PROGRAM}")
        SOURCE_GROUP("" FILES ${PROGRAM}.cpp)

        ADD_EXECUTABLE(${PROGRAM} ${PROGRAM}.cpp)
        SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES 
                              FOLDER demos
                              COMPILE_DEFINITIONS "${CH_CVXOPT_DEFS}"
                              COMPILE_FLAGS "${COMPILER_FLAGS} ${CH_IRRLICHT_CXX_FLAGS}"
                              LINK_FLAGS ${LINKER_FLAGS})
        TARGET_INCLUDE_DIRECTORIES(${PROGRAM} PUBLIC ${CH_IRRLICHTINC})
        TARGET_LINK_LIBRARIES(${PROGRAM} "${LIBRARIES};ChronoEngine_irrlicht;${CH_IRRLICHTLIB}")
        INSTALL(TARGETS ${PROGRAM} DESTINATION ${CH_INSTALL_DEMO})
        
    ENDFOREACH()
ENDIF (ENABLE_MODULE_IRRLICHT)



# ---------------------------------------------------------------------


