#--------------------------------------------------------------
# Add executables

SET(TESTS
  test_MR_system_reduction
  test_MR_beam_reduction
)

INCLUDE_DIRECTORIES(${CH_MODELREDUCTION_INCLUDES})

if(ENABLE_MODULE_IRRLICHT)
    include_directories(${CH_IRRLICHTINC})
    set(COMPILER_FLAGS "${COMPILER_FLAGS} ${CH_IRRLICHT_CXX_FLAGS}")
    list(APPEND LIBS "ChronoEngine_irrlicht")
    list(APPEND LIBS "${CH_IRRLICHTLIB}")
endif()

MESSAGE(STATUS "Additional test programs for Chrono_ModelReduction core...")

FOREACH(PROGRAM ${TESTS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER tests
        COMPILE_FLAGS "${CH_CXX_FLAGS}"
        LINK_FLAGS "${CH_LINKERFLAG_EXE}"
    )

    TARGET_LINK_LIBRARIES(${PROGRAM} ChronoEngine ChronoEngine_modelreduction ChronoEngine_fea ${LIBS})

    INSTALL(TARGETS ${PROGRAM} DESTINATION bin)

ENDFOREACH(PROGRAM)
