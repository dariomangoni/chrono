# ------------------------------------------------------------------------------
# Additional include paths and libraries
# ------------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${CH_MULTICORE_INCLUDES})

SET(LIBRARIES
    ChronoEngine
    ChronoEngine_multicore
)

IF(ENABLE_MODULE_OPENGL)
    INCLUDE_DIRECTORIES(${CH_OPENGL_INCLUDES})
    SET(LIBRARIES ${LIBRARIES} ChronoEngine_opengl)
ENDIF()

# ------------------------------------------------------------------------------
# List of all executables
# ------------------------------------------------------------------------------

# Programs that can work with or without OpenGL
SET(DEMOS
    demo_MCORE_ballsSMC
    demo_MCORE_ballsNSC
    demo_MCORE_mixerSMC
    demo_MCORE_mixerNSC
    demo_MCORE_fluidNSC
    demo_MCORE_snowMPM
    demo_MCORE_particlesNSC
    demo_MCORE_cratering
)

# Add programs that require OpenGL
if(ENABLE_MODULE_OPENGL)
    set(DEMOS
        ${DEMOS} 
        demo_MCORE_motors
        demo_MCORE_friction
        demo_MCORE_mesh_contact
        demo_MCORE_callbackSMC
        demo_MCORE_callbackNSC
    )
endif()

# Add programs that require Irrlicht
if(ENABLE_MODULE_IRRLICHT)
    set(LIBRARIES ${LIBRARIES} ChronoEngine_irrlicht)

	set(DEMOS
        ${DEMOS}
		demo_MCORE_custom_contact
        demo_MCORE_collision_object)
endif()

# ------------------------------------------------------------------------------
# Add all executables
# ------------------------------------------------------------------------------

MESSAGE(STATUS "Demo programs for MULTICORE module...")

FOREACH(PROGRAM ${DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_PROPERTY(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")
    TARGET_LINK_LIBRARIES(${PROGRAM} ${LIBRARIES})

    INSTALL(TARGETS ${PROGRAM} DESTINATION ${CH_INSTALL_DEMO})
ENDFOREACH(PROGRAM)
