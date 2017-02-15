# Demo for InteriorPoint module

SET(IP_DEMOS
    demo_IP_ballinbucket
    demo_IP_bricks
    demo_IP_convergence
	demo_IP_trivial
	demo_IP_FEA
)

INCLUDE_DIRECTORIES( ${CH_IRRLICHTINC} )

# # Check which solver is being used to solve the kernel
# IF (ENABLE_MODULE_MUMPS)
	# # SET(KERNEL_SOLVER_MODULE_LIB "ChronoEngine_MUMPS")
	# INCLUDE_DIRECTORIES( ${CH_MUMPS_INCLUDES} )
# ELSEIF (ENABLE_MODULE_MKL)
	# # SET(KERNEL_SOLVER_MODULE_LIB "ChronoEngine_MKL")
	# INCLUDE_DIRECTORIES( ${CH_MKL_INCLUDES} )
	# MESSAGE(SEND_ERROR "InteriorPoint with MKL not implemented yet.")
	# RETURN()
# ELSE()
	# MESSAGE(SEND_ERROR "ChronoEngine InteriorPoint module needs MKL or MUMPS modules to be built.")
	# RETURN()
# ENDIF()
    
FOREACH(PROGRAM ${IP_DEMOS})
	MESSAGE(STATUS "...add ${PROGRAM}")

	ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
	SOURCE_GROUP(""  FILES  "${PROGRAM}.cpp")
	
	SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES 
		FOLDER demos
		COMPILE_FLAGS "${CH_CXX_FLAGS} ${CH_IRRLICHT_CXX_FLAGS}"
		LINK_FLAGS "${CH_LINKERFLAG_EXE}")
		
	TARGET_LINK_LIBRARIES(${PROGRAM}
		${CH_IRRLICHTLIB}
		ChronoEngine
		ChronoEngine_fea
		ChronoEngine_irrlicht
		ChronoEngine_interiorpoint
		# ${KERNEL_SOLVER_MODULE_LIB}
		)
		
	ADD_DEPENDENCIES(${PROGRAM}
		ChronoEngine
		ChronoEngine_fea
		ChronoEngine_irrlicht
		ChronoEngine_interiorpoint
		# ${KERNEL_SOLVER_MODULE_LIB}
		)
	INSTALL(TARGETS ${PROGRAM} DESTINATION ${CH_INSTALL_DEMO})

ENDFOREACH()