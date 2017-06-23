# Demo for InteriorPoint module

SET(IP_DEMOS
    demo_IP_ballinbucket
    demo_IP_bricks
    demo_IP_convergence
	demo_IP_trivial
	demo_IP_FEA
)

INCLUDE_DIRECTORIES( ${CH_IRRLICHTINC} )

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
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
		${VTK_LIBRARIES}
		)
		
	ADD_DEPENDENCIES(${PROGRAM}
		ChronoEngine
		ChronoEngine_fea
		ChronoEngine_irrlicht
		ChronoEngine_interiorpoint
		)
	INSTALL(TARGETS ${PROGRAM} DESTINATION ${CH_INSTALL_DEMO})

ENDFOREACH()