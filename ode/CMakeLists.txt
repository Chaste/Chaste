project(ode)
include_directories("${CHASTE_ODE_INCLUDE_DIRS}" "${CHASTE_LINALG_INCLUDE_DIRS}" 
	"${CHASTE_IO_INCLUDE_DIRS}"  "${CHASTE_GLOBAL_INCLUDE_DIRS}")
file(GLOB_RECURSE CHASTE_ODE_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
add_library(ode STATIC ${CHASTE_ODE_SOURCES})
install(TARGETS ode DESTINATION lib)

#Component testing for ode
if(ENABLE_ODE_TESTING)
	add_subdirectory(test)
endif(ENABLE_ODE_TESTING)