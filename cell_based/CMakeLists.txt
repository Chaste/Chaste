project(cell_based)
#The compiler flag /bigobj deals with the "fatal error C1128": 'number of sections exceeded object file format limit'
add_definitions(-bigobj) 
include_directories("${CHASTE_CELL_BASED_INCLUDE_DIRS}" "${CHASTE_PDE_INCLUDE_DIRS}" 
	"${CHASTE_ODE_INCLUDE_DIRS}" "${CHASTE_MESH_INCLUDE_DIRS}" "${CHASTE_LINALG_INCLUDE_DIRS}" 
	"${CHASTE_IO_INCLUDE_DIRS}"  "${CHASTE_GLOBAL_INCLUDE_DIRS}")
file(GLOB_RECURSE CHASTE_CELL_BASED_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
add_library(cell_based STATIC ${CHASTE_CELL_BASED_SOURCES})
install(TARGETS cell_based DESTINATION lib)

#Component testing for cell_based
if(ENABLE_CELL_BASED_TESTING)
	add_subdirectory(test)
endif(ENABLE_CELL_BASED_TESTING)