project(mesh)

#Define TETLIBRARY: prevents Tetgen from including main() method
add_definitions(-DTETLIBRARY)
add_definitions(-DANSI_DECLARATORS -DTRILIBRARY)

include_directories("${CHASTE_MESH_INCLUDE_DIRS}" "${CHASTE_LINALG_INCLUDE_DIRS}"
    "${CHASTE_GLOBAL_INCLUDE_DIRS}")
file(GLOB_RECURSE CHASTE_MESH_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
add_library(mesh STATIC ${CHASTE_MESH_SOURCES})
install(TARGETS mesh DESTINATION lib)

#Component testing for mesh
if(ENABLE_MESH_TESTING)
    add_subdirectory(test)
endif(ENABLE_MESH_TESTING)