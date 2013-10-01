project(linalg)
include_directories("${CHASTE_LINALG_INCLUDE_DIRS}" "${CHASTE_GLOBAL_INCLUDE_DIRS}")
file(GLOB_RECURSE CHASTE_LINALG_SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.cpp *.hpp)
add_library(linalg STATIC ${CHASTE_LINALG_SOURCES})
install(TARGETS linalg DESTINATION lib)

#Component testing for linalg
if(ENABLE_LINALG_TESTING)
    add_subdirectory(test)
endif(ENABLE_LINALG_TESTING)