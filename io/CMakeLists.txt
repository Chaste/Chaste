project(io)
include_directories("${CHASTE_IO_INCLUDE_DIRS}" "${CHASTE_GLOBAL_INCLUDE_DIRS}")
file(GLOB_RECURSE CHASTE_IO_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
add_library(io STATIC ${CHASTE_IO_SOURCES})
install(TARGETS io DESTINATION lib)

#Component testing for io
if(ENABLE_IO_TESTING)
    add_subdirectory(test)
endif()