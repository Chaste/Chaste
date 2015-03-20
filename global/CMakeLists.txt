project(global)
chaste_do_component(global ${config_dir}/Version.cpp ${config_dir}/ChasteBuildInfo.cpp)
if (MSCV)
	set_property(TARGET global APPEND PROPERTY COMPILE_DEFINITIONS CHASTE_BUILD_DIR_NAME="")
else()
	set_property(TARGET global APPEND PROPERTY COMPILE_DEFINITIONS CHASTE_BUILD_DIR_NAME="$<CONFIGURATION>")
endif()
add_dependencies(global generateTimestamp)

#message("Configuring component global")

#Build the global library
#include_directories("${CHASTE_global_INCLUDE_DIRS}")
#file(GLOB_RECURSE CHASTE_GLOBAL_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
#set(CHASTE_BUILD_DIR_NAME "")
#set(COMPILE_DEFINTIONS ${COMPILE_DEFINITONS};CHASTE_BUILD_DIR_NAME=${CHASTE_BUILD_DIR_NAME})
#add_library(global ${CHASTE_GLOBAL_SOURCES})
#target_link_libraries(global ${CHASTE_LINK_LIBRARIES})
#add_dependencies(${component} generateTimestamp)
#install(TARGETS global DESTINATION lib COMPONENT global_libraries)
#Component testing for global
#if(ENABLE_global_TESTING)
#    add_subdirectory(test)
#endif()
