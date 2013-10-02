project(testglobal)

if(ENABLE_CHASTE_TESTING)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}" "${CHASTE_GLOBAL_INCLUDE_DIRS}")

#Create testglobal library to link global's tests against
file(GLOB_RECURSE CHASTE_TEST_GLOBAL_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
add_library(testglobal STATIC ${CHASTE_TEST_GLOBAL_SOURCES})

#Auto-generate test suites
set(COMPONENT_LIBRARIES global testglobal)
file(GLOB_RECURSE TEST_SUITES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Test*.hpp)
chaste_generate_test_suites()

endif(ENABLE_CHASTE_TESTING)