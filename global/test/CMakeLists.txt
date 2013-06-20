project(testglobal)

#copy test data
file(COPY data DESTINATION "${CHASTE_BUILD_ROOT_DIR}/global/test")

include_directories(#"${CHASTE_GLOBAL_INCLUDE_DIRS}" 
"${CMAKE_CURRENT_SOURCE_DIR}" "${CXXTEST_INCLUDES}" "${CHASTE_LINALG_INCLUDE_DIRS}" )

#Create testglobal library to link global's tests against
file(GLOB_RECURSE CHASTE_TEST_GLOBAL_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
add_library(testglobal STATIC ${CHASTE_TEST_GLOBAL_SOURCES})

set(TEST_MPIEXEC_ENABLED TRUE)

#Auto-generate test suites
set(COMPONENT_LIBRARIES global testglobal)
file(GLOB_RECURSE TEST_SUITES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Test*.hpp)

foreach(test ${TEST_SUITES})
 string(REGEX REPLACE "([a-zA-Z0-9_/]+)[.]hpp" "\\1" testName "${test}")
 chaste_add_test("${testName}Runner" "${testName}.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/${testName}.hpp")
 target_link_libraries("${testName}Runner" ${COMPONENT_LIBRARIES} ${CHASTE_LINK_LIBRARIES})
 set_target_properties("${testName}Runner" PROPERTIES LINK_FLAGS /NODEFAULTLIB:LIBCMT)
endforeach(test)

