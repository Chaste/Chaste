project(testio)

if(ENABLE_CHASTE_TESTING)
#copy test data
file(COPY data DESTINATION "${CHASTE_BUILD_ROOT_DIR}/io/test")

include_directories("${CHASTE_IO_INCLUDE_DIRS}" 
"${CMAKE_CURRENT_SOURCE_DIR}" "${CXXTEST_INCLUDES}" "${CHASTE_GLOBAL_INCLUDE_DIRS}" )

#Auto-generate test suites
set(COMPONENT_LIBRARIES io global )
file(GLOB_RECURSE TEST_SUITES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Test*.hpp)

foreach(test ${TEST_SUITES})
 string(REGEX REPLACE "([a-zA-Z0-9_/]+)[.]hpp" "\\1" testName "${test}")
 chaste_add_test("${testName}Runner" "${testName}.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/${testName}.hpp")
 target_link_libraries("${testName}Runner" ${COMPONENT_LIBRARIES} ${CHASTE_LINK_LIBRARIES})
 set_target_properties("${testName}Runner" PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:MSVCRTD")
endforeach(test)
endif(ENABLE_CHASTE_TESTING)