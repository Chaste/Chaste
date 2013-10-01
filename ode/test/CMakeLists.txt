project(testode)

if(ENABLE_CHASTE_TESTING)

add_definitions(-D_USE_MATH_DEFINES)

include_directories("${CHASTE_ODE_INCLUDE_DIRS}" 
    "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/odes" "${CXXTEST_INCLUDES}" "${CHASTE_GLOBAL_INCLUDE_DIRS}" 
    "${CHASTE_LINALG_INCLUDE_DIRS}" "${CHASTE_IO_INCLUDE_DIRS}" )

#Auto-generate test suites
set(COMPONENT_LIBRARIES ode linalg io global)
file(GLOB_RECURSE TEST_SUITES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Test*.hpp)

chaste_generate_test_suites()

endif(ENABLE_CHASTE_TESTING)