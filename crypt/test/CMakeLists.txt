project(testcrypt)

if(ENABLE_CHASTE_TESTING)

#Create testcrypt library to link crypt's tests against
file(GLOB_RECURSE CHASTE_TEST_CRYPT_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
add_definitions(-D_USE_MATH_DEFINES)#cmath and M_PI

include_directories("${CHASTE_CRYPT_INCLUDE_DIRS}" 
"${CMAKE_CURRENT_SOURCE_DIR}" "${CXXTEST_INCLUDES}" "${CHASTE_GLOBAL_INCLUDE_DIRS}" 
"${CHASTE_LINALG_INCLUDE_DIRS}" "${CHASTE_IO_INCLUDE_DIRS}" )

#Auto-generate test suites
set(COMPONENT_LIBRARIES crypt cell_based pde ode mesh linalg io global)
file(GLOB_RECURSE TEST_SUITES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Test*.hpp)

set(HAS_OWN_LINKER_FLAGS TRUE)
set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:MSVCRT /NODEFAULTLIB:MSVCRTD")

chaste_generate_test_suites()

endif(ENABLE_CHASTE_TESTING)