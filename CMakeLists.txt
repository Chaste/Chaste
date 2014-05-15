cmake_minimum_required(VERSION 2.8.10)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/overrides.cmake)

project(Chaste)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/definitions.cmake)

option(RUN_TESTS OFF "This option simply runs Chaste tests. You should also set the test family.")
set(TEST_FAMILY "Continuous" CACHE STRING "The name of the test family, e.g, Continuous, Failing, Nightly, Parallel etc.")

option(TEST_MPIEXEC_ENABLED OFF "Whether to run tests in parallel using mpiexec, or sequentially as normal executables.")

if(RUN_TESTS)
set(TestPackTypes "Continuous;Failing;Nightly;Parallel;Production;ProfileAssembly;Profile")
foreach(type ${TestPackTypes})
    set(result "")
    file(GLOB_RECURSE TEST_PACKS "${CMAKE_CURRENT_SOURCE_DIR}" ${type}TestPack.txt)
        foreach(testp ${TEST_PACKS})
            file(STRINGS "${testp}" testpack)
                foreach(s ${testpack})
                    # Note that we append an underscore to test names, to ensure that ctest's regexp
                    # match doesn't grab more tests than we want!  This needs to match what we do in
                    # the chaste_add_test macro too.
                    string(REGEX REPLACE "(.*/)?([a-zA-Z0-9_]+)[.]hpp" "\\2_" s2 "${s}")
                    string(REGEX MATCH ".*[.]py" match "${s2}")
                    if(NOT match)
                        set(result "${result}" "${s2}")
                    endif(NOT match)
                endforeach(s ${testpack})
        endforeach(testp ${TEST_PACKS})
        list(REMOVE_AT result 0)#remove the first empty string.
        list(SORT result)
        string(REPLACE ";" "\n" result "${result}")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/test_runner/${type}TestsToRun.txt" "${result}")
endforeach(type ${TestPackTypes})

list(FIND TestPackTypes ${TEST_FAMILY} found)
if(found EQUAL -1)
    message(FATAL_ERROR "Test family ${TEST_FAMILY} does not exist. Must be one of ${TestPackTypes}. Aborting.")
else(found EQUAL -1)
    file(STRINGS "${CMAKE_CURRENT_BINARY_DIR}/test_runner/${TEST_FAMILY}TestsToRun.txt" tests)
    string(REPLACE ";" "|" tests "${tests}")
    #get date and time, to append to test result filename
    execute_process(COMMAND cmd /c echo %DATE% %TIME%
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    OUTPUT_VARIABLE date_time
                    )
    string(REGEX REPLACE "[:/. \n]" "_" date_time "${date_time}")
    # Note: set 6 minute (360s) timeout for each test
    execute_process(COMMAND ctest -C Debug --output-on-failure -O ${TEST_FAMILY}TestOutputs_${date_time}.txt --timeout 360 -R ${tests}
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    OUTPUT_VARIABLE t_out
                    RESULT_VARIABLE t_res
                    ERROR_VARIABLE t_err
                    )
    message("STDOUT______________\n${t_out}")
    message("STDERR______________\n${t_err}")
endif(found EQUAL -1)

endif(RUN_TESTS)

include(AdditionalIncludes.cmake)

#Chaste source root
set(CHASTE_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}")

#Chaste build directory
set(CHASTE_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}" CACHE PATH "The root directory of Chaste build.")

#Chaste build directory name
set(CHASTE_BUILD_DIR_NAME "Debug" CACHE STRING "The name of Chaste build directory.")

#Chaste Build Type
set(CHASTE_BUILD_TYPE "Windows_Debug" CACHE PATH
    "The build type. Options are Windows_Debug and Windows_Release")

#Set whether this is a statically or dynamically-linked build
option(CHASTE_SHARED_LIBRARY OFF
    "Set whether we are set whether this is a statically or dynamically-linked build. OFF by default")

#Some Chaste-specific #defines
add_definitions(-DCHASTE_CVODE -DCHASTE_SUNDIALS_VERSION=20500)
add_definitions(-DCHASTE_VTK)

#Third-party libraries
set(THIRD_PARTY_LIBS "${CMAKE_CURRENT_SOURCE_DIR}/../external_libs" CACHE PATH "Location of external libraries used by Chaste.")

#Perhaps use finders to locate these libraries?
#set(BOOST_INCLUDES "${THIRD_PARTY_LIBS}/boost_1_53_0/include/boost-1_53")
set(PETSC_ARCH "WINDOWS_BUILD" CACHE STRING
    "PETSc Architecture. Usually determines which build of PETSc to use")
#set(PETSC_INCLUDES "${THIRD_PARTY_LIBS}/petsc-3.3-p6/include"
#"${THIRD_PARTY_LIBS}/petsc-3.3-p6/${PETSC_ARCH}/include")
#set(HDF5_INCLUDES "${THIRD_PARTY_LIBS}/hdf5_1_8_10_p1/include")
#set(PARMETIS_INCLUDES "${THIRD_PARTY_LIBS}/parmetis_4_0_2_p3/include")
#set(METIS_INCLUDES "${THIRD_PARTY_LIBS}/parmetis_4_0_2_p3/metis/include")
set(CXXTEST_INCLUDES "${CMAKE_SOURCE_DIR}/cxxtest")

#MS Includes
set(MS_MPI_INCLUDES "C:/MS_HPC_PACK_2012/Inc" CACHE PATH "Path to MS HPC Pack header files.")
set(WINDOWS_SDK "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.1A/Include" CACHE PATH "Path to Windows SDK headers.")
set(WINDOWS_KITS "C:/Program Files (x86)/Windows Kits/8.0/Include" CACHE PATH "Path to Windows kits headers.")

if(MSVC11)
  set(VS_11_INCLUDES "C:/Program Files (x86)/Microsoft Visual Studio 11.0/VC/include" CACHE PATH "You are compiling with MSVC 2012. Set Visual Studio 11 header files.")
  set(VS_INCLUDES "${VS_11_INCLUDES}")
endif(MSVC11)

if(MSVC10)
  set(VS_10_INCLUDES "C:/Program Files (x86)/Microsoft Visual Studio 10.0/VC/include" CACHE PATH "You are compiling with MSVC 2010. Set Visual Studio 10 header files.")
  set(VS_INCLUDES "${VS_10_INCLUDES}")
endif(MSVC10)

set(CHASTE_INCLUDES #These are now mainly set externally in AdditionalIncludes.cmake
    "${WINDOWS_SDK}" "${CXXTEST_INCLUDES}" "${CHASTE_INCLUDES}" "${CMAKE_CURRENT_SOURCE_DIR}"
    "${VS_INCLUDES}" "${MS_MPI_INCLUDES}"
    #"${BOOST_INCLUDES}" "${PETSC_INCLUDES}" "${HDF5_INCLUDES}"
    #"${PARMETIS_INCLUDES}" "${METIS_INCLUDES}"
)
include_directories("${CHASTE_INCLUDES}")

# Where generated source/header files are put
include_directories("${CHASTE_BUILD_ROOT}/src")

# A macro to recursively find Chaste header locations
macro(HEADER_DIRS base_dir return_list)
    file(GLOB_RECURSE new_list base_dir *.hpp *.h)
    set(dir_list "")
    foreach(file_path ${new_list})
        get_filename_component(dir_path ${file_path} PATH)
        set(dir_list ${dir_list} ${dir_path})
    endforeach()
    list(REMOVE_DUPLICATES dir_list)
    if(NOT ${PROJECT_NAME} MATCHES "^test")
        list(APPEND dir_list "${CHASTE_INCLUDES}")
    endif()
    set(${return_list} ${dir_list})
endmacro()

#Third-party libraries linked against chaste test executables
if(MSVC11)
    set(CHASTE_LINK_LIBRARIES
        GKlib
        metis
        parmetis
        petsc
        sundials_cvodes
        sundials_nvecserial
        f2cblas
        f2clapack
        libboost_filesystem-vc110-sgd-1_53
        libboost_serialization-vc110-sgd-1_53
        libboost_system-vc110-sgd-1_53
        msmpi
        hdf5d)
endif(MSVC11)

if(MSVC10)
    set(CHASTE_LINK_LIBRARIES
        GKlib
        metis
        parmetis
        petsc
        sundials_cvodes
        sundials_nvecserial
        f2cblas
        f2clapack
        libboost_filesystem-vc100-sgd-1_53
        libboost_serialization-vc100-sgd-1_53
        libboost_system-vc100-sgd-1_53
        msmpi
        hdf5d)
endif(MSVC10)

set(CHASTE_LINK_LIBRARIES ${CHASTE_LINK_LIBRARIES} vtkFiltering vtkIO vtkCommon vtksys vtkzlib vtkexpat vtkGraphics)

#This is for a 64-bit build. 32-bit may need to be manually set. Would be nice to automate
#These switches were required due to import of time.h, Windows.h and Winsock2.h in "mesh"
add_definitions(-D_WIN64 -D_AMD64_)

link_directories("${THIRD_PARTY_LIBS}/all_lib")

if(NOT CHASTE_SHARED_LIBRARY)
    add_definitions(-MTd)
else()
    add_definitions(-MDd)
endif()

# Ensure M_PI is always defined in cmath
add_definitions(-D_USE_MATH_DEFINES)

#Locate MPI
find_package(MPI)

#Build the component Test infrastructure
option (ENABLE_CHASTE_TESTING "Set whether to generate the CHASTE testing infrastructure." ON)
set(CMAKE_INCLUDE_PATH "${CMAKE_INCLUDE_PATH}" "${CMAKE_SOURCE_DIR}/cxxtest")

if(ENABLE_CHASTE_TESTING)
    enable_testing()
    find_package(CxxTest)
    #Selective configuration of which components to test
    include(cmake/ConfigureComponentTesting.cmake)

    #Chaste Testing Macro. The predefined cxxtest_add_test is not suitable because of little control over
    #the test's working directory
    macro(CHASTE_ADD_TEST _testname _test_output_filename)
        set(_test_real_output_filename "${CMAKE_CURRENT_BINARY_DIR}/${_test_output_filename}")
        add_custom_command(
            OUTPUT "${_test_real_output_filename}"
            DEPENDS ${ARGN}
            COMMAND ${CXXTEST_TESTGEN_INTERPRETER} ${CXXTEST_TESTGEN_EXECUTABLE} ${CXXTEST_TESTGEN_ARGS} -o "${_test_real_output_filename}" ${ARGN}
        )

        set_source_files_properties("${_test_real_output_filename}" PROPERTIES GENERATED true)
        add_executable(${_testname} "${_test_real_output_filename}" ${ARGN})
        if("${CHASTE_BUILD_TYPE}" STREQUAL "Windows_Debug")
            get_property(test_location TARGET ${_testname} PROPERTY LOCATION_Debug)
        else()
            get_property(test_location TARGET ${_testname} PROPERTY LOCATION_Release)
        endif()

        if(TEST_MPIEXEC_ENABLED)
            #Note: "${MPIEXEC} /np 1 master : subordinate" means that we run one master process and n subordinate processes
            # on the local host with n+1 cores.
            # Here we are using the form ${MPIEXEC} /np 2 ${test}.
            # A figure-it-out-yourself form would be ${MPIEXEC} /np * ${test} which runs on all available cores
            # See http://technet.microsoft.com/en-us/library/cc947675%28v=ws.10%29.aspx
            # Note the underscore appended to the test name, to match with the RUN_TESTS block above, and ensure we don't
            # run more tests than intended!
            add_test(NAME "${_testname}_" WORKING_DIRECTORY "${CHASTE_SOURCE_ROOT}/" COMMAND "${MPIEXEC}" /np 2 ${test_location})
        else()
            add_test(NAME "${_testname}_" WORKING_DIRECTORY "${CHASTE_SOURCE_ROOT}/" COMMAND ${test_location})
        endif(TEST_MPIEXEC_ENABLED)
    endmacro(CHASTE_ADD_TEST)

  macro(CHASTE_GENERATE_TEST_SUITES)
    if(NOT HAS_OWN_LINKER_FLAGS)
        set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT")
    endif(NOT HAS_OWN_LINKER_FLAGS)

    #disable linker warnings 4217, 4049: locally-defined symbol imported in function ...
    set(LINKER_FLAGS "${LINKER_FLAGS} /IGNORE:4217 /IGNORE:4049")
    #message("Linker flags for project ${PROJECT_NAME} = ${LINKER_FLAGS}")

    foreach(test ${TEST_SUITES})
      string(REGEX REPLACE "([a-zA-Z0-9_/]+)[.]hpp" "\\1" testName "${test}")
      string(REPLACE "/" ";" testPath "${testName}")
      list(LENGTH testPath pathLength)
      if(${pathLength} EQUAL 1)
        set(testName ${testPath})
        set(testPath "")
        chaste_add_test("${testName}Runner" "${testName}.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/${testName}.hpp")
        target_link_libraries("${testName}Runner" ${COMPONENT_LIBRARIES} ${CHASTE_LINK_LIBRARIES})
        set_target_properties("${testName}Runner" PROPERTIES LINK_FLAGS "${LINKER_FLAGS}")
      else()
        math(EXPR index "${pathLength} - 1")
        list(GET testPath ${index} testName)
        list(REMOVE_AT testPath ${index})
        string(REPLACE ";" "_" _testPath_ "${testPath}")
        string(REPLACE ";" "/" testPath "${testPath}")
        chaste_add_test("${testName}_${_testPath_}_Runner" "${testName}_${_testPath_}_.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/${testPath}/${testName}.hpp")
        target_link_libraries("${testName}_${_testPath_}_Runner" ${COMPONENT_LIBRARIES} ${CHASTE_LINK_LIBRARIES})
        set_target_properties("${testName}_${_testPath_}_Runner" PROPERTIES LINK_FLAGS "${LINKER_FLAGS}")
      endif()
    endforeach(test)
  endmacro(CHASTE_GENERATE_TEST_SUITES)

  macro(CHASTE_DO_COMPONENT component)
    message("Configuring component ${component}")
    # Figure out include path
    include_directories("${CHASTE_${component}_INCLUDE_DIRS}")
    foreach(dep ${CHASTE_DEPENDS_${component}})
        include_directories("${CHASTE_${dep}_INCLUDE_DIRS}")
    endforeach(dep)
    # Make & install component library
    file(GLOB_RECURSE CHASTE_${component}_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
    add_library(${component} STATIC ${CHASTE_${component}_SOURCES})
    install(TARGETS ${component} DESTINATION lib)
    # Do testing if requested
    if(ENABLE_${component}_TESTING)
        add_subdirectory(test)
    endif()
	# Build applications if present
	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/apps")
		message("Building apps in ${CMAKE_CURRENT_SOURCE_DIR} comp ${component}")
		file(GLOB CHASTE_${component}_APPS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/apps/src/*.cpp)
		foreach(app ${CHASTE_${component}_APPS})
			string(REGEX REPLACE ".*/([a-zA-Z0-9_]+)[.]cpp" "\\1" appName "${app}")
			message("Building app ${app} named ${appName}.")
			message("Libraries=${component} ${CHASTE_DEPENDS_${component}} ${CHASTE_LINK_LIBRARIES}")
			add_executable(${appName} ${app})
			target_link_libraries(${appName} ${component} ${CHASTE_DEPENDS_${component}} ${CHASTE_LINK_LIBRARIES})
			set_target_properties(${appName} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT /IGNORE:4217 /IGNORE:4049")
			get_target_property(incDirs ${appName} INCLUDE_DIRECTORIES)
			list(APPEND incDirs "${CXXTEST_INCLUDES}")
			foreach(dep ${CHASTE_DEPENDS_${component}})
				list(APPEND incDirs "${CHASTE_${dep}_INCLUDE_DIRS}")
			endforeach(dep)
			set_target_properties(${appName} PROPERTIES INCLUDE_DIRECTORIES "${incDirs}")
		endforeach(app)
	endif()
  endmacro(CHASTE_DO_COMPONENT)

  macro(CHASTE_DO_TEST_COMPONENT component)
    if(ENABLE_CHASTE_TESTING)
        message("Configuring tests for ${component}")
        # Figure out include path for tests
        header_dirs("${CMAKE_CURRENT_SOURCE_DIR}" CHASTE_${component}_TEST_DIRS)
        include_directories("${CHASTE_${component}_TEST_DIRS}" "${CHASTE_${component}_INCLUDE_DIRS}" "${CXXTEST_INCLUDES}")
        foreach(dep ${CHASTE_DEPENDS_${component}})
            include_directories("${CHASTE_${dep}_INCLUDE_DIRS}")
        endforeach(dep)
        # Make test library if sources exist
        set(COMPONENT_LIBRARIES ${component})
        file(GLOB_RECURSE test_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
        if(test_sources)
            add_library(test${component} STATIC ${test_sources})
            set(COMPONENT_LIBRARIES ${COMPONENT_LIBRARIES} test${component})
        endif()
        set(COMPONENT_LIBRARIES ${COMPONENT_LIBRARIES} ${CHASTE_DEPENDS_${component}})
        # Generate test suites
        file(GLOB_RECURSE TEST_SUITES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} Test*.hpp)
        chaste_generate_test_suites()
    endif(ENABLE_CHASTE_TESTING)
  endmacro(CHASTE_DO_TEST_COMPONENT)

endif(ENABLE_CHASTE_TESTING)


#The generated timekeeper.cpp code below keeps track of build timestamp.
#It is built and executed prior to starting a build and prints the timestamp
#in a given format. This timestamp is used by Version.cpp, which is also auto-generated.
file(WRITE "${CMAKE_BINARY_DIR}/timekeeper.cpp"
"#include <iostream>
#include <ctime>
int main( )
{
   time_t now = time(0);
   tm* loc_time = localtime(&now);
   char buffer[80];
   strftime(buffer, 80, \"%a, %d %b %Y %H:%M:%S +0000\", loc_time);
   std::cout << buffer;
   return 0;
}
")

add_executable(timekeeper "${CMAKE_BINARY_DIR}/timekeeper.cpp")
execute_process(COMMAND "${CMAKE_CURRENT_BINARY_DIR}/${CHASTE_BUILD_DIR_NAME}/timekeeper"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    OUTPUT_FILE build_timestamp)

# List the available Chaste components
set(CHASTE_COMPONENTS global io linalg mesh ode pde continuum_mechanics cell_based crypt NottinghamSTI)

# Specify which other components each depends on.
# This information is used to set up CMake dependencies, include search paths and libraries to link against.
set(CHASTE_DEPENDS_global "")
set(CHASTE_DEPENDS_io global)
set(CHASTE_DEPENDS_linalg global)
set(CHASTE_DEPENDS_mesh linalg global)
set(CHASTE_DEPENDS_ode linalg io global)
set(CHASTE_DEPENDS_pde ode mesh linalg io global)
set(CHASTE_DEPENDS_cell_based pde ode mesh linalg io global)
set(CHASTE_DEPENDS_crypt cell_based pde ode mesh linalg io global)
set(CHASTE_DEPENDS_continuum_mechanics pde ode mesh linalg io global)
set(CHASTE_DEPENDS_NottinghamSTI continuum_mechanics pde ode mesh linalg io global)

foreach(component ${CHASTE_COMPONENTS})
    # Set up the component source & include directories
    #string(TOUPPER ${component} c_upper)
    set(CHASTE_${component}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${component}/src")
    header_dirs(${CHASTE_${component}_SOURCE_DIR} CHASTE_${component}_INCLUDE_DIRS)
endforeach(component)

foreach(component ${CHASTE_COMPONENTS})
    # Build each component as a project
    add_subdirectory(${component})
    if(CHASTE_DEPENDS_${component})
        add_dependencies(${component} ${CHASTE_DEPENDS_${component}})
    endif()

    # Install Chaste headers for Chaste developers
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${component}/src/"
            DESTINATION include/${component}
            FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
endforeach(component)
