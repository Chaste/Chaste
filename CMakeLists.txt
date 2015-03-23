cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(CheckIncludeFile)
include(ChasteMacros)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/overrides.cmake)

project(Chaste)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/definitions.cmake)


option(TEST_MPIEXEC_ENABLED OFF "Whether to run tests in parallel using mpiexec, or sequentially as normal executables.")

option(RUN_TESTS OFF "This option simply runs Chaste tests. You should also set the test family.")
set(TEST_FAMILY "Continuous" CACHE STRING "The name of the test family, e.g, Continuous, Failing, Nightly, Parallel etc.")
set(TestPackTypes "Continuous;Failing;Nightly;Parallel;Production;ProfileAssembly;Profile")

if(RUN_TESTS)
	list(FIND TestPackTypes ${TEST_FAMILY} found)
	if(found EQUAL -1)
		message(FATAL_ERROR "Test family ${TEST_FAMILY} does not exist. Must be one of ${TestPackTypes}. Aborting.")
	else(found EQUAL -1)
		#get date and time, to append to test result filename
		execute_process(COMMAND cmd /c echo %DATE% %TIME%
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
			OUTPUT_VARIABLE date_time
		)
		string(REGEX REPLACE "[:/. \n]" "_" date_time "${date_time}")
		# Note: set 6 minute (360s) timeout for each test
		execute_process(COMMAND ctest -C Debug --output-on-failure -O ${TEST_FAMILY}TestOutputs_${date_time}.txt --timeout 360 -L ${TEST_FAMILY}
		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
		OUTPUT_VARIABLE t_out
		RESULT_VARIABLE t_res
		ERROR_VARIABLE t_err
		)
		message("STDOUT______________\n${t_out}")
		message("STDERR______________\n${t_err}")
	endif(found EQUAL -1)
endif(RUN_TESTS)

if (WIN32 OR CYGWIN)
	option(CHASTE_AUTO_INSTALL_DEPS
		"Set whether we will automatically download and install Chaste dependences (windows-only option). ON by default"
		ON)
endif()


#Set whether this is a statically or dynamically-linked build
if (WIN32 OR CYGWIN)
option(BUILD_SHARED_LIBS
    "Set whether we are set whether to generate dynamic-linked libraries. OFF by default"
    OFF)
else()
option(BUILD_SHARED_LIBS
    "Set whether we are set whether to generate dynamic-linked libraries. ON by default"
    ON)
endif()

option(ENABLE_CHASTE_TESTING "Enable Chaste Testing" ON)

#Some Chaste-specific #defines
add_definitions(-DCHASTE_CVODE -DCHASTE_SUNDIALS_VERSION=20500)
add_definitions(-DCHASTE_VTK)

option(CMAKE_LINK_DEPENDS_NO_SHARED
    "Whether to skip link dependencies on shared library files. ON by default"
    ON)


if(WIN32 OR CYGWIN)
    add_definitions(-D_WIN64 -D_AMD64_)
    # Ensure M_PI is always defined in cmath
    add_definitions(-D_USE_MATH_DEFINES)
endif(WIN32 OR CYGWIN)


set(CHASTE_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
set(CHASTE_LINK_LIBRARIES "")

#Locate Python
find_package(PythonInterp REQUIRED)

#Locate MPI
find_package(MPI REQUIRED)

if(CHASTE_AUTO_INSTALL_DEPS)
	set(CHASTE_DEPS_ROOT_DIR "${Chaste_BINARY_DIR}/../install/third_party_libs" CACHE PATH "Root directory for installed third party libraries")
	  FILE(GLOB children RELATIVE ${CHASTE_DEPS_ROOT_DIR} ${CHASTE_DEPS_ROOT_DIR}/*)
	foreach(subdir ${children})
		if(IS_DIRECTORY ${CHASTE_DEPS_ROOT_DIR}/${subdir})
			if (${subdir} MATCHES ".*boost.*" AND NOT BOOST_ROOT)
				set(BOOST_ROOT "${CHASTE_DEPS_ROOT_DIR}/${subdir}")
			#elseif (${subdir} MATCHES ".*petsc.*" AND NOT ENV{PETSC_DIR})
			#	set(ENV{PETSC_DIR} "${CHASTE_DEPS_ROOT_DIR}/${subdir}")
			#	set(PETSC_ARCH "")
			elseif (${subdir} MATCHES ".*vtk.*" AND NOT VTK_DIR)
				set(VTK_DIR "${CHASTE_DEPS_ROOT_DIR}/${subdir}/lib/vtk-5.8")
			elseif (${subdir} MATCHES ".*sundials.*" AND NOT ENV{SUNDIALS_ROOT})
				set(ENV{SUNDIALS_ROOT} "${CHASTE_DEPS_ROOT_DIR}/${subdir}")
			elseif (${subdir} MATCHES ".*hdf5.*" AND NOT ENV{HDF5_ROOT})
				set(ENV{HDF5_ROOT} "${CHASTE_DEPS_ROOT_DIR}/${subdir}")
			endif()
		endif()
	endforeach()
endif()

#Locate VTK
#set(VTK_FIND_QUIETLY ON)
find_package(VTK REQUIRED)
list(APPEND CHASTE_INCLUDES "${VTK_INCLUDE_DIRS}")
list(APPEND CHASTE_LINK_LIBRARIES "${VTK_LIBRARIES}")

#find Boost
add_definitions( -DBOOST_ALL_NO_LIB )
set( Boost_USE_STATIC_RUNTIME ON) 
if(BUILD_SHARED_LIBS)
	set( Boost_USE_STATIC_LIBS OFF) 
else()
	set( Boost_USE_STATIC_LIBS ON) 
endif()
find_package(Boost COMPONENTS filesystem system serialization REQUIRED)
list(APPEND CHASTE_INCLUDES "${Boost_INCLUDE_DIR}")
list(APPEND CHASTE_LINK_LIBRARIES "${Boost_LIBRARIES}")

#find HDF5
find_package(HDF5 REQUIRED)
list(APPEND CHASTE_INCLUDES "${HDF5_C_INCLUDE_DIR}")

#find PETSc
find_package(PETSc REQUIRED)
list(APPEND CHASTE_INCLUDES "${PETSC_INCLUDES}")
list(APPEND CHASTE_LINK_LIBRARIES 
    "${PETSC_LIBRARIES}"
    "${HDF5_LIBRARIES}"
)

if (NOT ("${PETSC_LIBRARIES}" MATCHES ".*parmetis.*"))
    #find ParMETIS
    find_package(ParMETIS REQUIRED)
    list(APPEND CHASTE_INCLUDES "${PARMETIS_INCLUDE_DIR}")
    list(APPEND CHASTE_LINK_LIBRARIES "${PARMETIS_LIBRARIES}")

    if (NOT ("${PETSC_LIBRARIES}" MATCHES ".*metis.*"))
        #find METIS
        find_package(METIS REQUIRED)
        list(APPEND CHASTE_INCLUDES "${METIS_INCLUDE_DIR}")
        list(APPEND CHASTE_LINK_LIBRARIES "${METIS_LIBRARIES}")
    endif()
endif()



#Locate Sundials
if(BUILD_SHARED_LIBS)
	set(SUNDIALS_USE_STATIC_LIBRARIES OFF)
else()
	set(SUNDIALS_USE_STATIC_LIBRARIES ON)
endif()
find_package(SUNDIALS REQUIRED)
list(APPEND CHASTE_INCLUDES "${SUNDIALS_INCLUDE_DIRS}")
list(APPEND CHASTE_LINK_LIBRARIES "${SUNDIALS_LIBRARIES}")


#Locate Xerces and XSD
if (NOT (WIN32 OR CYGWIN))
	find_package(Xerces REQUIRED)
	find_package(XSD REQUIRED)
	list(APPEND CHASTE_INCLUDES "${XERCESC_INCLUDE}" "${XSD_INCLUDE_DIRS}")
	list(APPEND CHASTE_LINK_LIBRARIES "${XERCESC_LIBRARY}")
endif()

set(CXXTEST_INCLUDES "${CMAKE_SOURCE_DIR}/cxxtest")

if (WIN32 OR CYGWIN)
	#MS Includes
	set(MS_MPI_INCLUDES "C:/MS_HPC_PACK_2012/Inc" CACHE PATH "Path to MS HPC Pack header files.")
	set(WINDOWS_SDK "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.1A/Include" CACHE PATH "Path to Windows SDK headers.")
	set(WINDOWS_KITS "C:/Program Files (x86)/Windows Kits/8.0/Include" CACHE PATH "Path to Windows kits headers.")

	if(MSVC11)
		set(VS_11_INCLUDES "C:/Program Files (x86)/Microsoft Visual Studio 11.0/VC/include" CACHE PATH "You are compiling with MSVC 2012. Set Visual Studio 11 header files.")
		set(VS_INCLUDES "${VS_11_INCLUnstalled libvtk-java and libvtk5-qt4-devDES}")
	endif(MSVC11)

	if(MSVC10)
		set(VS_10_INCLUDES "C:/Program Files (x86)/Microsoft Visual Studio 10.0/VC/include" CACHE PATH "You are compiling with MSVC 2010. Set Visual Studio 10 header files.")
	set(VS_INCLUDES "${VS_10_INCLUDES}")
endif(MSVC10)

list(APPEND CHASTE_INCLUDES "${WINDOWS_SDK}" "${VS_INCLUDES}" "${MS_MPI_INCLUDES}")

endif(WIN32 OR CYGWIN)

# List the available Chaste components
set(CHASTE_COMPONENTS global io linalg mesh ode pde continuum_mechanics cell_based crypt)
if (NOT (WIN32 OR CYGWIN))
    set(CHASTE_COMPONENTS ${CHASTE_COMPONENTS} lung heart)
endif()

if(ENABLE_CHASTE_TESTING)
    enable_testing()
    list(APPEND CMAKE_INCLUDE_PATH "${Chaste_SOURCE_DIR}/cxxtest")
    find_package(CxxTest)
    foreach(component ${CHASTE_COMPONENTS})
        option(ENABLE_${component}_TESTING "Generate the test infrastructure for the ${component} component" ON)
        mark_as_advanced(ENABLE_${component}_TESTING)
    endforeach()
endif()

include(ChasteBuildInfo)

include_directories("${CHASTE_INCLUDES}")

# Where generated source/header files are put
include_directories("${config_dir}")

# A macro to recursively find Chaste header locations
macro(HEADER_DIRS base_dir return_list)
    set(new_list "")
    set(dir_list "")
    #message("base dir = ${base_dir}")

    file(GLOB_RECURSE new_list ${base_dir}/*.hpp ${base_dir}/*.h)
    #message("new list = ${new_list}")
    foreach(file_path ${new_list})
        get_filename_component(dir_path ${file_path} PATH)
        set(dir_list ${dir_list} ${dir_path})
    endforeach()
    
    list(REMOVE_DUPLICATES dir_list)
    
    if(NOT ${PROJECT_NAME} MATCHES "^test")
        list(APPEND dir_list "${CHASTE_INCLUDES}")
    endif()
    

    #message("return list = ${return_list}")
    #message("dir list = ${dir_list}")

    set(${return_list} ${dir_list})
endmacro()


#Build the component Test infrastructure
option (ENABLE_CHASTE_TESTING "Set whether to generate the CHASTE testing infrastructure." ON)
set(CMAKE_INCLUDE_PATH "${CMAKE_INCLUDE_PATH}" "${CMAKE_SOURCE_DIR}/cxxtest")


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
        add_executable("${_testname}" "${_test_real_output_filename}" ${ARGN})

        if(TEST_MPIEXEC_ENABLED)
            #Note: "${MPIEXEC} /np 1 master : subordinate" means that we run one master process and n subordinate processes
            # on the local host with n+1 cores.
            # Here we are using the form ${MPIEXEC} /np 2 ${test}.
            # A figure-it-out-yourselfnstalled libvtk-java and libvtk5-qt4-dev form would be ${MPIEXEC} /np * ${test} which runs on all available cores
            # See http://technet.microsoft.com/en-us/library/cc947675%28v=ws.10%29.aspx
            # Note the underscore appended to the test name, to match with the RUN_TESTS block above, and ensure we don't
            # run more tests than intended!

            add_test(NAME "${_testname}" WORKING_DIRECTORY "${Chaste_SOURCE_DIR}/" COMMAND "${MPIEXEC}" /np 2 $<TARGET_FILE:${_testname}>)
        else()
            add_test(NAME "${_testname}" WORKING_DIRECTORY "${Chaste_SOURCE_DIR}/" COMMAND $<TARGET_FILE:${_testname}>)
        endif(TEST_MPIEXEC_ENABLED)
    endmacro(CHASTE_ADD_TEST)

  macro(CHASTE_GENERATE_TEST_NAME test outTestName)
      string(REGEX REPLACE "([a-zA-Z0-9_/]+)[.]hpp" "\\1" testName "${test}")
      string(REPLACE "/" ";" testPath "${testName}")
      list(LENGTH testPath pathLength)
      if(${pathLength} EQUAL 1)
        set(testName ${testPath})
        set(testPath "")
        set(${outTestName} ${testName})
      else()
        math(EXPR index "${pathLength} - 1")
        list(GET testPath ${index} testName)
        list(REMOVE_AT testPath ${index})
        string(REPLACE ";" "_" _testPath_ "${testPath}")
        string(REPLACE ";" "/" testPath "${testPath}")
        set(${outTestName} "${testName}_${_testPath_}_")
      endif()
  endmacro(CHASTE_GENERATE_TEST_NAME test outTestName)

  macro(CHASTE_DO_COMPONENT component)
    message("Configuring component ${component}")
    # Figure out include path
    include_directories("${CHASTE_${component}_INCLUDE_DIRS}")
    foreach(dep ${CHASTE_DEPENDS_${component}})
        include_directories("${CHASTE_${dep}_INCLUDE_DIRS}")
    endforeach(dep)
    # Make & install component library
    file(GLOB_RECURSE CHASTE_${component}_SOURCES 
			RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
			${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp 
			${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)

    add_library(${component} ${CHASTE_${component}_SOURCES} ${ARGN})
    #if (BUILD_SHARED_LIBS)
    	target_link_libraries(${component} ${CHASTE_LINK_LIBRARIES})
    #endif()
    install(TARGETS ${component} DESTINATION lib COMPONENT ${component}_libraries)

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
            if(MSVC)
			    set_target_properties(${appName} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT /IGNORE:4217 /IGNORE:4049")
            endif()
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

        if(MSVC)
            if(NOT HAS_OWN_LINKER_FLAGS)
                set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT")
            endif(NOT HAS_OWN_LINKER_FLAGS)

            #disable linker warnings 4217, 4049: locally-defined symbol imported in function ...
            set(LINKER_FLAGS "${LINKER_FLAGS} /IGNORE:4217 /IGNORE:4049")
            #message("Linker flags for project ${PROJECT_NAME} = ${LINKER_FLAGS}")
        endif(MSVC)

        foreach(type ${TestPackTypes})
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${type}TestPack.txt")
            file(STRINGS "${type}TestPack.txt" testpack)
            foreach(filename ${testpack})
                string(STRIP ${filename} filename)
                chaste_generate_test_name(${filename} "testName")
                set(targetName "${testName}Runner")
                if (NOT TARGET ${targetName})
                    chaste_add_test(${targetName} "${testName}.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/${filename}")
                    target_link_libraries(${targetName} ${COMPONENT_LIBRARIES} ${CHASTE_LINK_LIBRARIES})
                    set_target_properties(${targetName} PROPERTIES LINK_FLAGS "${LINKER_FLAGS}")
                    set_property(TEST ${targetName} PROPERTY LABELS ${component} ${type})
                    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${testName}.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/${filename}"
                            DESTINATION test/${component} COMPONENT  ${component}_tests)
                else()
                    get_property(myLabels TEST ${targetName} PROPERTY LABELS)
                    list(APPEND myLabels ${type})
                    set_property(TEST ${targetName} PROPERTY LABELS ${myLabels})
                endif()
            endforeach(filename ${testpack})
            endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${type}TestPack.txt")
        endforeach(type ${TestPackTypes})
    endif(ENABLE_CHASTE_TESTING)
  endmacro(CHASTE_DO_TEST_COMPONENT)


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
set(CHASTE_DEPENDS_heart ${CHASTE_DEPENDS_continuum_mechanics} continuum_mechanics)
set(CHASTE_DEPENDS_lung ${CHASTE_DEPENDS_continuum_mechanics} continuum_mechanics)

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
            COMPONENT ${component}_headers
            FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
            )
endforeach(component)

set(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
set(CPACK_PACKAGE_VENDOR "Computational Biology Group - Computer Science - University of Oxford")
set(CPACK_PACKAGE_CONTACT "Chaste Team <chaste-admin@maillist.ox.ac.uk>")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Chaste (Cancer, Heart and Soft Tissue Environment).")
set(CPACK_PACKAGE_DESCRIPTION  "
 Chaste is a general purpose simulation package aimed at multi-scale, 
 computationally demanding problems arising in biology and physiology. 
 Current functionality includes tissue and cell level electrophysiology, 
 discrete tissue modelling, and soft tissue modelling. The package is 
 being developed by a team mainly based in the Computational Biology Group 
 at Oxford University Computing Laboratory, and development draws on expertise
 from software engineering, high performance computing, mathematical modelling
 and scientific computing.
 .
 The main website for Chaste can be found at
 http://www.cs.ox.ac.uk/chaste
")
SET(CPACK_PACKAGE_VERSION_MAJOR "${Chaste_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${Chaste_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${chaste_revision}")


foreach(component ${CHASTE_COMPONENTS})
    set(CPACK_COMPONENT_${component}_libraries_GROUP "libraries")
    set(CPACK_COMPONENT_${component}_headers_GROUP "headers")
    set(CPACK_COMPONENT_${component}_tests_GROUP "tests")


    set(CPACK_COMPONENT_${component}_headers_GROUP "headers")
    set(CPACK_COMPONENT_${component}_tests_GROUP "tests")

    set(CPACK_COMPONENT_${component}_headers_DESCRIPTION "C++ header files (.hpp) for Chaste component ${component}")
    set(CPACK_COMPONENT_${component}_tests_DESCRIPTION "Test Suite (C++ headers and source files) for Chaste component ${component}")

    if (BUILD_SHARED_LIBRARIES)
        set(CPACK_COMPONENT_${component}_libraries_DISPLAY_NAME "Dynamic Libraries")
        set(CPACK_COMPONENT_${component}_libraries_DESCRIPTION "Dynamic Libraries for Chaste component ${component}")
    else(BUILD_SHARED_LIBRARIES)
        set(CPACK_COMPONENT_${component}_libraries_DISPLAY_NAME "Static Libraries")
        set(CPACK_COMPONENT_${component}_libraries_DESCRIPTION "Static Libraries for Chaste component ${component}")
    endif(BUILD_SHARED_LIBRARIES)

    set(CPACK_COMPONENT_${component}_headers_DISPLAY_NAME "C++ Headers")
    set(CPACK_COMPONENT_${component}_tests_DISPLAY_NAME "C++ Test Suite")

    foreach(depend_component CHASTE_DEPENDS_${component})
        set(CPACK_COMPONENT_${component}_headers_DEPENDS ${CPACK_COMPONENT_${component}_headers_DEPENDS} ${depend_component}_headers)
        set(CPACK_COMPONENT_${component}_libraries_DEPENDS ${CPACK_COMPONENT_${component}_libraries_DEPENDS} ${depend_component}_libraries)
        #set(CPACK_COMPONENT_${component}_tests_DEPENDS ${CPACK_COMPONENT_${component}_tests_DEPENDS} ${depend_component}_tests)
    endforeach(depend_component CHASTE_DEPENDS_${component})

endforeach(component ${CHASTE_COMPONENTS})



SET(CPACK_DEBIAN_PACKAGE_DEPENDS "cmake, g++, libopenmpi-dev, petsc-dev (>= 2.3.3-14), libhdf5-openmpi-dev, xsdcxx, libboost-serialization-dev, libboost-filesystem-dev, libparmetis-dev, libxerces-c2-dev, libsundials-serial-dev, libvtk5-dev, python-lxml, python-amara, python-rdflib")
set(CPACK_DEBIAN_PACKAGE_RECOMMENDS "valgrind, libfltk1.1")
set(CPACK_DEBIAN_PACKAGE_SUGGESTS "libgoogle-perftools-dev, doxygen, graphviz, eclipse-cdt, gnuplot")
set(CPACK_DEBIAN_PACKAGE_PRIORITY extra)
set(CPACK_DEBIAN_PACKAGE_SECTION science)



# This must always be last!
include(CPack)
