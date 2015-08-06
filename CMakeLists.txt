cmake_minimum_required(VERSION 2.8.10)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include(EnsureOutOfSourceBuild)
Ensure_Out_Of_Source_Build()
include(CheckIncludeFile)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/overrides.cmake)

project(Chaste)

include(ChasteBuildTypes)
include(ChasteCompilerFlags)
include(ChasteMacros)

set (Chaste_VERSION_MAJOR 3)
set (Chaste_VERSION_MINOR 3)

set(Chaste_NUM_CPUS_TEST 1 CACHE STRING  "Number of cpus to use when running tests.")

option(Chaste_VERBOSE OFF "Provide extra information when building Chaste")
if (Chaste_VERBOSE AND MSVC)
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "/VERBOSE")
    set(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS} "/VERBOSE")
endif()

option(Chaste_USE_VTK "Compile Chaste with VTK support" ON)
option(Chaste_USE_CVODE "Compile Chaste with CVODE support" ON)

if (NOT (WIN32 OR CYGWIN))
    option(Chaste_USE_XERCES "Compile Chaste with XERCES and XSD support" ON)
else()
    option(Chaste_USE_XERCES "Compile Chaste with XERCES and XSD support" OFF)
endif()

option(Chaste_USE_PETSC_PARMETIS "Prefer to compile Chaste with PARMETIS library used by PETSc" ON)
option(Chaste_USE_PETSC_HDF5 "Prefer to compile Chaste with HDF5 library used by PETSc" ON)
option(Chaste_LINK_LIBRARIES_WITH_UNRESOLVED_SYMBOLS "All Chaste libraries will have unresolved symbols. User must link in all dependencies (third party and Chaste internal)" OFF)
option(Chaste_COVERAGE "Build Chaste with coverage information" OFF)
if (Chaste_COVERAGE AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(SEND_ERROR "CMAKE_BUILD_TYPE must be \"Debug\" to build Chaste with coverage (Chaste_COVERAGE=ON)")
endif()

option(Chaste_UPDATE_PROVENANCE "Update build timestamp. Disable to prevent re-linking of all Chaste libraries" ON)

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
option(ENABLE_CHASTE_PROJECT_TESTING "Enable Chaste Project Testing" ON)

#Some Chaste-specific #defines
add_definitions(-DCHASTE_CMAKE)

if(WIN32 OR CYGWIN)
    add_definitions(-D_WIN64 -D_AMD64_)
    # Ensure M_PI is always defined in cmath
    add_definitions(-D_USE_MATH_DEFINES)
endif(WIN32 OR CYGWIN)


################################
#  FIND THIRD PARTY LIBRARIES  #
################################

set(CHASTE_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
set(CHASTE_LINK_LIBRARIES "")

#Locate Python
find_package(PythonInterp REQUIRED)

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
if (Chaste_USE_VTK)
    #set(VTK_FIND_QUIETLY ON)
    find_package(VTK REQUIRED)
    list(APPEND CHASTE_INCLUDES "${VTK_INCLUDE_DIRS}")
    list(APPEND CHASTE_LINK_LIBRARIES "${VTK_LIBRARIES}")
    add_definitions(-DCHASTE_VTK)
endif()

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


#find PETSc
find_package(PETSc REQUIRED COMPONENTS CXX)

list(APPEND CHASTE_INCLUDES "${PETSC_INCLUDES}")
list(APPEND CHASTE_LINK_LIBRARIES 
    "${PETSC_LIBRARIES}"
)

#find HDF5
if (Chaste_USE_PETSC_HDF5)
    set(ENV{HDF5_ROOT} "${PETSC_DIR}/${PETSC_ARCH}" "${PETSC_DIR}/externalpackages")
    find_program( HDF5_C_COMPILER_EXECUTABLE
        NAMES h5cc h5pcc
        PATHS ENV HDF5_ROOT
        PATH_SUFFIXES bin Bin
        NO_DEFAULT_PATH
        DOC "HDF5 Wrapper compiler.  Used only to detect HDF5 compile flags." )
    find_program( HDF5_CXX_COMPILER_EXECUTABLE
        NAMES h5c++ h5pc++
        PATHS ENV HDF5_ROOT
        PATH_SUFFIXES bin Bin
        NO_DEFAULT_PATH
        DOC "HDF5 C++ Wrapper compiler.  Used only to detect HDF5 compile flags." )
endif()
find_package(HDF5 MODULE REQUIRED)
list(APPEND CHASTE_INCLUDES "${HDF5_INCLUDE_DIRS}")
list(APPEND CHASTE_LINK_LIBRARIES 
    "${HDF5_LIBRARIES}"
)

#Locate MPI 
set(MPI_CXX_COMPILER ${PETSC_COMPILER})
find_program(MPIEXEC
  NAMES mpiexec mpirun lamexec srun
  PATHS ${PETSC_DIR}/${PETSC_ARCH}
  PATH_SUFFIXES bin
  DOC "Executable for running MPI programs."
  NO_DEFAULT_PATH)
find_package(MPI REQUIRED)
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}" )
SET( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}" )
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}" )
list(APPEND CHASTE_INCLUDES "${MPI_INCLUDES}")
list(APPEND CHASTE_LINK_LIBRARIES "${MPI_CXX_LIBRARIES}")

#find ParMETIS and METIS
if (Chaste_USE_PETSC_PARMETIS)
    set(PARMETIS_ROOT "${PETSC_DIR}/${PETSC_ARCH}" 
                      "${PETSC_DIR}/externalpackages"
                      )
endif()
find_package(ParMETIS REQUIRED)
list(APPEND CHASTE_LINK_LIBRARIES "${PARMETIS_LIBRARIES}")

#Locate Sundials
if (Chaste_USE_CVODE)
    if(BUILD_SHARED_LIBS)
	    set(SUNDIALS_USE_STATIC_LIBRARIES OFF)
    else()
	    set(SUNDIALS_USE_STATIC_LIBRARIES ON)
    endif()
    find_package(SUNDIALS REQUIRED)
    list(APPEND CHASTE_INCLUDES "${SUNDIALS_INCLUDE_DIRS}")
    list(APPEND CHASTE_LINK_LIBRARIES "${SUNDIALS_LIBRARIES}")
    add_definitions(-DCHASTE_CVODE)
    math(EXPR Chaste_SUNDIALS_VERSION "${SUNDIALS_VERSION_MAJOR}*10000 + ${SUNDIALS_VERSION_MINOR}*100 + ${SUNDIALS_VERSION_SUBMINOR}")
    add_definitions(-DCHASTE_SUNDIALS_VERSION=${Chaste_SUNDIALS_VERSION})
endif()


#Locate Xerces and XSD
if (Chaste_USE_XERCES)
	find_package(Xerces REQUIRED)
	find_package(XSD REQUIRED)
	list(APPEND CHASTE_INCLUDES "${XERCESC_INCLUDE}" "${XSD_INCLUDE_DIRS}")
	list(APPEND CHASTE_LINK_LIBRARIES "${XERCESC_LIBRARY}")
    add_definitions(-DCHASTE_XERCES)
endif()

set(CXXTEST_INCLUDES "${CMAKE_SOURCE_DIR}/cxxtest")

add_definitions(-DTRILIBRARY -DTETLIBRARY -DANSI_DECLARATORS)

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

if(ENABLE_CHASTE_TESTING)
    enable_testing()
    list(APPEND CMAKE_INCLUDE_PATH "${Chaste_SOURCE_DIR}/cxxtest")
    find_package(CxxTest)
endif()

###########################################
# SETUP AVAILABLE COMPONENTS AND PROJECTS #
###########################################

# List the available Chaste components
set(CHASTE_COMPONENTS global io linalg mesh ode pde continuum_mechanics cell_based crypt)
if (NOT (WIN32 OR CYGWIN))
    set(CHASTE_COMPONENTS ${CHASTE_COMPONENTS} lung heart)
endif()

# Find any projects
file(GLOB potential_dirs RELATIVE "${Chaste_SOURCE_DIR}/projects" "${Chaste_SOURCE_DIR}/projects/*")
set(Chaste_PROJECTS "")
foreach(potential_dir ${potential_dirs})
    if (IS_DIRECTORY "${Chaste_SOURCE_DIR}/projects/${potential_dir}")
        # test for CMakeLists.txt file
        if (EXISTS "${Chaste_SOURCE_DIR}/projects/${potential_dir}/CMakeLists.txt")
            list(APPEND Chaste_PROJECTS "${potential_dir}")
        else()
            message(WARNING "No CMakeLists.txt file found in project directory ${Chaste_SOURCE_DIR}/projects/${potential_dir}. This project will not be built")
        endif()
    endif ()
endforeach(potential_dir ${potential_dirs})

#######################################################
# SETUP COMPONENT DEPENDANCIES AND HEADER DIRECTORIES #
#######################################################

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
set(CHASTE_DEPENDS_core global io linalg mesh ode pde continuum_mechanics)

foreach(component ${CHASTE_COMPONENTS})
    set(CHASTE_${component}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${component}/src")
    header_dirs(${CHASTE_${component}_SOURCE_DIR} CHASTE_${component}_INCLUDE_DIRS)
endforeach(component)

set(Chaste_PYTHON_DIR "${Chaste_SOURCE_DIR}/python")

#####################################
#  SETUP CONFIG FOR IN-TREE BUILDS  #
#####################################

configure_file(${Chaste_SOURCE_DIR}/cmake/Config/ChasteConfig.cmake.in
    "${Chaste_BINARY_DIR}/ChasteConfig.cmake" @ONLY)

set(Chaste_DIR ${Chaste_BINARY_DIR})

####################
# BUILD COMPONENTS #
####################

foreach(component ${CHASTE_COMPONENTS})
    # Build each component as a project
    add_subdirectory(${component})
endforeach(component)

####################
#  BUILD MAIN APPS (TODO) #
####################


####################
#  BUILD PROJECTS  #
####################

foreach(projectName ${Chaste_PROJECTS})
    add_subdirectory(projects/${projectName})
endforeach(projectName)

########################################
# EXPORT CONFIG FOR OUT-OF-TREE BUILDS #
########################################

export(PACKAGE Chaste)

export(TARGETS ${CHASTE_COMPONENTS}
    FILE "${Chaste_BINARY_DIR}/ChasteTargets.cmake"
    #added in 2.8.12, need this?    EXPORT_LINK_INTERFACE_LIBRARIES 
    )

# Configure file for install dir
set(EXPORT_Chaste_PYTHON_DIR "\${Chaste_CMAKE_DIR}/python")
foreach(component ${CHASTE_COMPONENTS})
    set(EXPORT_CHASTE_${component}_INCLUDE_DIRS "")
    foreach(dir ${CHASTE_${component}_INCLUDE_DIRS})
        file(RELATIVE_PATH rel_dir "${Chaste_SOURCE_DIR}" "${dir}")
        list(APPEND EXPORT_CHASTE_${component}_INCLUDE_DIRS "\${Chaste_CMAKE_DIR}/include/${rel_dir}")
    endforeach()
endforeach(component)

configure_file(${Chaste_SOURCE_DIR}/cmake/Config/ChasteConfig.cmake.in
    "${Chaste_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ChasteConfig.cmake" @ONLY)

# Configure file for build dir
set(EXPORT_Chaste_PYTHON_DIR "${Chaste_PYTHON_DIR}")
foreach(component ${CHASTE_COMPONENTS})
    set(EXPORT_CHASTE_${component}_INCLUDE_DIRS ${CHASTE_${component}_INCLUDE_DIRS})
endforeach(component)

configure_file(${Chaste_SOURCE_DIR}/cmake/Config/ChasteConfig.cmake.in
    "${Chaste_BINARY_DIR}/ChasteConfig.cmake" @ONLY)


file(COPY ${Chaste_SOURCE_DIR}/cmake/Modules/ChasteMacros.cmake 
    DESTINATION ${Chaste_BINARY_DIR}/cmake/Modules)

install(FILES
    "${Chaste_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ChasteConfig.cmake"
    "${Chaste_BINARY_DIR}/ChasteTargets.cmake"
    "${Chaste_SOURCE_DIR}/cmake/Modules/ChasteMacros.cmake"
    DESTINATION .
     COMPONENT Config)

# Python folder
install(DIRECTORY ${Chaste_SOURCE_DIR}/python
    DESTINATION .
    COMPONENT Python)

####################
#    COVERAGE      #
####################


if (Chaste_COVERAGE)
    include(CodeCoverage)
    setup_target_for_coverage(coverage "${Chaste_SOURCE_DIR}/cmake/coverage_script.sh" coverage)
endif()


####################
#    PACKAGING     #
####################

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
