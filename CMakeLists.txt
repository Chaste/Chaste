# Copyright (c) 2005-2019, University of Oxford.
# All rights reserved.
# 
# University of Oxford means the Chancellor, Masters and Scholars of the
# University of Oxford, having an administrative office at Wellington
# Square, Oxford OX1 2JD, UK.
# 
# This file is part of Chaste.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * Neither the name of the University of Oxford nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 

cmake_minimum_required (VERSION 3.5.1)

# Set standard to C++14. Note that CMake < 3.6 does not know how to set the standard flag for the Intel compiler, so
# the ChasteCompilerFlags.cmake also sets the --std=c++14 flag, too.
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
include (EnsureOutOfSourceBuild)
ensure_out_of_source_build ()
include (CheckIncludeFile)

include (${CMAKE_CURRENT_SOURCE_DIR}/cmake/overrides.cmake)

message ("\n####################################")
message ("# Setting compilers and build type")
message ("####################################\n")

# Since CMake 3.0, version can be supplied to the project() command.  We set it directly, for users on older CMake
project (Chaste VERSION 2019.1)
set (Chaste_VERSION_MAJOR 2019 CACHE STRING "Chaste major version number")
set (Chaste_VERSION_MINOR 1 CACHE STRING "Chaste minor version number")

include (ChasteHostOperatingSystem)
include (ChasteBuildTypes)

option (Chaste_ERROR_ON_WARNING "Error on compiler warnings" ON)
option (Chaste_VERBOSE "Provide extra information when building Chaste" OFF)
if (Chaste_VERBOSE AND MSVC)
    set (CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "/VERBOSE")
    set (CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS} "/VERBOSE")
endif ()
if (UNIX)
    # Memory testing
    option (Chaste_MEMORY_TESTING "Run tests using valgrind for memory testing" OFF)
    set (Chaste_MEMORY_TESTING_OUTPUT_DIR ${Chaste_BINARY_DIR}/memtest)
    set (Chaste_MEMORY_TESTING_CPUS 1 CACHE STRING "Number of CPUs for memory testing (default 1)")
    if (Chaste_MEMORY_TESTING_CPUS GREATER 1)
        set (Chaste_MEMORY_TESTING ON)
    endif ()
    if (Chaste_MEMORY_TESTING)
        file (MAKE_DIRECTORY ${Chaste_MEMORY_TESTING_OUTPUT_DIR})
    endif ()
    # Coverage
    option (Chaste_COVERAGE "Build Chaste with coverage information" OFF)
    set (Chaste_COVERAGE_CPUS 1 CACHE STRING "Number of CPUs for coverage testing (default 1)")
    if (Chaste_COVERAGE_CPUS GREATER 1)
        set (Chaste_COVERAGE ON)
    endif ()
    # Clang tidy
    option (Chaste_CLANG_TIDY "Configure Chaste for Clang tidy" OFF)
    if (Chaste_CLANG_TIDY)
        # This creates a compile_commands.json in the build directory
        set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
        set (Chaste_CLANG_TIDY_OUTPUT_DIR ${Chaste_BINARY_DIR}/clang_tidy)
        file (MAKE_DIRECTORY ${Chaste_CLANG_TIDY_OUTPUT_DIR})
    endif ()
    # Profiling
    option (Chaste_PROFILE_GPERFTOOLS "Compile Chaste with gperftools profiling support" OFF)
    set (Chaste_PROFILE_GPERFTOOLS_CPUS 1 CACHE STRING "Number of CPUs for GPERFTOOLS profiling (default 1)")
    if (Chaste_PROFILE_GPERFTOOLS_CPUS GREATER 1)
        set (Chaste_PROFILE_GPERFTOOLS ON)
    endif ()
    option (Chaste_PROFILE_GPROF "Compile Chaste with gprof profiling support" OFF)
    set (Chaste_PROFILE_GPROF_CPUS 1 CACHE STRING "Number of CPUs for GPROF profiling (default 1)")
    if (Chaste_PROFILE_GPROF_CPUS GREATER 1)
        set (Chaste_PROFILE_GPROF ON)
    endif ()
    if (Chaste_PROFILE_GPROF AND Chaste_PROFILE_GPERFTOOLS)
        message (WARNING "Both Chaste_PROFILE_GPROF and Chaste_PROFILE_GPERFTOOLS = ON, using latter by default")
        set (Chaste_PROFILE_GPROF OFF)
    endif ()
    if (Chaste_PROFILE_GPROF OR Chaste_PROFILE_GPERFTOOLS)
        set (Chaste_PROFILE_OUTPUT_DIR ${Chaste_BINARY_DIR}/profile)
        file (MAKE_DIRECTORY ${Chaste_PROFILE_OUTPUT_DIR})
    endif ()
    if (Chaste_COVERAGE OR Chaste_MEMORY_TESTING OR Chaste_PROFILE_GPROF OR Chaste_PROFILE_GPERFTOOLS)
        set (CMAKE_BUILD_TYPE "Debug")
    endif ()
endif ()
include (ChasteCompilerFlags)

include (ChasteMacros)

set (Chaste_NUM_CPUS_TEST 1 CACHE STRING "Number of cpus to use when running tests.")


option (Chaste_USE_VTK "Compile Chaste with VTK support" ON)
option (Chaste_USE_CVODE "Compile Chaste with CVODE support" ON)

if (NOT (WIN32 OR CYGWIN))
    option (Chaste_USE_XERCES "Compile Chaste with XERCES and XSD support" ON)
else ()
    option (Chaste_USE_XERCES "Compile Chaste with XERCES and XSD support" OFF)
endif ()

option (Chaste_USE_PETSC_PARMETIS "Prefer to compile Chaste with PARMETIS library used by PETSc" ON)
option (Chaste_USE_PETSC_HDF5 "Prefer to compile Chaste with HDF5 library used by PETSc" ON)

option (Chaste_UPDATE_PROVENANCE "Update build timestamp. Disable to prevent re-linking of all Chaste libraries" ON)

option (RUN_TESTS "This option simply runs Chaste tests. You should also set the test family." OFF)
set (TEST_FAMILY "Continuous" CACHE STRING "The name of the test family, e.g, Continuous, Failing, Nightly, Parallel etc.")
set (TestPackTypes "Continuous;Failing;Nightly;Parallel;Production;Weekly;Profile;ProfileAssembly;ExtraSimulations")

if (RUN_TESTS)
    list (FIND TestPackTypes ${TEST_FAMILY} found)
    if (found EQUAL -1)
        message (FATAL_ERROR "Test family ${TEST_FAMILY} does not exist. Must be one of ${TestPackTypes}. Aborting.")
    else (found EQUAL -1)
        #get date and time, to append to test result filename
        execute_process (COMMAND cmd /c echo %DATE% %TIME%
                         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                         OUTPUT_VARIABLE date_time
                         )
        string (REGEX REPLACE "[:/. \n]" "_" date_time "${date_time}")
        # Note: set 6 minute (360s) timeout for each test
        execute_process (COMMAND ctest -C Debug --output-on-failure -O ${TEST_FAMILY}TestOutputs_${date_time}.txt --timeout 360 -L ${TEST_FAMILY}
                         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                         OUTPUT_VARIABLE t_out
                         RESULT_VARIABLE t_res
                         ERROR_VARIABLE t_err
                         )
        message ("STDOUT______________\n${t_out}")
        message ("STDERR______________\n${t_err}")
    endif (found EQUAL -1)
endif (RUN_TESTS)

if (WIN32 OR CYGWIN)
    option (Chaste_AUTO_INSTALL_DEPS
            "Set whether we will automatically download and install Chaste dependences (windows-only option). ON by default"
            ON)
endif ()


#Set whether this is a statically or dynamically-linked build
if (WIN32 OR CYGWIN)
    option (BUILD_SHARED_LIBS
            "Set whether we are set whether to generate dynamic-linked libraries. OFF by default"
            OFF)
else ()
    option (BUILD_SHARED_LIBS
            "Set whether we are set whether to generate dynamic-linked libraries. ON by default"
            ON)
endif ()

option (Chaste_ENABLE_TESTING "Enable Chaste Testing" ON)
option (Chaste_INSTALL_TESTS "Install Chaste Testing infrastructure" OFF)

set (Chaste_PYCML_EXTRA_ARGS "" CACHE STRING "Add any extra arguments for PyCML here")

#Some Chaste-specific #defines
add_definitions (-DCHASTE_CMAKE)

if (WIN32 OR CYGWIN)
    add_definitions (-D_WIN64 -D_AMD64_)
    # Ensure M_PI is always defined in cmath
    add_definitions (-D_USE_MATH_DEFINES)
endif (WIN32 OR CYGWIN)


################################
#  FIND THIRD PARTY LIBRARIES  #
################################

message ("\n####################################")
message ("# Finding libraries")
message ("####################################\n")

set (Chaste_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
set (Chaste_LINK_LIBRARIES "")

# Check prereqs
if (Chaste_COVERAGE)
    find_program (GCOV_PATH gcov)
    find_program (LCOV_PATH lcov)
    find_program (GENHTML_PATH genhtml)
    find_program (GCOVR_PATH gcovr PATHS ${CMAKE_SOURCE_DIR}/tests)

    if (NOT GCOV_PATH)
        message (FATAL_ERROR "gcov not found! Aborting...")
    endif () # NOT GCOV_PATH

    if (NOT LCOV_PATH)
        message (FATAL_ERROR "lcov not found! Aborting...")
    endif () # NOT LCOV_PATH

    if (NOT GENHTML_PATH)
        message (FATAL_ERROR "genhtml not found! Aborting...")
    endif () # NOT GENHTML_PATH
endif ()


################################
####  Gperftools
################################
if (Chaste_PROFILE_GPERFTOOLS)
    find_package (Gperftools REQUIRED)
    list (APPEND Chaste_LINK_LIBRARIES "-Wl,--no-as-needed ${GPERFTOOLS_PROFILER} -Wl,--as-needed")
endif ()

#Gprof executable
if (Chaste_PROFILE_GPROF)
    find_file (GPROF_EXECUTABLE gprof DOC "Gprof code profiler")
    if (GPROF_EXECUTABLE)
        message (STATUS "Found GPROF_EXECUTABLE = ${GPROF_EXECUTABLE}")
    else ()
        message (FATAL_ERROR "Cannot find Gprof executable for profiling")
    endif ()
endif ()


################################
####  Find Valgrind
################################
if (Chaste_MEMORY_TESTING)
    find_package (Valgrind REQUIRED)
    if (NOT (${Chaste_NUM_CPUS_TEST} EQUAL 1))
        message (WARNING "Memory testing is turned on (Chaste_MEMORY_TESTING=ON) but you are trying to setup testing in parallel. Please set Chaste_NUM_CPUS_TEST to 1 for memory testing. Configuration and generation of tests will continue, but all tests will be run in serial")
    endif ()
endif ()


################################
####  Find Python
################################
find_package (PythonInterp REQUIRED)
if (NOT ${PYTHON_VERSION_MAJOR} VERSION_EQUAL 2)
    message (FATAL_ERROR "Python 2 is required; found v${PYTHON_VERSION_STRING} instead. Specify interpreter with -DPYTHON_EXECUTABLE=/path/to/python2")
endif ()

if (Chaste_ENABLE_TESTING AND Chaste_ACCEPTANCE)
    find_package (TextTest)
    if (NOT TEXTTEST_FOUND)
        message (WARNING "Texttest not found, turning off acceptance tests")
    endif ()
endif ()

if (Chaste_AUTO_INSTALL_DEPS)
    set (Chaste_DEPS_ROOT_DIR "${Chaste_BINARY_DIR}/../install/third_party_libs" CACHE PATH "Root directory for installed third party libraries")
    file (GLOB children RELATIVE ${Chaste_DEPS_ROOT_DIR} ${Chaste_DEPS_ROOT_DIR}/*)
    foreach (subdir ${children})
        if (IS_DIRECTORY ${Chaste_DEPS_ROOT_DIR}/${subdir})
            if (${subdir} MATCHES ".*boost.*" AND NOT BOOST_ROOT)
                set (BOOST_ROOT "${Chaste_DEPS_ROOT_DIR}/${subdir}")
                #elseif (${subdir} MATCHES ".*petsc.*" AND NOT ENV{PETSC_DIR})
                #    set(ENV{PETSC_DIR} "${Chaste_DEPS_ROOT_DIR}/${subdir}")
                #    set(PETSC_ARCH "")
            elseif (${subdir} MATCHES ".*vtk.*" AND NOT VTK_DIR)
                set (VTK_DIR "${Chaste_DEPS_ROOT_DIR}/${subdir}/lib/vtk-5.8")
            elseif (${subdir} MATCHES ".*sundials.*" AND NOT ENV{SUNDIALS_ROOT})
                set (ENV{SUNDIALS_ROOT} "${Chaste_DEPS_ROOT_DIR}/${subdir}")
            elseif (${subdir} MATCHES ".*hdf5.*" AND NOT ENV{HDF5_ROOT})
                set (ENV{HDF5_ROOT} "${Chaste_DEPS_ROOT_DIR}/${subdir}")
            endif ()
        endif ()
    endforeach ()
endif ()


################################
####  Find VTK
################################
if (Chaste_USE_VTK)
    # First pass will identify the version number
    find_package (VTK REQUIRED)
    message (STATUS "VTK version: ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")

    # Depending on the VTK version, we need different components. 5.x require the following
    if (VTK_MAJOR_VERSION EQUAL 5)
        find_package (VTK COMPONENTS vtkIO vtkCommon vtkGraphics z REQUIRED)
        # 6.0 and 6.1 require the following
    elseif ((VTK_MAJOR_VERSION EQUAL 6) AND (VTK_MINOR_VERSION LESS 2))
        find_package (VTK COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersCore vtkFiltersGeneral vtkFiltersGeneric
                      vtkFiltersGeometry vtkFiltersModeling vtkFiltersSources vtkIOCore vtkIOGeometry
                      vtkIOLegacy vtkIOXML REQUIRED)
        # 6.2 and up, up to and including 7.1, require the following
    else ()
        find_package (VTK COMPONENTS vtkCommonCore vtkCommonDataModel vtkFiltersCore vtkFiltersGeneral vtkFiltersGeneric
                      vtkFiltersGeometry vtkFiltersModeling vtkFiltersSources vtkIOCore vtkIOGeometry
                      vtkIOLegacy vtkIOParallelXML vtkIOXML REQUIRED)
    endif ()

    add_definitions (-DCHASTE_VTK)
endif ()


################################
####  Find Boost
################################
add_definitions (-DBOOST_ALL_NO_LIB)
set (Boost_USE_STATIC_RUNTIME ON)
if (BUILD_SHARED_LIBS)
    set (Boost_USE_STATIC_LIBS OFF)
else ()
    set (Boost_USE_STATIC_LIBS ON)
endif ()
find_package (Boost COMPONENTS filesystem system serialization program_options REQUIRED)
list (APPEND Chaste_INCLUDES "${Boost_INCLUDE_DIR}")
list (APPEND Chaste_LINK_LIBRARIES "${Boost_LIBRARIES}")


################################
####  Find PETSc
################################
find_package (PETSc REQUIRED)
list (APPEND Chaste_LINK_LIBRARIES "${PETSC_LIBRARIES}")


################################
####  Find HDF5
################################

# This alters the order that CMake looks for h5pcc h5cc (parallel h5pcc first)
set (HDF5_PREFER_PARALLEL TRUE)

# If we're using PETSc's HDF5 (ON by default, for more reliable MPI compatibility), we override HDF5_ROOT accordingly
if (Chaste_USE_PETSC_HDF5)
    if (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/bin/h5diff")
        set (ENV{HDF5_ROOT} "${PETSC_DIR}/${PETSC_ARCH}")
    else ()
        message (STATUS "No HDF5 found in ${PETSC_DIR}/${PETSC_ARCH}: looking for alternative HDF5 instead")
    endif ()
endif ()

# Find the HDF5_C_COMPILER_EXECUTABLE to help the find_package call (else the system HDF5 will always be found)
find_program (HDF5_C_COMPILER_EXECUTABLE
              NAMES h5pcc h5cc
              HINTS ENV HDF5_ROOT
              PATH_SUFFIXES bin Bin
              DOC "HDF5 Wrapper compiler.  Used only to detect HDF5 compile flags.")
find_package (HDF5 REQUIRED)

# Chaste requires a parallel HDF5 library
if (NOT HDF5_IS_PARALLEL)
    message (SEND_ERROR "Hdf5 library found was not build with --enable-parallel, use the HDF5_ROOT environment variable to specify a parallel hdf5 library. Note: include dirs are ${HDF5_INCLUDE_DIRS} and libraries are ${HDF5_LIBRARIES}")
endif ()
list (APPEND Chaste_INCLUDES "${HDF5_INCLUDE_DIRS}")

# for some reason the CMake in Homebrew returns a NOTFOUND entry in the Hdf5 libraries
foreach (dir ${HDF5_LIBRARIES})
    if (dir)
        list (APPEND Chaste_LINK_LIBRARIES "${dir}")
    endif ()
endforeach ()

# put petsc includes after hdf5 includes or else the hdf5 headers
# in the petsc include dir will clobber those chosen on FindHdf5.cmake
list (APPEND Chaste_INCLUDES "${PETSC_INCLUDES}")


################################
####  Find MPI
################################
if (PETSC_COMPILER)
    string (REPLACE mpicc mpicxx MPI_CXX_COMPILER ${PETSC_COMPILER})
endif ()

if (PETSC_MPIEXEC)
    if (IS_ABSOLUTE ${PETSC_MPIEXEC})
        set (MPIEXEC ${PETSC_MPIEXEC})  # MPIEXEC deprecated in newer FindMPI.cmake (#2967)
        set (MPIEXEC_EXECUTABLE ${PETSC_MPIEXEC})
    else ()
        find_program (MPIEXEC ${PETSC_MPIEXEC})
    endif ()
endif ()
find_package (MPI REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
list (APPEND Chaste_INCLUDES "${MPI_CXX_INCLUDE_PATH}")

if (Chaste_MEMORY_TESTING)
    get_filename_component (openmpi_supp_path ${MPIEXEC} PATH)
    get_filename_component (openmpi_supp_path ${openmpi_supp_path} PATH)
    set (openmpi_supp_path "${openmpi_supp_path}/share/openmpi/openmpi-valgrind.supp")
    set (Chaste_MEMORY_TESTING_SUPPS "--suppressions=${Chaste_SOURCE_DIR}/chaste.supp")
    set (Chaste_MEMORY_TESTING_SUPPS "${Chaste_MEMORY_TESTING_SUPPS} --suppressions=${Chaste_SOURCE_DIR}/chaste-lucid.supp")
    if (EXISTS ${openmpi_supp_path})
        set (Chaste_MEMORY_TESTING_SUPPS "${Chaste_MEMORY_TESTING_SUPPS} --suppressions=${openmpi_supp_path}")
    endif ()
endif ()

message (STATUS "Found MPIEXEC_EXECUTABLE: ${MPIEXEC_EXECUTABLE}")
message (STATUS "Found MPIEXEC: ${MPIEXEC}")
message (STATUS "Found MPI_CXX_HEADER_DIR: ${MPI_CXX_HEADER_DIR}")
message (STATUS "Found MPI_C_COMPILER: ${MPI_C_COMPILER}")


################################
####  Find ParMETIS and METIS
################################
# If we're using PETSc's ParMETIS (ON by default, for more reliable MPI compatibility), override PARMETIS_ROOT accordingly
if (Chaste_USE_PETSC_PARMETIS)
    if (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include/parmetis.h")
        set (PARMETIS_ROOT "${PETSC_DIR}/${PETSC_ARCH}")
    else ()
        message (STATUS "No ParMETIS found in ${PETSC_DIR}/${PETSC_ARCH}: looking for system ParMETIS instead")
    endif ()
endif ()
find_package (ParMETIS REQUIRED)
list (APPEND Chaste_LINK_LIBRARIES "${PARMETIS_LIBRARIES}")


################################
####  Find Sundials
################################
if (Chaste_USE_CVODE)
    if (BUILD_SHARED_LIBS)
        set (SUNDIALS_USE_STATIC_LIBRARIES OFF)
    else ()
        set (SUNDIALS_USE_STATIC_LIBRARIES ON)
    endif ()
    find_package (SUNDIALS COMPONENTS sundials_cvode sundials_nvecserial REQUIRED)
    list (APPEND Chaste_INCLUDES "${SUNDIALS_INCLUDE_DIRS}")
    #chaste_add_libraries(Chaste_LINK_LIBRARIES Chaste_THIRD_PARTY_STATIC_LIBRARIES Chaste_LINK_LIBRARIES)
    list (APPEND Chaste_LINK_LIBRARIES "${SUNDIALS_LIBRARIES}")
    add_definitions (-DCHASTE_CVODE)
    math (EXPR Chaste_SUNDIALS_VERSION "${SUNDIALS_VERSION_MAJOR}*10000 + ${SUNDIALS_VERSION_MINOR}*100 + ${SUNDIALS_VERSION_SUBMINOR}")
    add_definitions (-DCHASTE_SUNDIALS_VERSION=${Chaste_SUNDIALS_VERSION})
endif ()


# ParMETIS and Sundials might need MPI, so add MPI libraries after these
#chaste_add_libraries(MPI_CXX_LIBRARIES Chaste_THIRD_PARTY_STATIC_LIBRARIES Chaste_LINK_LIBRARIES)
list (APPEND Chaste_LINK_LIBRARIES "${MPI_CXX_LIBRARIES}")

# make sure VTK libraries added after HDF5 so VTK's link with HDF5 isn't used
if (Chaste_USE_VTK)
    list (APPEND Chaste_INCLUDES "${VTK_INCLUDE_DIRS}")
    list (APPEND Chaste_LINK_LIBRARIES "${VTK_LIBRARIES}")
endif ()


################################
####  Find Xerces and XSD
################################
if (Chaste_USE_XERCES)
    find_package (Xerces REQUIRED)
    find_package (XSD REQUIRED)
    list (APPEND Chaste_INCLUDES "${XERCESC_INCLUDE}" "${XSD_INCLUDE_DIRS}")
    list (APPEND Chaste_LINK_LIBRARIES "${XERCESC_LIBRARY}")
    add_definitions (-DCHASTE_XERCES)
endif ()

add_definitions (-DTRILIBRARY -DTETLIBRARY -DANSI_DECLARATORS)

set (CXXTEST_INCLUDES "${Chaste_BINARY_DIR}/cxxtest")
set (CXXTEST_PYTHON_TESTGEN_EXECUTABLE ${CXXTEST_INCLUDES}/cxxtestgen.py)

if (WIN32 OR CYGWIN)
    #MS Includes
    set (MS_MPI_INCLUDES "C:/MS_HPC_PACK_2012/Inc" CACHE PATH "Path to MS HPC Pack header files.")
    set (WINDOWS_SDK "C:/Program Files (x86)/Microsoft SDKs/Windows/v7.1A/Include" CACHE PATH "Path to Windows SDK headers.")
    set (WINDOWS_KITS "C:/Program Files (x86)/Windows Kits/8.0/Include" CACHE PATH "Path to Windows kits headers.")

    if (MSVC11)
        set (VS_11_INCLUDES "C:/Program Files (x86)/Microsoft Visual Studio 11.0/VC/include" CACHE PATH "You are compiling with MSVC 2012. Set Visual Studio 11 header files.")
        set (VS_INCLUDES "${VS_11_INCLUnstalled libvtk-java and libvtk5-qt4-devDES}")
    endif (MSVC11)

    if (MSVC10)
        set (VS_10_INCLUDES "C:/Program Files (x86)/Microsoft Visual Studio 10.0/VC/include" CACHE PATH "You are compiling with MSVC 2010. Set Visual Studio 10 header files.")
        set (VS_INCLUDES "${VS_10_INCLUDES}")
    endif (MSVC10)

    list (APPEND Chaste_INCLUDES "${WINDOWS_SDK}" "${VS_INCLUDES}" "${MS_MPI_INCLUDES}")

endif (WIN32 OR CYGWIN)

if (Chaste_ENABLE_TESTING)
    enable_testing ()
    list (APPEND CMAKE_INCLUDE_PATH "${Chaste_SOURCE_DIR}/cxxtest")
    find_package (CxxTest)
endif ()

###########################################
# SETUP AVAILABLE COMPONENTS AND PROJECTS #
###########################################

message ("\n####################################")
message ("# Configuring components")
message ("####################################\n")

# List the available Chaste components
set (Chaste_COMPONENTS global io linalg mesh ode pde continuum_mechanics cell_based crypt)
if (NOT (WIN32 OR CYGWIN))
    set (Chaste_COMPONENTS ${Chaste_COMPONENTS} lung heart)
endif ()

# Find any projects
file (GLOB potential_dirs RELATIVE "${Chaste_SOURCE_DIR}/projects" "${Chaste_SOURCE_DIR}/projects/*")
set (Chaste_PROJECTS "")
foreach (potential_dir ${potential_dirs})
    if (IS_DIRECTORY "${Chaste_SOURCE_DIR}/projects/${potential_dir}")
        # test for CMakeLists.txt file
        if (EXISTS "${Chaste_SOURCE_DIR}/projects/${potential_dir}/CMakeLists.txt")
            option (Chaste_ENABLE_project_${potential_dir} "Turn ${potential_dir} ON or OFF." ON)
            option (Chaste_ENABLE_project_${potential_dir}_INSTALL "Install ${potential_dir} along with standard Chaste libraries." OFF)
            file (READ "${Chaste_SOURCE_DIR}/projects/${potential_dir}/CMakeLists.txt" cmake_lists_file)
            if (Chaste_ENABLE_project_${potential_dir})
                if (cmake_lists_file MATCHES "\nfind_package\\(Chaste COMPONENTS.*project_.*\\)")
                    list (APPEND Chaste_PROJECTS "${potential_dir}")
                else ()
                    list (INSERT Chaste_PROJECTS 0 "${potential_dir}")
                endif ()
            endif ()
        else ()
            message (WARNING "No CMakeLists.txt file found in project directory ${Chaste_SOURCE_DIR}/projects/${potential_dir}. This project will not be built")
        endif ()
    endif ()
endforeach (potential_dir ${potential_dirs})

####################################
# setup tutorial generation target #
####################################

add_custom_target (tutorials)

#################################################
# setup test pack and component testing targets #
#################################################

set (Chaste_ALL_LIBRARIES "")
foreach (type ${TestPackTypes})
    add_custom_target (${type})
endforeach ()

foreach (component ${Chaste_COMPONENTS})
    # ${component} target will build library and tests
    add_custom_target (${component})
    # chaste_${component} is the actual library target
    list (APPEND Chaste_ALL_LIBRARIES chaste_${component})
endforeach ()

foreach (project ${Chaste_PROJECTS})
    # project_${project} target will build library and tests
    add_custom_target (project_${project})
    # chaste_project_${project} target is the actual library target
    if (Chaste_ENABLE_project_${project}_INSTALL)
        list (APPEND Chaste_ALL_LIBRARIES chaste_project_${project})
    endif ()
endforeach ()

# Targets for making core components and core component libraries, respectively
add_custom_target (core)
add_custom_target (chaste_core)

# Targets for making all components and all component libraries, respectively
add_custom_target (all_components)
add_custom_target (chaste_all_components)

#######################################################
# SETUP COMPONENT DEPENDANCIES AND HEADER DIRECTORIES #
#######################################################

# Specify which other components each depends on.
# This information is used to set up CMake dependencies, include search paths and libraries to link against.
set (Chaste_DEPENDS_global "")
set (Chaste_DEPENDS_io global)
set (Chaste_DEPENDS_linalg global)
set (Chaste_DEPENDS_mesh linalg global)
set (Chaste_DEPENDS_ode linalg io global)
set (Chaste_DEPENDS_pde ode mesh linalg io global)
set (Chaste_DEPENDS_cell_based pde ode mesh linalg io global)
set (Chaste_DEPENDS_crypt cell_based pde ode mesh linalg io global)
set (Chaste_DEPENDS_continuum_mechanics pde ode mesh linalg io global)
set (Chaste_DEPENDS_heart ${Chaste_DEPENDS_continuum_mechanics} continuum_mechanics)
set (Chaste_DEPENDS_lung ${Chaste_DEPENDS_continuum_mechanics} continuum_mechanics)

# These custom targets depend on other components but are not components themselves
set (Chaste_DEPENDS_core global io linalg mesh ode pde continuum_mechanics)
set (Chaste_DEPENDS_all_components ${Chaste_COMPONENTS})

# Set dependencies for custom targets
add_dependencies (core ${Chaste_DEPENDS_core})
add_dependencies (all_components ${Chaste_DEPENDS_all_components})

foreach (component ${Chaste_COMPONENTS})
    set (Chaste_${component}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${component}/src")
    header_dirs (${Chaste_${component}_SOURCE_DIR} Chaste_${component}_SOURCE_INCLUDE_DIRS)
    cellml_dirs (${Chaste_${component}_SOURCE_DIR} Chaste_${component}_CELLML_DIRS)
endforeach (component)

foreach (project ${Chaste_PROJECTS})
    set (Chaste_project_${project}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/projects/${project}/src")
    header_dirs (${Chaste_project_${project}_SOURCE_DIR} Chaste_project_${project}_SOURCE_INCLUDE_DIRS)
    cellml_dirs (${Chaste_project_${project}_SOURCE_DIR} Chaste_project_${project}_CELLML_DIRS)
endforeach (project)

# for now, assume only heart has xsd files
set (Chaste_heart_XSD_DIRS ${Chaste_heart_SOURCE_DIR}/io)

# generate include dirs
foreach (component ${Chaste_COMPONENTS})
    set (Chaste_${component}_INCLUDE_DIRS ${Chaste_${component}_SOURCE_INCLUDE_DIRS})
    foreach (dir ${Chaste_${component}_CELLML_DIRS} ${Chaste_${component}_XSD_DIRS})
        file (RELATIVE_PATH rel_dir "${Chaste_SOURCE_DIR}" "${dir}")
        list (APPEND Chaste_${component}_INCLUDE_DIRS ${Chaste_BINARY_DIR}/${rel_dir})
    endforeach ()
endforeach ()

foreach (project ${Chaste_PROJECTS})
    set (Chaste_project_${project}_INCLUDE_DIRS ${Chaste_project_${project}_SOURCE_INCLUDE_DIRS})
    foreach (dir ${Chaste_project_${project}_CELLML_DIRS})
        file (RELATIVE_PATH rel_dir "${Chaste_SOURCE_DIR}" "${dir}")
        list (APPEND Chaste_project_${project}_INCLUDE_DIRS ${Chaste_BINARY_DIR}/${rel_dir})
    endforeach ()
endforeach ()

set (Chaste_PYTHON_DIR "${Chaste_SOURCE_DIR}/python")

#####################################
#  SETUP CONFIG FOR IN-TREE BUILDS  #
#####################################

configure_file (${Chaste_SOURCE_DIR}/cmake/Config/ChasteConfig.cmake.in
                "${Chaste_BINARY_DIR}/ChasteConfig.cmake" @ONLY)

set (Chaste_DIR ${Chaste_BINARY_DIR})

####################
# BUILD COMPONENTS #
####################

foreach (component ${Chaste_COMPONENTS})
    # Build each component as a project
    add_subdirectory (${component})
endforeach (component)

# After components are built, add dependencies for the custom chaste_core and chaste_all_components libraries
foreach (component ${Chaste_DEPENDS_core})
    add_dependencies (chaste_core chaste_${component})
endforeach ()

foreach (component ${Chaste_DEPENDS_all_components})
    add_dependencies (chaste_all_components chaste_${component})
endforeach ()

####################
#  BUILD MAIN APPS #
####################

#note, heart not supported in windows
if (NOT (WIN32 OR CYGWIN))
    add_subdirectory (apps)
endif ()

####################
# RUN PYTHON TESTS #
####################

if (Chaste_ENABLE_TESTING)
    add_subdirectory (python/test)
endif ()

####################
#  BUILD PROJECTS  #
####################

foreach (projectName ${Chaste_PROJECTS})
    add_subdirectory (projects/${projectName})
endforeach (projectName)

########################################
# EXPORT CONFIG FOR OUT-OF-TREE BUILDS #
########################################

export (PACKAGE Chaste)

export (TARGETS ${Chaste_ALL_LIBRARIES}
        FILE "${Chaste_BINARY_DIR}/ChasteTargets.cmake"
        #added in 2.8.12, need this?    EXPORT_LINK_INTERFACE_LIBRARIES
        )

# Configure file for install dir
set (EXPORT_Chaste_PYTHON_DIR "\${Chaste_CMAKE_DIR}/python")
foreach (component ${Chaste_COMPONENTS})
    set (EXPORT_Chaste_${component}_INCLUDE_DIRS "")
    foreach (dir ${Chaste_${component}_SOURCE_INCLUDE_DIRS} ${Chaste_${component}_CELLML_DIRS} ${Chaste_${component}_XSD_DIRS})
        file (RELATIVE_PATH rel_dir "${Chaste_SOURCE_DIR}/${component}/src" "${dir}")
        list (APPEND EXPORT_Chaste_${component}_INCLUDE_DIRS "\${Chaste_CMAKE_DIR}/../../include/chaste/${component}/${rel_dir}")
    endforeach ()
endforeach (component)


set (EXPORT_Chaste_ADDITIONAL_PROJECTS "")
foreach (project ${Chaste_PROJECTS})
    if (Chaste_ENABLE_project_${project}_INSTALL)
        set (EXPORT_Chaste_ADDITIONAL_PROJECTS "${EXPORT_Chaste_ADDITIONAL_PROJECTS}\n\tset(Chaste_project_${project}_INCLUDE_DIRS ${Chaste_project_${project}_INCLUDE_DIRS})")
    endif ()
endforeach ()

configure_file (${Chaste_SOURCE_DIR}/cmake/Config/ChasteConfig.cmake.in
                "${Chaste_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ChasteConfig.cmake" @ONLY)

# Configure file for build dir
set (EXPORT_Chaste_PYTHON_DIR "${Chaste_PYTHON_DIR}")
foreach (component ${Chaste_COMPONENTS})
    set (EXPORT_Chaste_${component}_INCLUDE_DIRS ${Chaste_${component}_INCLUDE_DIRS})
endforeach (component)

configure_file (${Chaste_SOURCE_DIR}/cmake/Config/ChasteConfig.cmake.in
                "${Chaste_BINARY_DIR}/ChasteConfig.cmake" @ONLY)

file (COPY ${Chaste_SOURCE_DIR}/cmake/Modules/ChasteMacros.cmake
      DESTINATION ${Chaste_BINARY_DIR}/cmake/Modules)

file (COPY ${Chaste_SOURCE_DIR}/cmake/Modules/ChasteRunTestAndPostProcess.cmake
      DESTINATION ${Chaste_BINARY_DIR}/cmake/Modules)

file (COPY ${Chaste_SOURCE_DIR}/python
      DESTINATION ${Chaste_BINARY_DIR})

file (COPY ${Chaste_SOURCE_DIR}/cxxtest
      DESTINATION ${Chaste_BINARY_DIR})


install (FILES
         "${Chaste_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ChasteConfig.cmake"
         #    "${Chaste_BINARY_DIR}/ChasteTargets.cmake"
         DESTINATION lib/chaste
         COMPONENT Config)

install (FILES
         "${Chaste_BINARY_DIR}/cmake/Modules/ChasteMacros.cmake"
         "${Chaste_BINARY_DIR}/cmake/Modules/ChasteRunTestAndPostProcess.cmake"
         DESTINATION lib/chaste/cmake/Modules
         COMPONENT Config)

# CxxTest folder
install (DIRECTORY ${Chaste_BINARY_DIR}/cxxtest
         DESTINATION lib/chaste
         USE_SOURCE_PERMISSIONS
         COMPONENT CxxTest)

# Python folder
install (DIRECTORY ${Chaste_BINARY_DIR}/python
         DESTINATION lib/chaste
         USE_SOURCE_PERMISSIONS
         COMPONENT Python)

# libraries
install (EXPORT chaste-targets
         DESTINATION lib/chaste
         FILE ChasteTargets.cmake
         #EXPORT_LINK_INTERFACE_LIBRARIES
         COMPONENT ChasteTargets
         )


####################
#    COVERAGE      #
####################

if (Chaste_COVERAGE)

    # Run ctest with a low priority (+15)
    set (NICE_COMMAND nice)
    set (NICENESS -n15)
    set (CTEST_COMMAND ctest)
    set (CTEST_COMMAND ctest)
    set (_outputname coverage)
    add_custom_target (coverage

                       # Cleanup lcov
                       ${LCOV_PATH} --directory . --zerocounters

                       # Run tests
                       COMMAND ${NICE_COMMAND} ${NICENESS} ${CTEST_COMMAND} "-j${Chaste_COVERAGE_CPUS}" "-L" "Continuous" "--output-on-failure"
                       COMMAND ${NICE_COMMAND} ${NICENESS} ${CTEST_COMMAND} "-L" "Parallel" "--output-on-failure"

                       # Capturing lcov counters and generating report
                       COMMAND ${LCOV_PATH} --config-file ${Chaste_SOURCE_DIR}/cmake/Config/lcovrc --directory . --capture --output-file ${_outputname}.info
                       COMMAND ${LCOV_PATH} --config-file ${Chaste_SOURCE_DIR}/cmake/Config/lcovrc --remove ${_outputname}.info /home/bob/petsc* /usr/* */fortests/* */test/* */3rdparty/* */global/src/random/* Debug/* cxxtest/* --output-file ${_outputname}.info.cleaned
                       set (_page_title "\"Chaste Coverage Results for commit ${Chaste_REVISION}\"")
                       COMMAND ${GENHTML_PATH} --title "${_page_title}" --config-file ${Chaste_SOURCE_DIR}/cmake/Config/lcovrc --no-function-coverage -o ${_outputname} ${_outputname}.info.cleaned
                       COMMAND ${CMAKE_COMMAND} -E remove ${_outputname}.info ${_outputname}.info.cleaned
                       COMMAND ${PYTHON_EXECUTABLE} "${Chaste_SOURCE_DIR}/cmake/process_coverage_output.py" "${_outputname}"

                       DEPENDS Continuous Parallel
                       WORKING_DIRECTORY ${Chaste_BINARY_DIR}
                       COMMENT "Resetting code coverage counters to zero.\nProcessing code coverage counters and generating report."
                       VERBATIM
                       )

    # Show info where to find the report
    add_custom_command (TARGET coverage POST_BUILD
                        COMMAND ;
                        COMMENT "Open ./${_outputname}/index.html in your browser to view the coverage report."
                        )
endif ()


####################
#    CLANG TIDY    #
####################

if (Chaste_CLANG_TIDY)

    # Find the executable; highest version number for preference
    find_program (
            CLANG_TIDY
            NAMES clang-tidy-5.0 clang-tidy-4.0 clang-tidy-3.9 clang-tidy-3.8
            HINTS $ENV{CLANG_TIDY_DIR} /usr/bin
            DOC "Clang tidy executable."
    )

    # Get a list of all source files we wish to process with clang_tidy
    set (LIST_OF_TRANSLATION_UNITS "")
    foreach (component ${Chaste_COMPONENTS})
        # Every translation unit including HeartConfig.hpp will #include "ChasteParameters_3_4.hpp" which doesn't exist
        # until compilation. \todo: how can we get around this, ideally without having to compile!
        if (NOT (${component} MATCHES heart))
            # Get the cpp files recursively from this component
            file (GLOB_RECURSE potential_srcs RELATIVE ${Chaste_SOURCE_DIR} ${component}/src/*.cpp)
            # Check them for specific excludes
            foreach (potential_src ${potential_srcs})
                if (NOT (${potential_src} MATCHES 3rdparty))
                    list (APPEND LIST_OF_TRANSLATION_UNITS ${potential_src})
                endif ()
            endforeach ()
        endif (NOT (${component} MATCHES heart))
    endforeach (component)

    # Add custom commands: one for each src file to process
    foreach (src_file ${LIST_OF_TRANSLATION_UNITS})
        # Generate output file name from source file path
        string (REPLACE "/" "_" out_file_name ${src_file})
        set (out_file "${Chaste_BINARY_DIR}/clang_tidy/${out_file_name}")

        # Generate path to input file
        set (infile "${Chaste_SOURCE_DIR}/${src_file}")

        # Custom command to do the processing
        add_custom_command (
                OUTPUT "${out_file}"
                COMMAND ${CLANG_TIDY} -p=${Chaste_BINARY_DIR} ${Chaste_SOURCE_DIR}/${src_file} > ${out_file}
                COMMENT "Producing clang_tidy report: ${out_file}"
        )

        # Finally set out_files for dependencies
        set (out_files ${out_files} "${out_file}")
    endforeach (src_file)

    # Create a target that will run clang_tidy on each target
    add_custom_target (clang_tidy ALL DEPENDS ${out_files})
endif ()

####################
#     Doxygen      #
####################

add_custom_target (doxygen
                   COMMAND ${PYTHON_EXECUTABLE} "${Chaste_SOURCE_DIR}/cmake/run-doxygen.py" "${Chaste_SOURCE_DIR}"
                   "${Chaste_BINARY_DIR}/doxygen" "${Chaste_REVISION}"
                   WORKING_DIRECTORY ${Chaste_BINARY_DIR}
                   COMMENT "Generating Doxygen documentation"
                   VERBATIM)

add_custom_target (doxygen_coverage
                   COMMAND ${PYTHON_EXECUTABLE} "${Chaste_SOURCE_DIR}/cmake/run-doxygen.py" "${Chaste_SOURCE_DIR}"
                   "${Chaste_BINARY_DIR}/doxygen_coverage" "${Chaste_REVISION}" "True"
                   WORKING_DIRECTORY ${Chaste_BINARY_DIR}
                   COMMENT "Checking Doxygen coverage"
                   VERBATIM)

####################
#  MEMORY TESTING  #
####################

if (Chaste_MEMORY_TESTING)
    # Run ctest with a low priority (+15)
    set (NICE_COMMAND nice)
    set (NICENESS -n15)
    set (CTEST_COMMAND ctest)
    add_custom_target (memtest
                       COMMAND ${NICE_COMMAND} ${NICENESS} ${CTEST_COMMAND} "-j${Chaste_MEMORY_TESTING_CPUS}" "-L" Continuous "--output-on-failure"
                       COMMAND ${PYTHON_EXECUTABLE} "${Chaste_SOURCE_DIR}/cmake/process_valgrind_output.py" "${Chaste_MEMORY_TESTING_OUTPUT_DIR}"
                       DEPENDS Continuous
                       WORKING_DIRECTORY ${Chaste_BINARY_DIR}
                       VERBATIM)
endif ()

####################
#     Profiling    #
####################

if (Chaste_PROFILE_GPROF OR Chaste_PROFILE_GPERFTOOLS)
    if (Chaste_PROFILE_GPERFTOOLS)
        set (extension svg)
        set (NUM_CPUS ${Chaste_PROFILE_GPERFTOOLS_CPUS})
    else ()
        set (extension gmon)
        set (NUM_CPUS ${Chaste_PROFILE_GPROF_CPUS})
    endif ()

    # Run ctest with a neutral priority (0): buildbot doesn't have permission to set higher priorities
    set (NICE_COMMAND nice)
    set (NICENESS -n0)
    set (CTEST_COMMAND ctest)
    add_custom_target (profile
                       COMMAND ${NICE_COMMAND} ${NICENESS} ${CTEST_COMMAND} "-j${NUM_CPUS}" "-L" "^Profile_" "--output-on-failure"
                       COMMAND ${PYTHON_EXECUTABLE} "${Chaste_SOURCE_DIR}/cmake/process_profile.py" "${Chaste_PROFILE_OUTPUT_DIR}" ${extension}
                       DEPENDS Profile
                       WORKING_DIRECTORY ${Chaste_BINARY_DIR}
                       VERBATIM)
endif ()

##############################
#     Infrastructure Tests   #
##############################


add_custom_target (infrastructure
                   COMMAND ${PYTHON_EXECUTABLE} ${Chaste_SOURCE_DIR}/python/infra/CheckForCopyrights.py
                   COMMAND ${PYTHON_EXECUTABLE} ${Chaste_SOURCE_DIR}/python/infra/CheckForDuplicateFileNames.py
                   COMMAND ${PYTHON_EXECUTABLE} ${Chaste_SOURCE_DIR}/python/infra/CheckForOrphanedTests.py
                   COMMAND ${PYTHON_EXECUTABLE} ${Chaste_SOURCE_DIR}/python/infra/CheckSchemas.py

                   WORKING_DIRECTORY ${Chaste_SOURCE_DIR}
                   VERBATIM)


####################
#    PACKAGING     #
####################

set (CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
set (CPACK_PACKAGE_VENDOR "Computational Biology Group - Computer Science - University of Oxford")
set (CPACK_PACKAGE_CONTACT "Chaste Team <chaste-admin@maillist.ox.ac.uk>")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "Chaste (Cancer, Heart and Soft Tissue Environment).")
set (CPACK_PACKAGE_DESCRIPTION "
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
set (CPACK_PACKAGE_VERSION_MAJOR "${Chaste_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${Chaste_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${chaste_revision}")


set (CPACK_COMPONENT_Config_GROUP "config")
set (CPACK_COMPONENT_Config_DESCRIPTION "configuration files")
set (CPACK_COMPONENT_Config_DISPLAY_NAME "Config files")

set (CPACK_COMPONENT_Python_GROUP "python")
set (CPACK_COMPONENT_Python_DESCRIPTION "Support python files for Chaste project")
set (CPACK_COMPONENT_Python_DISPLAY_NAME "Python")

set (CPACK_COMPONENT_CxxTest_GROUP "cxxtest")
set (CPACK_COMPONENT_CxxTest_DESCRIPTION "Bundled CxxTest library (http://cxxtest.com/)")
set (CPACK_COMPONENT_CxxTest_DISPLAY_NAME "CxxTest")

foreach (component ${Chaste_COMPONENTS})
    set (CPACK_COMPONENT_${component}_libraries_GROUP "libraries")
    set (CPACK_COMPONENT_${component}_headers_GROUP "headers")
    set (CPACK_COMPONENT_${component}_headers_DESCRIPTION "C++ header files (.hpp) for Chaste component ${component}")

    if (BUILD_SHARED_LIBRARIES)
        set (CPACK_COMPONENT_${component}_libraries_DISPLAY_NAME "Dynamic Libraries")
        set (CPACK_COMPONENT_${component}_libraries_DESCRIPTION "Dynamic Libraries for Chaste component ${component}")
    else (BUILD_SHARED_LIBRARIES)
        set (CPACK_COMPONENT_${component}_libraries_DISPLAY_NAME "Static Libraries")
        set (CPACK_COMPONENT_${component}_libraries_DESCRIPTION "Static Libraries for Chaste component ${component}")
    endif (BUILD_SHARED_LIBRARIES)

    set (CPACK_COMPONENT_${component}_headers_DISPLAY_NAME "C++ Headers")

    if (Chaste_INSTALL_TESTS)
        set (CPACK_COMPONENT_${component}_tests_GROUP "tests")
        set (CPACK_COMPONENT_${component}_tests_DESCRIPTION "Test Suite (C++ headers and source files) for Chaste component ${component}")
        set (CPACK_COMPONENT_${component}_tests_DISPLAY_NAME "C++ Test Suite")
    endif ()

    foreach (depend_component Chaste_DEPENDS_${component})
        set (CPACK_COMPONENT_${component}_headers_DEPENDS ${CPACK_COMPONENT_${component}_headers_DEPENDS} ${depend_component}_headers)
        set (CPACK_COMPONENT_${component}_libraries_DEPENDS ${CPACK_COMPONENT_${component}_libraries_DEPENDS} ${depend_component}_libraries)
        #set(CPACK_COMPONENT_${component}_tests_DEPENDS ${CPACK_COMPONENT_${component}_tests_DEPENDS} ${depend_component}_tests)
    endforeach (depend_component Chaste_DEPENDS_${component})

endforeach (component ${Chaste_COMPONENTS})


set (CPACK_DEBIAN_PACKAGE_DEPENDS "cmake, g++, libopenmpi-dev, petsc-dev (>= 2.3.3-14), libhdf5-openmpi-dev, xsdcxx, libboost-serialization-dev, libboost-filesystem-dev, libparmetis-dev, libxerces-c3-dev, libsundials-serial-dev, libvtk5-dev, python-lxml, python-rdflib")
set (CPACK_DEBIAN_PACKAGE_RECOMMENDS "valgrind, libfltk1.1")
set (CPACK_DEBIAN_PACKAGE_SUGGESTS "libgoogle-perftools-dev, doxygen, graphviz, eclipse-cdt, gnuplot")
set (CPACK_DEBIAN_PACKAGE_PRIORITY extra)
set (CPACK_DEBIAN_PACKAGE_SECTION science)

# This must always be last!
include (CPack)
