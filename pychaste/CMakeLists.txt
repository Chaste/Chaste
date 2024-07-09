# Copyright (c) 2005-2024, University of Oxford.
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


find_package(Chaste COMPONENTS ${Chaste_DEPENDS_pychaste})

chaste_do_component(pychaste)

################################
# Compiler options
################################
add_compile_options(-Wno-unused-local-typedefs)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # https://stackoverflow.com/questions/25365160/boostmultiprecisionfloat128-and-c11
    add_compile_options(-fext-numeric-literals)
endif()

################################
# PyChaste non-wrapper C++ code
################################
set(PYCHASTE_NON_WRAPPER_INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cell_based/
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ode/
    ${CMAKE_CURRENT_SOURCE_DIR}/src/visualization/)

################################
# Wrapper generation targets
################################
add_custom_target(generate_pychaste_wrappers)
add_custom_command(TARGET generate_pychaste_wrappers
    COMMAND . ${chaste_python3_venv}/activate && 
        ${chaste_python3_venv}/cppwg ${CMAKE_SOURCE_DIR}
            -w ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/wrappers
            -p ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/package_info.yaml
            -i ${PYCHASTE_NON_WRAPPER_INCLUDE_DIRS}
               ${Chaste_INCLUDE_DIRS}
               ${Chaste_THIRD_PARTY_INCLUDE_DIRS}
            --std c++17 &&
        deactivate
)

add_custom_target(generate_pychaste_wrappers_diff)
add_custom_command(TARGET generate_pychaste_wrappers_diff
    COMMAND . ${chaste_python3_venv}/activate && 
        ${chaste_python3_venv}/cppwg ${CMAKE_SOURCE_DIR}
            -w ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/wrappers_diff
            -p ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/package_info.yaml
            -i ${PYCHASTE_NON_WRAPPER_INCLUDE_DIRS}
               ${Chaste_INCLUDE_DIRS}
               ${Chaste_THIRD_PARTY_INCLUDE_DIRS}
            --std c++17 &&
        deactivate
)

################################
# Build python modules from wrappers
################################
include_directories(${PYCHASTE_NON_WRAPPER_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dynamic/wrappers)
include_directories(${PYTHON3_INCLUDE_DIRS})

# Copy python package structure to build directory, ignoring existing shared libraries
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/python/
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/
     PATTERN "*${CMAKE_SHARED_LIBRARY_SUFFIX}" EXCLUDE)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test/python/
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/test/)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/doc/
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/python/doc/)

# List of pychaste modules to build into shared libraries
set (PYCHASTE_PYTHON_MODULES
        core
        ode
        pde
        mesh
        cell_based
        visualization)

# Locations to put each module library
set (PYCHASTE_PYTHON_MODULE_LOCATIONS
        ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/core/
        ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/ode/
        ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/pde/
        ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/mesh/
        ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/cell_based/
        ${CMAKE_CURRENT_BINARY_DIR}/python/chaste/visualization/)

# Target for building all module shared libraries
add_custom_target(chaste_pychaste_modules)

# The module library name must be the same as the pybind11 module name
# defined in the main wrapper e.g. `dynamic/wrappers/ode/ode.main.cpp`
# defines a pybind11 module `_pychaste_ode`. By convention,
# the name starts with an underscore. The usual 'lib' prefix is disabled.
set(PYCHASTE_PYTHON_MODULE_PREFIX "_pychaste_")

# Create a shared library target for each module
list(LENGTH PYCHASTE_PYTHON_MODULES max_idx)
math(EXPR max_idx "${max_idx} - 1")

foreach(idx RANGE ${max_idx})
    list(GET PYCHASTE_PYTHON_MODULES ${idx} module)
    list(GET PYCHASTE_PYTHON_MODULE_LOCATIONS ${idx} module_dir)

    # Glob the module's wrapper code from the `dynamic` directory
    file(GLOB module_sources ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/wrappers/${module}/*.cpp)

    set(module_library_name ${PYCHASTE_PYTHON_MODULE_PREFIX}${module})
    add_library(${module_library_name} SHARED ${module_sources})

    set_target_properties(${module_library_name}
        PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${module_dir}
            PREFIX "${PYTHON_MODULE_PREFIX}"
            SUFFIX "${CMAKE_SHARED_LIBRARY_SUFFIX}")

    # The order here is important - pybind11 and python come first
    target_link_libraries(${module_library_name} 
        pybind11::module
        ${Python3_LIBRARIES}
        ${Chaste_THIRD_PARTY_LIBRARIES}
        Chaste_COMMON_DEPS
        ${Chaste_LIBRARIES}
        chaste_pychaste)

    add_dependencies(chaste_pychaste_modules ${module_library_name})
endforeach()
