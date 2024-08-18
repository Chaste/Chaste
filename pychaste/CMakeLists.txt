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

if(NOT Chaste_ENABLE_PYCHASTE)
    return()
endif()

find_package(Chaste COMPONENTS ${Chaste_DEPENDS_pychaste})
chaste_do_component(pychaste)

set(PYCHASTE_MODULES cell_based core mesh ode pde visualization)

##################################################
# Python wrapper generation
##################################################
set(PYCHASTE_WRAPPER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/wrappers)

# Add a target for re-generating wrappers manually
set(PYCHASTE_CPPWG_COMMAND
    . ${chaste_python3_venv}/activate &&
    ${chaste_python3_venv}/cppwg ${CMAKE_SOURCE_DIR}
        -w ${PYCHASTE_WRAPPER_DIR}
        -p ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/config.yaml
        -i ${Chaste_INCLUDE_DIRS}
           ${Chaste_pychaste_INCLUDE_DIRS}
           ${Chaste_THIRD_PARTY_INCLUDE_DIRS}
        --std c++17 &&
    deactivate
)
add_custom_target(pychaste_wrappers)
add_custom_command(TARGET pychaste_wrappers COMMAND ${PYCHASTE_CPPWG_COMMAND})

# Generate wrappers if they are missing
foreach(module ${PYCHASTE_MODULES})
    if(NOT EXISTS ${PYCHASTE_WRAPPER_DIR}/${module}/${module}.main.cpp)
        message(STATUS "Generating PyChaste wrappers")
        string(JOIN " " PYCHASTE_CPPWG_CMDSTR ${PYCHASTE_CPPWG_COMMAND})
        execute_process(COMMAND sh -c "${PYCHASTE_CPPWG_CMDSTR}")
        break()
    endif()
endforeach()

header_dirs(${PYCHASTE_WRAPPER_DIR} PYCHASTE_WRAPPER_INCLUDE_DIRS)

##################################################
# Copy Python package structure to build directory
##################################################
set(PYCHASTE_PACKAGE_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/package)
file(
    COPY ${CMAKE_CURRENT_SOURCE_DIR}/py-src/
    DESTINATION ${PYCHASTE_PACKAGE_BINARY_DIR}/
    PATTERN "*${CMAKE_SHARED_LIBRARY_SUFFIX}" EXCLUDE # Skip shared libraries
)

##################################################
# Build Python modules from wrappers
##################################################

# Create a shared library for each pychaste module. The library name is 
# the same as the module name in the main wrapper (e.g. in `ode.main.cpp`).

foreach(module ${PYCHASTE_MODULES})
    file(
        GLOB_RECURSE
        module_SOURCES
        RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
        ${PYCHASTE_WRAPPER_DIR}/${module}/*.cpp
        ${PYCHASTE_WRAPPER_DIR}/${module}/*.hpp
    )

    set(lib_name _pychaste_${module})

    add_library(
        ${lib_name}
        SHARED
        ${module_SOURCES}
        ${PYCHASTE_WRAPPER_DIR}/${module}/${module}.main.cpp
    )

    set_target_properties(
        ${lib_name}
        PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${PYCHASTE_PACKAGE_BINARY_DIR}/chaste/${module}
            PREFIX "${PYTHON_MODULE_PREFIX}"
            SUFFIX "${CMAKE_SHARED_LIBRARY_SUFFIX}"
    )

    target_link_libraries(
        ${lib_name}
        PRIVATE
        ${Python3_LIBRARIES}  # Python and pybind11 come first
        pybind11::module
        Chaste_COMMON_DEPS
        chaste_pychaste  # Library built from non-wrapper code in src/
    )

    target_include_directories(
        ${lib_name}
        PRIVATE
        ${PYTHON3_INCLUDE_DIRS}
        ${Chaste_pychaste_INCLUDE_DIRS}  # Non-wrapper code in src/
        ${PYCHASTE_WRAPPER_INCLUDE_DIRS}
    )

    target_compile_options(${lib_name} PRIVATE -Wno-unused-local-typedefs)

    add_dependencies(pychaste ${lib_name})
endforeach()
