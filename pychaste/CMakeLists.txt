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

# Create the command for generating wrappers.
# The resulting command will generate wrappers to OUTPUT_DIR.
# If QUIET is set, wrapper generation will run in quiet mode.
# The command will be stored as a list in CPPWG_COMMAND_LST.
# It will also be stored as a string in CPPWG_COMMAND_STR.
function(cppwg_command OUTPUT_DIR)
    cmake_parse_arguments(ARG "QUIET" "" "" ${ARGN})
    if(ARG_QUIET)
        set(QUIET_FLAG "-q")
    else()
        set(QUIET_FLAG "")
    endif()
    set(COMMAND_LST
        . ${chaste_python3_venv}/activate &&
        ${chaste_python3_venv}/cppwg ${CMAKE_SOURCE_DIR}
            -w ${OUTPUT_DIR}
            -p ${CMAKE_SOURCE_DIR}/pychaste/dynamic/config.yaml
            -i ${Chaste_INCLUDE_DIRS}
               ${Chaste_pychaste_INCLUDE_DIRS}
               ${Chaste_THIRD_PARTY_INCLUDE_DIRS}
            -l ${CMAKE_CURRENT_BINARY_DIR}/cppwg.log
            --std c++17 ${QUIET_FLAG} &&
        deactivate
    )
    string(JOIN " " COMMAND_STR ${COMMAND_LST})
    set(CPPWG_COMMAND_LST "${COMMAND_LST}" PARENT_SCOPE)
    set(CPPWG_COMMAND_STR "${COMMAND_STR}" PARENT_SCOPE)
endfunction()

chaste_do_component(pychaste)

set(PYCHASTE_WRAPPER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dynamic/wrappers")
set(PYCHASTE_MODULES cell_based core mesh ode pde visualization)

##################################################
# Python wrapper generation
##################################################
# Add `pychaste_wrappers` target for re-generating wrappers manually
add_custom_target(pychaste_wrappers)
cppwg_command("${PYCHASTE_WRAPPER_DIR}")
add_custom_command(TARGET pychaste_wrappers COMMAND ${CPPWG_COMMAND_LST})

# Generate wrappers if they are missing
foreach(module ${PYCHASTE_MODULES})
    if(NOT EXISTS ${PYCHASTE_WRAPPER_DIR}/${module}/${module}.main.cpp)
        message(STATUS "Generating PyChaste wrappers")
        cppwg_command("${PYCHASTE_WRAPPER_DIR}")
        execute_process(COMMAND sh -c "${CPPWG_COMMAND_STR}")
        break()
    endif()
endforeach()

# Find wrapper header directories
header_dirs(${PYCHASTE_WRAPPER_DIR} PYCHASTE_WRAPPER_INCLUDE_DIRS)

##################################################
# Copy Python package to build directory
##################################################
set(PYCHASTE_PACKAGE_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/package)

set(file_paths "")
file(
    GLOB_RECURSE
    file_paths
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/src/py
    ${CMAKE_CURRENT_SOURCE_DIR}/src/py/*
)
list(FILTER file_paths EXCLUDE REGEX ".*\\.(dll|dylib|pyc|so)$") # Remove compiled files

foreach(file_path ${file_paths})
    if(file_path MATCHES ".*\\.(cfg|in|js|py|toml)$") # Files to track changes on
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/src/py/${file_path}
            ${PYCHASTE_PACKAGE_BINARY_DIR}/${file_path}
            COPYONLY
        )
    else()
        file(
            COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/py/${file_path}
            DESTINATION ${PYCHASTE_PACKAGE_BINARY_DIR}/
        )
    endif()
endforeach()

##################################################
# Build Python modules from wrappers
##################################################
# Create a shared library for each pychaste module. The library name is
# the same as the module name in the main wrapper e.g. in `ode.main.cpp`,
# the module name is `_pychaste_ode`, which is imported in Python code as
# `from chaste.ode._pychaste_ode import *`.

foreach(module ${PYCHASTE_MODULES})
    file(
        GLOB_RECURSE
        module_sources
        RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
        ${PYCHASTE_WRAPPER_DIR}/${module}/*.cpp
        ${PYCHASTE_WRAPPER_DIR}/${module}/*.hpp
    )

    set(module_lib _pychaste_${module})

    add_library(${module_lib} SHARED ${module_sources})

    set_target_properties(
        ${module_lib}
        PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${PYCHASTE_PACKAGE_BINARY_DIR}/chaste/${module}
            PREFIX "${PYTHON_MODULE_PREFIX}"
            SUFFIX "${CMAKE_SHARED_LIBRARY_SUFFIX}"
    )

    target_link_libraries(
        ${module_lib}
        PRIVATE
        ${Python3_LIBRARIES}  # Python and pybind11 come first
        pybind11::module
        Chaste_COMMON_DEPS
        chaste_pychaste  # Library built from non-wrapper cpp code in src/
    )

    target_include_directories(
        ${module_lib}
        PRIVATE
        ${PYTHON3_INCLUDE_DIRS}
        ${Chaste_pychaste_INCLUDE_DIRS}  # Non-wrapper cpp code in src/
        ${PYCHASTE_WRAPPER_INCLUDE_DIRS}
    )

    target_compile_options(${module_lib} PRIVATE -Wno-unused-local-typedefs)

    add_dependencies(pychaste ${module_lib})
endforeach()
