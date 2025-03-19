# Copyright (c) 2005-2025, University of Oxford.
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

################################
####  Generate Python wrappers
################################

# Get includes from `chaste_pychaste` (target for non-wrapper code in pychaste/src/cpp)
get_target_property(PYCHASTE_INCLUDE_DIRS chaste_pychaste INCLUDE_DIRECTORIES)

# Add includes for VTK 9+ (includes for VTK < 9 will already be in the list)
if(TARGET VTK::CommonCore)
    get_target_property(_include_dirs VTK::CommonCore INTERFACE_INCLUDE_DIRECTORIES)
    list(APPEND PYCHASTE_INCLUDE_DIRS ${_include_dirs})
endif()

# Add PETSC4Py includes
list(APPEND PYCHASTE_INCLUDE_DIRS ${PETSC4PY_INCLUDES})

# Add includes for pybind11 typecasters
header_dirs("${CMAKE_CURRENT_SOURCE_DIR}/dynamic" _include_dirs)
list(APPEND PYCHASTE_INCLUDE_DIRS ${_include_dirs})

# Command for generating wrappers
set(PYCHASTE_WRAPPER_GENERATION_COMMAND
    ${chaste_python3_venv}/cppwg ${CMAKE_SOURCE_DIR}
        -w ${CMAKE_CURRENT_BINARY_DIR}/wrappers
        -p ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/config.yaml
        -i ${PYCHASTE_INCLUDE_DIRS}
        -l ${CMAKE_CURRENT_BINARY_DIR}/cppwg.log
        --std c++17
)

# Generate wrappers if not already created
if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/wrappers")
    execute_process(COMMAND ${PYCHASTE_WRAPPER_GENERATION_COMMAND})
endif()

# Target for manually re-generating the wrappers
add_custom_target(pychaste_wrappers
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/wrappers
    COMMAND ${PYCHASTE_WRAPPER_GENERATION_COMMAND}
)

# Add generated wrappers to includes
header_dirs("${CMAKE_CURRENT_BINARY_DIR}/wrappers" _include_dirs)
list(APPEND PYCHASTE_INCLUDE_DIRS ${_include_dirs})

################################
####  Build Python module
################################
# Creates a `_pychaste_all` shared library from the Python wrappers.
# The library name `_pychaste_all` is the same as the module name in the main 
# module wrapper source file `wrappers/all/_pychaste_all.main.cppwg.cpp`.
# When the shared library has been built, it can be used in Python code 
# as `from chaste._pychaste_all import *`.

file(GLOB_RECURSE PYCHASTE_WRAPPER_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/wrappers/all/*.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/wrappers/all/*.hpp
)

pybind11_add_module(_pychaste_all SHARED ${PYCHASTE_WRAPPER_SOURCES})

target_link_libraries(_pychaste_all PRIVATE
    ${Python3_LIBRARIES}
    Chaste_COMMON_DEPS
    chaste_pychaste
)

target_include_directories(_pychaste_all PRIVATE ${PYCHASTE_INCLUDE_DIRS})

target_compile_options(_pychaste_all PRIVATE -flto=auto -Wno-unused-local-typedefs)

set_target_properties(_pychaste_all PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/package/chaste
)

add_dependencies(pychaste _pychaste_all)

################################
####  Copy Python package
################################

# Get list of files in the Python package directory: pychaste/src/py/
set(package_file_paths "")
file(
    GLOB_RECURSE
    package_file_paths
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/src/py
    ${CMAKE_CURRENT_SOURCE_DIR}/src/py/*
)

# Exclude compiled files
list(FILTER package_file_paths EXCLUDE REGEX ".*\\.(dll|dylib|pyc|so)$")

# Copy files to the build directory
foreach(file_path ${package_file_paths})
    if(file_path MATCHES ".*\\.(cfg|in|js|py|toml)$")
        # Copy file and track changes
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/src/py/${file_path}
            ${CMAKE_CURRENT_BINARY_DIR}/package/${file_path}
            COPYONLY
        )
    else()
        # Copy file, but don't track changes
        file(
            COPY ${CMAKE_CURRENT_SOURCE_DIR}/src/py/${file_path}
            DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/package/
        )
    endif()
endforeach()
