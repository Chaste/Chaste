# Copyright (c) 2005-2026, University of Oxford.
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
    foreach(_dir IN LISTS _include_dirs)
        if(_dir MATCHES "^\\$<BUILD_INTERFACE:(.+)>$")
            # VTK 9.5+ have $<BUILD_INTERFACE:path> generator expressions in their includes
            list(APPEND PYCHASTE_INCLUDE_DIRS "${CMAKE_MATCH_1}")
        elseif(NOT _dir MATCHES "^\\$<")
            # VTK 9.0 to 9.4 have plain include paths, no generator expressions
            list(APPEND PYCHASTE_INCLUDE_DIRS "${_dir}")
        endif()
    endforeach()
endif()

# Add PETSC4Py includes
list(APPEND PYCHASTE_INCLUDE_DIRS ${PETSC4PY_INCLUDES})

# Add includes for pybind11 typecasters
header_dirs("${CMAKE_CURRENT_SOURCE_DIR}/dynamic" _include_dirs)
list(APPEND PYCHASTE_INCLUDE_DIRS ${_include_dirs})

# Expose the includes for user projects
foreach(_dir IN LISTS PYCHASTE_INCLUDE_DIRS)
    target_include_directories(chaste_pychaste INTERFACE "$<BUILD_INTERFACE:${_dir}>")
endforeach()

# Command for generating wrappers.
#  - gnu++17 (rather than c++17) and -Wno-missing-template-arg-list-after-template-kw
#    silence a castxml/LLVM error on use of the template keyword without a
#    following argument list, which appears in some boost headers.
#  - -U__SIZEOF_FLOAT128__ stops boost.multiprecision from including quadmath.h,
#    a GCC-only header that castxml's bundled Clang cannot find.
set(PYCHASTE_WRAPPER_GENERATION_COMMAND
    ${chaste_python3_venv}/cppwg ${CMAKE_SOURCE_DIR}
        -w ${CMAKE_CURRENT_BINARY_DIR}/wrappers
        -p ${CMAKE_CURRENT_SOURCE_DIR}/dynamic/config.yaml
        -i ${PYCHASTE_INCLUDE_DIRS}
        -l ${CMAKE_CURRENT_BINARY_DIR}/cppwg.log
        --std "gnu++17"
        --castxml_cflags="-Wno-missing-template-arg-list-after-template-kw -U__SIZEOF_FLOAT128__"
        --castxml_compiler ${CMAKE_CXX_COMPILER}
)

# Generate wrappers if not already created
if(NOT EXISTS "${CMAKE_CURRENT_BINARY_DIR}/wrappers")
    execute_process(COMMAND ${PYCHASTE_WRAPPER_GENERATION_COMMAND} COMMAND_ERROR_IS_FATAL ANY)
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
# Creates a `_pychaste_all` library from the Python wrappers.
# The library name `_pychaste_all` is the same as the module name in the main
# module wrapper source file `wrappers/all/_pychaste_all.main.cppwg.cpp`.
# When the library has been built, it can be used in Python code
# as `from chaste._pychaste_all import *`.

file(GLOB_RECURSE PYCHASTE_WRAPPER_SOURCES
    ${CMAKE_CURRENT_BINARY_DIR}/wrappers/all/*.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/wrappers/all/*.hpp
)

pybind11_add_module(_pychaste_all OPT_SIZE ${PYCHASTE_WRAPPER_SOURCES})

target_link_libraries(_pychaste_all PRIVATE
    Chaste_COMMON_DEPS
    chaste_pychaste
)

target_include_directories(_pychaste_all PRIVATE ${PYCHASTE_INCLUDE_DIRS})

# -g0: Disable debug symbols in the Python modules as they are not needed.
target_compile_options(_pychaste_all PRIVATE -flto=auto -Wno-unused-local-typedefs -g0)

set_target_properties(_pychaste_all PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/package/chaste
)

add_dependencies(pychaste _pychaste_all)

################################
####  Copy Python package
################################

file(GLOB_RECURSE package_file_paths
    CONFIGURE_DEPENDS
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/src/py
    ${CMAKE_CURRENT_SOURCE_DIR}/src/py/*
)
list(FILTER package_file_paths EXCLUDE REGEX ".*\\.(dll|dylib|pyc|so)$")
foreach(file_path ${package_file_paths})
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/src/py/${file_path}
        ${CMAKE_CURRENT_BINARY_DIR}/package/${file_path}
        COPYONLY
    )
endforeach()
