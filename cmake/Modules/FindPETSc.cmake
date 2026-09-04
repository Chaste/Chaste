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
# * Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# * Neither the name of the University of Oxford nor the names of its
# contributors may be used to endorse or promote products derived from this
# software without specific prior written permission.
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

# This file is based on work by Frédéric Simonis @fsimonis, and is used and
# redistributed with permission.

# FindPETSc
# ---------
#
# Locates the PETSc library using pkg-config module PETSc
#
# Imported Targets
# ^^^^^^^^^^^^^^^^
#
# This module defines the following IMPORTED target:
#
# PETSc::PETSc        - the PETSc library
#
# Result Variables
# ^^^^^^^^^^^^^^^^
#
# This module will set the following variables in your project:
#
# PETSc_FOUND          - if false, do not try to link to PETSc
# PETSc_LIBRARIES      - a list of the full paths to all libraries
# PETSc_INCLUDE_DIRS   - a list of all include directories
# PETSc_VERSION        - the full version of PETSc MAJOR.MINOR.PATCH
# PETSc_VERSION_MAJOR  - the MAJOR part of PETSc_VERSION
# PETSc_VERSION_MINOR  - the MINOR part of PETSc_VERSION
# PETSc_VERSION_PATCH  - the PATCH part of PETSc_VERSION

cmake_policy(VERSION 3.10)

# Generate a argument for cmake pkg-config call
if(PETSc_FIND_QUIETLY)
    find_package(PkgConfig QUIET REQUIRED)
else()
    find_package(PkgConfig REQUIRED)
endif()

# Collect candidate PKG_CONFIG_PATH entries in preference order
set(_custom_pkg_paths "")

# ... from PETSC_DIR and PETSC_ARCH (user-defined)
if(DEFINED ENV{PETSC_DIR} AND DEFINED ENV{PETSC_ARCH})
    set(_petsc_env_path "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig")

    if(EXISTS "${_petsc_env_path}")
        list(APPEND _custom_pkg_paths "${_petsc_env_path}")
    endif()

    unset(_petsc_env_path)
endif()

# Prepend to the current PKG_CONFIG_PATH
if(_custom_pkg_paths)
    list(JOIN _custom_pkg_paths ":" _joined_paths)
    set(ENV{PKG_CONFIG_PATH} "${_joined_paths}:$ENV{PKG_CONFIG_PATH}")
endif()

# Build the pkg-config version spec
set(_pkg_version_spec "")

if(DEFINED PETSc_FIND_VERSION)
    if(PETSc_FIND_VERSION_EXACT)
        set(_pkg_version_spec "=${PETSc_FIND_VERSION}")
    else()
        set(_pkg_version_spec ">=${PETSc_FIND_VERSION}")
    endif()
endif()

# Allow system flags
set(ENV{PKG_CONFIG_ALLOW_SYSTEM_CFLAGS} 1)
set(ENV{PKG_CONFIG_ALLOW_SYSTEM_LIBS} 1)

# Use pkg-config to find PETSc
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH "YES")

if(PETSc_FIND_QUIETLY)
    pkg_check_modules(PETSc QUIET IMPORTED_TARGET GLOBAL "PETSc${_pkg_version_spec}")
else()
    pkg_check_modules(PETSc IMPORTED_TARGET GLOBAL "PETSc${_pkg_version_spec}")
endif()

unset(_pkg_version_spec)

# Extract version parts from the version information
if(PC_PETSc_VERSION)
    set(_petsc_versions "")
    string(REGEX MATCHALL "[0-9]+" _petsc_versions ${PETSc_VERSION})
    list(GET _petsc_versions 0 _petsc_version_major)
    list(GET _petsc_versions 1 _petsc_version_minor)
    list(GET _petsc_versions 2 _petsc_version_patch)

    set(PETSc_VERSION ${PC_PETSc_VERSION} CACHE STRING "Full version of PETSc")
    set(PETSc_VERSION_MAJOR ${_petsc_version_major} CACHE INTERNAL "Major version of PETSc")
    set(PETSc_VERSION_MINOR ${_petsc_version_minor} CACHE INTERNAL "Minor version of PETSc")
    set(PETSc_VERSION_PATCH ${_petsc_version_patch} CACHE INTERNAL "Patch version of PETSc")

    unset(_petsc_versions)
    unset(_petsc_version_major)
    unset(_petsc_version_minor)
    unset(_petsc_version_patch)
endif()

# Detect whether the found PETSc was built with CUDA support, by inspecting petscconf.h
set(PETSc_HAVE_CUDA FALSE)
foreach(_petsc_include_dir ${PETSc_INCLUDE_DIRS})
    if(EXISTS "${_petsc_include_dir}/petscconf.h")
        file(STRINGS "${_petsc_include_dir}/petscconf.h" _petsc_have_cuda_line REGEX "^#define PETSC_HAVE_CUDA 1")
        if(_petsc_have_cuda_line)
            set(PETSc_HAVE_CUDA TRUE)
        endif()
        unset(_petsc_have_cuda_line)
    endif()
endforeach()
unset(_petsc_include_dir)
set(PETSc_HAVE_CUDA ${PETSc_HAVE_CUDA} CACHE BOOL "Whether the found PETSc was built with CUDA support")
mark_as_advanced(PETSc_HAVE_CUDA)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
    REQUIRED_VARS PETSc_FOUND PETSc_INCLUDE_DIRS PETSc_LIBRARIES
    VERSION_VAR PETSc_VERSION
)

if(NOT TARGET PETSc::PETSc)
    add_library(PETSc::PETSc ALIAS PkgConfig::PETSc)
endif()

mark_as_advanced(PETSc_INCLUDE_DIRS PETSc_LIBRARIES PETSc_VERSION_MAJOR PETSc_VERSION_MINOR PETSc_VERSION_PATCH VERSION_VAR PETSc_VERSION)
