# Acknowledgement: This file is based on the FindGperftools.cmake file from the GitHub repo https://github.com/mavam/vast (License below)
# 
# Copyright (c) 2014, Matthias Vallentin
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
# 
# Tries to find Gperftools.
#
# Usage of this module as follows:
#
#     find_package(Gperftools)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  Gperftools_ROOT_DIR  Set this variable to the root installation of
#                       Gperftools if the module has problems finding
#                       the proper installation path.
#
# Variables defined by this module:
#
#  GPERFTOOLS_FOUND              System has Gperftools libs/headers
#  GPERFTOOLS_LIBRARIES          The Gperftools libraries (tcmalloc & profiler)
#  GPERFTOOLS_INCLUDE_DIR        The location of Gperftools headers

find_library(GPERFTOOLS_TCMALLOC
  NAMES tcmalloc
  HINTS ${Gperftools_ROOT_DIR}/lib)

find_library(GPERFTOOLS_PROFILER
  NAMES profiler
  HINTS ${Gperftools_ROOT_DIR}/lib)

find_library(GPERFTOOLS_TCMALLOC_AND_PROFILER
  NAMES tcmalloc_and_profiler
  HINTS ${Gperftools_ROOT_DIR}/lib)

find_path(GPERFTOOLS_INCLUDE_DIR
  NAMES gperftools/heap-profiler.h google/heap-profiler.h
  HINTS ${Gperftools_ROOT_DIR}/include)

find_file(GPERFTOOLS_PPROF_EXE
  NAMES pprof google-pprof
  HINTS ${Gperftools_ROOT_DIR}/bin)

if (GPERFTOOLS_TCMALLOC_AND_PROFILER)
    message("-- Gperftools: Found combined tcmalloc and profiler")
    message("-- \t ${GPERFTOOLS_TCMALLOC_AND_PROFILER}")
    set(GPERFTOOLS_LIBRARIES ${GPERFTOOLS_TCMALLOC_AND_PROFILER})
elseif (GPERFTOOLS_PROFILER AND GPERFTOOLS_TCMALLOC)
    message("-- Gperftools: Found separate tcmalloc and profiler")
    message("-- \t ${GPERFTOOLS_TCMALLOC} ${GPERFTOOLS_PROFILER}")
    set(GPERFTOOLS_LIBRARIES ${GPERFTOOLS_TCMALLOC} ${GPERFTOOLS_PROFILER})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Gperftools
  DEFAULT_MSG
  GPERFTOOLS_LIBRARIES
  GPERFTOOLS_PPROF_EXE
  GPERFTOOLS_INCLUDE_DIR)

mark_as_advanced(
  Gperftools_ROOT_DIR
  GPERFTOOLS_TCMALLOC
  GPERFTOOLS_PROFILER
  GPERFTOOLS_PPROF_EXE
  GPERFTOOLS_TCMALLOC_AND_PROFILER
  GPERFTOOLS_LIBRARIES
  GPERFTOOLS_INCLUDE_DIR)
