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

message(STATUS "Adding build types...")


##################################
#  COVERAGE for coverage testing #
##################################


SET(CMAKE_CXX_FLAGS_COVERAGE
        "-g -fprofile-arcs -ftest-coverage"
        CACHE STRING "Flags used by the compiler during coverage builds."
        FORCE)
SET(CMAKE_C_FLAGS_COVERAGE
        "-g -fprofile-arcs -ftest-coverage"
        CACHE STRING "Flags used by the compiler during coverage builds."
        FORCE)
SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE
        ""
        CACHE STRING "Flags used for linking binaries during coverage builds."
        FORCE)
SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE
        ""
        CACHE STRING "Flags used by the shared libraries linker during coverage builds."
        FORCE)
MARK_AS_ADVANCED(
        CMAKE_CXX_FLAGS_COVERAGE
        CMAKE_C_FLAGS_COVERAGE
        CMAKE_EXE_LINKER_FLAGS_COVERAGE
        CMAKE_SHARED_LINKER_FLAGS_COVERAGE)


########################################################################
#  AGGRESSIVEOPT for release builds with more aggressive optimisation  #
########################################################################

if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
    set(AGGRESSIVE_OPT_FLAGS "-Ofast -march=native -DNDEBUG")
else ()
    set(AGGRESSIVE_OPT_FLAGS "${CMAKE_CXX_FLAGS_RELEASE}")
endif ()

SET(CMAKE_CXX_FLAGS_AGGRESSIVEOPT
        "${AGGRESSIVE_OPT_FLAGS}" CACHE STRING
        "Flags used by the compiler during release builds with more aggressive optimisation."
        FORCE)
SET(CMAKE_C_FLAGS_AGGRESSIVEOPT
        "${AGGRESSIVE_OPT_FLAGS}" CACHE STRING
        "Flags used by the compiler during release builds with more aggressive optimisation."
        FORCE)
SET(CMAKE_EXE_LINKER_FLAGS_AGGRESSIVEOPT
        "${CMAKE_EXE_LINKER_FLAGS_RELEASE}" CACHE STRING
        "Flags used for linking binaries during release builds with more aggressive optimisation."
        FORCE)
SET(CMAKE_SHARED_LINKER_FLAGS_AGGRESSIVEOPT
        "${CMAKE_EXE_LINKER_FLAGS_RELEASE}" CACHE STRING
        "Flags used by the shared libraries linker during release builds with more aggressive optimisation."
        FORCE)
MARK_AS_ADVANCED(
        CMAKE_CXX_FLAGS_AGGRESSIVEOPT
        CMAKE_C_FLAGS_AGGRESSIVEOPT
        CMAKE_EXE_LINKER_FLAGS_AGGRESSIVEOPT
        CMAKE_SHARED_LINKER_FLAGS_AGGRESSIVEOPT)


#############################################
#  DEBUGOPT for edit-compile-debug workflow #
#############################################

if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") OR (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang"))
    set(DEBUGOPT_OPT_FLAGS "-Og -g")
else ()
    set(DEBUGOPT_OPT_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
endif ()

SET(CMAKE_CXX_FLAGS_DEBUGOPT
        "${DEBUGOPT_OPT_FLAGS}" CACHE STRING
        "Flags used by the compiler during debug builds in an edit-compile-debug workflow."
        FORCE)
SET(CMAKE_C_FLAGS_DEBUGOPT
        "${DEBUGOPT_OPT_FLAGS}" CACHE STRING
        "Flags used by the compiler during debug builds in an edit-compile-debug workflow."
        FORCE)
SET(CMAKE_EXE_LINKER_FLAGS_DEBUGOPT
        "${CMAKE_EXE_LINKER_FLAGS_RELEASE}" CACHE STRING
        "Flags used for linking binaries during debug builds in an edit-compile-debug workflow."
        FORCE)
SET(CMAKE_SHARED_LINKER_FLAGS_DEBUGOPT
        "${CMAKE_EXE_LINKER_FLAGS_RELEASE}" CACHE STRING
        "Flags used by the shared libraries linker debug builds in an edit-compile-debug workflow."
        FORCE)
MARK_AS_ADVANCED(
        CMAKE_CXX_FLAGS_DEBUGOPT
        CMAKE_C_FLAGS_DEBUGOPT
        CMAKE_EXE_LINKER_FLAGS_DEBUGOPT
        CMAKE_SHARED_LINKER_FLAGS_DEBUGOPT)


################################################################
# Update the documentation string of CMAKE_BUILD_TYPE for GUIs #
################################################################

SET(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
        "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel AggressiveOpt Coverage DebugOpt."
        FORCE)


############################################
# Select default if no build type selected #
############################################

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING
            "Choose the type of build: None Debug Release RelWithDebInfo MinSizeRel AggressiveOpt Coverage DebugOpt."
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

message(STATUS "Current build type is : ${CMAKE_BUILD_TYPE}")

