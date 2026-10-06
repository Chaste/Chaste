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
#

# Preprocessor definitions that Chaste itself no longer uses, but which are deliberately still defined.
#
# VTK, CVODE and Xerces are required dependencies of Chaste (see issues #578 and #580), so the code
# that used to be guarded by `#ifdef CHASTE_VTK`, `#ifdef CHASTE_CVODE` and `#ifdef CHASTE_XERCES` has
# been made unconditional and nothing inside Chaste tests these macros any more.
#
# Code outside this repository (user projects, forks, downstream applications) may still contain
# `#ifdef` guards on them. If the definitions were removed, that code would silently compile out.
# They are therefore kept defined, and must remain so until a deliberate deprecation cycle says otherwise.
# global/test/TestChasteBuildInfo.hpp checks that they are still defined.
#
# This module is included both from the top-level CMakeLists.txt and from ChasteConfig.cmake.in, so that
# projects using find_package(Chaste) receive the same definitions.

add_definitions (-DCHASTE_VTK)
add_definitions (-DCHASTE_CVODE)
add_definitions (-DCHASTE_XERCES)
