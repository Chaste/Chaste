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

# Options that used to control optional Chaste dependencies but no longer have any effect.
# If a user or a cached build directory still has one of these set, warn and unset it rather
# than silently ignoring it.

# A plain WARNING (rather than DEPRECATION) is used below so the message is always shown,
# regardless of CMAKE_WARN_DEPRECATED / CMAKE_ERROR_DEPRECATED, which are meant to govern
# CMake's own API deprecations rather than project-specific ones like this.
#
# ${Chaste_USE_VTK} etc. are referenced in the messages below (not just DEFINED) so that CMake
# does not also report them as unused manually-specified variables.

# VTK is now a required dependency of Chaste (see https://github.com/Chaste/Chaste/issues/578).
if (DEFINED Chaste_USE_VTK)
    message (WARNING "Chaste_USE_VTK (set to '${Chaste_USE_VTK}') is deprecated and is now ignored: VTK is a required dependency of Chaste.")
    unset (Chaste_USE_VTK CACHE)
endif ()

# CVODE is now a required dependency of Chaste (see https://github.com/Chaste/Chaste/issues/578).
if (DEFINED Chaste_USE_CVODE)
    message (WARNING "Chaste_USE_CVODE (set to '${Chaste_USE_CVODE}') is deprecated and is now ignored: CVODE is a required dependency of Chaste.")
    unset (Chaste_USE_CVODE CACHE)
endif ()

# Xerces is now a required dependency of Chaste (see https://github.com/Chaste/Chaste/issues/580).
if (DEFINED Chaste_USE_XERCES)
    message (WARNING "Chaste_USE_XERCES (set to '${Chaste_USE_XERCES}') is deprecated and is now ignored: Xerces is a required dependency of Chaste.")
    unset (Chaste_USE_XERCES CACHE)
endif ()
