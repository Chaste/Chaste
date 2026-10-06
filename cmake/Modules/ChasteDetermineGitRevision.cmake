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

# chaste_determine_git_revision(<source dir> <out revision var> <out modified var>)
#
# Determines a single revision token for a Chaste-style checkout: the last dot-separated
# component of ReleaseVersion.txt's version string for an official release tarball, the
# current git commit hash for an ordinary checked out working copy, or "0" if neither is
# available. Also reports ("true"/"false") whether a git working copy has uncommitted changes;
# this is always "false" for the release-tarball and unknown cases.
#
# Shared between the top-level CMakeLists.txt (which only needs a configure-time snapshot, for
# labelling the doxygen and coverage targets) and cmake/Modules/ChasteGenerateProvenance.cmake
# (which needs a fresh value on every build, since a commit can happen between two builds with
# no cmake reconfigure in between).
#
# Status messages are only printed at configure time: in script mode (cmake -P, as used by the
# generate_provenance target on every build) this runs silently, to avoid noise in build logs.
macro(CHASTE_DETERMINE_GIT_REVISION source_dir out_revision out_modified)
    set(${out_modified} "false")
    set(_cdgr_release_version_file "${CMAKE_CURRENT_BINARY_DIR}/ReleaseVersion.txt")
    if (EXISTS "${_cdgr_release_version_file}")
        file(STRINGS "${_cdgr_release_version_file}" _cdgr_release_data)
        list(GET _cdgr_release_data 0 _cdgr_full_version)
        string(REPLACE "." ";" _cdgr_full_version_list "${_cdgr_full_version}")
        list(LENGTH _cdgr_full_version_list _cdgr_len)
        math(EXPR _cdgr_len "${_cdgr_len}-1")
        list(GET _cdgr_full_version_list ${_cdgr_len} ${out_revision})
        set(_cdgr_message "Chaste Release Full Version = ${_cdgr_full_version}, Revision = ${${out_revision}}")
    elseif (EXISTS "${source_dir}/.git")
        find_package(Git REQUIRED QUIET)
        Git_WC_INFO("${source_dir}" _cdgr)
        set(${out_revision} "${_cdgr_WC_REVISION}")
        set(${out_modified} "${_cdgr_WC_MODIFIED}")
        set(_cdgr_message "Current Chaste Git Revision = ${_cdgr_WC_REVISION}. Chaste Modified = ${_cdgr_WC_MODIFIED}")
    else ()
        set(${out_revision} "0")
        set(_cdgr_message "Cannot find ReleaseVersion.txt or Git revision")
    endif ()
    if (NOT CMAKE_SCRIPT_MODE_FILE)
        message(STATUS "${_cdgr_message}")
    endif ()
endmacro()
