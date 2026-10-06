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

# Writes the plain files under provenance_dir that must reflect the current build, not just the
# last cmake reconfigure: the build timestamp, Chaste's own git revision and modified flag, and
# the same for each checked out project under projects/. A project's own commit state, like
# Chaste's, can change between two builds with no reconfigure in between.
#
# Invoked by the generate_provenance custom target in global/CMakeLists.txt, which has no
# tracked output and so runs every time chaste_global is built. Nothing here is compiled, so
# writing a new provenance file has no effect on the build graph: it never marks chaste_global,
# or anything that depends on it, as needing to recompile or relink.
#
# Everything else that used to be baked into Version.cpp (version numbers, compiler info, xsd
# and chaste_codegen versions, and so on) is fixed for the life of a cmake configuration, and is
# generated separately, directly in global/CMakeLists.txt, by a plain configure_file() call.

# Run as a cmake -P script, which starts with no policies set, so match the main project's.
cmake_minimum_required(VERSION 3.22.1)

include(ChasteDetermineGitRevision)

# Writes each file to a temporary name and then renames it into place, which replaces the old file
# atomically: a test that is still running from this build tree, and that reads its provenance at
# just that moment, sees either the old content or the new, never a partly written file.
function(write_provenance_file path content)
    file(WRITE "${path}.tmp" "${content}")
    file(RENAME "${path}.tmp" "${path}")
endfunction()

chaste_determine_git_revision("${Chaste_SOURCE_DIR}" Chaste_revision Chaste_modified)
write_provenance_file("${provenance_dir}/git_revision" "${Chaste_revision}")
write_provenance_file("${provenance_dir}/git_modified" "${Chaste_modified}")

# Local time with its true offset from UTC. The %z specifier needs CMake 3.26, so older versions
# fall back to UTC, labelled accordingly: either way the offset shown is correct.
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.26)
    string(TIMESTAMP build_time "%a, %d %b %Y %H:%M:%S %z")
else ()
    string(TIMESTAMP build_time "%a, %d %b %Y %H:%M:%S +0000" UTC)
endif ()
write_provenance_file("${provenance_dir}/timestamp" "${build_time}")

# Determine project versions (either the git hash or svn revision number), and whether there
# are uncommitted revisions, for each checked out project under projects/. A subdirectory per
# project, rather than one file naming them all, means the reader (ChasteBuildInfo, via
# FileFinder::FindMatches) never needs a separate manifest of which projects exist. Directories
# for projects that are no longer checked out are removed, so they are not reported.
file(GLOB existing_projects LIST_DIRECTORIES true RELATIVE "${provenance_dir}/projects" "${provenance_dir}/projects/*")
foreach (existing_project ${existing_projects})
    if (NOT existing_project IN_LIST Chaste_PROJECTS)
        file(REMOVE_RECURSE "${provenance_dir}/projects/${existing_project}")
    endif ()
endforeach ()
foreach (project ${Chaste_PROJECTS})
    if (IS_DIRECTORY "${Chaste_SOURCE_DIR}/projects/${project}/.git")
        execute_process(
                COMMAND git -C "${Chaste_SOURCE_DIR}/projects/${project}" rev-parse HEAD
                OUTPUT_VARIABLE this_project_version
        )
        execute_process(
                COMMAND git -C "${Chaste_SOURCE_DIR}/projects/${project}" diff-index HEAD --
                OUTPUT_VARIABLE diff_index_result
        )
        if (diff_index_result STREQUAL "")
            set(this_project_modified "False")
        else ()
            set(this_project_modified "True")
        endif ()
    elseif (IS_DIRECTORY "${Chaste_SOURCE_DIR}/projects/${project}/.svn")
        execute_process(
                COMMAND svnversion "${Chaste_SOURCE_DIR}/projects/${project}"
                OUTPUT_VARIABLE this_project_version
        )
        if (this_project_version MATCHES "M")
            set(this_project_modified "True")
        else ()
            set(this_project_modified "False")
        endif ()
    else ()
        set(this_project_version "Unknown")
        set(this_project_modified "False")
    endif ()

    string(STRIP "${this_project_version}" this_project_version)

    file(MAKE_DIRECTORY "${provenance_dir}/projects/${project}")
    write_provenance_file("${provenance_dir}/projects/${project}/version" "${this_project_version}")
    write_provenance_file("${provenance_dir}/projects/${project}/modified" "${this_project_modified}")
endforeach ()
