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

# chaste_get_source_revision(<source tree> <out commit var> <out modified var>)
#
# Determines which commit a source tree (Chaste itself, or a checked out project) is at, and
# whether it has uncommitted changes ("true"/"false"):
#
#  * A git working copy (<source tree>/.git exists, as a directory, or as a file for a worktree):
#    the full commit hash, and modified if git reports any change to tracked files.
#  * An archive created by git, such as the source tarball GitHub attaches to each release: the
#    full commit hash git substitutes into <source tree>/.git_archival.txt (see .gitattributes).
#    Whether files were edited after extraction cannot be known, so this reports unmodified.
#  * Otherwise, "unknown" and unmodified.
#
# If git itself fails (for example refusing a repository owned by another user), the commit is
# "unknown" rather than stopping the build. That is reported by a warning at configure time; the
# refresh on every build (cmake -P script mode) stays silent.
function(chaste_get_source_revision source_tree out_commit out_modified)
    set(commit "unknown")
    set(modified "false")
    if (EXISTS "${source_tree}/.git")
        find_package(Git QUIET)
        if (Git_FOUND)
            execute_process(
                COMMAND "${GIT_EXECUTABLE}" -C "${source_tree}" rev-parse --verify HEAD
                RESULT_VARIABLE head_result
                OUTPUT_VARIABLE head
                ERROR_VARIABLE git_error
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_STRIP_TRAILING_WHITESPACE)
            # --no-optional-locks stops this taking git's index lock, which runs on every build and
            # would otherwise clash with any git command being run in the repository at the time.
            execute_process(
                COMMAND "${GIT_EXECUTABLE}" --no-optional-locks -C "${source_tree}" status --porcelain --untracked-files=no
                RESULT_VARIABLE status_result
                OUTPUT_VARIABLE status
                ERROR_VARIABLE status_error
                ERROR_STRIP_TRAILING_WHITESPACE)
            if (head_result EQUAL 0 AND status_result EQUAL 0)
                set(commit "${head}")
                if (NOT status STREQUAL "")
                    set(modified "true")
                endif ()
            elseif (NOT CMAKE_SCRIPT_MODE_FILE)
                if (git_error STREQUAL "")
                    set(git_error "${status_error}")
                endif ()
                message(WARNING "Could not determine the git revision of ${source_tree}, so it will be reported as unknown: ${git_error}")
            endif ()
        elseif (NOT CMAKE_SCRIPT_MODE_FILE)
            message(WARNING "git was not found, so the revision of ${source_tree} will be reported as unknown")
        endif ()
    elseif (EXISTS "${source_tree}/.git_archival.txt")
        # Outside an archive the file holds an unsubstituted placeholder, which this doesn't match
        file(STRINGS "${source_tree}/.git_archival.txt" node REGEX "^node: [0-9a-f]+$")
        if (node MATCHES "^node: ([0-9a-f]+)$")
            set(commit "${CMAKE_MATCH_1}")
        endif ()
    endif ()
    set(${out_commit} "${commit}" PARENT_SCOPE)
    set(${out_modified} "${modified}" PARENT_SCOPE)
endfunction()
