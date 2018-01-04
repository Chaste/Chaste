# Copyright (c) 2005-2018, University of Oxford.
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

###########################
#   Generate Version.cpp  #
###########################

set(Chaste_WC_MODIFIED "false")
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/ReleaseVersion.txt")
if(EXISTS "${version_file}")
    file(STRINGS "${version_file}" v_data)
    list(GET v_data 0 full_version)
    string(REPLACE "." ";" full_version_list "${full_version}")
    list(LENGTH full_version_list len)
    math(EXPR length ${len}-1)
    list(GET full_version_list ${length} Chaste_revision)
    message(STATUS "Chaste Release Full Version = ${full_version}, Revision = ${Chaste_revision}")
elseif(EXISTS "${Chaste_SOURCE_DIR}/.git")
    # If ReleaseVersion.txt not found, obtain revision information from Git
    find_package(Git REQUIRED)
    Git_WC_INFO("${Chaste_SOURCE_DIR}" Chaste)
    set(Chaste_revision "${Chaste_WC_REVISION}")
    message(STATUS "Current Chaste Git Revision = ${Chaste_WC_REVISION}. Chaste Modified = ${Chaste_WC_MODIFIED}")
else()
    set(Chaste_revision "0")
    message(STATUS "Cannot find ReleaseVersion.txt or Git revision")
endif()

if (Chaste_UPDATE_PROVENANCE)
    set(Chaste_REVISION ${Chaste_revision} CACHE STRING "Current Chaste Git Revision" FORCE)
else()
    set(Chaste_REVISION ${Chaste_revision} CACHE STRING "Current Chaste Git Revision")
endif()


#string(TIMESTAMP build_time)
if (NOT EXISTS build_timestamp OR Chaste_UPDATE_PROVENANCE)
    message(STATUS "updating buildtime...")
    execute_process(COMMAND ${timekeeper_exe})
endif()
file(READ build_timestamp build_time)
execute_process(COMMAND ${XSD_EXECUTABLE} "--version" ERROR_VARIABLE xsd_version_full)
string(REGEX MATCH "^XML Schema Definition Compiler" xsd_version_2 "${xsd_version_full}")
string(REGEX MATCH "^CodeSynthesis XSD XML Schema to C\\+\\+ compiler" xsd_version_3 "${xsd_version_full}")
if (xsd_version_2) 
    set(xsd_version "2")
elseif (xsd_version_3)
    set(xsd_version "3")
else()
    set(xsd_version "undertermined")
endif()

set(time_size 80)
set(time_format "%a, %d %b %Y %H:%M:%S +0000")

# Determine project versions (either the git hash or svn revision number), and whether there are uncommited revisions
foreach(project ${Chaste_PROJECTS})
    # Project is a git repo
    if (IS_DIRECTORY "${Chaste_SOURCE_DIR}/projects/${project}/.git")
        # Determine the git hash as a string
        execute_process(
                COMMAND git -C ${Chaste_SOURCE_DIR}/projects/${project} rev-parse --short HEAD
                OUTPUT_VARIABLE this_project_version
        )
        # Determine whether there are uncommitted revisions
        execute_process(
                COMMAND git -C ${Chaste_SOURCE_DIR}/projects/${project} diff-index HEAD --
                OUTPUT_VARIABLE diff_index_result
        )
        if (diff_index_result STREQUAL "")
            set(this_project_modified "False")
        else()
            set(this_project_modified "True")
        endif()
    # Project is an svn repo
    elseif(IS_DIRECTORY "${Chaste_SOURCE_DIR}/projects/${project}/.svn")
        # Determine the svn revision number as a string
        execute_process(
                COMMAND svnversion ${Chaste_SOURCE_DIR}/projects/${project}
                OUTPUT_VARIABLE this_project_version
        )
        # Determine whether there are uncommitted revisions
        if (${this_project_version} MATCHES "M")
            set(this_project_modified "True")
        else()
            set(this_project_modified "False")
        endif()
    # If it's not git or svn, we pass out some default values to indicate unknown version
    else()
        set(this_project_version "Unknown")
        set(this_project_modified "False")
    endif()

    # Strip trailing whitespace from project version
    string(STRIP ${this_project_version} this_project_version)

    set(project_versions "${project_versions} versions[\"${project}\"] = \"${this_project_version}\";\n")
    set(projects_modified "${projects_modified} modified[\"${project}\"] = \"${this_project_modified}\";\n")
endforeach()

find_package(PythonInterp QUIET)
execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" "from CheckForCopyrights import current_notice; print current_notice"
 WORKING_DIRECTORY "${Chaste_SOURCE_DIR}/python/infra"
 OUTPUT_VARIABLE licence)
string(REPLACE "\nThis file is part of Chaste.\n" "" licence "${licence}")
string(REPLACE "\n" "\\n" licence "${licence}")
set(quote "\"")
string(REPLACE ${quote} "\\${quote}" licence "${licence}")
string(REGEX REPLACE "\\\\n$" "" licence "${licence}")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${Chaste_SOURCE_DIR}/global/src/Version_cmake.cpp.in"
  ${generate_dir}/Version.cpp
  )

##################################
#  Generate ChasteBuildInfo.cpp  #
##################################

set(additional "")
if(MSVC)
    set(additional
    "
    //This has been put here to satisfy an MSVC linker issue 
    int __cdecl _purecall(void){return 0;}
    ")
endif()

configure_file (
  "${Chaste_SOURCE_DIR}/global/src/ChasteBuildInfo_cmake.cpp.in"
  ${generate_dir}/ChasteBuildInfo.cpp
  )

