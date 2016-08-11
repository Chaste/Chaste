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
    message("Chaste Release Full Version = ${full_version}, Revision = ${Chaste_revision}")
else()
    # ReleaseVersion file not found, obtain revision information from SVN
	# The following requires a proper command-line svn client to be installed, not
	# just an ordinary shell extension like TortoiseSVN'
	# Install SlikSVN or the distribution from Collabnet (if you don't mind registering)
	find_package(Subversion QUIET)
	if(SUBVERSION_FOUND)
       execute_process(COMMAND ${Subversion_SVN_EXECUTABLE} info "${Chaste_SOURCE_DIR}"
            OUTPUT_VARIABLE dummy
            ERROR_VARIABLE Subversion_svn_info_error
            RESULT_VARIABLE Subversion_svn_info_result
            OUTPUT_STRIP_TRAILING_WHITESPACE)

       if(${Subversion_svn_info_result} EQUAL 0)
          Subversion_WC_INFO("${Chaste_SOURCE_DIR}" Chaste)
	      set(Chaste_revision "${Chaste_WC_REVISION}")
	      message("Current Chaste SVN Revision = ${Chaste_WC_REVISION}. Chaste Last Changed Revision = ${Chaste_WC_LAST_CHANGED_REV}")
	      if(${Chaste_WC_REVISION} EQUAL ${Chaste_WC_LAST_CHANGED_REV})
	         set(Chaste_WC_MODIFIED "false")
	      else()
	         set(Chaste_WC_MODIFIED "true")
	      endif()
       endif()
    endif(SUBVERSION_FOUND)
    if (NOT (SUBVERSION_FOUND AND (${Subversion_svn_info_result} EQUAL 0)))
       # assume its a git repo 
       find_package(Git REQUIRED)
       Git_WC_INFO("${Chaste_SOURCE_DIR}" Chaste)
       set(Chaste_revision "${Chaste_WC_REVISION}")
       message("Current Chaste Git Revision = ${Chaste_WC_REVISION}. Chaste Modified = ${Chaste_WC_MODIFIED}")
    endif()
endif()

if (Chaste_UPDATE_PROVENANCE)
    set(Chaste_REVISION ${Chaste_revision} CACHE STRING "Current Chaste Git Revision" FORCE)
else()
    set(Chaste_REVISION ${Chaste_revision} CACHE STRING "Current Chaste Git Revision")
endif()


#string(TIMESTAMP build_time)
if (NOT EXISTS build_timestamp OR Chaste_UPDATE_PROVENANCE)
    message("updating buildtime...")
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
set(project_versions ";")
#TODO: update project versions
foreach(project ${Chaste_PROJECTS})
    set(project_versions "${project_versions} versions[\"${project}\"] = 1.0;\n")
endforeach()

list(APPEND Chaste_PROJECTS "${potential_dir}")


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

