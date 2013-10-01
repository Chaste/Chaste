project(global)

# The auto-generated source files
set(version_cpp "${CHASTE_BUILD_ROOT}/src/Version.cpp")
set(build_root_cpp "${CHASTE_BUILD_ROOT}/src/GetChasteBuildRoot.cpp")

#Generate Version.cpp
set(CHASTE_WC_MODIFIED "false")
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/ReleaseVersion.txt")
if(EXISTS "${version_file}")
    file(STRINGS "${version_file}" v_data)
    list(GET v_data 0 full_version)
    string(REPLACE "." ";" full_version_list "${full_version}")
    list(LENGTH full_version_list len)
    math(EXPR length ${len}-1)
    list(GET full_version_list ${length} chaste_revision)
    message("Chaste Release Full Version = ${full_version}, Revision = ${chaste_revision}")
else()
    # ReleaseVersion file not found, obtain revision information from SVN

    # The following requires a proper command-line svn client to be installed, not
    # just an ordinary shell extension like TortoiseSVN'
    # Install SlikSVN or the distribution from Collabnet (if you don't mind registering)
    find_package(Subversion)
    if(SUBVERSION_FOUND)
        Subversion_WC_INFO("${CHASTE_SOURCE_ROOT}" chaste)
        set(chaste_revision "${chaste_WC_REVISION}")
        message("Current Chaste SVN Revision = ${chaste_WC_REVISION}. Chaste Last Changed Revision = ${chaste_WC_LAST_CHANGED_REV}")
    if(${chaste_WC_REVISION} EQUAL ${chaste_WC_LAST_CHANGED_REV})
        set(CHASTE_WC_MODIFIED "false")
    else()
        set(CHASTE_WC_MODIFIED "true")
    endif()
    endif(SUBVERSION_FOUND)
endif()
file(STRINGS "${CHASTE_BUILD_ROOT}/build_timestamp" time_stamp)
set(compiler_type "${CMAKE_CXX_COMPILER}")
set(compiler_version "${CMAKE_CXX_COMPILER_VERSION}")
set(xsd_version "undertermined")

find_package(PythonInterp)
execute_process(COMMAND "${PYTHON_EXECUTABLE}" "-c" "from CheckForCopyrights import current_notice; print current_notice"
 WORKING_DIRECTORY "${CHASTE_SOURCE_ROOT}/python/infra"
 OUTPUT_VARIABLE licence)
string(REPLACE "\nThis file is part of Chaste.\n" "" licence "${licence}")
string(REPLACE "\n" "\\n" licence "${licence}")
set(quote "\"")
string(REPLACE ${quote} "\\${quote}" licence "${licence}")

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/src/Version.cpp.in" template)
string(REPLACE "%(licence)s" "${licence}" template "${template}")
string(REPLACE "%(chaste_root)s" "${CHASTE_SOURCE_ROOT}/" template "${template}")
string(REPLACE "%(revision)s" "${chaste_revision}" template "${template}")
string(REPLACE "%(wc_modified)s" "${CHASTE_WC_MODIFIED}" template "${template}")
string(REPLACE "%(build_time)s" "${time_stamp}" template "${template}")
string(REPLACE "%(time_size)d" "80" template "${template}")
string(REPLACE "%(time_format)s" "%a, %d %b %Y %H:%M:%S +0000" template "${template}")
string(REPLACE "%(uname)s" "Microsoft Windows" template "${template}")
string(REPLACE "%(build_info)s" "Windows Build" template "${template}")
string(REPLACE "%(compiler_type)s" "${compiler_type}" template "${template}")
string(REPLACE "%(compiler_version)s" "${compiler_version}" template "${template}")
string(REPLACE "%(cc_flags)s" "${CMAKE_CXX_FLAGS}" template "${template}")
string(REPLACE "%(project_versions)s" "\"TODO\"" template "${template}")
string(REPLACE "%(xsd_version)s" "${xsd_version}" template "${template}")
file(WRITE "${version_cpp}" "${template}")


#generate GetChasteBuildRoot.cpp
file(WRITE "${build_root_cpp}" 
"#include \"ChasteBuildRoot.hpp\" 

const char* ChasteBuildRootDir() 
{ 
    return \"${CHASTE_SOURCE_ROOT}/\"; 
}

std::string ChasteComponentBuildDir(const std::string& rComponent)
{
    return std::string(ChasteBuildRootDir()) + rComponent + \"/${CHASTE_BUILD_DIR_NAME}/\" ;
}

std::string ChasteBuildDirName()
{
    return \"${CHASTE_BUILD_DIR_NAME}\";
}

std::string ChasteBuildType()
{
    return \"${CHASTE_BUILD_TYPE}\";
}

//This has been put here to satisfy an MSVC linker issue 
int __cdecl _purecall(void){return 0;}

")

#Build the global library
include_directories("${CHASTE_GLOBAL_INCLUDE_DIRS}")
file(GLOB_RECURSE CHASTE_GLOBAL_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp *.hpp)
add_library(global STATIC ${CHASTE_GLOBAL_SOURCES} ${version_cpp} ${build_root_cpp})
install(TARGETS global DESTINATION lib)

#Component testing for global
if(ENABLE_GLOBAL_TESTING)
    add_subdirectory(test)
endif()