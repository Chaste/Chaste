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


find_package(Chaste COMPONENTS none)

set(generate_dir ${CMAKE_CURRENT_BINARY_DIR}/src)

#The generated timekeeper.cpp code below keeps track of build timestamp.
#It is built and executed prior to starting a build and prints the timestamp
#in a given format. This timestamp is used by Version.cpp, which is also auto-generated.
file(WRITE ${generate_dir}/timekeeper.cpp
"#include <iostream>
#include <fstream>
#include <ctime>
int main( )
{
   time_t now = time(0);
   tm* loc_time = localtime(&now);
   char buffer[80];
   strftime(buffer, 80, \"%a, %d %b %Y %H:%M:%S +0000\", loc_time);

   std::ofstream timestampFile;
   timestampFile.open (\"build_timestamp\");
   timestampFile << buffer;
   timestampFile.close();
   return 0;
}
")

add_executable(timekeeper "${generate_dir}/timekeeper.cpp")

###################################################################
#  Setup command to generate Version and BuildInfo at build-time  #
###################################################################

# Combine the CMAKE_CXX_FLAGS and CMAKE_CXX_FLAGS_<build_type> into one string with all the flags
string(TOUPPER ${CMAKE_BUILD_TYPE} UPPERCASE_BUILD_TYPE)
string(STRIP ${CMAKE_CXX_FLAGS} STRIPPED_CXX_FLAGS)
string(STRIP ${CMAKE_CXX_FLAGS_${UPPERCASE_BUILD_TYPE}} STRIPPED_BUILD_TYPE_FLAGS)
set(ALL_CXX_FLAGS "${STRIPPED_CXX_FLAGS} ${STRIPPED_BUILD_TYPE_FLAGS}")

add_custom_command(OUTPUT 
    ${generate_dir}/Version.cpp
    ${generate_dir}/ChasteBuildInfo.cpp
    COMMAND ${CMAKE_COMMAND} 
        "-DChaste_SOURCE_DIR=${Chaste_SOURCE_DIR}"
        "-DChaste_BINARY_DIR=${Chaste_BINARY_DIR}"
        "-DChaste_BUILD_DIR_NAME=${CMAKE_CFG_INTDIR}"
        "-DChaste_BUILD_TYPE=$<CONFIGURATION>"
        "-Dtimekeeper_exe=$<TARGET_FILE:timekeeper>"
        "-Dxsd_version=${xsd_version}"
        "-Dgenerate_dir=${generate_dir}"
        "-DChaste_UPDATE_PROVENANCE=${Chaste_UPDATE_PROVENANCE}"
        "-DChaste_PROJECTS=${Chaste_PROJECTS}"
        "-DChaste_VERSION_MAJOR=${Chaste_VERSION_MAJOR}"
        "-DChaste_VERSION_MINOR=${Chaste_VERSION_MINOR}"
        "-DCMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME}"
        "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
        "-DCMAKE_CXX_COMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}"
        "-DALL_CXX_FLAGS=${ALL_CXX_FLAGS}"
        "-DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}"
        "-DMSVC=${MSVC}"
        "-P" "${Chaste_SOURCE_DIR}/cmake/Modules/ChasteGenerateVersionAndBuildInfo.cmake"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS timekeeper
    COMMENT "Generating Build Config and Version Info"
    VERBATIM 
    )


chaste_do_component(global ${generate_dir}/Version.cpp ${generate_dir}/ChasteBuildInfo.cpp)
