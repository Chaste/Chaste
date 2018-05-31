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
        "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
        "-DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}"
        "-DMSVC=${MSVC}"
        "-P" "${Chaste_SOURCE_DIR}/cmake/Modules/ChasteGenerateVersionAndBuildInfo.cmake"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    DEPENDS timekeeper
    COMMENT "Generating Build Config and Version Info"
    VERBATIM 
    )


chaste_do_component(global ${generate_dir}/Version.cpp ${generate_dir}/ChasteBuildInfo.cpp)
