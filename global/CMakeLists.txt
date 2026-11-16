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


find_package(Chaste COMPONENTS none)

set(generate_dir ${CMAKE_CURRENT_BINARY_DIR}/src)
set(provenance_dir ${CMAKE_CURRENT_BINARY_DIR}/provenance)

##########################################################################
#  Generate Version.cpp and ChasteBuildInfo.cpp, fixed for the life of  #
#  this cmake configuration                                             #
##########################################################################

# With a multi-config generator (such as Ninja Multi-Config) the configuration is only chosen at
# build time, so anything depending on it is written as a generator expression below, and each
# configuration gets its own copy of the two generated files. A single-config generator (Unix
# Makefiles, plain Ninja) gets exactly the values CMAKE_BUILD_TYPE implies.
get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)

# Combine the CMAKE_CXX_FLAGS and CMAKE_CXX_FLAGS_<build_type> into one string with all the flags
string(STRIP "${CMAKE_CXX_FLAGS}" STRIPPED_CXX_FLAGS)
if (is_multi_config)
    set(ALL_CXX_FLAGS "${STRIPPED_CXX_FLAGS}")
    foreach (config ${CMAKE_CONFIGURATION_TYPES})
        string(TOUPPER "${config}" upper_config)
        string(STRIP "${CMAKE_CXX_FLAGS_${upper_config}}" config_flags)
        # Escape the characters that would otherwise end or split the generator expression
        string(REPLACE ">" "$<ANGLE-R>" config_flags "${config_flags}")
        string(REPLACE "," "$<COMMA>" config_flags "${config_flags}")
        string(APPEND ALL_CXX_FLAGS "$<$<CONFIG:${config}>: ${config_flags}>")
    endforeach ()
    set(Chaste_BUILD_DIR_NAME "$<CONFIG>")
else ()
    string(TOUPPER ${CMAKE_BUILD_TYPE} UPPERCASE_BUILD_TYPE)
    string(STRIP "${CMAKE_CXX_FLAGS_${UPPERCASE_BUILD_TYPE}}" STRIPPED_BUILD_TYPE_FLAGS)
    set(ALL_CXX_FLAGS "${STRIPPED_CXX_FLAGS} ${STRIPPED_BUILD_TYPE_FLAGS}")
    set(Chaste_BUILD_DIR_NAME "${CMAKE_CFG_INTDIR}")
endif ()
set(Chaste_BUILD_TYPE "$<CONFIG>")
set(Chaste_DEFAULT_TEST_OUTPUT "${Chaste_BINARY_DIR}/testoutput")
set(time_size 80)
set(time_format "%a, %d %b %Y %H:%M:%S %z")

# The licence text is compiled into Version.cpp as a raw string literal, so the banner printed
# by Chaste executables never depends on finding a file at runtime. The "BSD 3-Clause License."
# heading and the "This file is part of Chaste." line (which belongs in source file headers) are
# dropped, so that banner is unchanged. LICENSE is a configure dependency, so editing it (such as
# the annual copyright update) reconfigures and recompiles Version.cpp automatically.
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${Chaste_SOURCE_DIR}/LICENSE")
file(READ "${Chaste_SOURCE_DIR}/LICENSE" licence_text)
string(REPLACE "BSD 3-Clause License.\n\n" "" licence_text "${licence_text}")
string(REPLACE "This file is part of Chaste.\n\n" "" licence_text "${licence_text}")
string(STRIP "${licence_text}" licence_text)

# Every value substituted into these two files (version numbers, compiler info, build type,
# xsd and chaste_codegen versions, and so on) is fixed for the life of a configuration, so
# neither has any custom target or tie to a build-time step. configure_file() substitutes the
# cmake variables, then file(GENERATE) evaluates the generator expressions, writing one output
# per configuration. Both only touch an output whose content actually changed, so an ordinary
# build, or a reconfigure that changes nothing, never recompiles either file.
foreach (template Version.cpp ChasteBuildInfo_cmake.cpp)
    string(REPLACE "_cmake" "" generated "${template}")
    configure_file(
        "${Chaste_SOURCE_DIR}/global/src/${template}.in"
        "${generate_dir}/${generated}.genex"
        )
    file(GENERATE
        OUTPUT "${generate_dir}/$<CONFIG>/${generated}"
        INPUT "${generate_dir}/${generated}.genex"
        )
endforeach ()

##########################################################################
#  Generate build-time provenance: build timestamp and git state        #
##########################################################################

# generate_provenance has no tracked OUTPUT, so cmake always considers it out of date: its
# recipe (cmake/Modules/ChasteGenerateProvenance.cmake) runs every time chaste_global is built,
# giving the build timestamp and git state (Chaste's own, and that of each checked out project)
# genuine per-build freshness. Nothing here is compiled, so refreshing these files never marks
# chaste_global, or anything that depends on it, as needing to recompile or relink. Each
# configuration writes to its own directory, named as ChasteBuildDirName() reports it at runtime:
# "." for a single-config generator, so just provenance/ itself, or the configuration name for a
# multi-config generator, so that building one configuration never overwrites the build time or
# git state reported by another.
add_custom_target(generate_provenance
    COMMAND ${CMAKE_COMMAND}
        "-DChaste_SOURCE_DIR=${Chaste_SOURCE_DIR}"
        "-DChaste_PROJECTS=${Chaste_PROJECTS}"
        "-Dprovenance_dir=${provenance_dir}/${Chaste_BUILD_DIR_NAME}"
        "-DCMAKE_MODULE_PATH=${CMAKE_MODULE_PATH}"
        "-P" "${Chaste_SOURCE_DIR}/cmake/Modules/ChasteGenerateProvenance.cmake"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Refreshing build provenance"
    VERBATIM
    )

chaste_do_component(global ${generate_dir}/$<CONFIG>/Version.cpp ${generate_dir}/$<CONFIG>/ChasteBuildInfo.cpp)

add_dependencies(chaste_global generate_provenance)

# Installed alongside the libraries and headers (matching the existing cxxtest/ and python/
# folders below) so that an installed, out-of-tree Chaste (via ChasteConfig.cmake) still has
# provenance to read, resolved at runtime via ChasteInstallRootDir().
install(DIRECTORY "${provenance_dir}/${Chaste_BUILD_DIR_NAME}/"
        DESTINATION lib/chaste/provenance
        COMPONENT Config)
