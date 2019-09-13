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


find_package(Chaste COMPONENTS ${Chaste_DEPENDS_heart})

set(Chaste_heart_SOURCES "")

file(GLOB xsd_config_files src/io/*.xsd)
foreach(xsd_config ${xsd_config_files})
    xsd_schema(Chaste_heart_SOURCES ${xsd_config} 
                "--generate-serialization"
                "--hxx-suffix" ".hpp" 
                "--cxx-suffix" ".cpp"
                "--prologue-file" "${CMAKE_CURRENT_SOURCE_DIR}/src/io/XsdPrologue.txt"
                "--epilogue-file" "${CMAKE_CURRENT_SOURCE_DIR}/src/io/XsdEpilogue.txt"
                "--namespace-regex" "X.* $Xchaste::parametersX"
                "--namespace-regex" "X.* https://chaste.comlab.ox.ac.uk/nss/parameters/(.+)Xchaste::parameters::v$1X"
                )

endforeach()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/io)
set(Chaste_heart_INCLUDE_DIRS ${Chaste_heart_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/src/io)

# Make dynamic directory for dynamic cellml compilation
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dynamic)

chaste_do_component(heart ${Chaste_heart_SOURCES})

# make sure to link dl libraries
target_link_libraries(chaste_heart LINK_PUBLIC ${CMAKE_DL_LIBS})

# Compile dynamic models
set(dyn_sources dynamic/DynamicallyLoadableLr91.cpp dynamic/NotACellModel.cpp)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dynamic)
foreach(dyn_source ${dyn_sources})
    get_filename_component(dyn_source_name ${dyn_source} NAME_WE)
    add_library(${dyn_source_name} SHARED ${dyn_source})
    set_property(TARGET ${dyn_source_name} PROPERTY CXX_STANDARD 14)
    set_property(TARGET ${dyn_source_name} PROPERTY CXX_STANDARD_REQUIRED ON)
    set_property(TARGET ${dyn_source_name} PROPERTY CMAKE_CXX_EXTENSIONS OFF)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        target_link_libraries(${dyn_source_name} LINK_PRIVATE "-Wl,-undefined,dynamic_lookup")
        #set_target_properties(${dyn_source_name} PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
    else()
        #target_link_libraries(${dyn_source_name} LINK_PRIVATE "-Wl,--undefined")
        set_target_properties(${dyn_source_name} PROPERTIES LINK_FLAGS "-Wl,--undefined")
    endif()

    # add dynamic models to the heart_tests and Continuous targets
    add_dependencies(heart ${dyn_source_name})
    add_dependencies(Continuous ${dyn_source_name})
endforeach(dyn_source)
