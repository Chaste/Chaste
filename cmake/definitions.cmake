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

# Compiler flag definitions common to Chaste & third party library builds


if(MSVC) # Applies only to MSVC
#For GUI configs. Change C, and CXX compiler flags dynamically to static, debug build.
#The overrides.cmake include takes care of non GUI builds.
foreach(flag_var
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
   if(${flag_var} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
   endif(${flag_var} MATCHES "/MD")
endforeach(flag_var)

# A bunch of compiler flags. Since we now propagate the compiler and linker flags to all external projects
# this is one place to set all the flags.
add_definitions(
    -Z7     # => embed debugging info in library as opposed to using an external .pdb database
    -wd4996 # => disable insecure api warnings
    -wd4267 # => disable "possible loss of data due to 'narrowing' conversion, e.g. size_t to int"
    -wd4290 # => disable warning "C++ exception specification ignored except to indicate a function is not __declspec(nothrow)"
    -wd4005 # => 'identifier' : macro redefinition
    -wd4018 # => 'expression' : signed/unsigned mismatch
    -wd4244 # => 'argument' : conversion from 'type1' to 'type2', possible loss of data
    -wd4101 # => 'identifier' : unreferenced local variable
    -wd4661 # => 'identifier' : no suitable definition provided for explicit template instantiation request
)

if(NOT "${PROJECT_NAME}" STREQUAL "Chaste")
    # Some third party libaries need further suppressions
    message("Adding further compiler flags for third-party libraries")
    add_definitions(
    #    -MTd    # => static debug build
        -Yu     # => use precompiled headers: this prevents a linker warning: "H5FDdirect.obj : warning LNK4221: This object
                # file does not define any previously undefined public symbols, so it will not be used by any link operation
                # that consumes this library".
        -wd4554 # => possible operator precedence error warning
        -wd4305 # => truncation from type1 to type2
        -wd4133 # => 'type' : incompatible types - from 'type1' to 'type2'
    )
endif()

#message("Set standard compiler flags.")
#message("  C: ${CMAKE_C_FLAGS}; debug: ${CMAKE_C_FLAGS_DEBUG}")
#message("  CXX: ${CMAKE_CXX_FLAGS}; debug: ${CMAKE_CXX_FLAGS_DEBUG}")
#message("  LINK: ${CMAKE_EXE_LINKER_FLAGS}; debug: ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
#get_directory_property(DirDefs DEFINITIONS)
#message("  DEFS: ${DirDefs}")
endif(MSVC)
