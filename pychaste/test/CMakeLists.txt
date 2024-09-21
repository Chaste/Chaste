# Copyright (c) 2005-2024, University of Oxford.
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

if(NOT Chaste_ENABLE_PYCHASTE)
    return()
endif()

project(testpychaste)
chaste_do_test_component(pychaste)

# Add TestWrapperGeneration
cppwg_command("${CMAKE_CURRENT_BINARY_DIR}/wrappers.gen")
add_test(
    NAME TestPyWrapperGeneration
    COMMAND sh -c
        "${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/wrappers &&
         ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/wrappers.gen &&
         ${CMAKE_COMMAND} -E copy_directory \
             ${CMAKE_SOURCE_DIR}/pychaste/dynamic/wrappers \
             ${CMAKE_CURRENT_BINARY_DIR}/wrappers &&
         ${CPPWG_COMMAND_STR}"
)
set_tests_properties(TestPyWrapperGeneration PROPERTIES FIXTURES_SETUP GEN_WRAPPERS)

# TestPyWrapperChanges depends on output from TestWrapperGeneration 
set_tests_properties(TestPyWrapperChanges PROPERTIES FIXTURES_REQUIRED GEN_WRAPPERS)
set_tests_properties(TestPyWrapperChanges PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
