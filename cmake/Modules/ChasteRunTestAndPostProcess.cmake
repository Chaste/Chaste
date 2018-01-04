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

if( NOT test_cmd )
   message( FATAL_ERROR "Variable test_cmd not defined" )
endif( NOT test_cmd )

if( NOT post_cmd )
   message( FATAL_ERROR "Variable post_cmd not defined" )
endif( NOT post_cmd )

# output_file contains the name of the output file the post_cmd will produce 
if( NOT output_file )
    message( FATAL_ERROR "Variable output_file not defined" )
endif( NOT output_file )

# convert the space-separated string to a list
separate_arguments( test_args )
separate_arguments( post_args )

if (env_var)
    set(ENV{${env_var}} ${env_var_value})
endif()

message("executing command:")
message("\t${test_cmd} ${test_args}")
execute_process(
    COMMAND ${test_cmd} ${test_args}
    RESULT_VARIABLE test_not_successful
    OUTPUT_VARIABLE std_out
    ERROR_VARIABLE std_err
    )

if( test_not_successful )
    message( SEND_ERROR "Test failed. \n
                        test command was ${test_cmd} \n
                        test arguements was ${test_args} \n
                        standard output was ${std_out} \n
                        standard error was ${std_err}")
else()
    message("executing post-processing command:")
    message("\t${post_cmd} ${post_args}")
    execute_process(
        COMMAND ${post_cmd} ${post_args}
        OUTPUT_FILE ${output_file}
        RESULT_VARIABLE post_not_successful
        )
    if (post_not_successful)
        message( SEND_ERROR "post_processing ${post_cmd} failed")
    endif()
endif()
