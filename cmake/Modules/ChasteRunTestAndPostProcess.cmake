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
