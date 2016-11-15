execute_process(
        COMMAND ${CTEST_COMMAND} "-L" "Continuous" "--output-on-failure" 
        COMMAND ${CTEST_COMMAND} "-L" "Parallel" "--output-on-failure" 
        WORKING_DIRECTORY ${Chaste_BINARY_DIR}
        ERROR_QUIET
        )
