execute_process(
        COMMAND ${CTEST_COMMAND} "-L" "Continuous|Parallel" "--output-on-failure" 
        WORKING_DIRECTORY ${Chaste_BINARY_DIR}
        ERROR_QUIET
        )
