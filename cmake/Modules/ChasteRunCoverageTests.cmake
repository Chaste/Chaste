execute_process(
        COMMAND nice "-n15" ${CTEST_COMMAND} "-j${Chaste_COVERAGE_CPUS}" "-L" "Continuous|Parallel" "--output-on-failure"
        WORKING_DIRECTORY ${Chaste_BINARY_DIR}
        ERROR_QUIET
        )
