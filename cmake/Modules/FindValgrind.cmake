# Find Valgrind
#
# Sets the following cmake variables
#  VALGRIND_COMMAND, the valgrind executable


find_program(VALGRIND_COMMAND valgrind
            PATH /usr/bin /usr/local/bin)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(Valgrind DEFAULT_MSG
    VALGRIND_COMMAND)
