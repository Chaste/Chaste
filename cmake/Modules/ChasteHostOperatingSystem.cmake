# Determine the OS infomration and print it

set(DIST_NAME "OS name not identified")

# Determine the "name" of the OS (CMake gives us the kernel version for free)
if (UNIX AND NOT (APPLE OR CYGWIN))
    if (EXISTS "/etc/os-release")
        execute_process(
                COMMAND cat /etc/os-release
                COMMAND grep -i PRETTY_NAME
                COMMAND cut -d = -f2
                COMMAND tr -d \"
                OUTPUT_VARIABLE DIST_NAME
        )
    elseif (EXISTS "/etc/lsb-release")
        execute_process(
                COMMAND cat /etc/lsb-release
                COMMAND grep -i DISTRIB_DESCRIPTION
                COMMAND cut -d = -f2
                COMMAND tr -d \"
                OUTPUT_VARIABLE DIST_NAME
        )
    endif()
elseif (APPLE)
    execute_process(COMMAND sw_vers -productVersion OUTPUT_VARIABLE DIST_NAME)
elseif (WIN32)
    # Annoying to do on windows: systeminfo takes AGES to run
endif()

string(STRIP ${DIST_NAME} DIST_NAME)

message(STATUS "Operating system detected as...")
message(STATUS "\t... ${DIST_NAME}: ${CMAKE_HOST_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_VERSION}")
