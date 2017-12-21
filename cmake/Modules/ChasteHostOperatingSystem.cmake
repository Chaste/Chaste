# Determine the OS infomration and print it

execute_process(
        COMMAND ${PYTHON_EXECUTABLE} "${Chaste_SOURCE_DIR}/cmake/Modules/ChasteHostOperatingSystem.py"
        OUTPUT_VARIABLE DIST_NAME
)

message(STATUS "Operating system detected as...")
message(STATUS "\t... ${DIST_NAME}: ${CMAKE_HOST_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_VERSION}")