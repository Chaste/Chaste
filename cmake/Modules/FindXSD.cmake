# This file is based on the FindXSD.cmake file:
#     http://wiki.codesynthesis.com/uploads/8/86/FindXSD.cmake.gz
# 
# Attempt to find the xsd application in various places. If found, the full
# path will be in XSD_EXECUTABLE. Look in the usual locations, as well as in
# the 'bin' directory in the path given in the XSD_ROOT environment variable.

# Mono also has a utility called xsd, so try to find xsdcxx first


# Prefer $ENV{XSD_ROOT}
find_program (XSD_EXECUTABLE NAMES xsdcxx xsd
              HINTS $ENV{XSD_ROOT}/bin
              NO_DEFAULT_PATH)


# Now try other system paths
find_program (XSD_EXECUTABLE NAMES xsdcxx xsd
              PATHS /usr/local/xsd-3.2.0-i686-macosx/bin
              /usr/local/xsd-3.2.0-x86_64-linux-gnu/bin
              /usr/local/bin
              /usr/bin
              /opt/xsd-3.2.0-i686-macosx/bin
              /opt/xsd-3.2.0-x86_64-linux-gnu/bin
              /usr/bin
              ENV PATH)

if (XSD_EXECUTABLE)
    get_filename_component (XSD_BIN_DIR "${XSD_EXECUTABLE}" PATH)
    get_filename_component (XSD_ROOT_DIR "${XSD_BIN_DIR}" PATH)

    # Obtain the include directory that one can use with INCLUDE_DIRECTORIES() to
    # access the xsd include files.
    find_path (XSD_INCLUDE_DIR xsd/cxx/version.hxx
               HINTS "$ENV{XSD_ROOT}/libxsd"
               PATHS /usr/include
               /usr/local/opt/libxsd)

    # Idenfity the XSD version from the executable
    # Look for X.Y.Z in the result of xsd --version
    execute_process (COMMAND ${XSD_EXECUTABLE} "--version" OUTPUT_VARIABLE xsd_version_full)

    string (REGEX REPLACE
            ".*([0-9]+\\.[0-9]+\\.[0-9]+).*"
            "\\1"
            xsd_version
            "${xsd_version_full}")

    if (NOT xsd_version MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
        set (xsd_version "undertermined")
        message (WARNING "XSD found, but version undetermined")
    endif ()

endif (XSD_EXECUTABLE)


# General CMake package configuration.
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (XSD DEFAULT_MSG XSD_EXECUTABLE XSD_INCLUDE_DIR)

if (XSD_FOUND)
    set (XSD_INCLUDE_DIRS ${XSD_INCLUDE_DIR})
else (XSD_FOUND)
    set (XSD_INCLUDE_DIRS)
endif (XSD_FOUND)

mark_as_advanced (XSD_INCLUDE_DIR XSD_EXECUTABLE)


# Macro that attempts to generate C++ files from an XML schema. The NAME
# argument is the name of the CMake variable to use to store paths to the
# derived C++ source file. The FILE argument is the path of the schema file to
# process. Additional arguments should be XSD command-line options.
#
# Example:
#
# XSD_SCHEMA( FOO_SRCS Foo.xsd --root-element-first --generate-serialization )
#
# On return, FOO_SRCS will contain Foo.cxx.
#
macro (xsd_schema NAME FILE)

    # Make a full path from the soource directory
    set (xs_SRC "${FILE}")

    # XSD will generate two or three C++ files (*.cxx,*.hxx,*.ixx). Get the
    # destination file path sans any extension and then build paths to the
    # generated files.
    get_filename_component (xs_FILE "${FILE}" NAME_WE)
    file (RELATIVE_PATH xs_FILE_REL "${CMAKE_SOURCE_DIR}" "${FILE}")
    set (xs_OUT_TMP "${CMAKE_BINARY_DIR}/${xs_FILE_REL}")
    get_filename_component (xs_OUT_DIR "${xs_OUT_TMP}" PATH)
    set (xs_CPP "${xs_OUT_DIR}/${xs_FILE}.cpp")
    set (xs_HPP "${xs_OUT_DIR}/${xs_FILE}.hpp")
    #SET( xs_IPP "${xs_FILE_DIR}/${xs_FILE}.ipp" )

    # Add the source files to the NAME variable, which presumably will be used to
    # define the source of another target.
    list (APPEND ${NAME} ${xs_CPP} ${xs_HPP})

    # Set up a generator for the output files from the given schema file using
    # the XSD cxx-tree command.
    add_custom_command (OUTPUT "${xs_CPP}" "${xs_HPP}"
                        COMMAND ${XSD_EXECUTABLE}
                        ARGS "cxx-tree" "--output-dir" "${xs_OUT_DIR}" "--std" "c++11" ${ARGN} ${xs_SRC}
                        DEPENDS ${xs_SRC}
                        COMMENT "Processing XML schema ${xs_FILE_REL}"
                        VERBATIM)

    # Don't fail if a generated file does not exist.
    set_source_files_properties ("${xs_CPP}" "${xs_HPP}" PROPERTIES GENERATED TRUE)

endmacro (xsd_schema)
