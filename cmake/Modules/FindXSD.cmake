# This file is based on the FindXSD.cmake file:
#     http://wiki.codesynthesis.com/uploads/8/86/FindXSD.cmake.gz
# 
# Attempt to find the xsd application in various places. If found, the full
# path will be in XSD_EXECUTABLE. Look in the usual locations, as well as in
# the 'bin' directory in the path given in the XSD_ROOT environment variable.

# Mono also has a utility called xsd, so try to find xsdcxx first


# Prefer $ENV{XSD_ROOT}
FIND_PROGRAM(XSD_EXECUTABLE NAMES xsdcxx xsd
             HINTS $ENV{XSD_ROOT}/bin
             NO_DEFAULT_PATH)


# Now try other system paths
FIND_PROGRAM(XSD_EXECUTABLE NAMES xsdcxx xsd
             PATHS /usr/local/xsd-3.2.0-i686-macosx/bin
             /usr/local/xsd-3.2.0-x86_64-linux-gnu/bin
             /usr/local/bin
             /usr/bin
             /opt/xsd-3.2.0-i686-macosx/bin
             /opt/xsd-3.2.0-x86_64-linux-gnu/bin
             /usr/bin
             ENV PATH)

if (XSD_EXECUTABLE)
    GET_FILENAME_COMPONENT(XSD_BIN_DIR "${XSD_EXECUTABLE}" PATH)
    GET_FILENAME_COMPONENT(XSD_ROOT_DIR "${XSD_BIN_DIR}" PATH)

    # Obtain the include directory that one can use with INCLUDE_DIRECTORIES() to
    # access the xsd include files.
    FIND_PATH(XSD_INCLUDE_DIR xsd/cxx/version.hxx
              HINTS "$ENV{XSD_ROOT}/libxsd"
              PATHS /usr/include
              /usr/local/opt/libxsd)

endif (XSD_EXECUTABLE)


# General CMake package configuration.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XSD DEFAULT_MSG XSD_EXECUTABLE XSD_INCLUDE_DIR)

if (XSD_FOUND)
    SET(XSD_INCLUDE_DIRS ${XSD_INCLUDE_DIR})
else (XSD_FOUND)
    SET(XSD_INCLUDE_DIRS)
endif (XSD_FOUND)

MARK_AS_ADVANCED(XSD_INCLUDE_DIR XSD_EXECUTABLE)


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
macro (XSD_SCHEMA NAME FILE)

    # Make a full path from the soource directory
    SET(xs_SRC "${FILE}")

    # XSD will generate two or three C++ files (*.cxx,*.hxx,*.ixx). Get the
    # destination file path sans any extension and then build paths to the
    # generated files.
    GET_FILENAME_COMPONENT(xs_FILE "${FILE}" NAME_WE)
    file(RELATIVE_PATH xs_FILE_REL "${CMAKE_SOURCE_DIR}" "${FILE}")
    set(xs_OUT_TMP "${CMAKE_BINARY_DIR}/${xs_FILE_REL}")
    GET_FILENAME_COMPONENT(xs_OUT_DIR "${xs_OUT_TMP}" PATH)
    SET(xs_CPP "${xs_OUT_DIR}/${xs_FILE}.cpp")
    SET(xs_HPP "${xs_OUT_DIR}/${xs_FILE}.hpp")
    #SET( xs_IPP "${xs_FILE_DIR}/${xs_FILE}.ipp" )

    # Add the source files to the NAME variable, which presumably will be used to
    # define the source of another target.
    LIST(APPEND ${NAME} ${xs_CPP} ${xs_HPP})

    # Set up a generator for the output files from the given schema file using
    # the XSD cxx-tree command.
    ADD_CUSTOM_COMMAND(OUTPUT "${xs_CPP}" "${xs_HPP}"
                       COMMAND ${XSD_EXECUTABLE}
                       ARGS "cxx-tree" "--output-dir" "${xs_OUT_DIR}" ${ARGN} ${xs_SRC}
                       DEPENDS ${xs_SRC}
                       COMMENT "Processing XML schema ${xs_FILE_REL}"
                       VERBATIM)

    # Don't fail if a generated file does not exist.
    SET_SOURCE_FILES_PROPERTIES("${xs_CPP}" "${xs_HPP}" PROPERTIES GENERATED TRUE)

endmacro (XSD_SCHEMA)
