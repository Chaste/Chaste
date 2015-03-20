# 
# Attempt to find the xsd application in various places. If found, the full
# path will be in XSD_EXECUTABLE. Look in the usual locations, as well as in
# the 'bin' directory in the path given in the XSD_ROOT environment variable.
#
FIND_PROGRAM( XSD_EXECUTABLE xsd
		   	  HINTS ${RWSL_DEPS}/xsd/bin $ENV{XSD_ROOT}/bin
			  PATHS /usr/local/xsd-3.2.0-i686-macosx/bin
			  		/usr/local/xsd-3.2.0-x86_64-linux-gnu/bin
			  		/usr/local/bin
					/opt/xsd-3.2.0-i686-macosx/bin
			  		/opt/xsd-3.2.0-x86_64-linux-gnu/bin
			  		/usr/bin
					ENV PATH )

IF( XSD_EXECUTABLE )

  # 
  # Obtain the include directory that one can use with INCLUDE_DIRECTORIES() to
  # access the xsd include files.
  #
  GET_FILENAME_COMPONENT( XSD_BIN_DIR "${XSD_EXECUTABLE}" PATH )
  GET_FILENAME_COMPONENT( XSD_ROOT_DIR "${XSD_BIN_DIR}" PATH )
  SET( XSD_INCLUDE_DIR "${XSD_ROOT_DIR}/libxsd" )
ENDIF( XSD_EXECUTABLE )

#
# General CMake package configuration.
#
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( XSD DEFAULT_MSG XSD_EXECUTABLE
								   XSD_INCLUDE_DIR )
IF( XSD_FOUND )
  SET( XSD_INCLUDE_DIRS ${XSD_INCLUDE_DIR} )
ELSE( XSD_FOUND )
  SET( XSD_INCLUDE_DIRS )
ENDIF( XSD_FOUND )

MARK_AS_ADVANCED( XSD_INCLUDE_DIR XSD_EXECUTABLE )

# 
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
MACRO( XSD_SCHEMA NAME FILE )

  #
  # Make a full path from the soource directory
  #
  SET( xs_SRC "${FILE}" )

  # 
  # XSD will generate two or three C++ files (*.cxx,*.hxx,*.ixx). Get the
  # destination file path sans any extension and then build paths to the
  # generated files.
  #
  GET_FILENAME_COMPONENT( xs_FILE "${FILE}" NAME_WE )
  SET( xs_CXX "${CMAKE_CURRENT_BINARY_DIR}/${xs_FILE}.cxx" )
  SET( xs_HXX "${CMAKE_CURRENT_BINARY_DIR}/${xs_FILE}.hxx" )
  SET( xs_IXX "${CMAKE_CURRENT_BINARY_DIR}/${xs_FILE}.ixx" )

  #
  # Add the source files to the NAME variable, which presumably will be used to
  # define the source of another target.
  #
  LIST( APPEND ${NAME} ${xs_CXX} )

  #
  # Set up a generator for the output files from the given schema file using
  # the XSD cxx-tree command.
  #
  ADD_CUSTOM_COMMAND( OUTPUT "${xs_CXX}" "${xs_HXX}" "${xs_IXX}"
  					  COMMAND ${XSD_EXECUTABLE}
					  ARGS "cxx-tree" ${ARGN} ${xs_SRC}
					  DEPENDS ${xs_SRC} )

  #
  # Don't fail if a generated file does not exist.
  #
  SET_SOURCE_FILES_PROPERTIES( "${xs_CXX}" "${xs_HXX}" "${xs_IXX}"
  							   PROPERTIES GENERATED TRUE )

ENDMACRO( XSD_SCHEMA )
