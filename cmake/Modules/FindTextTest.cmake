# 
# Attempt to find the texttest application in various places. If found, the full
# path will be in TESTTEST_PY. Look in the usual locations, as well as in
# the 'bin' directory in the path given in the TEXTTEST_ROOT environment variable.
#

FIND_PROGRAM( TEXTTEST_PY NAMES texttest.py texttest_release.py texttest
              HINTS $ENV{TEXTTEST_ROOT}/bin
			  PATHS /usr/local/texttest-3.19/bin
			  		/usr/local/bin
			  		/usr/bin
					/opt/texttest-3.19/bin
					ENV PATH )

IF( TEXTTEST_PY )
  GET_FILENAME_COMPONENT( TEXTTEST_BIN_DIR "${TEXTTEST_PY}" PATH )
  GET_FILENAME_COMPONENT( TEXTTEST_ROOT_DIR "${TEXTTEST_BIN_DIR}" PATH )
ENDIF( TEXTTEST_PY )

#
# General CMake package configuration.
#
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( TextTest DEFAULT_MSG TEXTTEST_PY)

MARK_AS_ADVANCED( TEXTTEST_PY )
