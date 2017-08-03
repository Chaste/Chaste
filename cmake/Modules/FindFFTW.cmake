# Acknowledgement: This file is based on the FindParMETIS.cmake file from the GitHub repo
# https://github.com/jedbrown/cmake-modules (BSD-2-Clause)
#
# Find FFTW3 includes and libs
#
# This module sets the following variables:
#
#  FFTW_INCLUDES    - where to find fftw3.h
#  FFTW_LIBRARIES   - List of libraries when using FFTW.
#  FFTW_FOUND       - True if FFTW found.

if (FFTW_INCLUDES)
  # Already in cache, be silent
  set (FFTW_FIND_QUIETLY TRUE)
endif (FFTW_INCLUDES)

find_path (FFTW_INCLUDES fftw3.h HINTS $ENV{FFTW_ROOT} PATH_SUFFIXES include)

find_library (FFTW_LIBRARIES NAMES fftw3 HINTS $ENV{FFTW_ROOT} PATH_SUFFIXES lib)

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARIES FFTW_INCLUDES)

mark_as_advanced (FFTW_LIBRARIES FFTW_INCLUDES)
