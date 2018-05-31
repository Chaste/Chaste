# Find Sprctra includes.  Spectra extends Eigen3 to solve dense and sparse eigenvalue problems.
#
# https://github.com/yixuan/spectra
#
# This module sets the following variables:
#
#  SPECTRA_INCLUDES   - where to find GenEigsComplexShiftSolver.h
#  SPECTRA_FOUND      - True if Spectra is found

if (SPECTRA_INCLUDES)
    # Already in cache, be silent
    set (SPECTRA_FIND_QUIETLY TRUE)
endif (SPECTRA_INCLUDES)

find_path (SPECTRA_INCLUDES GenEigsComplexShiftSolver.h HINTS $ENV{SPECTRA_ROOT} $ENV{SPECTRA_ROOT}/include PATH_SUFFIXES)

# handle the QUIETLY and REQUIRED arguments and set SPECTRA_FOUND to TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SPECTRA DEFAULT_MSG SPECTRA_INCLUDES)

mark_as_advanced (SPECTRA_INCLUDES)
