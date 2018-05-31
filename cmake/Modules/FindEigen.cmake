# Find Eigen3 includes
#
# This module sets the following variables:
#
#  EIGEN_INCLUDES   - where to find signature_of_eigen3_matrix_library
#  EIGEN_FOUND      - True if Eigen3 is found

if (EIGEN_INCLUDES)
    # Already in cache, be silent
    set (EIGEN_FIND_QUIETLY TRUE)
endif (EIGEN_INCLUDES)

find_path (EIGEN_INCLUDES signature_of_eigen3_matrix_library
           HINTS $ENV{EIGEN3_ROOT} /usr/include/eigen3
           PATH_SUFFIXES Eigen)

# handle the QUIETLY and REQUIRED arguments and set EIGEN_FOUND to TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (EIGEN DEFAULT_MSG EIGEN_INCLUDES)

mark_as_advanced (EIGEN_INCLUDES)
