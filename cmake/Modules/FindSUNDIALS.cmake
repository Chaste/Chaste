# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson,
#                         Andreas Puettmann¹, Sebastian Schnittert¹,
#                         Samuel Leweke¹
#                                      
#    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Find SUNDIALS, the SUite of Nonlinear and DIfferential/ALgebraic equation Solvers.
#
# The module will optionally accept the COMPONENTS argument.  If no COMPONENTS
# are specified, then the find module will default to find all the SUNDIALS
# libraries.  If one or more COMPONENTS are specified, the module will attempt to
# find the specified components.
#
# On UNIX systems, this module will read the variable SUNDIALS_USE_STATIC_LIBRARIES
# to determine whether or not to prefer a static link to a dynamic link for SUNDIALS
# and all of it's dependencies.  To use this feature, make sure that the
# SUNDIALS_USE_STATIC_LIBRARIES variable is set before the call to find_package.
#
# To provide the module with a hint about where to find your SUNDIALS installation,
# you can set the environment variable SUNDIALS_ROOT. The FindSUNDIALS module will
# then look in this path when searching for SUNDIALS paths and libraries.
#
# This module will define the following variables:
#  SUNDIALS_INCLUDE_DIRS - Location of the SUNDIALS includes
#  SUNDIALS_VERSION_MAJOR 
#  SUNDIALS_VERSION_MINOR
#  SUNDIALS_VERSION_SUBMINOR
#  SUNDIALS_FOUND - true if SUNDIALS was found on the system
#  SUNDIALS_LIBRARIES - Required libraries for all requested components

include(FindPackageHandleStandardArgs)


# Option that allows users to specify custom SUNDIALS path
if (NOT "$ENV{SUNDIALS_ROOT}" STREQUAL "")
    list (APPEND CMAKE_INCLUDE_PATH "$ENV{SUNDIALS_ROOT}")
    list (APPEND CMAKE_LIBRARY_PATH "$ENV{SUNDIALS_ROOT}")
#else ()
#    message (STATUS "No environment variable 'SUNDIALS_ROOT' found,\n\tyou may specify custom path in 'PATH_SUNDIALS_ROOT' cache variable")
endif ()


#if (NOT PATH_SUNDIALS_ROOT)
#    set (PATH_SUNDIALS_ROOT "" CACHE STRING "Optional path to the SUNDIALS lib and include direrctory" FORCE)
#else ()
#    list (APPEND CMAKE_INCLUDE_PATH "${PATH_SUNDIALS_ROOT}")
#    list (APPEND CMAKE_LIBRARY_PATH "${PATH_SUNDIALS_ROOT}")
#endif ()


# List of user definable search paths
set (SUNDIALS_USER_PATHS
#    ${PATH_SUNDIALS_ROOT}
    $ENV{HOME}/.local
    $ENV{HOME}
    $ENV{HOME}/libs/sundials-2.5.0
)


# List of the valid SUNDIALS components
set( SUNDIALS_VALID_COMPONENTS
    sundials_cvode
    sundials_cvodes
    sundials_ida
    sundials_idas
    sundials_kinsol
    sundials_nvecserial
)


if( NOT SUNDIALS_FIND_COMPONENTS )
    set( SUNDIALS_WANT_COMPONENTS ${SUNDIALS_VALID_COMPONENTS} )
else()
    # add the extra specified components, ensuring that they are valid.
    foreach( _COMPONENT ${SUNDIALS_FIND_COMPONENTS} )
        string (TOLOWER ${_COMPONENT} _COMPONENT_LOWER)
        list( FIND SUNDIALS_VALID_COMPONENTS ${_COMPONENT_LOWER} COMPONENT_LOCATION )
        if( ${COMPONENT_LOCATION} EQUAL -1 )
            message( FATAL_ERROR
                "\"${_COMPONENT_LOWER}\" is not a valid SUNDIALS component." )
        else()
            list( APPEND SUNDIALS_WANT_COMPONENTS ${_COMPONENT_LOWER} )
        endif()
    endforeach()
endif()


# find the SUNDIALS include directories
find_path( SUNDIALS_INCLUDE_DIR sundials_types.h
    ENV
        SUNDIALS_ROOT
    PATHS
        ${SUNDIALS_USER_PATHS}
    PATH_SUFFIXES
        include
        include/sundials
)
if (SUNDIALS_INCLUDE_DIR)
    GET_FILENAME_COMPONENT( root_dir "${SUNDIALS_INCLUDE_DIR}" PATH)
    set( SUNDIALS_INCLUDE_DIRS
        "${root_dir}"
        "${SUNDIALS_INCLUDE_DIR}"
    )
endif()

# get SUNDIALS version
if (SUNDIALS_INCLUDE_DIR)
    file(STRINGS "${SUNDIALS_INCLUDE_DIR}/sundials_config.h" _SUNDIALS_VERSION_HPP_CONTENTS REGEX "#define SUNDIALS_(PACKAGE_|)VERSION ")
    string(REGEX MATCH "\"([0-9]*)\\." _DUMMY "${_SUNDIALS_VERSION_HPP_CONTENTS}")
    set(SUNDIALS_VERSION_MAJOR ${CMAKE_MATCH_1})
    string(REGEX MATCH "\\.([0-9]*)\\." _DUMMY "${_SUNDIALS_VERSION_HPP_CONTENTS}")
    set(SUNDIALS_VERSION_MINOR ${CMAKE_MATCH_1})
    string(REGEX MATCH "\\.([0-9]*)\"" _DUMMY "${_SUNDIALS_VERSION_HPP_CONTENTS}")
    set(SUNDIALS_VERSION_SUBMINOR ${CMAKE_MATCH_1})
endif(SUNDIALS_INCLUDE_DIR)


# find the SUNDIALS libraries
foreach( LIB ${SUNDIALS_WANT_COMPONENTS} )
    if( UNIX AND SUNDIALS_USE_STATIC_LIBRARIES )
        # According to bug 1643 on the CMake bug tracker, this is the
        # preferred method for searching for a static library.
        # See http://www.cmake.org/Bug/view.php?id=1643.  We search
        # first for the full static library name, but fall back to a
        # generic search on the name if the static search fails.
        set( THIS_LIBRARY_SEARCH lib${LIB}.a ${LIB} )
    else()
        set( THIS_LIBRARY_SEARCH ${LIB} )
    endif()

    find_library( SUNDIALS_${LIB}_LIBRARY
        NAMES ${THIS_LIBRARY_SEARCH}
        ENV
            SUNDIALS_ROOT
        PATHS
            ${SUNDIALS_USER_PATHS}
        PATH_SUFFIXES
            lib
            Lib
    )

    list( APPEND SUNDIALS_LIBRARIES ${SUNDIALS_${LIB}_LIBRARY} )
    mark_as_advanced(SUNDIALS_${LIB}_LIBRARY)
endforeach()

# find blas and lapack
# (sundials on ubuntu 12.04 need to explicitly link against lapack)
if (BUILD_SHARED_LIBS)
    find_package(LAPACK)
else()
    set(BLA_STATIC ON)
    find_package(LAPACK)
endif()
    

list(APPEND SUNDIALS_LIBRARIES ${LAPACK_LIBRARIES})

find_package_handle_standard_args( SUNDIALS DEFAULT_MSG
    SUNDIALS_LIBRARIES
    SUNDIALS_INCLUDE_DIRS
)

mark_as_advanced(
    SUNDIALS_LIBRARIES
    SUNDIALS_INCLUDE_DIR
    SUNDIALS_INCLUDE_DIRS
)
