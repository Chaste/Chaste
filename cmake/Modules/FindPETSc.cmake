# - Try to find PETSc
# Acknowledgement: This file is based on the FindPETSc.cmake file from the GitHub repo https://github.com/jedbrown/cmake-modules.git (BSD-2 License)
# Once done this will define
#
#  PETSC_FOUND        - system has PETSc
#  PETSC_INCLUDES     - the PETSc include directories
#  PETSC_LIBRARIES    - Link these to use PETSc
#  PETSC_COMPILER     - Compiler used by PETSc, helpful to find a compatible MPI
#  PETSC_DEFINITIONS  - Compiler switches for using PETSc
#  PETSC_MPIEXEC      - Executable for running MPI programs
#  PETSC_VERSION      - Version string (MAJOR.MINOR.SUBMINOR)
#
#  Usage:
#  find_package(PETSc COMPONENTS CXX)  - required if build --with-clanguage=C++ --with-c-support=0
#  find_package(PETSc COMPONENTS C)    - standard behavior of checking build using a C compiler
#  find_package(PETSc)                 - same as above
#
# Setting these changes the behavior of the search
#  PETSC_DIR - directory in which PETSc resides
#  PETSC_ARCH - build architecture
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#

set(PETSC_VALID_COMPONENTS
  C
  CXX)

if(NOT PETSc_FIND_COMPONENTS)
  set(PETSC_LANGUAGE_BINDINGS "C")
else()
  # Right now, this is designed for compatability with the --with-clanguage option, so
  # only allow one item in the components list.
  list(LENGTH ${PETSc_FIND_COMPONENTS} components_length)
  if(${components_length} GREATER 1)
    message(FATAL_ERROR "Only one component for PETSc is allowed to be specified")
  endif()
  # This is a stub for allowing multiple components should that time ever come. Perhaps
  # to also test Fortran bindings?
  foreach(component ${PETSc_FIND_COMPONENTS})
    list(FIND PETSC_VALID_COMPONENTS ${component} component_location)
    if(${component_location} EQUAL -1)
      message(FATAL_ERROR "\"${component}\" is not a valid PETSc component.")
    else()
      list(APPEND PETSC_LANGUAGE_BINDINGS ${component})
    endif()
  endforeach()
endif()

function (petsc_get_version)

  # Find the file petscversion.h
  find_path(VERSION_HEADER_PATH include/petscversion.h
            HINTS ${PETSC_DIR} ${PETSC_DIR}/${PETSC_ARCH}
            DOC "PETSc Directory")
  set(PETSC_VERSION_HEADER ${VERSION_HEADER_PATH}/include/petscversion.h)

  if (EXISTS "${PETSC_VERSION_HEADER}")
    file (STRINGS "${PETSC_VERSION_HEADER}" vstrings REGEX "#define PETSC_VERSION_(RELEASE|MAJOR|MINOR|SUBMINOR|PATCH) ")
    foreach (line ${vstrings})
      string (REGEX REPLACE " +" ";" fields ${line}) # break line into three fields (the first is always "#define")
      list (GET fields 1 var)
      list (GET fields 2 val)
      set (${var} ${val} PARENT_SCOPE)
      set (${var} ${val})         # Also in local scope so we have access below
    endforeach ()
    if (PETSC_VERSION_RELEASE)
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}p${PETSC_VERSION_PATCH}" PARENT_SCOPE)
    else ()
      # make dev version compare higher than any patch level of a released version
      set (PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.99" PARENT_SCOPE)
    endif ()
  else ()
    message (SEND_ERROR "PETSC_DIR can not be used, ${VERSION_HEADER_PATH}/include/petscversion.h does not exist")
  endif ()
endfunction ()

###############################################
# Begin determining PETSC_DIR and PETSC_ARCH
###############################################

# This will define PETSC_DIR and PETSC_ARCH if and only if they are defined in the corresponding environment variables
set(PETSC_DIR $ENV{PETSC_DIR})
set(PETSC_ARCH $ENV{PETSC_ARCH})

# If both variables are set, use them
if(DEFINED PETSC_DIR AND DEFINED PETSC_ARCH)
  if(EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include/petsc.h" OR EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h")
    set(PETSC_DIR "${PETSC_DIR}" CACHE FILEPATH "PETSc install directory")
    set(PETSC_ARCH "${PETSC_ARCH}" CACHE STRING "PETSc build architecture")
    message(STATUS "Determined PETSc install from the pair PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH}")
  else()
    message(STATUS "The pair PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} does not specify a valid PETSc installation")
    message(STATUS "Attempting to look in default locations instead...")
    unset(PETSC_DIR)
    unset(PETSC_ARCH)
  endif()
endif()

# Look in the following hard-coded PETSc locations first
if(NOT DEFINED PETSC_DIR)
  # If it's Linux...
  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(debian_dir /usr/lib/petscdir)
    # Hardcode default package locations for...
    # ... Ubuntu 14.04
    if(IS_DIRECTORY "/usr/lib/petscdir/3.4.2/linux-gnu-c-debug")
      set(PETSC_DIR "/usr/lib/petscdir/3.4.2" CACHE FILEPATH "PETSc install directory")
      set(PETSC_ARCH "linux-gnu-c-debug" CACHE STRING "PETSc build architecture")
      message(STATUS "Found candidate PETSc in default Ubuntu location: ${PETSC_DIR}/${PETSC_ARCH}")
    # ... Ubuntu 16.04
    elseif(IS_DIRECTORY "/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real")
      set(PETSC_DIR "/usr/lib/petscdir/3.6.2" CACHE FILEPATH "PETSc install directory")
      set(PETSC_ARCH "x86_64-linux-gnu-real" CACHE STRING "PETSc build architecture")
      message(STATUS "Found candidate PETSc in default Ubuntu location: ${PETSC_DIR}/${PETSC_ARCH}")
    # ... Ubuntu 16.10
    elseif(IS_DIRECTORY "/usr/lib/petscdir/3.7.3/x86_64-linux-gnu-real")
      set(PETSC_DIR "/usr/lib/petscdir/3.7.3" CACHE FILEPATH "PETSc install directory")
      set(PETSC_ARCH "x86_64-linux-gnu-real" CACHE STRING "PETSc build architecture")
      message(STATUS "Found candidate PETSc in default Ubuntu location: ${PETSC_DIR}/${PETSC_ARCH}")
    # ... Ubuntu 17.04
    elseif(IS_DIRECTORY "/usr/lib/petscdir/3.7.5/x86_64-linux-gnu-real")
      set(PETSC_DIR "/usr/lib/petscdir/3.7.5" CACHE FILEPATH "PETSc install directory")
      set(PETSC_ARCH "x86_64-linux-gnu-real" CACHE STRING "PETSc build architecture")
      message(STATUS "Found candidate PETSc in default Ubuntu location: ${PETSC_DIR}/${PETSC_ARCH}")
    # ... Ubuntu 17.10
    elseif(IS_DIRECTORY "/usr/lib/petscdir/3.7.6/x86_64-linux-gnu-real")
      set(PETSC_DIR "/usr/lib/petscdir/3.7.6" CACHE FILEPATH "PETSc install directory")
      set(PETSC_ARCH "x86_64-linux-gnu-real" CACHE STRING "PETSc build architecture")
      message(STATUS "Found candidate PETSc in default Ubuntu location: ${PETSC_DIR}/${PETSC_ARCH}")
    # ... Ubuntu 18.04
    elseif(IS_DIRECTORY "/usr/lib/petscdir/3.7.7/x86_64-linux-gnu-real")
      set(PETSC_DIR "/usr/lib/petscdir/3.7.7" CACHE FILEPATH "PETSc install directory")
      set(PETSC_ARCH "x86_64-linux-gnu-real" CACHE STRING "PETSc build architecture")
      message(STATUS "Found candidate PETSc in default Ubuntu location: ${PETSC_DIR}/${PETSC_ARCH}")
    # ... Ubuntu 18.10
    elseif(IS_DIRECTORY "/usr/lib/petscdir/petsc3.9/x86_64-linux-gnu-real")
      set(PETSC_DIR "/usr/lib/petscdir/petsc3.9" CACHE FILEPATH "PETSc install directory")
      set(PETSC_ARCH "x86_64-linux-gnu-real" CACHE STRING "PETSc build architecture")
      message(STATUS "Found candidate PETSc in default Ubuntu location: ${PETSC_DIR}/${PETSC_ARCH}")
    # ... Ubuntu 19.04
    elseif(IS_DIRECTORY "/usr/lib/petscdir/petsc3.10/x86_64-linux-gnu-real")
      set(PETSC_DIR "/usr/lib/petscdir/petsc3.10" CACHE FILEPATH "PETSc install directory")
      set(PETSC_ARCH "x86_64-linux-gnu-real" CACHE STRING "PETSc build architecture")
      message(STATUS "Found candidate PETSc in default Ubuntu location: ${PETSC_DIR}/${PETSC_ARCH}")
    # ... Anything else with dir /usr/lib/petscdir
    elseif(IS_DIRECTORY "/usr/lib/petscdir")
      # Find a petsc dir that ends in a.b.c
      file(GLOB potential_dirs "/usr/lib/petscdir/?.?.?")
      foreach(potential_dir ${potential_dirs})
        if(NOT DEFINED PETSC_DIR)
          set(PETSC_DIR ${potential_dir} CACHE FILEPATH "PETSc install directory")
        endif()
      endforeach()
      if(DEFINED PETSC_DIR)
        file(GLOB potential_archs RELATIVE ${PETSC_DIR} "${PETSC_DIR}/*")
        foreach(potential_arch ${potential_archs})
          if( (NOT DEFINED PETSC_ARCH) AND (${potential_arch} MATCHES "linux-gnu") )
            set(PETSC_ARCH ${potential_arch} CACHE STRING "PETSc build architecture")
            message(STATUS "Found candidate PETSc in default Ubuntu location: ${PETSC_DIR}/${PETSC_ARCH}")
          endif()
        endforeach()
      endif()
    # ... generic case.  Other Linux, perhaps
    else()
    endif()
  # If it's Apple OS X or macOS
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if(IS_DIRECTORY "/usr/local/Cellar/petsc")
      file(GLOB potential_dirs "/usr/local/Cellar/petsc/?.?.?")
      foreach(potential_dir ${potential_dirs})
        if(NOT DEFINED PETSC_DIR)
          if(IS_DIRECTORY "${potential_dir}/real")
            set(PETSC_DIR "${potential_dir}/real" CACHE FILEPATH "PETSc install directory")
            message(STATUS "Found candidate PETSc in default homebrew location: ${PETSC_DIR}")
          endif()
        endif()
      endforeach()
    endif()
  # If it's Windows
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if(IS_DIRECTORY "C:/Program\ Files/PETSc\ for\ Windows/PETSc/c-opt_icl_mkl")
      set(PETSC_DIR "C:/Program\ Files/PETSc\ for\ Windows/PETSc" CACHE FILEPATH "PETSc install directory")
      set(PETSC_ARCH "c-opt_icl_mkl" CACHE STRING "PETSc build architecture")
    endif()
  # In case we need a default case
  else()
    message(STATUS "CMAKE_SYSTEM_NAME is ${CMAKE_SYSTEM_NAME} is not Windows, Darwin or Linux")
  endif()
endif(NOT DEFINED PETSC_DIR)

# If it's still not found, use the previous infrastructure to have a stab at it
if(NOT DEFINED PETSC_DIR)
  set(homebrew_dir /usr/local/Cellar/petsc)
  file(GLOB homebrew_dirs ${homebrew_dir}/*/real)

  set(debian_dir /usr/lib/petscdir)
  file(GLOB ubuntu_dirs ${debian_dir}/*/*)
  file(GLOB debian_dirs ${debian_dir}/*)

  find_path (PETSC_DIR include/petsc.h
    HINTS
    $ENV{PETSC_DIR}
    # Homebrew
    ${homebrew_dirs}
    PATHS
    # Debian paths
    ${debian_dirs}
    # Ubuntu paths
    ${ubuntu_dirs}
    # MacPorts path
    /opt/local/lib/petsc
    # PETSc for Windows
    C:/Program\ Files/PETSc\ for\ Windows/PETSc
    $ENV{HOME}/petsc
    DOC "PETSc Directory")

endif(NOT DEFINED PETSC_DIR)

find_program (MAKE_EXECUTABLE NAMES make gmake)

if (PETSC_DIR AND NOT PETSC_ARCH)
    set (_petsc_arches
    $ENV{PETSC_ARCH}                   # If set, use environment variable first
    c-debug_icl_mkl c-opt_icl_mkl #PETSc for Windows
    x86_64-linux-gnu-real x86_64-linux-gnu-real-debug 
    linux-gnu-c-debug linux-gnu-c-opt  # Debian defaults
    linux-gnu-debug linux-gnu linux-gnu-opt linux-gnu-profile
    linux-intel-debug linux-intel-opt linux-intel-opt-mkl
    x86_64-unknown-linux-gnu i386-unknown-linux-gnu)
  set (petscconf "NOTFOUND" CACHE FILEPATH "Cleared" FORCE)
  foreach (arch ${_petsc_arches})
    if (NOT PETSC_ARCH)
      find_path (petscconf petscconf.h
        HINTS ${PETSC_DIR}
        PATH_SUFFIXES ${arch}/include bmake/${arch}
        NO_DEFAULT_PATH)
      if (petscconf)
        set (PETSC_ARCH "${arch}" CACHE STRING "PETSc build architecture")
      endif (petscconf)
    endif (NOT PETSC_ARCH)
  endforeach (arch)
  set (petscconf "NOTFOUND" CACHE INTERNAL "Scratch variable" FORCE)
endif (PETSC_DIR AND NOT PETSC_ARCH)


###############################################
# End determining PETSC_DIR and PETSC_ARCH
###############################################


set (petsc_slaves LIBRARIES_SYS LIBRARIES_VEC LIBRARIES_MAT LIBRARIES_DM LIBRARIES_KSP LIBRARIES_SNES LIBRARIES_TS
  INCLUDE_DIR INCLUDE_CONF)
include (FindPackageMultipass)
find_package_multipass (PETSc petsc_config_current
  STATES DIR ARCH
  DEPENDENTS INCLUDES LIBRARIES COMPILER MPIEXEC ${petsc_slaves})

if (PETSC_DIR)
    petsc_get_version()
	#search for config file in ${PETSC_DIR}/${PETSC_ARCH}. if exists use this
	if (EXISTS "${PETSC_DIR}/${PETSC_ARCH}/PETScConfig.cmake")
		find_package(PETSc NO_MODULE PATHS ${PETSC_DIR}/${PETSC_ARCH} NO_DEFAULT_PATH)
	endif()
	if (PETSC_FOUND)
		return()
	endif()
endif()

# Find the `conf/rules` and `conf/variables` files, that will help define the installation
find_path(petsc_conf_dir "conf/rules"
          HINTS
          ${PETSC_DIR} ${PETSC_DIR}/lib/petsc
          ${PETSC_DIR}/${PETSC_ARCH} ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc
          NO_DEFAULT_PATH)

if (EXISTS "${petsc_conf_dir}/conf/rules" AND EXISTS "${petsc_conf_dir}/conf/variables")
  set (petsc_conf_rules "${petsc_conf_dir}/conf/rules")
  set (petsc_conf_variables "${petsc_conf_dir}/conf/variables")
else ()
  message (SEND_ERROR "The pair PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} do not specify a valid PETSc installation")
endif ()

if (petsc_conf_rules AND petsc_conf_variables AND NOT petsc_config_current)
	
  # Put variables into environment since they are needed to get
  # configuration (petscvariables) in the PETSc makefile
  set (ENV{PETSC_DIR} "${PETSC_DIR}")
  set (ENV{PETSC_ARCH} "${PETSC_ARCH}")

  # A temporary makefile to probe the PETSc configuration
  set (petsc_config_makefile "${PROJECT_BINARY_DIR}/Makefile.petsc")
  file (WRITE "${petsc_config_makefile}"
"## This file was autogenerated by FindPETSc.cmake
# PETSC_DIR  = ${PETSC_DIR}
# PETSC_ARCH = ${PETSC_ARCH}
include ${petsc_conf_rules}
include ${petsc_conf_variables}
show :
\t-@echo -n \${\${VARIABLE}}
")


  macro (PETSC_GET_VARIABLE name var)
    set (${var} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
    execute_process (COMMAND ${MAKE_EXECUTABLE} --no-print-directory -f ${petsc_config_makefile} show VARIABLE=${name}
      OUTPUT_VARIABLE ${var}
      RESULT_VARIABLE petsc_return)
  endmacro (PETSC_GET_VARIABLE)
  

  petsc_get_variable (PETSC_LIB_DIR            petsc_lib_dir)
  petsc_get_variable (PETSC_EXTERNAL_LIB_BASIC petsc_libs_external)
  petsc_get_variable (PETSC_CCPPFLAGS          petsc_cpp_line)
  petsc_get_variable (PETSC_INCLUDE            petsc_include)
  petsc_get_variable (PCC                      petsc_cc)
  petsc_get_variable (PCC_FLAGS                petsc_cc_flags)
  petsc_get_variable (MPIEXEC                  petsc_mpiexec)
  # We are done with the temporary Makefile, calling PETSC_GET_VARIABLE after this point is invalid!
  file (REMOVE ${petsc_config_makefile})

  include (ResolveCompilerPaths)
  # Extract include paths and libraries from compile command line
  resolve_includes (petsc_includes_all "${petsc_cpp_line}")

  #on windows we need to make sure we're linking against the right
  #runtime library
  if (WIN32)
    if (petsc_cc_flags MATCHES "-MT*")
      set(using_md False)
      foreach(flag_var
          CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
          CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
          CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
          CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
        if(${flag_var} MATCHES "/MD")
          set(using_md True)
        endif(${flag_var} MATCHES "/MD")
      endforeach(flag_var)
      if(${using_md} MATCHES "True")
        message(WARNING "PETSc was built with /MT, but /MD is currently set.
 See http://www.cmake.org/Wiki/CMake_FAQ#How_can_I_build_my_MSVC_application_with_a_static_runtime.3F")
      endif(${using_md} MATCHES "True")
    endif (petsc_cc_flags MATCHES "-MT*")
  endif (WIN32)

  include (CorrectWindowsPaths)
  convert_cygwin_path(petsc_lib_dir)

  macro (PETSC_FIND_LIBRARY suffix name)
    set (PETSC_LIBRARY_${suffix} "NOTFOUND" CACHE INTERNAL "Cleared" FORCE) # Clear any stale value, if we got here, we need to find it again
    if (WIN32)
      set (libname lib${name}) #windows expects "libfoo", linux expects "foo"
    else (WIN32)
      set (libname ${name})
    endif (WIN32)
    find_library (PETSC_LIBRARY_${suffix} NAMES ${libname} ${libname}_real ${libname}_complex cray${libname}_cray_real  HINTS ${petsc_lib_dir} NO_DEFAULT_PATH)
    if (NOT PETSC_LIBRARY_${suffix})
        set (libname ${name})
		find_library (PETSC_LIBRARY_${suffix} NAMES ${libname} HINTS ${petsc_lib_dir} NO_DEFAULT_PATH)
    endif()
    set (PETSC_LIBRARIES_${suffix} "${PETSC_LIBRARY_${suffix}}")
    mark_as_advanced (PETSC_LIBRARY_${suffix})
  endmacro (PETSC_FIND_LIBRARY suffix name)

  # Look for petscvec first, if it doesn't exist, we must be using single-library
  petsc_find_library (VEC petscvec)
  if (PETSC_LIBRARY_VEC)
    petsc_find_library (SYS  "petscsys;petsc") # libpetscsys is called libpetsc prior to 3.1 (when single-library was introduced)
    petsc_find_library (MAT  petscmat)
    petsc_find_library (DM   petscdm)
    petsc_find_library (KSP  petscksp)
    petsc_find_library (SNES petscsnes)
    petsc_find_library (TS   petscts)
    macro (PETSC_JOIN libs deps)
      list (APPEND PETSC_LIBRARIES_${libs} ${PETSC_LIBRARIES_${deps}})
    endmacro (PETSC_JOIN libs deps)
    petsc_join (VEC  SYS)
    petsc_join (MAT  VEC)
    petsc_join (DM   MAT)
    petsc_join (KSP  DM)
    petsc_join (SNES KSP)
    petsc_join (TS   SNES)
    petsc_join (ALL  TS)
  else ()
    set (PETSC_LIBRARY_VEC "NOTFOUND" CACHE INTERNAL "Cleared" FORCE) # There is no libpetscvec
    petsc_find_library (SINGLE petsc)
    foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
      set (PETSC_LIBRARIES_${pkg} "${PETSC_LIBRARY_SINGLE}")
    endforeach ()
  endif ()
  if (PETSC_LIBRARY_TS)
    message (STATUS "Recognized PETSc install with separate libraries for each package")
  else ()
    message (STATUS "Recognized PETSc install with single library for all packages")
  endif ()

  include(Check${PETSC_LANGUAGE_BINDINGS}SourceRuns)

  find_path (PETSC_INCLUDE_DIR petscts.h HINTS "${PETSC_DIR}" PATH_SUFFIXES include NO_DEFAULT_PATH)
  find_path (PETSC_INCLUDE_CONF petscconf.h HINTS "${PETSC_DIR}" PATH_SUFFIXES "${PETSC_ARCH}/include" "bmake/${PETSC_ARCH}" NO_DEFAULT_PATH)
  mark_as_advanced (PETSC_INCLUDE_DIR PETSC_INCLUDE_CONF)


  if (petsc_libs_external AND NOT petsc_libraries_external)
      resolve_libraries (petsc_libraries_external "${petsc_libs_external}")
  endif()
  foreach (pkg SYS VEC MAT DM KSP SNES TS ALL)
      list (APPEND PETSC_LIBRARIES_${pkg}  ${petsc_libraries_external})
  endforeach (pkg)
  set (petsc_includes_needed ${petsc_includes_all})

  # We do an out-of-source build so __FILE__ will be an absolute path, hence __INSDIR__ is superfluous
  if (${PETSC_VERSION} VERSION_LESS 3.1)
    set (PETSC_DEFINITIONS "-D__SDIR__=\"\"" CACHE STRING "PETSc definitions" FORCE)
  else ()
    set (PETSC_DEFINITIONS "-D__INSDIR__=" CACHE STRING "PETSc definitions" FORCE)
  endif ()
  # Sometimes this can be used to assist FindMPI.cmake
  set (PETSC_MPIEXEC ${petsc_mpiexec} CACHE FILEPATH "Executable for running PETSc MPI programs" FORCE)
  set (PETSC_INCLUDES ${petsc_includes_needed} CACHE STRING "PETSc include path" FORCE)
  set (PETSC_LIBRARIES ${PETSC_LIBRARIES_ALL} CACHE STRING "PETSc libraries" FORCE)
  set (PETSC_COMPILER ${petsc_cc} CACHE FILEPATH "PETSc compiler" FORCE)
  # Note that we have forced values for all these choices.  If you
  # change these, you are telling the system to trust you that they
  # work.  It is likely that you will end up with a broken build.
  mark_as_advanced (PETSC_INCLUDES PETSC_LIBRARIES PETSC_COMPILER PETSC_DEFINITIONS PETSC_MPIEXEC)
endif (petsc_conf_rules AND petsc_conf_variables AND NOT petsc_config_current)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSc
  "PETSc could not be found.  Be sure to set PETSC_DIR and PETSC_ARCH."
  PETSC_INCLUDES PETSC_LIBRARIES)

