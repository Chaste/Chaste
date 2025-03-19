# - Try to find petsc4py
# Once done this will define
#
#  PETSC4PY_FOUND    - system has petsc4py
#  PETSC4PY_INCLUDES - the petsc4py include directories
#  PETSC4PY_VERSION  - version string (MAJOR.MINOR.SUBMINOR)
#
# Usage: find_package(PETSc4py)

execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import petsc4py; print(petsc4py.get_include(), end='')"
    OUTPUT_VARIABLE PETSC4PY_INCLUDES
)

if(PETSC4PY_INCLUDES)
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -c "import petsc4py; print(petsc4py.__version__, end='')"
        OUTPUT_VARIABLE PETSC4PY_VERSION
    )
endif()

mark_as_advanced(PETSC4PY_INCLUDES, PETSC4PY_VERSION)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    PETSc4py
    REQUIRED_VARS PETSC4PY_INCLUDES
    VERSION_VAR PETSC4PY_VERSION
    FAIL_MESSAGE "PETSc4py could not be found."
)
