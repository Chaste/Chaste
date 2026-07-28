message(STATUS "Adding compiler flags...")

# Detect whether the standard library actually implements floating-point std::to_chars.
# This was added to libstdc++ in GCC 11; some officially-supported compilers (e.g. GCC 10)
# compile as C++17 but only implement the integer overloads of <charconv>, so checking
# __cplusplus/CMAKE_CXX_STANDARD alone is not sufficient here.
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
    #include <charconv>
    int main()
    {
        char buffer[32];
        double value = 3.14;
        std::to_chars_result result = std::to_chars(buffer, buffer + sizeof(buffer), value, std::chars_format::general, 6);
        return result.ec != std::errc();
    }
" Chaste_HAS_FLOAT_TO_CHARS)
if (Chaste_HAS_FLOAT_TO_CHARS)
    message(STATUS "Standard library supports floating-point std::to_chars")
    add_definitions(-DCHASTE_HAS_FLOAT_TO_CHARS)
else()
    message(STATUS "Standard library does not support floating-point std::to_chars; falling back to iostream formatting")
endif()

# default flags added to all compilers
set(default_flags "-Wall")
if (Chaste_ERROR_ON_WARNING)
    set(default_flags "${default_flags} -Werror")
endif()

set(default_exe_linker_flags "")
if (UNIX)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set(default_shared_link_flags "-Wl,-undefined,error")
    else()
        set(default_shared_link_flags "-Wl,--no-undefined")
    endif()
endif()

# Define a custom project-wide variable for linker flags
set(Chaste_SHARED_LINKER_FLAGS "${default_shared_link_flags}")

if (Chaste_COVERAGE)
    message(STATUS "adding --coverage to CXX flags for coverage checking")
    #--coverage seems to be the preferred flag
    #set(default_flags "${default_flags} -fprofile-arcs -ftest-coverage")
    #set(default_exe_linker_flags "${default_exe_linker_flags} -fprofile-arcs -ftest-coverage")
    set(default_flags "${default_flags} --coverage")
    set(default_shared_link_flags "${default_shared_link_flags} --coverage")
    set(default_exe_linker_flags "${default_exe_linker_flags} --coverage")
endif()

if (Chaste_PROFILE_GPROF)
    message(STATUS "adding -Og and -pg to CXX flags for Gprof profiling")
    set(default_flags "${default_flags} -Og -pg -Wno-array-bounds")
    set(default_shared_link_flags "${default_shared_link_flags} -pg")
    set(default_exe_linker_flags "${default_exe_linker_flags} -pg")
endif()


if (Chaste_PROFILE_GPERFTOOLS)
    message(STATUS "adding -O3 to CXX flags for Gperftools profiling")
    set(default_flags "${default_flags} -O3")
endif()

# Allow easier checks in source files for specific compilers. See here for possible
# values: https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_ID.html
add_compile_definitions(Chaste_COMPILER_IS_${CMAKE_CXX_COMPILER_ID})

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray")
    message(STATUS "\t...for Cray compiler, version ${CMAKE_CXX_COMPILER_VERSION}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${default_flags} -Wnon-virtual-dtor -Woverloaded-virtual -Wextra -Wno-unused-parameter -Wvla")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${default_flags} -Wextra -Wno-unused-parameter -Wvla")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    message(STATUS "\t...for GNU compiler, version ${CMAKE_CXX_COMPILER_VERSION}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${default_flags} -Wnon-virtual-dtor -Woverloaded-virtual -Wextra -Wno-attributes -Wno-unused-parameter -Wvla")
    if (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7))
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wimplicit-fallthrough=2")  # See https://developers.redhat.com/blog/2017/03/10/wimplicit-fallthrough-in-gcc-7/
    endif (NOT (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 7))

    # The following flags allow us to ignore (in a judicious manner) some spurious GCC warnings about potentially
    # uninitialised c_vectors in Release mode. See https://github.com/Chaste/Chaste/issues/231 for full details.
    if (CMAKE_BUILD_TYPE MATCHES "Release|RelWithDebInfo|MinSizeRel"
            AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9.4
            AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16)

        if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11.1)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized -Wno-array-bounds -Wno-stringop-overflow")
        else ()
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-maybe-uninitialized -Wno-array-bounds -Wno-stringop-overflow -Wno-stringop-overread")
        endif ()
    endif ()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${default_flags}  -Wextra -Wno-unused-parameter -Wvla")
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "AppleClang" OR ${CMAKE_CXX_COMPILER_ID} STREQUAL "IntelLLVM")
    message(STATUS "\t... for ${CMAKE_CXX_COMPILER_ID} compiler, version ${CMAKE_CXX_COMPILER_VERSION}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${default_flags} -Wnon-virtual-dtor -Woverloaded-virtual -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-undefined-var-template -Wno-unknown-warning-option -Wno-tautological-constant-compare -ftemplate-depth-512")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${default_flags}  -Wextra -Wno-unused-parameter -Wno-unused-variable -ftemplate-depth-512")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${default_exe_linker_flags}")
    if (${CMAKE_CXX_COMPILER_ID} STREQUAL "IntelLLVM")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Rno-debug-disables-optimization -Wno-unused-but-set-variable")
        option(Chaste_USE_INTEL_PRECISE_FP_MODEL "Use Intel -fp-model precise" ON)
        if (Chaste_USE_INTEL_PRECISE_FP_MODEL)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model precise")
            message(STATUS
                    "Using Intel floating point model \"precise\". "
                    "To use the default fast floating point model, set -DChaste_USE_INTEL_PRECISE_FP_MODEL=OFF. "
                    "Note that this may reduce the precision of calculations and change other behaviour.")
        endif ()
    endif ()
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel")
    message(STATUS "\t... for Intel compiler, version ${CMAKE_CXX_COMPILER_VERSION}")
    set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ")
    set(Intel_flags
        # Not available on 10.0
        -Wnon-virtual-dtor -Woverloaded-virtual -Wno-unused-parameter

        -wr2304 #2304: non-explicit constructor with single argument

        # Switch these ones on for compatibility
        -wr271 #271: trailing comma is nonstandard

        #Following doesn't seem to play
        -wr810 #810: conversion from "double" to "unsigned int" may lose significant bits

        # This is where the statement is unreachable in a particular instatiation of the template.  e.g. "if (SPACE_DIM<3){return;}" will complain that the SPACE_DIM=3 specific code is unreachable.
       -wr111 #111: statement is unreachable (DUE TO INSTANTIATED TEMPLATES)

        # This is where the statement is unreachable in a particular instatiation of the template.  e.g. "if (ELEMENT_DIM<SPACE_DIM){return;}" will complain that the ELEMENT_DIM == SPACE_DIM dynamic initialization is unreachable.
        -wr185 #185: dynamic initialization in unreachable code (DUE TO INSTANTIATED TEMPLATES)

        # This happens when a switch is based on an unsigned template parameter
        -wr280 #280: selector expression is constant

        # This is seen when used templates to access the is_abstract base class definition
        -wr304 #304: access control not specified ("public" by default)

        # This is when we pass an explict string to a std::string reference: e.g. FileFinder save_bidomain_dir("some_directory", RelativeTo::ChasteSourceRoot);
        -wr383 #383: value copied to temporary, reference to temporary used

        # Noncopyable doesn't have a virtual destructor.  The derived class should not have access to it either
        -wr444 #444: destructor for base class "boost::noncopyable_::noncopyable" ... is not virtual

        # Most commonly seen in archiving where the "version" variable is often redundant
        -wr869 #869: parameter "..." was never referenced

        # Triggered by macros such as TS_ASSERT_EQUALS(a,b)
        -wr981 #981: operands are evaluated in unspecified order

        # We do this when we need to define templated functions in the header file
        -wr1418 #1418: external function definition with no prior declaration

        # There are times when we want a local helper function (RecursiveCopy in FileFinder) or when we need to refer to KSPConvergedReasons
        -wr1419 #1419: external declaration in primary source file

        # This one is potentially useful for telling us where we might want to use CompareDoubles::WithinRelativeTolerance, but in our core code (TimeStepper) the tests should ensure we aren't doing anything silly
        -wr1572 #1572: floating-point equality and inequality comparisons are unreliable

        #2289: proper signature for "auto_ptr" is "Type(const Type&)"

        #2026: Effective C++ Item 14 Make sure base classes have virtual destructors

         #2305: declaration of 'explicit' constructor without a single argument is redundant
        )
    string (REPLACE ";" " " Intel_flags_str "${Intel_flags}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${default_flags} ${Intel_flags_str}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${default_flags} ${Intel_flags_str}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${default_exe_linker_flags}")
else()
    message(WARNING "Unknown CXX compiler type ${CMAKE_CXX_COMPILER_ID}")
endif()

set(Chaste_SHARED_LINKER_FLAGS "${Chaste_SHARED_LINKER_FLAGS}" CACHE STRING "Project-wide shared linker flags" FORCE)
