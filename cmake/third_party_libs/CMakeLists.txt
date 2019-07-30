# Copyright (c) 2005-2019, University of Oxford.
# All rights reserved.
# 
# University of Oxford means the Chancellor, Masters and Scholars of the
# University of Oxford, having an administrative office at Wellington
# Square, Oxford OX1 2JD, UK.
# 
# This file is part of Chaste.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * Neither the name of the University of Oxford nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#This CMake build file automatically builds Chaste's third-party libraries.
cmake_minimum_required (VERSION 2.8.10)

include(${CMAKE_CURRENT_SOURCE_DIR}/../overrides.cmake)

project(ChasteThirdPartyLibraries)
include(ExternalProject)

include(${CMAKE_CURRENT_SOURCE_DIR}/../definitions.cmake)


if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Default install directory. Change as required." FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(MS_HPC_PACK_DIR "C:/MS_HPC_PACK_2012" CACHE PATH
    "The root directory of Microsoft's HPC Pack. Make sure that the path has no space in it (create a symbolic link if necessary to ensure this).")
set(OUTPUT_LIB_DIR "${MS_HPC_PACK_DIR}/Lib/amd64")
set(OUTPUT_INCLUDE_DIR "${MS_HPC_PACK_DIR}/Inc")

option(BUILD_BOOST "Do you want to build Boost?" ON)
option(BUILD_PETSC "Do you want to build PETSc?" OFF)

set(DOWNLOAD_DIR "${CMAKE_BINARY_DIR}/downloads")
string(REGEX REPLACE "([a-zA-Z]):" "/cygdrive/\\1" CYG_DOWNLOAD_DIR "${DOWNLOAD_DIR}")

set(patch_message "#Patched by Chaste\n") # used to prevent re-patching what we've already patched

#Automatic install of Cygwin and the required software, if the user requests it.
option(AUTO_INSTALL_CYGWIN "Install Cygwin and required packages automatically" OFF)
if(AUTO_INSTALL_CYGWIN)
    if(NOT EXISTS "${DOWNLOAD_DIR}/cygwin_installer")
        file(MAKE_DIRECTORY "${DOWNLOAD_DIR}/cygwin_installer")
    endif()

    if(NOT CYGWIN_ROOT_DIR)
        set(CYGWIN_ROOT_DIR "${DOWNLOAD_DIR}/cygwin_installer/cygwin" CACHE PATH "Cygwin root directory.")
        message(FATAL_ERROR "Please set a root directory for Cygwin. One has been generated in the variable CYGWIN_ROOT_DIR above. Run 'Configure' again to accept it, or you may change it.")
    endif(NOT CYGWIN_ROOT_DIR)

    message(STATUS "Attempting to install Cygwin, and software needed for library build.
            This may take a while ...")

    #Download Cygwin (an improvement would be to verify the binary signature before running it!)
    if(NOT EXISTS "${DOWNLOAD_DIR}/cygwin_installer/setup.exe")
        file(DOWNLOAD "http://cygwin.com/setup.exe" "${DOWNLOAD_DIR}/cygwin_installer/setup.exe"
                SHOW_PROGRESS
                STATUS status
                  LOG log)
    endif()

    message("Download log: ${log}")
    message("Download status: ${status}")

   if(CYGWIN_ROOT_DIR AND NOT EXISTS "${CYGWIN_ROOT_DIR}")
    set(C_COMMAND ${DOWNLOAD_DIR}/cygwin_installer/setup.exe --root ${CYGWIN_ROOT_DIR}
        --site http://ftp.heanet.ie/mirrors/cygwin --no-shortcuts --quiet-mode --disable-buggy-antivirus
        --packages mingw64-i686-gcc-core,gendef,python,cmake,make,openssh)
    #gcc4-core and zlib seem to be optional Cygwin packages. I added them while hunting down the MS MPI PETSc integration failure
    #In the end, they don't seem to contribute to the solution. Keeping a record here for reference.
    #Openssh is needed for the automated testing platform for Chaste. It is not be needed if the intention is to just build PETSc.

    message(STATUS "Attempting to run downloaded Cygwin installer.
    This may take a while depending on your connection and machine speed ...
    Not to worry, the installation should not happen often :-).
    A timeout of 20 minutes has been set. Please WAIT ...")

      execute_process(
        COMMAND ${C_COMMAND}
        WORKING_DIRECTORY "${DOWNLOAD_DIR}/cygwin_installer"
        OUTPUT_VARIABLE cyg_log_out
        ERROR_VARIABLE cyg_log_err
        RESULT_VARIABLE cyg_result
        TIMEOUT 1200 #20 minutes max!
    )
    message("Outputs and result of command ${DOWNLOAD_DIR}/cygwin_installer/setup.exe ${C_COMMAND}:\n
Standard error
==============
${cyg_log_err}
Standard output
===============
${cyg_log_out}\n
Result = ${cyg_result}\n
=======End of Outputs for\n${DOWNLOAD_DIR}/cygwin_installer/setup.exe ${C_COMMAND}\n
====================================================")
    else()
        message(STATUS "It seems Cygwin is already installed. If not, delete the directory ${CYGWIN_ROOT_DIR}")
    endif(CYGWIN_ROOT_DIR AND NOT EXISTS "${CYGWIN_ROOT_DIR}")
endif(AUTO_INSTALL_CYGWIN)


#============================= Locate Cygwin (Needed by PETSc)

if(CYGWIN_ROOT_DIR)
    set(CYG_HINTS "${CYGWIN_ROOT_DIR}/bin" "C:/Cygwin/bin" "D:/Cygwin/bin" "E:/Cygwin/bin" "F:/Cygwin/bin")
else()
    set(CYG_HINTS "C:/Cygwin/bin" "D:/Cygwin/bin" "E:/Cygwin/bin" "F:/Cygwin/bin")
endif()

#Find Cygwin bash
find_program(CYGWINBASH bash HINTS ${CYG_HINTS})

if(CYGWINBASH STREQUAL "CYGWINBASH-NOTFOUND")
    option(AUTO_INSTALL_CYGWIN OFF "Should I attempt to automatically install Cygwin and the required software?")
    message(FATAL_ERROR "Cygwin is required to build PETSc. I can auto-install it for you if you enable the option AUTO_INSTALL_CYGWIN above.")
else()
    unset(AUTO_INSTALL_CYGWIN CACHE)
endif()


string(REPLACE "/" ";" CYG_BIN "${CYGWINBASH}")
list(LENGTH CYG_BIN length)
math(EXPR index "${length} - 1")
list(REMOVE_AT CYG_BIN ${index})
string(REPLACE ";" "/" CYG_BIN "${CYG_BIN}")
string(REGEX REPLACE "([a-zA-Z]):" "/cygdrive/\\1" CYG_BIN "${CYG_BIN}")
string(REGEX REPLACE "([a-zA-Z0-9_/:.-]+cygwin.*)/bin/bash.exe$" "\\1" CYG_PATH "${CYGWINBASH}")

string(REGEX REPLACE "([a-zA-Z]):" "/cygdrive/\\1" MS_HPC_PACK_DIR_CYG "${MS_HPC_PACK_DIR}")

#fix an issue with python thinking it is not properly configured
if(NOT EXISTS "${CYG_PATH}/include/python2.7")
    file(MAKE_DIRECTORY "${CYG_PATH}/include/python2.7")
endif()

if(NOT EXISTS "${CYG_PATH}/include/python2.7/pyconfig.h")
string(REPLACE "/" "\\" SYMLINK "${CYG_PATH}/include/python2.7/pyconfig.h")
string(REPLACE "/" "\\" SYMLINK_TARGET "${CYG_PATH}/usr/include/python2.7/pyconfig.h")

execute_process(COMMAND cmd /c mklink ${SYMLINK} ${SYMLINK_TARGET})
endif()


include(ChasteThirdPartyLibs.cmake)

macro(download_from_url url package)
    #remove leading and trailing whitespaces and newlines from url
    string(STRIP "${url}" url)
    string(REPLACE "/" ";" urllist "${url}")
    list(LENGTH urllist length)
    math(EXPR index "${length} - 1")
    list(GET urllist ${index} filename)

    if (NOT EXISTS "${DOWNLOAD_DIR}/${package}/${filename}")
        message(STATUS "Downloading ${filename}")
        file(DOWNLOAD "${url}" "${DOWNLOAD_DIR}/${package}/${filename}"
            SHOW_PROGRESS
            STATUS status
              LOG log)
    endif()
    string(REPLACE ".tar" "" bname "${filename}")
    string(REPLACE ".gz" "" bname "${bname}")
    string(REPLACE ".zip" "" bname "${bname}")
    string(REPLACE ".tgz" "" bname "${bname}")

    if (NOT EXISTS "${DOWNLOAD_DIR}/${package}/${bname}")
        # Unzip
        message(STATUS "Unzipping ${filename}")
        execute_process(
        COMMAND ${CMAKE_COMMAND} -E tar xzf "${DOWNLOAD_DIR}/${package}/${filename}"
        WORKING_DIRECTORY "${DOWNLOAD_DIR}/${package}")
    endif()
endmacro(download_from_url url package)



if(BUILD_BOOST)
#Build the Boosts
foreach(url ${BOOST_URLS})

    download_from_url(${url} "boost")

    set(boost_basicname "${bname}")

    #Boost's library naming convention
    string(REGEX REPLACE ".*_([0-9])+" "\\1" X "${boost_basicname}")
    if(X EQUAL 0)
      string(REGEX REPLACE "boost_(.*)_[0-9]+" "boost-\\1" SUFFIX "${boost_basicname}")
    else()
      string(REGEX REPLACE "boost_(.*)_([0-9]+)" "boost-\\1_\\2" SUFFIX "${boost_basicname}")
    endif()
    set(OUTPUT_LIB_DIR "${OUTPUT_LIB_DIR}" "${CMAKE_INSTALL_PREFIX}/boost_${boost_basicname}/lib")
    set(OUTPUT_INCLUDE_DIR "${OUTPUT_INCLUDE_DIR}" "${CMAKE_INSTALL_PREFIX}/boost_${boost_basicname}/include/${SUFFIX}")



    #Check whether Boost has been previously configured.
    find_program(BJAM_${boost_basicname} bjam HINTS "${DOWNLOAD_DIR}/boost/${boost_basicname}")
    if(BJAM_${boost_basicname} STREQUAL "BJAM_${boost_basicname}-NOTFOUND")
        #Configure Boost
        message(STATUS "Configuring Boost: ${boost_basicname}")
        set(C_COMMAND bootstrap.bat)
        execute_process(
                COMMAND ${C_COMMAND}
                WORKING_DIRECTORY "${DOWNLOAD_DIR}/boost/${boost_basicname}"
                OUTPUT_VARIABLE boost_log_out
                ERROR_VARIABLE boost_log_err
                RESULT_VARIABLE boost_result
            )
        message("Outputs and result of command ${C_COMMAND}:\n
Standard error\n
============\n
${boost_log_err}
Standard output\n
=============\n
${boost_log_out}\n
Result = ${boost_result}\n
=======End of Outputs for\n${C_COMMAND}\n
====================================================")

    else()
        #unset(BJAM_${boost_basicname})
        message(STATUS "Boost ${boost_basicname} is already configured. Reconfigure manually, or delete ${DOWNLOAD_DIR}/boost/${boost_basicname} for automatic configuration.")
    endif()

    #Patch Boost to enable MPI ...
    #Patch $BOOST_SRC/tools/build/v2/user-config.jam
        file(READ "${DOWNLOAD_DIR}/boost/${boost_basicname}/tools/build/v2/user-config.jam" boost_userconf)
          #check whether we have patched this already
        string(FIND "${boost_userconf}" "${patch_message}" patched)
        if(patched EQUAL -1)#not patched yet
            message(STATUS "Patching Boost ${boost_basicname} to enable MPI")
              file(WRITE "${DOWNLOAD_DIR}/boost/${boost_basicname}/tools/build/v2/user-config.jam"
                      "${patch_message}\n\n# Enable MPI\n   using mpi ;\n\n${boost_userconf}")
          endif()
      #Help boost find MS MPI
          file(READ "${DOWNLOAD_DIR}/boost/${boost_basicname}/tools/build/v2/tools/mpi.jam" boost_mpi)
          #check whether we have patched this already
        string(FIND "${boost_mpi}" "${patch_message}" patched)
        if(patched EQUAL -1)#not patched yet

            string(REGEX REPLACE "(local[ ]*cluster_pack_path_native[ ]*=[ ]*)[^;]+;" "\\1 \"${MS_HPC_PACK_DIR}\" ;"
                boost_mpi "${boost_mpi}")
            string(REGEX REPLACE "(if[ ]*[[][ ]*GLOB[ ]*[$][(][ ]*cluster_pack_path_native[ ]*[)])[^]]+[]]" "\\1\\\\\\\\Inc : mpi.h ]"
                boost_mpi "${boost_mpi}")
            string(REGEX REPLACE "(options[ ]*=[ ]*<include>[ ]*[$][ ]*[(][ ]*cluster_pack_path[ ]*[)])/Include" "\\1/Inc"
                boost_mpi "${boost_mpi}")
              file(WRITE "${DOWNLOAD_DIR}/boost/${boost_basicname}/tools/build/v2/tools/mpi.jam"
                      "${patch_message}${boost_mpi}")
          endif()

          #Build Boost (Boost does not like too deep build paths)
        message(STATUS "Building Boost: ${boost_basicname}. Note this can take very long! Please be patient.")
        set(C_COMMAND b2 --build-type=complete msvc stage address-model=64 --build-dir=${CMAKE_BINARY_DIR} --without-python --without-regex --stagedir=${CMAKE_INSTALL_PREFIX}/boost_${boost_basicname})
        execute_process(
                COMMAND ${C_COMMAND}
                WORKING_DIRECTORY "${DOWNLOAD_DIR}/boost/${boost_basicname}"
                OUTPUT_VARIABLE boost_log_out
                ERROR_VARIABLE boost_log_err
                RESULT_VARIABLE boost_result
            )
        message("Outputs and result of command ${C_COMMAND}:\n
Standard error\n
============\n
${boost_log_err}
Standard output\n
=============\n
${boost_log_out}\n
Result = ${boost_result}\n
=======End of Outputs for\n"${C_COMMAND}"\n
====================================================")

        #Install Boost (if build succeeds)
        if(NOT boost_result EQUAL 1)
            message(STATUS "Installing Boost: ${boost_basicname} at ${CMAKE_INSTALL_PREFIX}/boost_${boost_basicname}. This also takes very long.")
            set(C_COMMAND b2 install --prefix=${CMAKE_INSTALL_PREFIX}/boost_${boost_basicname})
            execute_process(
                    COMMAND ${C_COMMAND}
                    WORKING_DIRECTORY "${DOWNLOAD_DIR}/boost/${boost_basicname}"
                    OUTPUT_VARIABLE boost_log_out
                    ERROR_VARIABLE boost_log_err
                    RESULT_VARIABLE boost_result
                )
            message("Outputs and result of command ${C_COMMAND}:\n
Standard error\n
============\n
${boost_log_err}
Standard output\n
=============\n
${boost_log_out}\n
Result = ${boost_result}\n
=======End of Outputs for\n${C_COMMAND}\n
====================================================")
    endif(NOT boost_result EQUAL 1)

endforeach(url ${BOOST_URLS})
endif(BUILD_BOOST)


if (BUILD_PETSC)
foreach(url ${PETSC_URLS})
    #remove leading and trailing whitespaces and newlines from url
    string(STRIP "${url}" url)
    string(REPLACE "/" ";" urllist "${url}")
    list(LENGTH urllist length)
    math(EXPR index "${length} - 1")
    list(GET urllist ${index} filename)
    string(REPLACE ".tar" "" basicname "${filename}")
    string(REPLACE ".gz" "" basicname "${basicname}")
    string(REPLACE ".zip" "" basicname "${basicname}")
    string(REPLACE ".tgz" "" basicname "${basicname}")

    #Download PETSc
    download_from_url(${url} "petsc")

    #Download corresponding parmetis
    file(STRINGS "${DOWNLOAD_DIR}/petsc/${basicname}/config/PETSc/packages/parmetis.py" parmetispy)
    set(parmetis_url "")
    foreach(s ${parmetispy})
        string(REGEX REPLACE "self[.]download[\t ]*=[\t ]*[[](.*)[]]" "" match ${s})
        string(STRIP "${match}" match)
    if(match STREQUAL "")
        set(parmetis_url "${CMAKE_MATCH_1}")
        string(REPLACE "'" "" parmetis_url "${parmetis_url}" ) #remove single quotes
    endif()
    endforeach()

    if(parmetis_url STREQUAL "")
        message(FATAL_ERROR "Could not locate PARMETIS URL in ${DOWNLOAD_DIR}/petsc/${basicname}/config/PETSc/packages/parmetis.py")
      else()
          set(url "${parmetis_url}")
          download_from_url(${url} "parmetis")
          set(parmetis_basicname "${bname}")
      endif()


    file(STRINGS "${DOWNLOAD_DIR}/petsc/${basicname}/config/PETSc/packages/metis.py" metispy)
    set(metis_url "")
    foreach(s ${metispy})
        string(REGEX REPLACE "self[.]download[\t ]*=[\t ]*[[](.*)[]]" "" match ${s})
        string(STRIP "${match}" match)
    if(match STREQUAL "")
        set(metis_url "${CMAKE_MATCH_1}")
        string(REPLACE "'" "" metis_url "${metis_url}" ) #remove single quotes
    endif()
    endforeach()

    if(metis_url STREQUAL "")
        message(FATAL_ERROR "Could not locate METIS URL in ${DOWNLOAD_DIR}/petsc/${basicname}/config/PETSc/packages/metis.py")
    else()
          set(url "${metis_url}")
          download_from_url(${url} "metis")
          set(metis_basicname "${bname}")
      endif()

          #patch the PARMETIS CMakeLists.txt
          file(READ "${DOWNLOAD_DIR}/parmetis/${parmetis_basicname}/CMakeLists.txt" parcmake)

          #check whether we have patched this already
        string(FIND "${parcmake}" "${patch_message}" patched)
        if(patched EQUAL -1)#not patched yet
          string(REPLACE "project(ParMETIS)" "project(ParMETIS)

#Use gk_regex.h instead of regex.h and a bunch of other flags for GKlib and metis
add_definitions(-DUSE_GKREGEX -DWIN32 -DMSC -D_CRT_SECURE_NO_DEPRECATE)
add_definitions(-MTd)#static debug build
add_definitions(-Z7)#embed debugging info in library as opposed to using an external .pdb database

 find_package(MPI)
 if(NOT MPI_FOUND)
   message(FATAL_ERROR \"MPI is not found\")
 endif()
 set(CMAKE_C_FLAGS \"\${CMAKE_C_FLAGS} \${MPI_COMPILE_FLAGS}\")
" parcmake "${parcmake}")
          string(REGEX REPLACE "set[ ]*[(][ ]*GKLIB_PATH[^)]+[)][^\n]*"
              "set(GKLIB_PATH \"${DOWNLOAD_DIR}/metis/${metis_basicname}/GKlib\" CACHE PATH \"path to GKlib\")"
              parcmake "${parcmake}")
          string(REGEX REPLACE "set[ ]*[(][ ]*METIS_PATH[^)]+[)][^\n]*"
              "set(METIS_PATH \"${DOWNLOAD_DIR}/metis/${metis_basicname}\" CACHE PATH \"path to METIS\")

#make sure that METIS headers and generated GKLibs headers are found
include_directories(headers)
include_directories(\${CMAKE_INSTALL_PREFIX}/include)
include_directories(\${METIS_INSTALL_DIR}/include)
include_directories(\${METIS_BINARY_DIR}/include)
              " #End of replacement string
              parcmake "${parcmake}")
        string(REGEX REPLACE "link_directories[ ]*[(][ ]*[$]{METIS_PATH}/lib[ ]*[)][^\n]*" "" parcmake "${parcmake}")
        string(REGEX REPLACE
            "link_directories[ ]*[(][ ]*[$]{CMAKE_INSTALL_PREFIX}/lib[ ]*[)][^\n]*" "" parcmake "${parcmake}")

        file(WRITE "${DOWNLOAD_DIR}/parmetis/${parmetis_basicname}/CMakeLists.txt" "${patch_message}${parcmake}")
endif()#check whether patched

        #Patch METIS CMakeLists.txt
        file(READ "${DOWNLOAD_DIR}/metis/${metis_basicname}/CMakeLists.txt" metiscmake)
          #check whether we have patched this already
        string(FIND "${metiscmake}" "${patch_message}" patched)
        if(patched EQUAL -1)#not patched yet
              string(REPLACE "project(METIS)" "project(METIS)
set(METIS_INSTALL TRUE CACHE BOOL \"Enable the independent install of METIS\")
add_definitions(-MTd)#static debug build
add_definitions(-Z7)#embed debugging info in library as opposed to using an external .pdb database
include_directories(\"\${CMAKE_BINARY_DIR}/include\")"
                    metiscmake "${metiscmake}")
              string(REGEX REPLACE "if[ ]*[(][ ]*MSVC[ ]*[)].*endif[ ]*[(][^)]*[)]" "include_directories(\${GKLIB_PATH}/include)"
                  metiscmake "${metiscmake}")

        file(WRITE "${DOWNLOAD_DIR}/metis/${metis_basicname}/CMakeLists.txt"
            "${patch_message}${metiscmake}\nadd_subdirectory(\"GKlib\")")
          endif()#check whether patched


          #Build METIS as an enternal project dependency
        externalproject_add(METIS_${metis_basicname}
            SOURCE_DIR ${DOWNLOAD_DIR}/metis/${metis_basicname}
            CMAKE_GENERATOR ${CMAKE_GENERATOR}
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/metis_${metis_basicname}
                -DMETIS_INSTALL:BOOL=TRUE
                -DMETIS_USE_LONGINDEX:BOOL=TRUE
                -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
                -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
                -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
                -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
                -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
                -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
                -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
                -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
            BINARY_DIR ${CMAKE_BINARY_DIR}/metis_${metis_basicname}
            INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/metis_${metis_basicname}
            )
        set(OUTPUT_LIB_DIR "${OUTPUT_LIB_DIR}" "${CMAKE_INSTALL_PREFIX}/metis_${metis_basicname}/lib")
        set(OUTPUT_INCLUDE_DIR "${OUTPUT_INCLUDE_DIR}" "${CMAKE_INSTALL_PREFIX}/metis_${metis_basicname}/include")

        #Build ParMETIS as an external project
        externalproject_add(ParMETIS_${parmetis_basicname}
            SOURCE_DIR ${DOWNLOAD_DIR}/parmetis/${parmetis_basicname}
            CMAKE_GENERATOR ${CMAKE_GENERATOR}
            BINARY_DIR ${CMAKE_BINARY_DIR}/parmetis_${parmetis_basicname}
            INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/parmetis_${parmetis_basicname}
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/parmetis_${parmetis_basicname}
                -DMETIS_INSTALL_DIR=${CMAKE_INSTALL_PREFIX}/metis_${metis_basicname}
                -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
                -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
                -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
                -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
                -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
                -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
                -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
                -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
            DEPENDS METIS_${metis_basicname}
        )
        set(OUTPUT_LIB_DIR "${OUTPUT_LIB_DIR}" "${CMAKE_INSTALL_PREFIX}/parmetis_${parmetis_basicname}/lib" )
        set(OUTPUT_INCLUDE_DIR "${OUTPUT_INCLUDE_DIR}" "${CMAKE_INSTALL_PREFIX}/parmetis_${parmetis_basicname}/include")

    #Don't forget to build a native Windows version of F2CBlas and Lapack!
    file(STRINGS "${DOWNLOAD_DIR}/petsc/${basicname}/config/BuildSystem/config/packages/f2cblaslapack.py" f2cpy)
    set(f2c_url "")
    foreach(s ${f2cpy})
        string(REGEX REPLACE "self[.]download[\t ]*=[\t ]*[[](.*)[]]" "" match ${s})
        string(STRIP "${match}" match)
    if(match STREQUAL "")
        set(f2c_url "${CMAKE_MATCH_1}")
        string(REPLACE "'" "" f2c_url "${f2c_url}" ) #remove single quotes
    endif()
    endforeach()

    if(f2c_url STREQUAL "")
        message(FATAL_ERROR "Could not locate F2CBlasLapack URL in ${DOWNLOAD_DIR}/petsc/${basicname}/config/PETSc/packages/f2cblaslapack.py")
    else()
          set(url "${f2c_url}")
          download_from_url(${url} "f2cblaslapack")
          set(f2c_basicname "${bname}")
          # Generate CMake build scripts for the libraries
          # We get some help on what sources to add from the makefiles in the directories
          file(READ "${DOWNLOAD_DIR}/f2cblaslapack/${f2c_basicname}/blas/makefile" blasmake)
          string(REGEX MATCHALL "(AUXO|SINGLEO|DOUBLEO)[ ]*=[ ]*[^\n]+\n" blasmakematch "${blasmake}")
          string(REGEX MATCHALL "[a-zA-Z0-9_]+[.]o" blasmakematch "${blasmakematch}")
          string(REGEX REPLACE "([a-zA-Z0-9_]+)[.]o" "blas/\\1.c" blasmakematch "${blasmakematch}")

          file(READ "${DOWNLOAD_DIR}/f2cblaslapack/${f2c_basicname}/lapack/makefile" lapackmake)
          string(REGEX MATCHALL "(AUXO|SINGLEO|DOUBLEO)[ ]*=[ ]*[^\n]+\n" lapackmakematch "${lapackmake}")
          string(REGEX MATCHALL "[a-zA-Z0-9_]+[.]o" lapackmakematch "${lapackmakematch}")
          string(REGEX REPLACE "([a-zA-Z0-9_]+)[.]o" "lapack/\\1.c" lapackmakematch "${lapackmakematch}")
          string(REPLACE ";" " " blasmakematch "${blasmakematch}")
          string(REPLACE ";" " " lapackmakematch "${lapackmakematch}")


          set(f2cmake "#Auto-generated CMake build file for f2cblaslapack libraries
cmake_minimum_required(VERSION 2.8)
project(f2cblaslapack C)

add_definitions(-U__LAPACK_PRECISION_QUAD) #remove dependency on quadmath.h
#Note: -MTd => static link with debugging information, -wd4996 => disable insecure api warnings,
# -wd4244 => conversion from 'real' to 'integer', possible loss of data etc.
# -wd4554 => possible operator precedence error warning
#and -Z7 => embed debugging info in library as opposed to using an external .pdb database
add_definitions (-MTd -wd4996 -wd4244 -wd4554 -Z7)
set(CMAKE_INCLUDE_DIRECTORIES_BEFORE ON)#Allows us to change the default ordering of include directory searches
include_directories(\${CMAKE_CURRENT_BINARY_DIR}/blas)
set(BLAS_SOURCES ${blasmakematch})
set(LAPACK_SOURCES ${lapackmakematch})
add_library(f2cblas STATIC \${BLAS_SOURCES})
include_directories(BEFORE \${CMAKE_CURRENT_BINARY_DIR}/lapack) #prepend this directory so it is searched first
add_library(f2clapack STATIC \${LAPACK_SOURCES})
install(TARGETS f2cblas f2clapack DESTINATION lib)
")
          file(WRITE "${DOWNLOAD_DIR}/f2cblaslapack/${f2c_basicname}/CMakeLists.txt" "${f2cmake}")
      endif()


          #Build F2CBlasLapack as an external project
        externalproject_add(F2CBlasLapack_${f2c_basicname}
            SOURCE_DIR ${DOWNLOAD_DIR}/f2cblaslapack/${f2c_basicname}
            CMAKE_GENERATOR ${CMAKE_GENERATOR}
            BINARY_DIR ${CMAKE_BINARY_DIR}/f2cblaslapack_${f2c_basicname}
            INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/f2cblaslapack_${f2c_basicname}
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/f2cblaslapack_${f2c_basicname}
                -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
                -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
                -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
                -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
                -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
                -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
                -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
                -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
        )
        set(OUTPUT_LIB_DIR "${OUTPUT_LIB_DIR}" "${CMAKE_INSTALL_PREFIX}/f2cblaslapack_${f2c_basicname}/lib")
        set(OUTPUT_INCLUDE_DIR "${OUTPUT_INCLUDE_DIR}" "${CMAKE_INSTALL_PREFIX}/f2cblaslapack_${f2c_basicname}/include")

      #The original tmp directory on my platform was causing win32fe to fail
      #fix by explicitly setting the following environment variables.
      #Although, win32fe is no longer used, leave the following, as it is still useful to control this environment variable.
      if(NOT EXISTS "${DOWNLOAD_DIR}/petsc/${basicname}/tmp")
        file(MAKE_DIRECTORY "${DOWNLOAD_DIR}/petsc/${basicname}/tmp")
    endif()
      string(REGEX REPLACE "([a-zA-Z]):" "/cygdrive/\\1" TEMP_DIR_CYG "${DOWNLOAD_DIR}/petsc/${basicname}/tmp")
    set(ENV{TMPDIR} "${TEMP_DIR_CYG}")
    set(ENV{TMP} "${TEMP_DIR_CYG}")
    set(ENV{TEMP} "${TEMP_DIR_CYG}")

    #ensure that msmpi library is available to and usable by i686-w64-mingw32-gcc
    if(NOT EXISTS "${MS_HPC_PACK_DIR}/Lib/amd64/libmsmpi.a")
        message(STATUS "Generating ${MS_HPC_PACK_DIR}/Lib/amd64/libmsmpi.a")
        file(COPY "C:/Windows/System32/msmpi.dll" DESTINATION "${MS_HPC_PACK_DIR}/Lib/amd64")

        # Generate the msmpi.def file
        set(C_COMMAND ${CYGWINBASH} -c "
            export PATH=\"/usr/i686-w64-mingw32/bin:/usr/i686-w64-mingw32/sys-root/mingw/bin:/usr/local/bin:/usr/bin:${CYG_BIN}:${MS_HPC_PACK_DIR_CYG}/Bin:$PATH\"
            gendef msmpi.dll
        ")
        execute_process(
        COMMAND ${C_COMMAND}
        WORKING_DIRECTORY "${MS_HPC_PACK_DIR}/Lib/amd64"
        OUTPUT_VARIABLE cyg_log_out
        ERROR_VARIABLE cyg_log_err
        RESULT_VARIABLE cyg_result
    )

        # Fix an issue with calling convention mismatch
        message(STATUS "Fixing an issue with calling convention mismatch")
        file(READ "${MS_HPC_PACK_DIR}/Lib/amd64/msmpi.def" deffile)
        string(REGEX REPLACE "MPI_([a-zA-Z0-9_]+)@[0-9]+" "MPI_\\1" deffile_patch "${deffile}")
        file(WRITE "${MS_HPC_PACK_DIR}/Lib/amd64/msmpi.def" "${deffile_patch}")

        # Generate libmsmpi.a
        set(C_COMMAND ${CYGWINBASH} -c "
            export PATH=\"/usr/i686-w64-mingw32/bin:/usr/i686-w64-mingw32/sys-root/mingw/bin:/usr/local/bin:/usr/bin:${CYG_BIN}:${MS_HPC_PACK_DIR_CYG}/Bin:$PATH\"
            i686-w64-mingw32-dlltool -d msmpi.def -D msmpi.dll -l libmsmpi.a
        ")
        execute_process(
        COMMAND ${C_COMMAND}
        WORKING_DIRECTORY "${MS_HPC_PACK_DIR}/Lib/amd64"
        OUTPUT_VARIABLE cyg_log_out
        ERROR_VARIABLE cyg_log_err
        RESULT_VARIABLE cyg_result
    )
    message("Outputs and result of command ${C_COMMAND}:\n
Standard error
==============
${cyg_log_err}
Standard output
===============
${cyg_log_out}\n
Result = ${cyg_result}\n
=======End of Outputs for\n${C_COMMAND}\n
====================================================")
    endif()
    set(PETSC_ARCH "WINDOWS_BUILD" CACHE STRING
        "PETSc Architecture. Usually determines which build of PETSc to use")
    #The script that configures PETSc in Cygwin
      set(script "
        cd ${DOWNLOAD_DIR}/petsc/${basicname}
        export PETSC_DIR=`pwd`
        echo $PETSC_DIR
        export PETSC_ARCH=\"${PETSC_ARCH}\"
        export PATH=\"/usr/i686-w64-mingw32/bin:/usr/i686-w64-mingw32/sys-root/mingw/bin:/usr/local/bin:/usr/bin:${CYG_BIN}:${MS_HPC_PACK_DIR_CYG}/Bin:$PATH\"
        export TMPDIR=\"${TEMP_DIR_CYG}\"
        export TEMP=\"${TEMP_DIR_CYG}\"
        export TMP=\"${TEMP_DIR_CYG}\"
        config/configure.py --with-cc=i686-w64-mingw32-gcc --with-fc=0 --with-log=1 --with-info=1 --with-shared-libraries=0 --download-f2cblaslapack --useThreads=0 --with-mpi-include=${MS_HPC_PACK_DIR_CYG}/Inc --with-mpi-lib=${MS_HPC_PACK_DIR_CYG}/Lib/amd64/libmsmpi.a")

    if(EXISTS "${DOWNLOAD_DIR}/petsc/${basicname}/${PETSC_ARCH}")
        message(STATUS
        "It seems PETSc has already been configured. If you wish to create a fresh configure, you can either delete the folder ${DOWNLOAD_DIR}/petsc/${basicname}/${PETSC_ARCH} or change the variable PETSC_ARCH")
    else()
        message(STATUS "WAIT while PETSc is being configured. Time for a cuppa: this could take a while ...")
        set(C_COMMAND ${CYGWINBASH} -c "${script}")
          execute_process(
            COMMAND ${C_COMMAND}
            WORKING_DIRECTORY "${DOWNLOAD_DIR}/petsc/${basicname}"
            OUTPUT_VARIABLE cyg_log_out
            ERROR_VARIABLE cyg_log_err
            RESULT_VARIABLE cyg_result
        )
        message("Outputs and result of command ${C_COMMAND}:\n
Standard error\n
============\n
${cyg_log_err}
Standard output\n
=============\n
${cyg_log_out}\n
Result = ${cyg_result}\n
=======End of Outputs for\n${C_COMMAND}\n
====================================================")
            if(NOT cyg_result EQUAL 0)
                message("PETSc configuration has failed.")
            endif()
        endif(EXISTS "${DOWNLOAD_DIR}/petsc/${basicname}/${PETSC_ARCH}")

        #Delete CMake Cache generated by PETSc configuration.
        execute_process(
        COMMAND ${CMAKE_COMMAND} -E remove -f CMakeCache.txt
        WORKING_DIRECTORY "${DOWNLOAD_DIR}/petsc/${basicname}/${PETSC_ARCH}")

        set(petsc_basicname "${basicname}")

        #Patch and customise PETSc the configuration header file
        file(READ  "${CMAKE_CURRENT_SOURCE_DIR}/petscconf.h.in" petscconf)
        string(REGEX REPLACE "#[ ]*define[ ]+PETSC_DIR[ \t][^\n]+\n" "#define PETSC_DIR \"${DOWNLOAD_DIR}/petsc/${basicname}\" \n" petscconf "${petscconf}")
        string(REGEX REPLACE "#[ ]*define[ ]+PETSC_ARCH[ \t][^\n]+\n" "#define PETSC_ARCH \"${PETSC_ARCH}\" \n" petscconf "${petscconf}")
        string(REGEX REPLACE "#[ ]*define[ ]+PETSC_LIB_DIR[ \t][^\n]+\n" "#define PETSC_LIB_DIR \"${CMAKE_INSTALL_PREFIX}/petsc_${petsc_basicname}/lib\" \n" petscconf "${petscconf}")
        file(WRITE "${DOWNLOAD_DIR}/petsc/${basicname}/${PETSC_ARCH}/include/petscconf.h" "${petscconf}")

        # PATCH generated PETSc CMake files to enable PARMETIS and set up library locations properly for
        # native Windows build
        set(static_build "#Note: -MTd => static link with debugging information, -wd4996 => disable insecure api warnings,
#and -Z7 => embed debugging info in library as opposed to using an external .pdb database
#-wd4005 => disable macro redefinition
#-wd4305 => truncation from type1 to type2
#-wd4133 => 'function': incompatible types from type1 to type2
#-wd4267 => possible loss of data: conversion from type1 to type2
#-wd4244 => another possible loss of data
#-wd4101 => unreferenced local variable
add_definitions (-MTd -wd4996 -Z7 -wd4005 -wd4305 -wd4133 -wd4267 -wd4244 -wd4101)
include_directories(\"${CMAKE_INSTALL_PREFIX}/parmetis_${parmetis_basicname}/include\" \"${CMAKE_INSTALL_PREFIX}/metis_${metis_basicname}/include\")
")


        file(READ "${DOWNLOAD_DIR}/petsc/${basicname}/CMakeLists.txt" petscmake)
        #check whether we have patched this already
        string(FIND "${petscmake}" "${patch_message}" patched)
        if(patched EQUAL -1)#not patched yet
            string(REGEX REPLACE "(add_definitions[ ]*)[(](-D__INSDIR__[ ]*=)[ ]*[)](.*)$"
                "${static_build}\\1(\\2./)\\3" petscmake_patch "${petscmake}")
            file(WRITE "${DOWNLOAD_DIR}/petsc/${basicname}/CMakeLists.txt" "${patch_message}${petscmake_patch}\n
install(TARGETS petsc DESTINATION lib)
install(DIRECTORY \"${DOWNLOAD_DIR}/petsc/${basicname}/${PETSC_ARCH}/include/\"
DESTINATION include
FILES_MATCHING PATTERN \"*.h\")

install(DIRECTORY \"${DOWNLOAD_DIR}/petsc/${basicname}/include/\"
DESTINATION include
FILES_MATCHING PATTERN \"*.h\" PATTERN \"*.hh\")
")
        endif()


        # enable PARMETIS and other Windows-specific PETSc options
        file(READ "${CMAKE_CURRENT_SOURCE_DIR}/PETScConfig.cmake.in" petsconfig)
        string(REGEX REPLACE "set[ ]*[(][ ]*MS_HPC_PACK_LIB64[^\n]+\n"
            "
set(BLASLAPACK_HINT \"${CMAKE_INSTALL_PREFIX}/f2cblaslapack_${f2c_basicname}/lib\")
set(PARMETIS_HINT \"${CMAKE_INSTALL_PREFIX}/parmetis_${parmetis_basicname}/lib\")
set(MS_HPC_PACK_LIB64 \"${MS_HPC_PACK_DIR}/Lib/amd64\") \n"
            petsconfig "${petsconfig}")
        string(REGEX REPLACE "set[ ]*[(][ ]*MS_HPC_PACK_INCLUDES[^\n]+\n" "set(MS_HPC_PACK_INCLUDES \"${MS_HPC_PACK_DIR}/Inc\") \n" petsconfig "${petsconfig}")
        string(REGEX REPLACE "set[ ]*[(][ ]*PETSC_LIBRARIES[^\n]+\n" "set(PETSC_LIBRARIES \"${DOWNLOAD_DIR}/petsc/${basicname}/${PETSC_ARCH}/lib\") \n" petsconfig "${petsconfig}")
        file(WRITE "${DOWNLOAD_DIR}/petsc/${petsc_basicname}/${PETSC_ARCH}/conf/PETScConfig.cmake" "${petsconfig}")

        #Build PETSc itself!

        externalproject_add(PETSc_${petsc_basicname}
            SOURCE_DIR ${DOWNLOAD_DIR}/petsc/${petsc_basicname}
            CMAKE_GENERATOR ${CMAKE_GENERATOR}
            BINARY_DIR ${DOWNLOAD_DIR}/petsc/${petsc_basicname}/${PETSC_ARCH}
            INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/petsc_${petsc_basicname}
            CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/petsc_${petsc_basicname}
                -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
                -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
                -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
                -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
                -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
                -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
                -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
                -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
                -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
                -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
                -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
                -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
                -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
                -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
            DEPENDS ParMETIS_${parmetis_basicname} F2CBlasLapack_${f2c_basicname}
        )
        set(OUTPUT_LIB_DIR "${OUTPUT_LIB_DIR}" "${CMAKE_INSTALL_PREFIX}/petsc_${petsc_basicname}/lib")
        set(OUTPUT_INCLUDE_DIR "${OUTPUT_INCLUDE_DIR}" "${CMAKE_INSTALL_PREFIX}/petsc_${petsc_basicname}/include")

	#generate cmake config file for petsc
 	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/PETScConfigForChaste.cmake.in
  			"${CMAKE_INSTALL_PREFIX}/petsc_${petsc_basicname}/PETScConfig.cmake" @ONLY)

        #Delete an errant directory created by PETSc
        string(REGEX REPLACE "([a-zA-Z]:)[^\n]+" "\\1/" WORKING_DIR "${DOWNLOAD_DIR}")
        execute_process(
        COMMAND ${CMAKE_COMMAND} -E remove_directory cygdrive
        WORKING_DIRECTORY "${WORKING_DIR}")

        #Restore TEMP directory to native WINDOWS format
        set(TEMP_DIR "${DOWNLOAD_DIR}/petsc/${basicname}/tmp")
        set(ENV{TMPDIR} "${TEMP_DIR}")
        set(ENV{TMP} "${TEMP_DIR}")
        set(ENV{TEMP} "${TEMP_DIR}")

endforeach(url ${PETSC_URLS})
endif(BUILD_PETSC)

#Build the HDF5s
foreach(url ${HDF5_URLS})

 download_from_url(${url} "hdf5")
 set(hdf5_basicname "${bname}")

#Build HDF5 as an external project
externalproject_add(HDF5_${hdf5_basicname}
    SOURCE_DIR ${DOWNLOAD_DIR}/hdf5/${hdf5_basicname}
    CMAKE_GENERATOR ${CMAKE_GENERATOR}
    BINARY_DIR ${CMAKE_BINARY_DIR}/hdf5_${hdf5_basicname}
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/hdf5_${hdf5_basicname}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/hdf5_${hdf5_basicname}
        -DHDF5_ENABLE_PARALLEL:BOOL=ON
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
        -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
        -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
        -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
        -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
        -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
        -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}

        )
set(OUTPUT_LIB_DIR "${OUTPUT_LIB_DIR}" "${CMAKE_INSTALL_PREFIX}/hdf5_${hdf5_basicname}/lib")
set(OUTPUT_INCLUDE_DIR "${OUTPUT_INCLUDE_DIR}" "${CMAKE_INSTALL_PREFIX}/hdf5_${hdf5_basicname}/include")

endforeach(url ${HDF5_URLS})


#Build the Sundials
foreach(url ${SUNDIALS_URLS})

 download_from_url(${url} "sundials")
 set(sundials_basicname "${bname}")

#Build sundials as an external project
    message("trying to build from ${DOWNLOAD_DIR}/sundials/${sundials_basicname}")
    externalproject_add(Sundials_${sundials_basicname}
        SOURCE_DIR ${DOWNLOAD_DIR}/sundials/${sundials_basicname}
        CMAKE_GENERATOR ${CMAKE_GENERATOR}
        BINARY_DIR ${CMAKE_BINARY_DIR}/sundials_${sundials_basicname}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/sundials_${sundials_basicname}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/sundials_${sundials_basicname}
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
            -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
            -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
            -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
            -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
            -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
            -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
            -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}

        )

    file(GLOB_RECURSE SUNDIALS_INC ${CMAKE_INSTALL_PREFIX}/sundials_${sundials_basicname}/include/*)
    set(SUNDIALS_INC_DIRS "${CMAKE_INSTALL_PREFIX}/sundials_${sundials_basicname}/include")
    foreach(f ${SUNDIALS_INC})
        get_filename_component(p ${f} PATH)
        list(APPEND SUNDIALS_INC_DIRS ${p})
    endforeach()
    list(REMOVE_DUPLICATES SUNDIALS_INC_DIRS)
set(OUTPUT_LIB_DIR "${OUTPUT_LIB_DIR}" "${CMAKE_INSTALL_PREFIX}/sundials_${sundials_basicname}/lib")
set(OUTPUT_INCLUDE_DIR "${OUTPUT_INCLUDE_DIR}" "${SUNDIALS_INC_DIRS}")
endforeach(url ${SUNDIALS_URLS})


# Build VTK as an external project
foreach(url ${VTK_URLS})

    download_from_url(${url} "vtk")
    set(vtk_basicname "VTK")

    message("Building VTK with flags:")
    message("   C: ${CMAKE_C_FLAGS}; debug: ${CMAKE_C_FLAGS_DEBUG}")
    message("   CXX: ${CMAKE_CXX_FLAGS}; debug: ${CMAKE_CXX_FLAGS_DEBUG}")
    get_directory_property(DirDefs DEFINITIONS)
    string(REPLACE "-" "/" DirDefs "${DirDefs}")
    message("   DEFS: ${DirDefs}")

    externalproject_add(Vtk_${vtk_basicname}
        SOURCE_DIR ${DOWNLOAD_DIR}/vtk/${vtk_basicname}
        CMAKE_GENERATOR ${CMAKE_GENERATOR}
        BINARY_DIR ${CMAKE_BINARY_DIR}/vtk_${vtk_basicname}
        INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/vtk_${vtk_basicname}
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}/vtk_${vtk_basicname}
            -DBUILD_SHARED_LIBS:BOOL=OFF
            -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
            -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG} /Z7
            -DCMAKE_C_FLAGS_MINSIZEREL:STRING=${CMAKE_C_FLAGS_MINSIZEREL}
            -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
            -DCMAKE_C_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_C_FLAGS_RELWITHDEBINFO}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG} /Z7
            -DCMAKE_CXX_FLAGS_MINSIZEREL:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
            -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
            -DCMAKE_CXX_FLAGS_RELWITHDEBINFO:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
            -DCOMPILE_DEFINITIONS:STRING=${COMPILE_DEFINITIONS}
        )

#    file(GLOB_RECURSE VTK_INC ${CMAKE_INSTALL_PREFIX}/vtk_${vtk_basicname}/include/*)
    set(VTK_INC_DIRS "${CMAKE_INSTALL_PREFIX}/vtk_${vtk_basicname}/include/vtk-5.8")
#    foreach(f ${VTK_INC})
#        get_filename_component(p ${f} PATH)
#        list(APPEND VTK_INC_DIRS ${p})
#    endforeach()
    list(REMOVE_DUPLICATES VTK_INC_DIRS)
set(OUTPUT_LIB_DIR "${OUTPUT_LIB_DIR}" "${CMAKE_INSTALL_PREFIX}/vtk_${vtk_basicname}/lib/vtk-5.8")
set(OUTPUT_INCLUDE_DIR "${OUTPUT_INCLUDE_DIR}" "${VTK_INC_DIRS}")
endforeach(url ${VTK_URLS})




# Write a file in the Chaste source tree containing include and link paths for all the libraries installed

if(NOT CHASTE_SOURCE_DIR)
    set(CHASTE_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..)
endif()

file(WRITE "${CHASTE_SOURCE_DIR}/AdditionalIncludes.cmake"
"
include_directories(\"${OUTPUT_INCLUDE_DIR}\")

link_directories(\"${OUTPUT_LIB_DIR}\")
")
