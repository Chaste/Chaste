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


project(testmesh)
if (MSVC)
set(HAS_OWN_LINKER_FLAGS TRUE)
#GKLib uses/does C-style structured exception handling, using longjmp,
#which caused untold headaches -- with MSVC linker errors.
#In the end the solution was to comment out the line
#longjmp(gk_jbufs[gk_cur_jbufs], signum);
#in the functions "void gk_sigthrow(int signum)" and "void gk_NonLocalExit_Handler(int signum)", 
#both of which are in the file "error.c" under the GKLib directory under the metis source directory.

set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:MSVCRTD") # undefined __imp_longjmp - FIXED by description above
endif(MSVC)
#Various other combinations that didn't work
#set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:LIBCMTD")#libcpmtd, oldnames -- unresolved externals
#set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:MSVCRTD /NODEFAULTLIB:MSVCRT")# undefined __imp_longjmp
#set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:LIBCPMTD /NODEFAULTLIB:MSVCRT")#MSVCRTD clashing with LIBCMTD
#set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:LIBCPMTD /NODEFAULTLIB:MSVCRTD")# unresolved externals std::_Lockit
#set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT")#MSVCRTD clashing with LIBCMTD
#set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMTD")#LIBCMT clashing with other libs
#set(LINKER_FLAGS "/NODEFAULTLIB:MSVCRTD")#LIBCMT clashing with LIBCMTD
#set(LINKER_FLAGS "/NODEFAULTLIB:MSVCRT")#MSVCRTD clashing with LIBCMTD
#set(LINKER_FLAGS "/NODEFAULTLIB:MSVCRT /NODEFAULTLIB:LIBCMTD /NODEFAULTLIB:LIBCMT") unresolved externals atexit, operators new, delete etc.
#set(LINKER_FLAGS "/NODEFAULTLIB:MSVCRT /NODEFAULTLIB:MSVCRTD /NODEFAULTLIB:LIBCMTD")#LIBCMT clashes with other libs
#set(LINKER_FLAGS "/NODEFAULTLIB:MSVCRT /NODEFAULTLIB:MSVCRTD /NODEFAULTLIB:LIBCMTD /NODEFAULTLIB:LIBCMT")#unresolved externa; new
#set(LINKER_FLAGS "/NODEFAULTLIB:MSVCRTD /NODEFAULTLIB:LIBCMTD /NODEFAULTLIB:LIBCMT")# unresolved externals new etc

chaste_do_test_component(mesh)
