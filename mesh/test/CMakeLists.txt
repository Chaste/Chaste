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
