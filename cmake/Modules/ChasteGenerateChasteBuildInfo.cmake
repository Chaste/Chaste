if(MSVC)
    set(additional
    "
    //This has been put here to satisfy an MSVC linker issue 
    int __cdecl _purecall(void){return 0;}
    ")
endif()

configure_file (
  "${Chaste_SOURCE_DIR}/global/src/ChasteBuildInfo_cmake.cpp.in"
  ChasteBuildInfo.cpp
)
