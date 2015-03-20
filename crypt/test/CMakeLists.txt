project(testcrypt)

if(MSVC)
set(HAS_OWN_LINKER_FLAGS TRUE)
set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:MSVCRT /NODEFAULTLIB:MSVCRTD")
endif()

chaste_do_test_component(crypt)
