project(testlinalg)

if(MSVC)
set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:MSVCRTD")
set(HAS_OWN_LINKER_FLAGS TRUE)
endif(MSVC)

chaste_do_test_component(linalg)
