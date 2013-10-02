project(testio)

set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT /NODEFAULTLIB:MSVCRTD")
set(HAS_OWN_LINKER_FLAGS TRUE)

chaste_do_test_component(io)
