# Ensure consistent view of static/thread/debug settings.
# See e.g. http://msdn.microsoft.com/en-us/library/2kzt1wy3.aspx

set(CMAKE_USER_MAKE_RULES_OVERRIDE
   ${CMAKE_CURRENT_LIST_DIR}/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
   ${CMAKE_CURRENT_LIST_DIR}/cxx_flag_overrides.cmake)

