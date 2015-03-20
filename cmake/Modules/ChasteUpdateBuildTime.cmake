file(READ build_timestamp build_time)
file(READ Version.cpp version_cpp)
string(REGEX REPLACE "GetBuildTime\\(\\)[\r\n\t ]*{[^}]*}" 
"GetBuildTime()
{
    return \"${build_time}\";
}
" version_new_cpp "${version_cpp}")
file(WRITE Version.cpp "${version_new_cpp}")
