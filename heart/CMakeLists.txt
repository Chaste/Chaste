find_package(Chaste COMPONENTS ${Chaste_DEPENDS_heart})

set(Chaste_heart_SOURCES "")

file(GLOB xsd_config_files src/io/*.xsd)
foreach(xsd_config ${xsd_config_files})
    xsd_schema(Chaste_heart_SOURCES ${xsd_config} 
                "--generate-serialization"
                "--hxx-suffix" ".hpp" 
                "--cxx-suffix" ".cpp"
                "--prologue-file" "${CMAKE_CURRENT_SOURCE_DIR}/src/io/XsdPrologue.txt"
                "--epilogue-file" "${CMAKE_CURRENT_SOURCE_DIR}/src/io/XsdEpilogue.txt"
                "--namespace-regex" "X.* $Xchaste::parametersX"
                "--namespace-regex" "X.* https://chaste.comlab.ox.ac.uk/nss/parameters/(.+)Xchaste::parameters::v$1X"
                )

endforeach()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/io)
set(Chaste_heart_INCLUDE_DIRS ${Chaste_heart_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/src/io)

# Make dynamic directory for dynamic cellml compilation
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dynamic)

# XSD < 4.0.0 generates code with std::auto_ptr.  This is deprecated from C++11 (in favour of new standard smart
# pointers). This code turns off the deprecation warnings so that Chaste compiles.
# TODO: Once we drop Chaste support for XSD 3.3, this should be removed in favour of changing auto_ptr to unique_ptr
foreach(xsd_generated_cpp_file ${Chaste_heart_SOURCES})
    set_source_files_properties(${xsd_generated_cpp_file} PROPERTIES COMPILE_FLAGS "-Wno-deprecated-declarations")
endforeach(xsd_generated_cpp_file ${Chaste_heart_SOURCES})

chaste_do_component(heart ${Chaste_heart_SOURCES})

# make sure to link dl libraries
target_link_libraries(chaste_heart LINK_PUBLIC ${CMAKE_DL_LIBS})

# Compile dynamic models
set(dyn_sources dynamic/DynamicallyLoadableLr91.cpp dynamic/NotACellModel.cpp)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dynamic)
foreach(dyn_source ${dyn_sources})
    get_filename_component(dyn_source_name ${dyn_source} NAME_WE)
    add_library(${dyn_source_name} SHARED ${dyn_source})
    set_property(TARGET ${dyn_source_name} PROPERTY CXX_STANDARD 11)
    set_property(TARGET ${dyn_source_name} PROPERTY CXX_STANDARD_REQUIRED ON)
    if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        target_link_libraries(${dyn_source_name} LINK_PRIVATE "-Wl,-undefined,dynamic_lookup")
        #set_target_properties(${dyn_source_name} PROPERTIES LINK_FLAGS "-Wl,-undefined,dynamic_lookup")
    else()
        #target_link_libraries(${dyn_source_name} LINK_PRIVATE "-Wl,--undefined")
        set_target_properties(${dyn_source_name} PROPERTIES LINK_FLAGS "-Wl,--undefined")
    endif()

    # add dynamic models to the heart_tests and Continuous targets
    add_dependencies(heart ${dyn_source_name})
    add_dependencies(Continuous ${dyn_source_name})
endforeach(dyn_source)
