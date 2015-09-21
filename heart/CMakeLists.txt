find_package(Chaste COMPONENTS ${CHASTE_DEPENDS_heart})

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
set(CHASTE_heart_INCLUDE_DIRS ${CHASTE_heart_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/src/io PARENT_SCOPE)
set(CHASTE_heart_INCLUDE_DIRS ${CHASTE_heart_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/src/io)

file(GLOB cellml_files src/odes/cellml/*.cellml)
foreach(cellml_file ${cellml_files})
    chaste_do_cellml(Chaste_heart_SOURCES ${cellml_file} OFF  "--output-dir" "${CMAKE_CURRENT_BINARY_DIR}/src/odes/cellml")
endforeach()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src/odes/cellml)
set(CHASTE_heart_INCLUDE_DIRS ${CHASTE_heart_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/src/odes/cellml PARENT_SCOPE)
set(CHASTE_heart_INCLUDE_DIRS ${CHASTE_heart_INCLUDE_DIRS} ${CMAKE_CURRENT_BINARY_DIR}/src/odes/cellml)

# Make dynamic directory for dynamic cellml compilation
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dynamic)

chaste_do_component(heart ${Chaste_heart_SOURCES})

# make sure to link dl libraries
target_link_libraries(heart LINK_PUBLIC ${CMAKE_DL_LIBS})

# Compile dynamic models
set(dyn_sources dynamic/DynamicallyLoadableLr91.cpp dynamic/NotACellModel.cpp)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dynamic)
foreach(dyn_source ${dyn_sources})
    get_filename_component(dyn_source_name ${dyn_source} NAME_WE)
    add_library(${dyn_source_name} SHARED ${dyn_source})
    set_target_properties(${dyn_source_name} PROPERTIES LINK_FLAGS "-Wl,--undefined")
endforeach(dyn_source)
