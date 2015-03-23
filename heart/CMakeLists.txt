project(heart)


#Generate ChasteParameters_*.cpp and hpp
set(xsd_dir ${CMAKE_CURRENT_SOURCE_DIR}/src/io)
set(out_dir ${Chaste_BINARY_DIR}/src/heart)
set(CHASTE_heart_INCLUDE_DIRS ${CHASTE_heart_INCLUDE_DIRS} ${out_dir})
file(MAKE_DIRECTORY ${out_dir})

file(GLOB xsd_config_files RELATIVE ${xsd_dir} src/io/ChasteParameters_*.xsd)
set(CHASTE_heart_SOURCES "")
foreach(xsd_config ${xsd_config_files})
    get_filename_component(xsd_config_name ${xsd_config} NAME_WE)

    add_custom_command(OUTPUT ${out_dir}/${xsd_config_name}.cpp ${out_dir}/${xsd_config_name}.hpp
                        COMMAND "${XSD_EXECUTABLE}" "cxx-tree"
                                "--generate-serialization"
                                "--output-dir" "${out_dir}"
                                "--hxx-suffix" ".hpp" 
                                "--cxx-suffix" ".cpp"
                                "--prologue-file" "XsdPrologue.txt"
                                "--epilogue-file" "XsdEpilogue.txt"
                                "--namespace-regex" "X.* $Xchaste::parametersX"
                                "--namespace-regex" "X.* https://chaste.comlab.ox.ac.uk/nss/parameters/(.+)Xchaste::parameters::v$1X"
                                "${xsd_config}"
                    WORKING_DIRECTORY ${xsd_dir}
                    DEPENDS ${xsd_dir}/${xsd_config}
                    COMMENT "Generating ${xsd_config_name}.cpp and ${xsd_config_name}.hpp"
                    VERBATIM
                    )
    set(CHASTE_heart_SOURCES ${CHASTE_heart_SOURCES} "${out_dir}/${xsd_config_name}.cpp")
 
endforeach()

#Generate CellML models
set(cellml_dir ${CMAKE_CURRENT_SOURCE_DIR}/src/odes/cellml)
set(out_dir ${Chaste_BINARY_DIR}/src/heart/cellml)
set(CHASTE_heart_INCLUDE_DIRS ${CHASTE_heart_INCLUDE_DIRS} ${out_dir})
file(MAKE_DIRECTORY ${out_dir})
set(Chaste_python_dir ${Chaste_SOURCE_DIR}/python)

file(GLOB cellml_files RELATIVE ${Chaste_python_dir} src/odes/cellml/*.cellml)
foreach(cellml_file ${cellml_files})
    get_filename_component(cellml_file_name ${cellml_file} NAME_WE)
    set(pycml_args "-A" "-p" "--output-dir" "${out_dir}")
    #if(BUILD_SHARED_LIBS)
    #    set(pycml_args ${pycml_args} "-y")
    #else()
        set(pycml_args ${pycml_args} "--normal" "--opt" "--cvode")
    #endif()
    set(depends ${cellml_dir}/${cellml_file_name}.cellml)
    if(EXISTS ${cellml_dir}/${cellml_file_name}-conf.xml)
        set(depends ${depends} ${cellml_dir}/${cellml_file_name}-conf.xml)
        set(pycml_args ${pycml_args} "--conf=${cellml_dir}/${cellml_file_name}-conf.xml")
    endif()
    if(EXISTS ${cellml_dir}/${cellml_file_name}.out)
        set(depends ${depends} ${cellml_dir}/${cellml_file_name}.out)
        set(pycml_args ${pycml_args} "--backward-euler")
    endif()

    execute_process(COMMAND "${PYTHON_EXECUTABLE}" ConvertCellModel.py ${pycml_args} --show-outputs ${cellml_file}   
                        WORKING_DIRECTORY ${Chaste_python_dir}
                        OUTPUT_VARIABLE ConvertCellModelDepends
			OUTPUT_STRIP_TRAILING_WHITESPACE
                        )

    string(REGEX MATCHALL "[^\n]*\\.hpp" output_files_hpp "${ConvertCellModelDepends}")
    string(REGEX MATCHALL "[^\n]*\\.cpp" output_files_cpp "${ConvertCellModelDepends}")

    add_custom_command(OUTPUT ${output_files_hpp} ${output_files_cpp} 
                        COMMAND "${PYTHON_EXECUTABLE}" ConvertCellModel.py ${pycml_args} ${cellml_file}
                        WORKING_DIRECTORY ${Chaste_python_dir}
                        DEPENDS ${depends}
                        COMMENT "Generating ${cellml_file_name}*.cpp and ${cellml_file_name}*.hpp"
                        VERBATIM
                        )
    set(CHASTE_heart_SOURCES ${CHASTE_heart_SOURCES} ${output_files_cpp})
    
endforeach()

chaste_do_component(heart ${CHASTE_heart_SOURCES})
