# A macro to recursively find Chaste header locations
macro(HEADER_DIRS base_dir return_list)
    set(new_list "")
    set(dir_list "")
    #message("base dir = ${base_dir}")

    file(GLOB_RECURSE new_list ${base_dir}/*.hpp ${base_dir}/*.h)
    #message("new list = ${new_list}")
    foreach(file_path ${new_list})
        get_filename_component(dir_path ${file_path} PATH)
        set(dir_list ${dir_list} ${dir_path})
    endforeach()
    
    list(REMOVE_DUPLICATES dir_list)
    

    #message("return list = ${return_list}")
    #message("dir list = ${dir_list}")

    set(${return_list} ${dir_list})
endmacro()

macro(CHASTE_DO_CELLML output_sources cellml_file dynamic)
    get_filename_component(cellml_file_name ${cellml_file} NAME_WE)
    get_filename_component(cellml_dir ${cellml_file} DIRECTORY)
    file(RELATIVE_PATH cellml_file_rel "${CMAKE_SOURCE_DIR}" "${cellml_file}")
    set(pycml_args "-A" "-p")
    set(pycml_args ${pycml_args} ${ARGN})
    #if(BUILD_SHARED_LIBS)
    if (${dynamic})
        set(pycml_args ${pycml_args} "-y")
    else()
        set(pycml_args ${pycml_args} "--normal" "--opt" "--cvode")
        if(EXISTS ${cellml_dir}/${cellml_file_name}.out)
            set(depends ${depends} ${cellml_dir}/${cellml_file_name}.out)
            set(pycml_args ${pycml_args} "--backward-euler")
        endif()
    endif()
    set(depends ${cellml_dir}/${cellml_file_name}.cellml)
    if(EXISTS ${cellml_dir}/${cellml_file_name}-conf.xml)
        set(depends ${depends} ${cellml_dir}/${cellml_file_name}-conf.xml)
        set(pycml_args ${pycml_args} "--conf=${cellml_dir}/${cellml_file_name}-conf.xml")
    endif()
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" ${Chaste_PYTHON_DIR}/ConvertCellModel.py ${pycml_args} --show-outputs ${cellml_file}   
        OUTPUT_VARIABLE ConvertCellModelDepends
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    string(REGEX MATCHALL "[^\n]*\\.hpp" output_files_hpp "${ConvertCellModelDepends}")
    string(REGEX MATCHALL "[^\n]*\\.cpp" output_files_cpp "${ConvertCellModelDepends}")

    if (NOT Chaste_VERBOSE)
        set(pycml_args ${pycml_args} "--quiet")
    endif()

    add_custom_command(OUTPUT ${output_files_hpp} ${output_files_cpp} 
        COMMAND "${PYTHON_EXECUTABLE}" ${Chaste_PYTHON_DIR}/ConvertCellModel.py ${pycml_args} ${cellml_file}
        DEPENDS ${depends}
        COMMENT "Processing CellML file ${cellml_file_rel}" 
        VERBATIM
        )

    list(APPEND ${output_sources} ${output_files_cpp} ${output_files_hpp})
endmacro()
    


    #Chaste Testing Macro. The predefined cxxtest_add_test is not suitable because of little control over
    #the test's working directory
    macro(CHASTE_ADD_TEST _testTargetName )
        string(REGEX MATCH "^.*Parallel$" foundParallel ${_testTargetName})
        if (foundParallel)
            set(parallel ON)
            string(REGEX REPLACE "(Parallel)$" "" _testname "${_testTargetName}")
        else()
            set(parallel OFF)
            set(_testname ${_testTargetName})
        endif()

        set(_exeTargetName ${_testname}Runner)

        if (NOT TARGET ${exeTargetName})
            set(_test_real_output_filename "${CMAKE_CURRENT_BINARY_DIR}/${_testname}.cpp")
            add_custom_command(
                OUTPUT "${_test_real_output_filename}"
                DEPENDS ${ARGN}
                COMMAND ${CXXTEST_TESTGEN_INTERPRETER} ${CXXTEST_TESTGEN_EXECUTABLE} ${CXXTEST_TESTGEN_ARGS} -o "${_test_real_output_filename}" ${ARGN}
            )

            set_source_files_properties("${_test_real_output_filename}" PROPERTIES GENERATED true)

            add_executable(${exeTargetName} "${_test_real_output_filename}" ${ARGN})
        endif()


        if(${parallel} OR NOT (${Chaste_NUM_CPUS_TEST} EQUAL 1))
            #Note: "${MPIEXEC} /np 1 master : subordinate" means that we run one master process and n subordinate processes
            # on the local host with n+1 cores.
            # Here we are using the form ${MPIEXEC} /np 2 ${test}.
            # A figure-it-out-yourselfnstalled libvtk-java and libvtk5-qt4-dev form would be ${MPIEXEC} /np * ${test} which runs on all available cores
            # See http://technet.microsoft.com/en-us/library/cc947675%28v=ws.10%29.aspx
            # Note the underscore appended to the test name, to match with the RUN_TESTS block above, and ensure we don't
            # run more tests than intended!
            if (${Chaste_NUM_CPUS_TEST} EQUAL 1)
                set(num_cpus 2)
            else()
                set(num_cpus ${Chaste_NUM_CPUS_TEST})
            endif()
            add_test(NAME ${_testTargetName} WORKING_DIRECTORY "${Chaste_SOURCE_DIR}/" COMMAND  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${num_cpus} ${MPIEXEC_PREFLAGS}  $<TARGET_FILE:${exeTargetName}> ${MPIEXEC_POSTFLAGS})
            set_property(TEST ${testTargetName} PROPERTY PROCESSORS ${num_cpus})
        else()
            add_test(NAME ${_testTargetName} WORKING_DIRECTORY "${Chaste_SOURCE_DIR}/" COMMAND $<TARGET_FILE:${exeTargetName}>)
            set_property(TEST ${testTargetName} PROPERTY PROCESSORS 1)
        endif()
    endmacro(CHASTE_ADD_TEST)

  macro(CHASTE_GENERATE_TEST_NAME test outTestName)
      string(REGEX REPLACE "([a-zA-Z0-9_/]+)[.]hpp" "\\1" testName "${test}")
      string(REPLACE "/" ";" testPath "${testName}")
      list(LENGTH testPath pathLength)
      if(${pathLength} EQUAL 1)
        set(testName ${testPath})
        set(testPath "")
        set(${outTestName} ${testName})
      else()
        math(EXPR index "${pathLength} - 1")
        list(GET testPath ${index} testName)
        list(REMOVE_AT testPath ${index})
        string(REPLACE ";" "_" _testPath_ "${testPath}")
        string(REPLACE ";" "/" testPath "${testPath}")
        set(${outTestName} "${testName}_${_testPath_}_")
      endif()
  endmacro(CHASTE_GENERATE_TEST_NAME test outTestName)


  macro(CHASTE_DO_COMMON component)
    

    # Figure out include path
    #if (NOT CHASTE_${component}_INCLUDE_DIRS)
    #    set(CHASTE_${component}_INCLUDE_DIRS 
    #        "${CMAKE_CURRENT_SOURCE_DIR}/src" 
    #        "${CMAKE_CURRENT_SOURCE_DIR}/test" 
    #        "${CMAKE_CURRENT_SOURCE_DIR}/apps/src"
    #        )
    #endif (NOT CHASTE_${component}_INCLUDE_DIRS)

    if (NOT CHASTE_${component}_INCLUDE_DIRS)
        set(CHASTE_${component}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
        header_dirs(${CHASTE_${component}_SOURCE_DIR} CHASTE_${component}_INCLUDE_DIRS)
    endif()

    if (Chaste_THIRD_PARTY_INCLUDE_DIRS)
        include_directories(SYSTEM "${Chaste_THIRD_PARTY_INCLUDE_DIRS}")
    endif()
    if (CHASTE_${component}_INCLUDE_DIRS)
        include_directories("${CHASTE_${component}_INCLUDE_DIRS}")
    endif()
    if (Chaste_INCLUDE_DIRS)
        include_directories("${Chaste_INCLUDE_DIRS}")
    endif()
    
    # Make component library
    file(GLOB_RECURSE CHASTE_${component}_SOURCES 
			RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
			${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp 
			${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)

    add_library(${component} ${CHASTE_${component}_SOURCES} ${ARGN})
    if (BUILD_SHARED_LIBS)
        target_link_libraries(${component} LINK_PUBLIC ${Chaste_LIBRARIES})
        set(static_extension "a")
        foreach(library ${Chaste_THIRD_PARTY_LIBRARIES})
            if (library MATCHES ".*\\.${static_extension}")
                target_link_libraries(${component} LINK_PRIVATE ${library})
            else()
                target_link_libraries(${component} LINK_PUBLIC ${library})
            endif()
        endforeach()
    else()
        target_link_libraries(${component} LINK_PUBLIC ${Chaste_THIRD_PARTY_LIBRARIES})
    endif()


    if(NOT(${component} MATCHES "^project"))
        install(TARGETS ${component} 
            DESTINATION lib COMPONENT ${component}_libraries)

        install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/"
            DESTINATION include/${component}
            COMPONENT ${component}_headers
            FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
            )
    endif()

    if (ENABLE_CHASTE_TESTING) 
        if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")
            set(ENABLE_${component}_TESTING ON CACHE BOOL "Generate the test infrastructure for ${component} ")
            mark_as_advanced(ENABLE_project_${projectName}_TESTING)
        else()
            message(WARNING "No CMakeLists.txt file found in test directory ${CMAKE_CURRENT_SOURCE_DIR}/test. Tests for ${component} will not be built")
            set(ENABLE_project_${projectName}_TESTING OFF)
        endif()

        # Do testing if requested
        if(ENABLE_${component}_TESTING)
            add_subdirectory(test)
        endif()
    endif(ENABLE_CHASTE_TESTING)

	# Build applications if present
	if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/apps")
        add_subdirectory(apps)
    endif()
  endmacro(CHASTE_DO_COMMON)


  macro(CHASTE_DO_COMPONENT component)
    message("Configuring component ${component}")
    CHASTE_DO_COMMON(${component} ${ARGN})
  endmacro(CHASTE_DO_COMPONENT)

  macro(CHASTE_DO_PROJECT projectName)
    message("Configuring project ${projectName}")
    CHASTE_DO_COMMON(project_${projectName})
  endmacro(CHASTE_DO_PROJECT)

  macro(CHASTE_DO_APPS_COMMON component)
    include_directories("${Chaste_THIRD_PARTY_INCLUDE_DIRS}" "${Chaste_INCLUDE_DIRS}" "${CXXTEST_INCLUDES}")
    file(GLOB CHASTE_${component}_APPS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/*.cpp)
    foreach(app ${CHASTE_${component}_APPS})
        string(REGEX REPLACE ".*/([a-zA-Z0-9_]+)[.]cpp" "\\1" appName "${app}")
        if (component)
            message("Configuring ${appName} app for ${component}")
        else()
            message("Configuring ${appName} app")
        endif()
        add_executable(${appName} ${app})

        if (BUILD_SHARED_LIBS)
            target_link_libraries(${appName} LINK_PUBLIC ${component} ${Chaste_LIBRARIES})
        else()
            target_link_libraries(${appName} LINK_PUBLIC ${component}  ${Chaste_LIBRARIES} ${Chaste_THIRD_PARTY_LIBRARIES} )
        endif()
        if(MSVC)
            set_target_properties(${appName} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT /IGNORE:4217 /IGNORE:4049")
        endif()
    endforeach(app)
  endmacro(CHASTE_DO_APPS_COMMON)

  macro(CHASTE_DO_APPS_PROJECT projectName)
    message("Configuring apps for project ${projectName}")
    CHASTE_DO_APPS_COMMON(project_${projectName})
  endmacro(CHASTE_DO_APPS_PROJECT)

  macro(CHASTE_DO_APPS_MAIN)
    message("Configuring main Chaste apps")
    CHASTE_DO_APPS_COMMON("")
  endmacro(CHASTE_DO_APPS_MAIN)

  macro(CHASTE_DO_TEST_COMMON component)
        # make tutorial directories
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/UserTutorials)
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/PaperTutorials)

        # Figure out include path for tests
        header_dirs("${CMAKE_CURRENT_SOURCE_DIR}" CHASTE_${component}_TEST_DIRS)
        include_directories("${CHASTE_${component}_TEST_DIRS}" "${CXXTEST_INCLUDES}")

        # Make test library if sources exist
        set(COMPONENT_LIBRARIES ${component})
        file(GLOB_RECURSE test_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
        if(test_sources)
            add_library(test${component} STATIC ${test_sources})
            set(COMPONENT_LIBRARIES ${COMPONENT_LIBRARIES} test${component})
        endif()


        #set(COMPONENT_LIBRARIES ${COMPONENT_LIBRARIES} ${CHASTE_DEPENDS_${component}})
        # Generate test suites

        if(MSVC)
            if(NOT HAS_OWN_LINKER_FLAGS)
                set(LINKER_FLAGS "/NODEFAULTLIB:LIBCMT")
            endif(NOT HAS_OWN_LINKER_FLAGS)

            #disable linker warnings 4217, 4049: locally-defined symbol imported in function ...
            set(LINKER_FLAGS "${LINKER_FLAGS} /IGNORE:4217 /IGNORE:4049")
            #message("Linker flags for project ${PROJECT_NAME} = ${LINKER_FLAGS}")
        endif(MSVC)


        foreach(type ${TestPackTypes})
            if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${type}TestPack.txt")
            file(STRINGS "${type}TestPack.txt" testpack)
            foreach(filename ${testpack})
                string(STRIP ${filename} filename)
                chaste_generate_test_name(${filename} "testTargetName")
                set(parallel OFF)
                set(exeTargetName ${testTargetName}Runner)
                if (${type} STREQUAL "Parallel")
                    set(testTargetName ${testTargetName}Parallel)
                    set(parallel ON)
                endif()
                
                if (NOT DEFINED ${testTargetName})
                    set(${testTargetName} ON)
                    chaste_add_test(${testTargetName} "${CMAKE_CURRENT_SOURCE_DIR}/${filename}")
                    if (BUILD_SHARED_LIBS)
                        target_link_libraries(${exeTargetName} LINK_PUBLIC ${COMPONENT_LIBRARIES})
                    else()
                        target_link_libraries(${exeTargetName} LINK_PUBLIC ${COMPONENT_LIBRARIES} ${Chaste_LIBRARIES} ${Chaste_THIRD_PARTY_LIBRARIES} )
                    endif()
                    set_target_properties(${exeTargetName} PROPERTIES LINK_FLAGS "${LINKER_FLAGS}")
                    set_property(TEST ${testTargetName} PROPERTY LABELS ${component} ${type})
                    if(NOT(${component} MATCHES "^project"))
                        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${testName}.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/${filename}"
                            DESTINATION test/${component} COMPONENT  ${component}_tests)
                    endif(NOT(${component} MATCHES "^project"))

                    # filename is a user tutorial
                    if(filename MATCHES "Test(.*)Tutorial.(hpp|py)") 
                        set(out_filename  ${CMAKE_CURRENT_BINARY_DIR}/html/UserTutorials/${CMAKE_MATCH_1})
                        add_custom_command(OUTPUT ${out_filename}
                                           COMMAND ${PYTHON_EXECUTABLE} ARGS ${Chaste_SOURCE_DIR}/python/utils/CreateTutorial.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} ${out_filename}
                                           DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
                                           COMMENT "Generating user tutorial ${out_filename}" VERBATIM)
                        add_custom_target(${CMAKE_MATCH_1} DEPENDS ${out_filename})
                        add_dependencies(tutorials ${CMAKE_MATCH_1})
                    endif()

                    # filename is a paper tutorial
                    if(filename MATCHES "Test(.*)LiteratePaper.(hpp|py)") 
                        set(out_filename  ${CMAKE_CURRENT_BINARY_DIR}/html/PaperTutorials/${CMAKE_MATCH_1})
                        add_custom_command(OUTPUT ${out_filename}
                                           COMMAND ${PYTHON_EXECUTABLE} ARGS ${Chaste_SOURCE_DIR}/python/utils/CreateTutorial.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} ${out_filename} 
                                           DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
                                           COMMENT "Generating paper tutorial ${out_filename}" VERBATIM)
                        add_custom_target(${CMAKE_MATCH_1} DEPENDS ${out_filename})
                        add_dependencies(tutorials ${CMAKE_MATCH_1})
                    endif()

                else()
                    get_property(myLabels TEST ${testTargetName} PROPERTY LABELS)
                    list(APPEND myLabels ${type})
                    set_property(TEST ${testTargetName} PROPERTY LABELS ${myLabels})
                endif()
            endforeach(filename ${testpack})
            endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${type}TestPack.txt")
        endforeach(type ${TestPackTypes})
  endmacro(CHASTE_DO_TEST_COMMON)

  macro(CHASTE_DO_TEST_COMPONENT component)
    message("Configuring tests for ${component}")
    CHASTE_DO_TEST_COMMON(${component})
  endmacro(CHASTE_DO_TEST_COMPONENT)

  macro(CHASTE_DO_TEST_PROJECT projectName)
    message("Configuring tests for project ${projectName}")
    CHASTE_DO_TEST_COMMON(project_${projectName})
  endmacro(CHASTE_DO_TEST_PROJECT)
