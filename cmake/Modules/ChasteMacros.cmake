# ChasteMacros.cmake
#
# A collection of macros useful for the Chaste project
#
# header_dirs(base_dir return_list) - recursivly finds Chaste header locations
# cellml_dirs(base_dir return_list) - recursivly finds Chaste cellml locations
# chaste_do_cellml(output_sources cellml_file dynamic) - convert cellml file to source files
# chaste_add_test(_testTargetName) - add a new test
# chaste_generate_test_name(test outTestName) - generate test name from source file
# chaste_do_common(component) - process component, generating libraries, tests and apps
# chaste_do_component(component) - wraps chaste_do_common for main Chaste components
# chaste_do_project(component) - wraps chaste_do_common for Chaste projects
# chaste_do_apps_common(app) - process app folder
# chaste_do_apps_main() - wraps chaste_do_apps_common for main Chaste apps
# chaste_do_apps_project(projectName) - wraps chaste_do_apps_common for project apps
# chaste_do_test_common(component) - process test folder, generating test executables and targets
# chaste_do_test_component(component) - wraps chaste_do_test_common for main Chaste components
# chaste_do_test_project(projectName) - wraps chaste_do_test_common for Chaste projects

##########################################################
# header_dirs
# 
# A macro to recursively find Chaste header locations
##########################################################
macro(HEADER_DIRS base_dir return_list)
    set(new_list "")
    set(dir_list "")

    file(GLOB_RECURSE new_list ${base_dir}/*.hpp ${base_dir}/*.h)
    foreach(file_path ${new_list})
        get_filename_component(dir_path ${file_path} PATH)
        set(dir_list ${dir_list} ${dir_path})
    endforeach()

    list(REMOVE_DUPLICATES dir_list)

    set(${return_list} ${dir_list})
endmacro()

##########################################################
# cellml_dirs
# 
# A macro to recursively find Chaste cellml locations
##########################################################
macro(CELLML_DIRS base_dir return_list)
    set(new_list "")
    set(dir_list "")

    file(GLOB_RECURSE new_list ${base_dir}/*.cellml)
    foreach(file_path ${new_list})
        get_filename_component(dir_path ${file_path} PATH)
        set(dir_list ${dir_list} ${dir_path})
    endforeach()

    list(REMOVE_DUPLICATES dir_list)

    set(${return_list} ${dir_list})
endmacro()

##########################################################
# chaste_do_cellml 
# 
# convert cellml file to source files
##########################################################
macro(Chaste_DO_CELLML output_sources cellml_file dynamic)
    get_filename_component(cellml_file_name ${cellml_file} NAME_WE)
    get_filename_component(cellml_dir ${cellml_file} PATH)
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
    
    #set depends on everything in python/pycml/* except for *.pyc and protocol.py
    file(GLOB PyCML_SOURCES 
        ${Chaste_SOURCE_DIR}/python/pycml/* )
    file(GLOB PyCML_NOT_SOURCES 
        ${Chaste_SOURCE_DIR}/python/pycml/*.pyc )
    list(REMOVE_ITEM PyCML_SOURCES ${PYCML_NOT_SOURCES} ${Chaste_SOURCE_DIR}/python/pycml/protocol.py)

    set(depends ${depends} ${PyCML_SOURCES})

    if(EXISTS ${cellml_dir}/${cellml_file_name}-conf.xml)
        set(depends ${depends} ${cellml_dir}/${cellml_file_name}-conf.xml)
        set(pycml_args ${pycml_args} "--conf=${cellml_dir}/${cellml_file_name}-conf.xml")
    endif()
    execute_process(COMMAND "${PYTHON_EXECUTABLE}" ${Chaste_PYTHON_DIR}/ConvertCellModel.py ${pycml_args} ${Chaste_PYCML_EXTRA_ARGS} --show-outputs ${cellml_file}   
        OUTPUT_VARIABLE ConvertCellModelDepends
        OUTPUT_STRIP_TRAILING_WHITESPACE
        )

    string(REGEX MATCHALL "[^\n]*\\.hpp" output_files_hpp "${ConvertCellModelDepends}")
    string(REGEX MATCHALL "[^\n]*\\.cpp" output_files_cpp "${ConvertCellModelDepends}")

    if (NOT Chaste_VERBOSE)
        set(pycml_args ${pycml_args} "--quiet")
    endif()

    add_custom_command(OUTPUT ${output_files_hpp} ${output_files_cpp} 
        COMMAND "${PYTHON_EXECUTABLE}" ${Chaste_PYTHON_DIR}/ConvertCellModel.py ${pycml_args} ${Chaste_PYCML_EXTRA_ARGS} ${cellml_file}
        DEPENDS ${depends}
        COMMENT "Processing CellML file ${cellml_file_rel}" 
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        VERBATIM
        )

    list(APPEND ${output_sources} ${output_files_cpp} ${output_files_hpp})
endmacro()

##########################################################
# chaste_add_test
# 
# Chaste Testing Macro. The predefined cxxtest_add_test is 
# not suitable because of little control over the test's 
# working directory
##########################################################
macro(Chaste_ADD_TEST _testTargetName _filename)
    string(REGEX MATCH "^.*Parallel$" foundParallel ${_testTargetName})
    if (foundParallel)
        set(parallel ON)
        string(REGEX REPLACE "(Parallel)$" "" _testname "${_testTargetName}")
    else()
        set(parallel OFF)
        set(_testname ${_testTargetName})
    endif()

    if (${_filename} MATCHES ".py$")
        set(python ON)
    else()
        set(python OFF)
    endif()


    if (python)
        set(test_exe ${Chaste_BINARY_DIR}/python/infra/TestPythonCode.py ${_filename})
    else()
        set(_exeTargetName ${_testname})

        if (NOT TARGET ${exeTargetName})
            set(_test_real_output_filename "${CMAKE_CURRENT_BINARY_DIR}/${_testname}.cpp")
            add_custom_command(
                OUTPUT "${_test_real_output_filename}"
                DEPENDS ${_filename} ${ARGN}
                COMMAND ${PYTHON_EXECUTABLE} ${CXXTEST_PYTHON_TESTGEN_EXECUTABLE} --error-printer -o "${_test_real_output_filename}" ${_filename} ${ARGN}
                )

            set_source_files_properties("${_test_real_output_filename}" PROPERTIES GENERATED true)

            add_executable(${exeTargetName} "${_test_real_output_filename}" ${_filename} ${ARGN})
        endif()

        set(test_exe $<TARGET_FILE:${exeTargetName}>)
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
        if (python)
            set(test_command ${test_exe})
            set(test_args "--num-procs ${num_cpus}")
        else()
            set(test_command ${MPIEXEC})
            set(test_args "${MPIEXEC_NUMPROC_FLAG} ${num_cpus} ${MPIEXEC_PREFLAGS}  ${test_exe} ${MPIEXEC_POSTFLAGS}")
        endif()
    else()
        set(num_cpus 1)
        set(test_command ${test_exe})
        set(test_args "")
    endif()


    if (Chaste_MEMORY_TESTING AND NOT python)
        set(test_command ${VALGRIND_COMMAND})
        set(test_args "--tool=memcheck --log-file=${Chaste_MEMORY_TESTING_OUTPUT_DIR}/${_testname}_valgrind.out") 
        set(test_args "${test_args} --track-fds=yes --leak-check=yes --num-callers=50 ${Chaste_MEMORY_TESTING_SUPPS}")
        set(test_args "${test_args} --gen-suppressions=all $<TARGET_FILE:${exeTargetName}> -malloc_debug -malloc_dump -memory_info")
        set(num_cpus 1)
    elseif (Chaste_PROFILE_GPROF OR Chaste_PROFILE_GPERFTOOLS)
        if (python)
            set(test_args "${test_args} --profile")
        elseif (Chaste_PROFILE_GPERFTOOLS)
            set(profile_file ${Chaste_PROFILE_OUTPUT_DIR}/${_testname}.prof)
            set(post_command ${GPERFTOOLS_PPROF_EXE})
            set(post_args "--svg --nodefraction=0.0001 --edgefraction=0.0001 $<TARGET_FILE:${exeTargetName}> ${profile_file}")
            set(output_file ${Chaste_PROFILE_OUTPUT_DIR}/${_testname}.svg)
            set(env_var CPUPROFILE)
            set(env_var_value ${profile_file})
        else()
            set(output_file ${Chaste_PROFILE_OUTPUT_DIR}/${_testname}.gmon)
            set(post_command ${GPROF_EXECUTABLE})
            set(post_args $<TARGET_FILE:${exeTargetName}>)
        endif()
    endif()

    if (post_command)
        add_test(NAME ${_testTargetName} WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/" 
            COMMAND ${CMAKE_COMMAND}
            -Denv_var=${env_var}
            -Denv_var_value=${env_var_value}
            -Dtest_cmd=${test_command}
            -Dtest_args:string=${test_args}
            -Dpost_cmd=${post_command}
            -Dpost_args:string=${post_args}
            -Doutput_file=${output_file}
            -P ${Chaste_BINARY_DIR}/cmake/Modules/ChasteRunTestAndPostProcess.cmake)
    else()
        separate_arguments(test_args)
        add_test(NAME ${_testTargetName} WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/" 
            COMMAND ${test_command} ${test_args})
    endif()
    set_property(TEST ${_testTargetName} PROPERTY PROCESSORS ${num_cpus})
    if (python)
        set_property(TEST ${_testTargetName} PROPERTY
            ENVIRONMENT "PYTHONPATH=$ENV{PYTHONPATH}:${Chaste_BINARY_DIR}/python/pycml"
            )
    endif()


endmacro(Chaste_ADD_TEST)

##########################################################
# chaste_generate_test_name
# 
# returns the test name (in outTestName) from the source
# hpp file (in test)
##########################################################
macro(Chaste_GENERATE_TEST_NAME test outTestName)
    get_filename_component(${outTestName} ${test} NAME_WE)
#####################
# This functionality appends the path to each test: removed in #2906
#####################
#    string(REGEX REPLACE "([a-zA-Z0-9_/]+)[.](hpp|py)" "\\1" testName "${test}")
#    string(REPLACE "/" ";" testPath "${testName}")
#    list(LENGTH testPath pathLength)
#    if(${pathLength} EQUAL 1)
#        set(testName ${testPath})
#        set(testPath "")
#        set(${outTestName} ${testName})
#    else()
#        math(EXPR index "${pathLength} - 1")
#        list(GET testPath ${index} testName)
#        list(REMOVE_AT testPath ${index})
#        string(REPLACE ";" "_" _testPath_ "${testPath}")
#        string(REPLACE ";" "/" testPath "${testPath}")
#        set(${outTestName} "${testName}_${_testPath_}_")
#    endif()
endmacro(Chaste_GENERATE_TEST_NAME test outTestName)

##########################################################
# chaste_do_common
# 
# processes a component directory, generating libraries, 
# tests and apps according to the standard Chaste directory
# layout
##########################################################
macro(Chaste_DO_COMMON component)

    add_definitions(-DCOMPONENT_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")
    if (NOT TARGET ${component})
        add_custom_target(${component})
    endif()

    # check if include dirs exists yet, and generate it if not
    if (NOT Chaste_${component}_INCLUDE_DIRS)
        set(Chaste_${component}_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
        header_dirs(${Chaste_${component}_SOURCE_DIR} Chaste_${component}_SOURCE_INCLUDE_DIRS)
        cellml_dirs(${Chaste_${component}_SOURCE_DIR} Chaste_${component}_CELLML_DIRS)

        # generate include dirs
        set(Chaste_${component}_INCLUDE_DIRS ${Chaste_${component}_SOURCE_INCLUDE_DIRS})
        foreach(dir ${Chaste_${component}_CELLML_DIRS})
            file(RELATIVE_PATH rel_dir ${CMAKE_CURRENT_SOURCE_DIR} ${dir})
            list(APPEND Chaste_${component}_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/${rel_dir})
        endforeach()
    endif()

    # Find source files
    file(GLOB_RECURSE Chaste_${component}_SOURCES 
        RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} 
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)

    # Generate additional source files from cellml
    set(Chaste_${component}_SOURCES_CELLML)
    foreach(cellml_dir ${Chaste_${component}_CELLML_DIRS})
        file(RELATIVE_PATH cellml_rel_dir ${CMAKE_CURRENT_SOURCE_DIR} ${cellml_dir})
        set(cellml_output_dir ${CMAKE_CURRENT_BINARY_DIR}/${cellml_rel_dir})
        file(MAKE_DIRECTORY ${cellml_output_dir})
        file(GLOB cellml_files ${cellml_dir}/*.cellml)
        foreach(cellml_file ${cellml_files})
            chaste_do_cellml(Chaste_${component}_SOURCES ${cellml_file} OFF  "--output-dir" ${cellml_output_dir})
        endforeach()
    endforeach()

    # Add include directories
    if (Chaste_THIRD_PARTY_INCLUDE_DIRS)
        include_directories(SYSTEM "${Chaste_THIRD_PARTY_INCLUDE_DIRS}")
    endif()
    if (Chaste_${component}_INCLUDE_DIRS)
        include_directories("${Chaste_${component}_INCLUDE_DIRS}")
    endif()
    if (Chaste_INCLUDE_DIRS)
        include_directories("${Chaste_INCLUDE_DIRS}")
    endif()

    # Make component library, if component contains any source files
    if (NOT Chaste_${component}_SOURCES STREQUAL "")
        add_library(chaste_${component} ${Chaste_${component}_SOURCES} ${ARGN})
        if (BUILD_SHARED_LIBS)
            target_link_libraries(chaste_${component} LINK_PUBLIC ${Chaste_LIBRARIES})
            set(static_extension "a")
            set(keyword "")
            foreach(library ${Chaste_THIRD_PARTY_LIBRARIES})
                if (library STREQUAL debug OR library STREQUAL optimized OR library STREQUAL general)
                    set(keyword ${library})
                else()
                    if (library MATCHES ".*\\.${static_extension}")
                        target_link_libraries(chaste_${component} LINK_PRIVATE ${keyword} ${library})
                    else()
                        target_link_libraries(chaste_${component} LINK_PUBLIC ${keyword} ${library})
                    endif()
                    set(keyword "")
                endif()
            endforeach()
        else()
            target_link_libraries(chaste_${component} LINK_PUBLIC ${Chaste_THIRD_PARTY_LIBRARIES})
        endif()


        if(NOT(${component} MATCHES "^project"))
            # install component library
            install(TARGETS chaste_${component}
                EXPORT chaste-targets
                DESTINATION lib/chaste
                COMPONENT ${component}_libraries)

            # install component headers
            install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/"
                DESTINATION include/chaste/${component}
                COMPONENT ${component}_headers
                FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
                )

            # install generated headers
            install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/src/"
                DESTINATION include/chaste/${component}
                COMPONENT ${component}_headers
                FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
                )
        endif()
    endif()

    if (Chaste_ENABLE_TESTING) 
        if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")
            set(Chaste_ENABLE_${component}_TESTING ON CACHE BOOL "Generate the test infrastructure for ${component} ")
        else()
            message(WARNING "No CMakeLists.txt file found in test directory ${CMAKE_CURRENT_SOURCE_DIR}/test. Tests for ${component} will not be built")
            set(Chaste_ENABLE_${component}_TESTING OFF CACHE BOOL "Generate the test infrastructure for ${component} ")
        endif()

        # Do testing if requested
        if(Chaste_ENABLE_${component}_TESTING)
            add_subdirectory(test)
        endif()
    endif(Chaste_ENABLE_TESTING)

    # Build applications if present
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/apps")
        if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/apps/CMakeLists.txt")
            set(Chaste_ENABLE_${component}_APPS ON CACHE BOOL "Generate the applications infrastructure for ${component} ")
        else()
            message(WARNING "No CMakeLists.txt file found in test directory ${CMAKE_CURRENT_SOURCE_DIR}/apps. Applications for ${component} will not be built")
            set(Chaste_ENABLE_${component}_APPS OFF CACHE BOOL "Generate the applications infrastructure for ${component} ")
        endif()

        # Do apps if requested
        if(Chaste_ENABLE_${component}_APPS)
            add_subdirectory(apps)
        endif()
    endif()
endmacro(Chaste_DO_COMMON)

##########################################################
# chaste_do_component and chaste_do_project
# 
# wrapper for chaste_do_common for main Chaste components
# and other project directories
##########################################################
macro(Chaste_DO_COMPONENT component)
    message("Configuring component ${component}")
    Chaste_DO_COMMON(${component} ${ARGN})
endmacro(Chaste_DO_COMPONENT)

macro(Chaste_DO_PROJECT projectName)
    if (Chaste_ENABLE_project_${projectName})
        message("Configuring project ${projectName}")
        Chaste_DO_COMMON(project_${projectName})
    endif()
endmacro(Chaste_DO_PROJECT)

##########################################################
# chaste_do_apps_common
# 
# process the apps folder
##########################################################
macro(Chaste_DO_APPS_COMMON component)
    include_directories(SYSTEM "${Chaste_THIRD_PARTY_INCLUDE_DIRS}" "${Chaste_INCLUDE_DIRS}")
    include_directories(SYSTEM "${CXXTEST_INCLUDES}")
    file(GLOB_RECURSE Chaste_${component}_APPS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} src/*.cpp)
    foreach(app ${Chaste_${component}_APPS})
        string(REGEX REPLACE ".*/([a-zA-Z0-9_]+)[.]cpp" "\\1" appName "${app}")
        if (${component} MATCHES "project_")
            message("Configuring ${appName} app for ${component}")
            set(component_library chaste_${component})
        else()
            message("Configuring ${appName} app")
            set(component_library )
        endif()
        add_executable(${appName} ${app})
        #set_target_properties(${appName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/src)

        if (NOT ${component} STREQUAL "")
            add_dependencies(${component} ${appName})
        endif()

        if (BUILD_SHARED_LIBS)
            target_link_libraries(${appName} LINK_PUBLIC ${component_library} ${Chaste_LIBRARIES})
        else()
            target_link_libraries(${appName} LINK_PUBLIC ${component_library} ${Chaste_LIBRARIES} ${Chaste_THIRD_PARTY_LIBRARIES} )
        endif()
        if(MSVC)
            set_target_properties(${appName} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:LIBCMT /IGNORE:4217 /IGNORE:4049")
        endif()
    endforeach(app)
    if (Chaste_ENABLE_TESTING AND TEXTTEST_FOUND AND EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/texttest)
        configure_file(texttest/chaste/wrapper.cmake.in texttest/chaste/wrapper)
        file(GLOB test_directories texttest/*)
        foreach(tests_dir ${test_directories})
            file(RELATIVE_PATH acceptance_test ${CMAKE_CURRENT_SOURCE_DIR}/texttest ${tests_dir})
            set(output_dir ${CMAKE_CURRENT_BINARY_DIR}/texttest)
            set(texttest_report_dir ${output_dir}/texttest_reports/${acceptance_test})
            set(texttest_output_dir ${output_dir}/texttest_output/${acceptance_test})

            file(WRITE ${tests_dir}/run_acceptance.cmake
                "
                set(ENV{CHASTE_TEST_OUTPUT} ${output_dir})
                file(REMOVE_RECURSE ${texttest_report_dir})
                file(REMOVE_RECURSE ${texttest_output_dir})
                file(MAKE_DIRECTORY ${texttest_report_dir})
                file(MAKE_DIRECTORY ${texttest_output_dir})
                execute_process(COMMAND  ${PYTHON_EXECUTABLE} ${TEXTTEST_PY} -d ${tests_dir} -b default -c ${CMAKE_BINARY_DIR} 
                    RESULT_VARIABLE result)
                execute_process(COMMAND  ${PYTHON_EXECUTABLE} ${TEXTTEST_PY} -d ${tests_dir} -b default -c ${CMAKE_BINARY_DIR} -coll web)
                if (result)
                    message(SEND_ERROR \"Error running acceptance test\")
                endif()
                "
                )

            add_test(NAME acceptance_${acceptance_test} 
                COMMAND ${CMAKE_COMMAND} -P ${tests_dir}/run_acceptance.cmake
                )
            if (${acceptance_test} STREQUAL "weekly")
                set_property(TEST acceptance_${acceptance_test} PROPERTY LABELS ${component} Profile)
            else()
                set_property(TEST acceptance_${acceptance_test} PROPERTY LABELS ${component} Nightly)
            endif()

        endforeach()
    endif()
endmacro(Chaste_DO_APPS_COMMON)

##########################################################
# chaste_do_apps_main and chaste_do_apps_project
# 
# these wrap chaste_do_apps_common for both the main
# Chaste apps directory and external project apps dirs
##########################################################
macro(Chaste_DO_APPS_PROJECT projectName)
    message("Configuring apps for project ${projectName}")
    Chaste_DO_APPS_COMMON(project_${projectName})
endmacro(Chaste_DO_APPS_PROJECT)

macro(Chaste_DO_APPS_MAIN)
    message("Configuring main Chaste apps")
    Chaste_DO_APPS_COMMON("")
endmacro(Chaste_DO_APPS_MAIN)

##########################################################
# chaste_do_test_common
# 
# process the tests directory, generating tests for all
# enabled test packs.
##########################################################
macro(Chaste_DO_TEST_COMMON component)

    # Get the git revision (for tutorial tests)
    find_package(Git QUIET)

    # make tutorial directories
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tutorials)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tutorials/UserTutorials)
    file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tutorials/PaperTutorials)

    # Figure out include path for tests
    header_dirs("${CMAKE_CURRENT_SOURCE_DIR}" Chaste_${component}_TEST_DIRS)
    include_directories(${Chaste_${component}_TEST_DIRS})
    include_directories(SYSTEM "${CXXTEST_INCLUDES}")

    # Make test library if sources exist
    if (TARGET chaste_${component})
        set(COMPONENT_LIBRARIES chaste_${component})
    else()
        set(COMPONENT_LIBRARIES ${Chaste_LIBRARIES})
    endif()
    file(GLOB_RECURSE test_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
    if(test_sources)
        add_library(test${component} STATIC ${test_sources})
        set(COMPONENT_LIBRARIES ${COMPONENT_LIBRARIES} test${component})
    endif()


    #set(COMPONENT_LIBRARIES ${COMPONENT_LIBRARIES} ${Chaste_DEPENDS_${component}})
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

            # remove python tests from windows builds
            if (WIN32 OR CYGWIN) 
                set(testpack_new "")
                foreach(filename ${testpack})
                    if (NOT filename MATCHES ".py$")
                        list(APPEND testpack_new ${filename}) 
                    endif()
                endforeach()
                set(testpack ${testpack_new})
            endif(WIN32 OR CYGWIN)

            foreach(filename ${testpack})
                string(STRIP ${filename} filename)
                chaste_generate_test_name(${filename} "testTargetName")
                set(old_testTargetName ${testTargetName})
                set(parallel OFF)
                set(exeTargetName ${testTargetName})
                if (${type} STREQUAL "Parallel")
                    set(testTargetName ${testTargetName}Parallel)
                    set(parallel ON)
                endif()

                if (NOT DEFINED ${testTargetName})
                    set(${testTargetName} ON)
                    chaste_add_test(${testTargetName} "${CMAKE_CURRENT_SOURCE_DIR}/${filename}")

                    if (filename MATCHES ".hpp$")
                        if (BUILD_SHARED_LIBS)
                            target_link_libraries(${exeTargetName} LINK_PUBLIC ${COMPONENT_LIBRARIES})
                        else()
                            target_link_libraries(${exeTargetName} LINK_PUBLIC ${COMPONENT_LIBRARIES} ${Chaste_LIBRARIES} ${Chaste_THIRD_PARTY_LIBRARIES} )
                        endif()
                        set_target_properties(${exeTargetName} PROPERTIES LINK_FLAGS "${LINKER_FLAGS}")
                    endif()


                    set_property(TEST ${testTargetName} PROPERTY LABELS ${type}_${component})

                    if (Chaste_INSTALL_TESTS AND NOT(${component} MATCHES "^project")) 
                        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${old_testTargetName}.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/${filename}"
                            DESTINATION lib/chaste/tests/${component} COMPONENT  ${component}_tests)
                    endif()

                    # filename is a user tutorial
                    if(filename MATCHES "Test(.*)Tutorial.(hpp|py)")
                        # Get the git revision of last time this file was changed
                        if(DEFINED GIT_EXECUTABLE AND EXISTS "${Chaste_SOURCE_DIR}/.git")
                            execute_process(COMMAND git log -1 --format=%H --follow ${filename}
                                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                                    OUTPUT_VARIABLE Chaste_revision
                                    OUTPUT_STRIP_TRAILING_WHITESPACE)
                        endif()
                        if(DEFINED Chaste_revision)
                            set(revision_string "-r ${Chaste_revision}")
                        else()
                            set(revision_string "")
                        endif()
                        set(out_filename  ${CMAKE_BINARY_DIR}/tutorials/UserTutorials/${CMAKE_MATCH_1})
                        add_custom_command(OUTPUT ${out_filename}
                            COMMAND ${PYTHON_EXECUTABLE} ARGS ${Chaste_BINARY_DIR}/python/utils/CreateTutorial.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} ${out_filename} ${revision_string}
                            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
                            COMMENT "Generating user tutorial ${out_filename}" VERBATIM)
                        add_custom_target(${CMAKE_MATCH_1} DEPENDS ${out_filename})
                        add_dependencies(tutorials ${CMAKE_MATCH_1})
                    endif()

                    # filename is a paper tutorial
                    if(filename MATCHES "Test(.*)LiteratePaper.(hpp|py)") 
                        set(out_filename  ${CMAKE_BINARY_DIR}/tutorials/PaperTutorials/${CMAKE_MATCH_1})
                        add_custom_command(OUTPUT ${out_filename}
                            COMMAND ${PYTHON_EXECUTABLE} ARGS ${Chaste_BINARY_DIR}/python/utils/CreateTutorial.py ${CMAKE_CURRENT_SOURCE_DIR}/${filename} ${out_filename} 
                            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${filename}
                            COMMENT "Generating paper tutorial ${out_filename}" VERBATIM)
                        add_custom_target(${CMAKE_MATCH_1} DEPENDS ${out_filename})
                        add_dependencies(tutorials ${CMAKE_MATCH_1})
                    endif()
                endif()

                # add dependencies to component and type targets. Do not include the python component or tests in Python files
                if ((NOT ${component} STREQUAL python) AND (NOT (${filename} MATCHES ".py$")))
                    add_dependencies(${component} ${exeTargetName})
                    add_dependencies(${type} ${exeTargetName})
                endif()
            endforeach(filename ${testpack})
        endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${type}TestPack.txt")
    endforeach(type ${TestPackTypes})
endmacro(Chaste_DO_TEST_COMMON)

##########################################################
# chaste_do_test_component and chaste_do_test_project
# 
# wrapper for chaste_do_test_common for main Chaste tests
# and external project tests
##########################################################
macro(Chaste_DO_TEST_COMPONENT component)
    message("Configuring tests for ${component}")
    Chaste_DO_TEST_COMMON(${component})
endmacro(Chaste_DO_TEST_COMPONENT)

macro(Chaste_DO_TEST_PROJECT projectName)
    message("Configuring tests for project ${projectName}")
    Chaste_DO_TEST_COMMON(project_${projectName})
endmacro(Chaste_DO_TEST_PROJECT)
