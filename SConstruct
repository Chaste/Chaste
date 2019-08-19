"""Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

# Controlling SCons build script for Chaste.

# This script is executed within the root Chaste source directory.
# We need at least Python 2.3.
EnsurePythonVersion(2,3)

# We're also no longer compatible with SCons 0.96
EnsureSConsVersion(0,97)

# Avoid deprecation warnings by changing behaviour for new SCons versions
scons_version_two = True
try:
    import SCons
    if int(SCons.__version__[0]) < 2:
        scons_version_two = False
except:
    scons_version_two = False


Help("""
  Type: 'scons -c .' to remove all the compiled files (clean build),
        'scons' to do a default build,
        'scons <Path to TestWhatever.hpp>' to run a single test,
        'scons <component>' to build and test a single component.
  
  For other options, such as profiling, optimised builds and 
  memory testing please refer to:
  
  https://chaste.cs.ox.ac.uk/trac/wiki/SconsArchive/UserBuildGuide

""")

import sys
import os
import glob
import socket
import time


import SCons

sys.path[0:0] = ['python/infra', 'python/hostconfig', 'python']
import BuildTypes
import SConsTools
Export('SConsTools')

import hostconfig

# If building a loadable module at run-time
dyn_libs_only = int(ARGUMENTS.get('dyn_libs_only', 0))
Export('dyn_libs_only')
if dyn_libs_only:
    # Set some other options
    ARGUMENTS['test_summary'] = 0
    ARGUMENTS['do_inf_tests'] = 0
    # Note what folder is being built
    dyn_folder = os.path.join(Dir('#').abspath, COMMAND_LINE_TARGETS[0])
    Export('dyn_folder')

# Turn on some build-script debugging?
debug = int(ARGUMENTS.get('debug', 0))
Export('debug')

# The type of build to perform (see python/BuildTypes.py for options)
build_type = ARGUMENTS.get('build', ARGUMENTS.get('b', 'default'))
build = BuildTypes.GetBuildType(build_type)
build.SetRevision(ARGUMENTS.get('revision', ''))
build.debug = debug
build.quiet = GetOption('silent') or dyn_libs_only
Export('build')

# Whether to use static or shared libraries
static_libs = int(ARGUMENTS.get('static', 0))
if build.needs_static:
    static_libs = 1
Export('static_libs')

# Whether to build Chaste libraries, or link tests against object files directly
use_chaste_libs = int(ARGUMENTS.get('chaste_libs',  ARGUMENTS.get('cl', 0)))
Export('use_chaste_libs')

# Specify test_summary=0 to scons to *NOT* generate a summary html page
test_summary = int(ARGUMENTS.get('test_summary', 1))
Export('test_summary')

# Used by the automated build system
# Infrastructure tests default to only checking orphaned tests & duplicate file names, unless do_inf_tests=1 is given
run_infrastructure_tests = int(ARGUMENTS.get('do_inf_tests', getattr(hostconfig.conf, 'do_inf_tests', 2)))
check_failing_tests = int(ARGUMENTS.get('check_failing_tests', 0))

# Specifying extra run-time flags
run_time_flags = ARGUMENTS.get('run_time_flags', '')

# Specify all_tests=1 to select all tests for running (useful with
# compile_only=1)
all_tests = int(ARGUMENTS.get('all_tests', ARGUMENTS.get('at', 0)))
Export('all_tests')

# Specify compile_only=1 to not run any tests
compile_only = int(ARGUMENTS.get('compile_only', ARGUMENTS.get('co', 0)))
Export('compile_only')

# Special offline mode which doesn't run tests, but writes a script to run them separately
offline_mode = int(ARGUMENTS.get('offline_mode', ARGUMENTS.get('om', 0)))
Export('offline_mode')
if offline_mode:
    compile_only = 1
    Export('compile_only')
    offline_test_runners = []
    Export('offline_test_runners')

# To run a specific test suite only, give the path (relative to the Chaste root) as the test_suite=<path> argument.
# This will force the test suite to be run even if the source is unchanged.
# To run multiple tests in this way, separate the paths by commas.
# Also allow test .hpp files to be given as normal targets.
requested_tests = filter(None, ARGUMENTS.get('test_suite', ARGUMENTS.get('ts', '')).split(','))
for target in COMMAND_LINE_TARGETS[:]:
    if target.endswith('.hpp') or target.endswith('.py'):
        requested_tests.append(target)
        COMMAND_LINE_TARGETS.remove(target)
        BUILD_TARGETS.remove(target)
for test_idx, test_path in enumerate(requested_tests):
    parts = test_path.split(os.path.sep)
    if len(parts) < 3:
        raise ValueError('Path to test suite "%s" is too short' % test_path)
    for path_idx in range(-2, -len(parts), -1):
        if parts[path_idx] == 'test':
            requested_tests[test_idx] = (parts[path_idx-1], os.path.sep.join(parts[path_idx+1:]))
            break
    else:
        raise ValueError('Test suite "%s" is not in a test folder' % test_path)
    #print requested_tests
Export('requested_tests')

# Force re-running of all (selected) tests even if the source is unchanged.
force_test_runs = int(ARGUMENTS.get('force_test_runs', 0))
Export('force_test_runs')

# Don't update the provenance information (Version.cpp file).
update_provenance = int(ARGUMENTS.get('update_provenance', ARGUMENTS.get('up', 1)))

# Whether to kill the tests if they run too long (limit in seconds, 0 means don't kill)
test_time_limit = int(ARGUMENTS.get('test_time_limit', 0))

# Whether to build executables, or just tests
build_exes = int(ARGUMENTS.get('exe', 0))
Export('build_exes')
if build_exes:
    assert use_chaste_libs, "Cannot build executables unless building Chaste libraries"

# Whether to alter the path for the texttest acceptance suite.  This is used when
# exe=1, the target is "apps" and the "compile_only" flag is not set.
texttest_suite = ARGUMENTS.get('acceptance_suite', 'chaste')
Export('texttest_suite')

# Experimental support for installing Chaste as a normal collection of
# libraries and headers.
install_prefix = ARGUMENTS.get('install_prefix', '/usr/local')
Export('install_prefix')
install_files = 'install' in BUILD_TARGETS
if install_files:
    assert use_chaste_libs, "Cannot install unless building Chaste libraries"

# To run tests of only a single component, specify it with the
# test_component=<component> argument (deprecated).
test_component = ARGUMENTS.get('test_component', '')
Export('test_component')

# Special mode that doesn't store test results, to avoid the buffering that occurs when we do so
build.no_store_results = int(ARGUMENTS.get('no_store_results', 0))


# If building static libraries, get rid of any old shared libraries,
# in order to stop the automatic dependency algorithm getting confused.
if use_chaste_libs and static_libs:
    for lib in glob.glob('lib/lib*.so'):
        Execute(Delete(lib))


# Use a single file to store signatures.
# Forwards-compatible with SCons 0.97, and nicer for svn ignore.
if not dyn_libs_only:
    if sys.platform == 'cygwin':
        SConsignFile('.sconsign-cygwin')
    else:
        SConsignFile('.sconsign')
else:
    # Use a .sconsign file in the folder we're building to avoid conflicts.
    assert(len(COMMAND_LINE_TARGETS) == 1)
    SConsignFile(os.path.join(COMMAND_LINE_TARGETS[0], '.sconsign'))

# Chaste components (top level dirs).
# Direct dependencies of each component are given, and a topological sort of the
# full dependency graph done to figure out which components to link against in what order.
comp_deps = {'cell_based': ['pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'crypt': ['cell_based'],
             'notforrelease': ['heart'],
             'notforrelease_cell_based': ['crypt', 'cell_based'],
             'notforrelease_lung': ['lung', 'heart'],
             'lung': ['continuum_mechanics'],
             'heart': ['continuum_mechanics', 'pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'continuum_mechanics': ['pde', 'ode', 'mesh', 'linalg', 'io', 'global'],
             'pde': ['ode', 'mesh', 'linalg', 'io', 'global'],
             'mesh': ['linalg', 'global'],
             'linalg': ['global'],
             'ode': ['linalg', 'io', 'global'],
             'io': ['global'],
             'global': [],
             'python': []
             }
components = ['python', 'global', 'io', 'linalg', 'mesh', 'ode', 'pde', 'continuum_mechanics',
              'heart', 'cell_based', 'crypt', 'lung',
              'notforrelease', 'notforrelease_cell_based', 'notforrelease_lung']
# Ignore non-existent components, e.g. notforrelease wont appear in a release version
for comp in components[:]:
    if not os.path.isdir(comp):
        components.remove(comp)
        del comp_deps[comp]
Export('components', 'comp_deps')

# Virtual component aliasing all core components
comp_deps['core'] = Split('global io linalg mesh ode pde continuum_mechanics')
Alias('core', comp_deps['core'])

# Set extra paths to search for libraries and include files.
# Paths to PETSc, and any other external libraries, should be set here.
# The three variables exported are:
#   other_libs: names of libraries to link against.
#   other_libpaths: paths to search for libraries to link against.
#   other_includepaths: paths to search for header files.
# This is now done by the hostconfig subsystem.
hostconfig.Configure(build)
other_libs = hostconfig.libraries
other_libpaths = hostconfig.libpaths
other_includepaths = hostconfig.incpaths
if isinstance(build, BuildTypes.CovTool):
    build.UseCovTool(other_includepaths, other_libs)

Export("other_libpaths", "other_libs")


def RequestedProjects():
    """Return a list of projects explicitly mentioned on the command line."""
    projects = []
    for targ in COMMAND_LINE_TARGETS:
        if str(targ).startswith('projects'):
            projects.append(str(targ))
    for req_test in requested_tests:
        if req_test[0] in os.listdir('projects'):
            projects.append(os.path.join('projects', req_test[0]))
    return projects


# Set up the environment to use for building.
other_libpaths.append(os.path.abspath('lib'))
if os.environ.get('CHASTE_LOAD_ENV', ''):
    env = Environment(ENV = os.environ)
else:
    env = Environment(
        #tools = ['g++', 'gnulink', 'gas', 'ar', 'g77'],
        ENV={'PATH': '.:' + os.environ['PATH'],
             'PYTHONPATH': os.environ.get('PYTHONPATH', ''),
             'USER': os.environ['USER'],
             'INTEL_LICENSE_FILE': os.environ.get('INTEL_LICENSE_FILE', '.'),
             'CHASTE_TEST_OUTPUT':
                 os.environ.get('CHASTE_TEST_OUTPUT',
                                '/tmp/'+os.environ['USER']+'/testoutput/'),
             'CHASTE_DEBUG': str(debug),
             'CHASTE_LIBS': os.environ.get('CHASTE_LIBS', ''),
             'LD_LIBRARY_PATH': ':'.join(other_libpaths),
             'HOME': os.environ['HOME'],
            })
env.Append(BOPT = 'g_c++') # Needed for some versions of PETSc?
env.Replace(CXX = build.tools['mpicxx'])
env.Replace(AR = build.tools['ar'])
env.Replace(CXXFILESUFFIX = '.cpp')
env.Append(PYINCPATH=['#/python/pycml']) # Ensure Python tests can use PyCml easily
env['INSTALL_PREFIX'] = install_prefix
env['INSTALL_FILES'] = install_files

if int(ARGUMENTS.get('br', ARGUMENTS.get('brief', 0))):
    env.Replace(CXXCOMSTR = '$CXX -o $TARGET -c <flags etc. omitted> $SOURCES')
    env.Replace(SHCXXCOMSTR = '$SHCXX -o $TARGET -c <flags etc. omitted> $SOURCES')

# Any extra CCFLAGS and LINKFLAGS
extra_flags = build.CcFlags() + ' ' + hostconfig.CcFlags()
link_flags  = build.LinkFlags() + ' ' + hostconfig.LdFlags()
include_flag = ' ' + build.IncludeFlag() + ' '

env.Append(CCFLAGS = include_flag + include_flag.join(other_includepaths)
           + ' ' + extra_flags)
env.Append(LINKFLAGS = link_flags)
env.Append(CPPDEFINES = hostconfig.CppDefines() + ['TRILIBRARY', 'TETLIBRARY', 'ANSI_DECLARATORS'])
if debug:
    print "Core compiler flags used by Chaste:"
    print "    CCFLAGS =", env['CCFLAGS']
    print "    LDFLAGS =", env['LINKFLAGS']
    print "    Defines =", env['CPPDEFINES']

# Base search path for Chaste #includes, common to all components
cpppath = [Dir('.'), Dir('cxxtest')]
env.Replace(CPPPATH = cpppath)

# Some state needed by our build system
env['build'] = build
build.env = env
env['buildsig'] = build.GetSignature()
env['CHASTE_COMPONENTS'] = components + ['projects']
env['CHASTE_COMP_DEPS'] = comp_deps
env['CHASTE_LIBRARIES'] = {}
env['CHASTE_OBJECTS'] = {}
env['UPDATE_CHASTE_PROVENANCE'] = update_provenance
env['CHASTE_CPP_PATHS'] = {}
env['CHASTE_CPP_PATH'] = {}
env['CHASTE_PYINCPATHS'] = {}
env['CHASTE_PYINCPATH'] = {}

if not requested_tests:
    # Default is to build all components, but not user projects
    Default(components)


# Create Builders for generating test .cpp files, and running test executables
test = Builder(action='cxxtest/cxxtestgen.py --error-printer -o $TARGET $SOURCES')
import TestRunner
def TestDescription(target, source, env):
    return "Running '%s'" % (source[0])
test_action = Action(TestRunner.get_build_function(build, run_time_flags, test_time_limit),
                     TestDescription, varlist=['buildsig'])
env['BUILDERS']['Test'] = test
env['BUILDERS']['RunTest'] = Builder(action=test_action)
env['BUILDERS']['BuildTest'] = Builder(action=SConsTools.BuildTest)

# Test runner for Python tests needs to use a source scanner
# to find what Python files the test depends on.
env['BUILDERS']['PyRunTest'] = Builder(action=test_action,
                                       source_scanner=SConsTools.PyScanner())

# Faster builds of shared libraries
env['BUILDERS']['OriginalSharedLibrary'] = env['BUILDERS']['SharedLibrary']
env['BUILDERS']['SharedLibrary'] = SConsTools.FasterSharedLibrary

# Builder for generating C++ code from XML Schema files
SConsTools.CreateXsdBuilder(build, env, dyn_libs_only)

# Builder for generating C++ code from CellML files
SConsTools.CreatePyCmlBuilder(build, env)

# Find full path to valgrind, as parallel memory testing needs it to be
# given explicitly.
vg_path = env.WhereIs(build.tools['valgrind'])
if vg_path:
    build.tools['valgrind'] = vg_path
del vg_path

# Record key build info for the provenance system
SConsTools.RecordBuildInfo(env, build_type, static_libs, use_chaste_libs)

# Allow hostconfig scripts to modify the build environment
if hasattr(hostconfig.conf, 'ModifyEnv') and callable(hostconfig.conf.ModifyEnv):
    hostconfig.conf.ModifyEnv(env)

# Export the build environment to SConscript files
Export('env')

# Test log files to summarise
test_log_files = []

# 'Infrastructure' tests of the codebase layout etc.
def run_infra(test, out, run_time_flags=''):
    if run_time_flags:
        run_time_flags = ' "%s"' % run_time_flags
    os.system('python/infra/TestRunner.py python/infra/' + test + ' ' + str(out)
              + ' ' + build_type + ' --no-stdout' + run_time_flags)
if run_infrastructure_tests and not GetOption('clean'):
    if not os.path.exists(build.GetTestReportDir()):
        os.makedirs(build.GetTestReportDir())
    if not env.GetOption('silent'):
        print "Running infrastructure tests..."
    # Check for orphaned test files
    out = File(build.GetTestReportDir() + 'OrphanedTests.log')
    run_infra('CheckForOrphanedTests.py', out, ' '.join(RequestedProjects()))
    test_log_files.append(out)
    # Check for duplicate file names in multiple directories
    out = File(build.GetTestReportDir() + 'DuplicateFileNames.log')
    run_infra('CheckForDuplicateFileNames.py', out)
    test_log_files.append(out)
    if run_infrastructure_tests == 1:
        # These infrastructure tests don't run by default, just for core developer machines/people
        out = File(build.GetTestReportDir() + 'Copyrights.log')
        run_infra('CheckForCopyrights.py', out)
        test_log_files.append(out)
        ## Do not check for stale semaphores - it's only important on MPICH with Gnu Linux
        #out = File(build.GetTestReportDir() + 'Semaphores.log')
        #run_infra('CheckSemaphores.py', out)
        #test_log_files.append(out)
        # Check for stray schemas
        out = File(build.GetTestReportDir() + 'Schemas.log')
        run_infra('CheckSchemas.py', out)
        test_log_files.append(out)
if check_failing_tests:
    out = File(build.GetTestReportDir() + 'FailingTests.log')
    run_infra('CheckForFailingTests.py', out)
    test_log_files.append(out)

build_dir = build.build_dir
if not isinstance(build, BuildTypes.DoxygenCoverage):
    # Build each component.
    for toplevel_dir in components:
        bld_dir = os.path.join(toplevel_dir, 'build', build_dir)
        if not os.path.exists(bld_dir):
            os.mkdir(bld_dir)
        script = os.path.join(toplevel_dir, 'SConscript')
        if scons_version_two:
            test_log_files.append(SConscript(script, src_dir=toplevel_dir, variant_dir=bld_dir, duplicate=0))
        else:
            test_log_files.append(SConscript(script, src_dir=toplevel_dir, build_dir=bld_dir, duplicate=0))
    
    # Any user projects?
    project_path_list = set()
    for project in glob.glob('projects/*'):
        if not os.path.isdir(project):
            if debug:
                print "Found non-dir", project, "in projects folder"
            continue
        if not os.path.exists(os.path.join(project, 'SConscript')):
            print >>sys.stderr, "Unexpected folder", project, "in projects folder."
            continue
        real_path = os.path.realpath(project)
        if real_path in project_path_list:
            raise ValueError('Two project folders have the same real path (%s)! This would cause compilation failure; bailing out.' % real_path)
        else:
            project_path_list.add(real_path)
        if install_files and not (project in BUILD_TARGETS or project+os.sep in BUILD_TARGETS):
            # Only install projects if explicitly requested
            continue
        bld_dir = os.path.join(project, 'build', build_dir)
        if not os.path.exists(bld_dir):
            os.makedirs(bld_dir)
        script = os.path.join(project, 'SConscript')
        
        if scons_version_two:
            test_log_files.append(SConscript(script, src_dir=project, variant_dir=bld_dir, duplicate=0))
        else:
            test_log_files.append(SConscript(script, src_dir=project, build_dir=bld_dir, duplicate=0))
    del project_path_list
     
    # Calculate full library dependencies now we know what projects need
    # This also sets up the full include paths for each component & project
    SConsTools.DetermineLibraryDependencies(env, comp_deps)
    
    # Make sure test executables get built if compile_only=1
    env.Default(env.Alias('test_exes'))
    # Work around an error on some SCons versions, by ensuring the alias has a builder
    import SCons.Environment
    env.Alias('test_exes')[0].builder_set(SCons.Environment.AliasBuilder)


# Remove the contents of build.GetTestReportDir() on a clean build
test_output_files = glob.glob(build.GetTestReportDir() + '*')
Clean('.', test_output_files)
# Also remove the entire build.build_dir for each component, so we
# don't have stale tests, etc. still present
for toplevel_dir in components:
    Clean('.', os.path.join(toplevel_dir, 'build', build_dir))
# Also make sure we remove any libraries still hanging around, just in case
Clean('.', glob.glob('lib/*'))
Clean('.', glob.glob('linklib/*'))


# Test summary generation
if test_summary and not compile_only:
    # Copy the build env, since we change TargetSigs
    senv = SConsTools.CloneEnv(env)
    # Get the directory to put results & summary in
    output_dir = build.output_dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Remove old results. Note that this command gets run before anything is built.
    #for oldfile in os.listdir(output_dir):
    #    os.remove(os.path.join(output_dir, oldfile))
    # Record some additional info in an extra file there
    info_file = open(os.path.join(output_dir, 'info.log'), 'w')
    info_file.write('Targets:')
    for target in COMMAND_LINE_TARGETS:
        info_file.write(' ' + target)
    info_file.write('\nArguments: ' + str(ARGUMENTS) + '\n')
    info_file.close()
    # Add a summary generator to the list of things for scons to do
    show_output_folder = False
    if isinstance(build, BuildTypes.Coverage):
        # Remove old .gcda files before running more tests
        # First, find appropriate build directories
        build_dirs = glob.glob('*/build/' + build.build_dir)
        # Now find & remove .gcda files within there.
        # Also remove .log files so tests are re-run
        for build_dir in build_dirs:
            for dirpath, dirnames, filenames in os.walk(build_dir):
                for filename in filenames:
                    if filename[-5:] == '.gcda' or filename[-4:] == '.log':
                        os.remove(os.path.join(dirpath, filename))
        # For a Coverage build, run gcov & summarise instead
        summary_action = ('python python/infra/DisplayCoverage.py ' + output_dir + ' ' + build_type
                          + ' ' + ' '.join(RequestedProjects()))
    elif isinstance(build, BuildTypes.DoxygenCoverage):
        # Run Doxygen and parse the output
        doxy_conf = ['cat Doxyfile', 'echo "PROJECT_NUMBER=Build::r%s"' % build._revision]
        # Include projects?
        project_inputs = RequestedProjects()
        if project_inputs:
            doxy_conf.append('echo "INPUT += %s"' % (' '.join(map(lambda p: os.path.join(p, 'src'), project_inputs))))
        build.ExtendDoxygenConfig(doxy_conf)
        cmd = ('( ' + ' ; '.join(doxy_conf) + ' ) '
               + '| doxygen - 2>doxygen-error.log 1>doxygen-output.log')
        summary_action = cmd + '; python python/infra/ParseDoxygen.py doxygen-output.log doxygen-error.log ' + output_dir
    else:
        summary_action = 'python python/DisplayTests.py '+output_dir+' '+build_type
        show_output_folder = '@echo "Test output written to: ' + senv['ENV']['CHASTE_TEST_OUTPUT'] + '"'
  
    summary_index = os.path.join(output_dir, 'index.html')
    senv.Command(summary_index, Flatten(test_log_files), summary_action)
    if show_output_folder:
        senv.AddPostAction(summary_index, show_output_folder)
    # Avoid circular dependencies
    senv.Ignore(summary_index, summary_index)
    senv.Ignore(Dir(output_dir), summary_index)
    # Make sure the summary is always required by any build targets requested explicitly
    for targ in COMMAND_LINE_TARGETS:
        senv.Depends(targ, summary_index)
    senv.Default(summary_index)
    # Also allow other code to add dependencies, by making the summary depend on an Alias
    senv.Depends(summary_index, senv.Alias('test_summary_dependencies'))
    # Try to ensure it runs even if SCons thinks it's up-to-date, just to re-assure the user
    senv.AlwaysBuild(summary_index)


# Finally, tidy up the build targets if we modified the command-line list
# to do 'nice' test suite specification
if not BUILD_TARGETS:
    BUILD_TARGETS.extend(DEFAULT_TARGETS)


if offline_mode:
    import UserList
    ld_lib_paths = other_libpaths
    if use_chaste_libs and not static_libs:
        ld_lib_paths[0:0] = [Dir('lib').abspath]
    # \todo #136: Make the template file an argument 
    Execute(Copy('run-tests.sh','python/infra/offline_template.sh'))
    handle = open('run-tests.sh', 'a')
    print >> handle, ""
    print >> handle, "export LD_LIBRARY_PATH=" + ':'.join(ld_lib_paths)
    print >> handle, "export PATH=" + env['ENV']['PATH']
    print >> handle, ""
    for runner in offline_test_runners:
        if isinstance(runner, (list, UserList.UserList)):
            # For chaste_libs=1
            print >> handle, "$MPI_LAUNCH_COMMAND", runner[0]
        else:
            # For chaste_libs=0
            print >> handle, "$MPI_LAUNCH_COMMAND", runner
    print >> handle, ""
    handle.close()
    os.chmod('run-tests.sh', 0740)
