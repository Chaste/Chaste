
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


"""Useful functions for use by the build system."""

import glob
import os
import re
import subprocess
import sys
import time
import threading

from SCons.Script import Command, Dir, Value, Copy, Delete
import SCons
import SCons.Action
import SCons.Node
import SCons.Tool
import SCons.Script
import SCons.Scanner

# Other helper functions
our_path = os.path.dirname(os.path.realpath(__file__))
sys.path[0:0] = [our_path]
import BuildTools
relpath = BuildTools.relpath
set = BuildTools.set

def pns(nodes):
    """Pretty-print nodes for debugging."""
    return map(str, nodes)

#fp = open('tsp.txt', 'w')
def tsp(*args):
    """Thread-safer print, for debugging."""
    msg = ' '.join(map(str, args)) + '\n'
    try:
        print >>fp, msg,
    except:
        print msg,

# Should provide a whole-build global lock, if needed
_lock = threading.Lock()
# We also need per-object locks to reduce contention for parallel builds
_object_locks = {}

# Possible extensions for source files in Chaste
chaste_source_exts = ['.cpp', '.xsd', '.cellml']

def FindSourceFiles(env, rootDir, ignoreDirs=[], dirsOnly=False, includeRoot=False,
                    sourceExts=None, checkSourcesExist=False, otherVars={}):
    """Look for source files under rootDir.

    Returns 2 lists: the first of source (.cpp, .xsd, .cellml) files, and the second
    of the directories in which they may be found.

    Optionally:
     * specify ignoreDirs to not search within particular folder names
     * set dirsOnly to True to only find source directories.  In this case
       only a single list is returned
     * set includeRoot to True to include the rootDir in the returned folder list
     * set checkSourcesExist to True to only return directories that actually contain
       source files, not all non-ignored folders found
    """
    source_files = []
    source_dirs = []
    source_exts = sourceExts or chaste_source_exts
    ignoreDirs.append('.svn')
    for dirpath, dirnames, filenames in os.walk(rootDir):
        for dirname in dirnames[:]:
            if dirname in ignoreDirs:
                dirnames.remove(dirname)
        has_sources = not checkSourcesExist
        if not dirsOnly or checkSourcesExist:
            for filename in filenames:
                if os.path.splitext(filename)[1] in source_exts:
                    has_sources = True
                    if not dirsOnly:
                        filepath = os.path.join(dirpath, filename)
                        source_files.append(filepath)
                    if dirsOnly and source_exts == ['.py'] and filename == '__init__.py':
                        # This is a Python package root, so don't recurse further, and use its parent instead
                        dirnames[:] = []
                        has_sources = False
                        parent_dir = os.path.dirname(dirpath)
                        if not parent_dir in source_dirs:
                            source_dirs.append(parent_dir)
                        break
        if source_exts == ['.py'] and 'setup.py' in filenames and not otherVars.get('dyn_libs_only', False):
            # Special case hack for building Python package extension modules
            # We set environment variables (CFLAGS & LDFLAGS) so libraries like SUNDIALS can be found if not in standard locations
            setup_py_path = os.path.join(os.getcwd(), dirpath, 'setup.py')
            e = env.Clone(LIBPATH=otherVars['other_libpaths'])
            cflags = e.Split(e['CCFLAGS'])
            incflags = []
            for i, f in enumerate(cflags):
                if f == otherVars['build'].IncludeFlag():
                    incflags.append(f)
                    incflags.append(cflags[i+1])
            e.Replace(CCFLAGS=incflags)
            os.environ['CHASTE_CYTHON_CFLAGS'] = e.subst("$_CPPINCFLAGS $CCFLAGS")
            os.environ['CHASTE_CYTHON_LDFLAGS'] = e.subst("$_LIBDIRFLAGS $LINKFLAGS")
            e.Execute(SCons.Action.Action('CFLAGS="$_CPPINCFLAGS $CCFLAGS" LDFLAGS="$_LIBDIRFLAGS $LINKFLAGS"'
                                          ' python %s build_ext --inplace' % setup_py_path,
                                          chdir=os.path.dirname(setup_py_path)))
        if has_sources and (includeRoot or dirpath != rootDir):
            source_dirs.append(dirpath)
    if dirsOnly:
        return source_dirs
    elif '.cpp' in source_exts:
        component = os.path.basename(os.path.dirname(os.path.abspath(rootDir)))
        if component == 'global' and rootDir == 'src':
            # Special-case the version info files.
            file_name = os.path.join('src', 'Version.cpp')
            file_node = env.File(file_name)
            if not env['UPDATE_CHASTE_PROVENANCE'] and os.path.exists(file_node.abspath):
                # Don't update provenance info - just use the existing file
                version_value = Value(open(file_node.abspath).read())
            else:
                version_value = Value(GetVersionCpp(file_name + '.in', env))
            file_node = env.Command(file_name, [version_value], GenerateCppFromValue)[0]
            source_files.append(file_node)
            # This one just contains the path to Chaste
            source_files.append(env.Command(os.path.join('src', 'ChasteBuildRoot.cpp'),
                                            [Value(GetChasteBuildRootCpp(env))],
                                            GenerateCppFromValue)[0])
    return source_files, source_dirs


def FasterSharedLibrary(env, library, sources, **args):
    """Override SharedLibrary to update the libs in '#linklib' only if the symbols change.
    This is to avoid rebuilding binaries when a shared library has changes.

    Use it like:
    env['BUILDERS']['OriginalSharedLibrary'] = env['BUILDERS']['SharedLibrary']
    env['BUILDERS']['SharedLibrary'] = FasterSharedLibrary

    This would be MUCH simpler to implement if we could override the
    signature to be used for the SharedLibrary node itself directly. That is
    certainly possible, but would rely on the internal structure of SCons.

    Based on http://www.scons.org/wiki/SharedLibrarySignatureOverride.
    """
    # SCons version compatibility
    if type(library) != type([]):
        library = [library]
    # Use the 'quicker' shallow copy method!
    envContentSig = env.Clone()
    scons_ver = env._get_major_minor_revision(SCons.__version__)
    if scons_ver < (1,0,0):
        envContentSig.TargetSignatures('content')
    else:
        envContentSig.Decider('MD5-timestamp')

    cat = env.OriginalSharedLibrary(library, sources)

    # Copy all the latest libraries to ONE directory for our convenience.
    # Could modify the above to build directly to this dir instead.
    catLib = env.Install('#lib', cat)

    # Now generate the 'interface' file, using the content signature for its target
    catIF = envContentSig.Command(
        '%s.if' % library[0],
        catLib,
        # -g has the synonym --extern-only on GNU systems (but not on MAC OSX)
        'nm -g $SOURCES | cut -c 12- | sort > $TARGET')

    # Install command to copy lib to #linklib, where the link actually occurs.
    # Explicitly make this depend only on the catIF file, which has a target content signature.
    # Thus only if the global symbol list changes is the library copied, and any programs re-linked.
    catLink = env.Command(
        '#linklib/${SHLIBPREFIX}%s${SHLIBSUFFIX}' % library[0],
        '',
        Copy('$TARGET', str(catLib[0])))
    envContentSig.Depends(catLink, catIF)
    # Record the library to link against
    prefix = args.get('libIdPrefix', '')
    env['CHASTE_LIBRARIES'][prefix+library[0]] = catLink[0]

    return cat


def BuildTest(target, source, env):
    """A builder for test executables.

    Takes as a single source the object file compiled for the test,
    and uses its implicit dependencies to work out which other object
    files must be linked against.  For each header file, if it is a
    Chaste header and has a corresponding source file, then we
    should link with the associated object file.  This analysis is
    recursive - we analyse each object file in the same way.

    It requires a few attributes of the environment:
     * env['CHASTE_COMPONENTS']
       A list of the Chaste components.
     * env['CHASTE_OBJECTS']
       A dictionary mapping source file paths (relative to the Chaste root)
       to object file nodes.
     * env['RUNNER_EXE']
       The test runner executable (SCons File node)
     * env['TestBuilder']
       A callable for use in building the test runner executable.
       Should take keyword parameters target (filled by RUNNER_EXE)
       and source (will be given the required object files).
    """
    header_files = set()
    objects = []
    #import thread
    #tid = thread.get_ident()
    #tsp("\n", tid, source[0], 'BuildTest')#, sorted(env['CHASTE_OBJECTS'].keys()))
    #tsp(tid, source[0], sorted(env['CHASTE_COMPONENTS']))
    #S = SCons.Node.StateString

    def process(o):
        """Process an object file as described in BuildTest.__doc__"""
        #tsp(tid, source[0], "process", o, S[o.state], os.path.exists(str(o)), map(str, o.sources), map(str, o.depends), o.implicit is None, S[o.sources[0].state])
        while o.sources[0].state is SCons.Node.executing:
            # Wait for the sources to be generated, so scanning the object does something useful!
            time.sleep(0.01) # 10 ms
        # Ensure scons' dependencies are set up, in a thread-safe fashion
        _lock.acquire()
        obj_lock = _object_locks.setdefault(id(o), threading.Lock())
        _lock.release()
        obj_lock.acquire()
        o.scan()
        obj_lock.release()
        #tsp(tid, source[0], "state post scan", o, S[o.state], os.path.exists(str(o)), map(str, o.sources), map(str, o.depends), map(str, o.implicit))
        # Now process the dependencies
        for d in o.implicit:
            hdr = str(d)
            if hdr not in header_files:
                #tsp(tid, source[0], 'a dep of', o, 'is', hdr, S[o.state])
                header_files.add(hdr)
                # Is this a Chaste header?
                parts = hdr.split(os.path.sep)
                component = parts[0]
                if component in env['CHASTE_COMPONENTS']:
                    # Does it have a source file?
                    base, ext = os.path.splitext(hdr)
                    if base in ['global/src/Version', 'global/src/ChasteBuildRoot']:
                        # Special cases
                        has_source = True
                        source_filename = base + '.cpp'
                    else:
                        for ext in chaste_source_exts:
                            source_filename = base + ext
                            has_source = source_filename in env['CHASTE_OBJECTS']
                            if has_source:
                                break
                    #tsp(tid, source[0], 'has chaste dep', source_filename if has_source else base, has_source)
                    if has_source:
                        # Find the object file(s) and analyse it/them
                        objs = env['CHASTE_OBJECTS'][source_filename]
                        #tsp(tid, source[0], 'src objs', map(str, objs))
                        objects.extend(objs)
                        for obj in objs:
                            process(obj)

    for o in source:
        process(o)
    # Build the test itself
    runner = env['RUNNER_EXE']
    #tsp("Building", runner, "from", pns(source+objects), "by", tid, "\n")
    actual_runner = env['TestBuilder'](target=runner, source=source+objects)
    env.Alias('test_exes', actual_runner)
    assert actual_runner[0] is runner # Just in case
    return None

def RegisterObjects(env, key, objs):
    """Record how objects get built, for the benefit of BuildTest."""
    env['CHASTE_OBJECTS'][key] = objs
    # If the source is something from which C++ is generated, then we need to add objects
    # under other keys, too, to make sure they are found.
    for obj in objs:
        src = obj.sources[0]
        if src.is_derived():
            #tsp('RegisterObjects', key, map(str, objs), obj, src, src.path)
            env['CHASTE_OBJECTS'][src.path] = [obj]


def CloneEnv(env):
    """Clone a construction environment, but don't copy some objects.

    There are a few special dictionaries which need to be the same object in all
    environments, to ensure that updates are shared.  So after calling the normal Clone
    we need to refer back to the original.
    """
    newenv = env.Clone()
    newenv['CHASTE_OBJECTS'] = env['CHASTE_OBJECTS']
    newenv['CHASTE_COMP_DEPS'] = env['CHASTE_COMP_DEPS']
    newenv['CHASTE_LIBRARIES'] = env['CHASTE_LIBRARIES']
    newenv['CHASTE_CPP_PATH'] = env['CHASTE_CPP_PATH']
    newenv['CHASTE_CPP_PATHS'] = env['CHASTE_CPP_PATHS']
    newenv['CHASTE_PYINCPATH'] = env['CHASTE_PYINCPATH']
    newenv['CHASTE_PYINCPATHS'] = env['CHASTE_PYINCPATHS']
    return newenv


def DetermineLibraryDependencies(env, partialGraph):
    """Determine which Chaste libraries each component/project needs to link against.

    The supplied partial dependency graph maps each component/project to a list of its
    direct dependencies.  On exit from this method each will be mapped instead to an
    ordered list of all its direct and indirect dependencies.  The ordering is based
    on a topological sort, with each component listed prior to its dependencies, and
    so determines the linker command line.
    """
    library_mapping = env['CHASTE_LIBRARIES']
    if env['build'].debug:
        print "Initial component dependencies:", partialGraph
        print "Chaste libraries:", library_mapping
    WHITE, GRAY, BLACK = 0, 1, 2
    PROJECT_PREFIX = 'projects/'
    full_graph = {}
    # As well as projects/<name> mapping to the project library, add a mapping for plain <name>,
    # since that's what will be used for linking.
    for name in library_mapping.keys():
        if name.startswith(PROJECT_PREFIX):
            library_mapping[name[len(PROJECT_PREFIX):]] = library_mapping[name]
    def get_lib(comp, linkArg=False):
        """Get the library node for a component."""
        if not isinstance(comp, type('')):
            lib = comp # It's already a library node
        elif comp == 'core':
            lib = None # Not a real library
        elif linkArg:
            # Link with -lNAME, not linklib/libNAME.so
            if comp.startswith(PROJECT_PREFIX):
                lib = comp[len(PROJECT_PREFIX):]
            else:
                lib = comp
        else:
            lib = library_mapping[comp]
        return lib
    def process_comp(node, root, colours={}, gray_stack=[]):
        """Sort the dependencies for a single component."""
        colours[node] = GRAY
        gray_stack.append(node)
        for dep in partialGraph[node]:
            if dep not in partialGraph:
                raise ValueError("Chaste library %s is listed as a dependency of %s but doesn't exist"
                                 % (dep, node))
            dep_col = colours.get(dep, WHITE)
            if dep_col is GRAY:
                i = gray_stack.index(dep)
                raise ValueError("Chaste library dependency cycle found, with the following libraries: "
                                 + ', '.join(gray_stack[i:]))
            elif dep_col is WHITE:
                process_comp(dep, root, colours, gray_stack)
        colours[node] = BLACK
        gray_stack.pop()
        if node is not root:
            node_lib = get_lib(node, linkArg=True)
            if node_lib:
                full_graph[root].append(node_lib)
        else:
            # Prior to this we have each node listed *after* all its dependencies
            full_graph[root].reverse()
            colours.clear() # So it's ready for the next root
    for comp in partialGraph:
        if comp != 'core':
            full_graph[comp] = []
            process_comp(comp, comp)
            # If we're building libraries, make each depend on libraries from dependent components/projects
            if library_mapping:
                comp_lib = get_lib(comp)
                if comp_lib:
                    deps = map(get_lib, full_graph[comp])
                    if deps:
                        env.Depends(comp_lib, deps)
    # Set up construction variables that can be used in CPPPATH and PYINCPATH
    for comp in full_graph:
        comp_name = comp
        if comp_name.startswith(PROJECT_PREFIX):
            comp_name = comp_name[len(PROJECT_PREFIX):]
        env['CHASTE_CPP_PATH'][comp_name] = env.Flatten([env['CHASTE_CPP_PATHS'][c] for c in full_graph[comp]])
        env['CHASTE_PYINCPATH'][comp_name] = env.Flatten([env['CHASTE_PYINCPATHS'][c] for c in full_graph[comp]])
    # Early versions of SCons aren't as nice
    scons_ver = env._get_major_minor_revision(SCons.__version__)
    if scons_ver < (1,0,0):
        for comp in full_graph:
            for i, dep in enumerate(full_graph[comp]):
                assert isinstance(dep, type(''))
                full_graph[comp][i] = '-l' + dep
    if env['build'].debug:
        print "Complete component dependencies:", full_graph
    # Transfer results to partialGraph
    partialGraph.clear()
    partialGraph.update(full_graph)


def CompleteFolderPaths(folders):
    """
    Given a list of folders in the current component, return a list suitable for use in the CPPPATH
    for other components.  It must therefore contain full path information, and include the
    copies within the build folder.  Using SCons' Dir objects makes this easier.
    """
    paths = []
    for f in folders:
        paths.extend([Dir(f), Dir(f).srcnode()])
    return paths


def FindTestsToRun(env, build, BUILD_TARGETS, otherVars,
                   component=None, project=None):
    """Find header files defining tests to run.

    One of component or project must be specified; this says which Chaste
    component or user project to hunt for tests in.

    If otherVars['requested_tests'] is set, then just run those requested test(s)
    that are in this component (or project).
    If instead otherVars['all_tests'] is True, then find all tests listed in test packs
    in this component/project.
    Otherwise, if this component (or project) is being built (determined by
    checking BUILD_TARGETS) then search for all tests listed in the test packs
    specified by build.TestPacks().

    Returns an iterable of header file leaf paths (relative to the test folder).
    """
    testfiles = set()
    # Check arguments
    assert component or project
    if component:
        assert project is None
    else:
        component = project
    # Check for a single test
    if otherVars['requested_tests']:
        for test_comp, test_path in otherVars['requested_tests']:
            if test_comp == component:
                testfiles.add(test_path)
                # Remove any old test output file to force a re-run
                try:
                    base = os.path.splitext(test_path)[0]
                    os.remove(base + '.log')
                except OSError:
                    pass
    else:
        # Are we building (some of) this component/project?
        test_this_comp = False
        test_subset_folders = []
        root_dir = Dir('#').abspath
        this_comp_targets = set(['.', root_dir])
        if not project and component in env['CHASTE_COMP_DEPS']['core']:
            this_comp_targets.add('core')
        if project:
            comp_path = os.path.join('projects', project)
        else:
            comp_path = component
        this_comp_targets.add(comp_path)
        this_comp_targets.add(os.path.join(root_dir, comp_path))
        test_root_path = os.path.join(comp_path, 'test') + os.sep
#        print map(str, BUILD_TARGETS)
#        print component, project, this_comp_targets
        for targ in BUILD_TARGETS:
            targ = str(targ)
            if targ.endswith(os.sep):
                # Allow users to specify (e.g.) "global/" as a target (handy for use with tab completion)
                targ = targ[:-len(os.sep)]
            if targ in this_comp_targets:
                test_this_comp = True
                break
            if targ.startswith(test_root_path):
                # A particular test folder has been given
                test_subset_folders.append(targ[len(test_root_path):])
        if test_this_comp or test_subset_folders:
            if otherVars['all_tests']:
                pack_names = set()
            else:
                pack_names = set(build.TestPacks())
            if test_this_comp:
                testfiles.update(BuildTools.GetTestsInTestPacks('../../test', pack_names))
            else:
                for folder in test_subset_folders:
                    testfiles.update(BuildTools.GetTestsInTestPacks('../../test', pack_names, subfolder=folder))
#    print component, project, testfiles
    return testfiles


def ExeName(env, exePath):
    """Figure out the real name of an executable.

    Given the Linux-style path, this works out what an executable is actually
    called, so that we can run on cygwin.
    """
    pre = env.subst('$PROGPREFIX')
    suf = env.subst('$PROGSUFFIX')
    dirpath = os.path.dirname(exePath)
    name = os.path.basename(exePath)
    return os.path.join(dirpath, pre+name+suf)


def GetPathRevision(path):
    """Use svnversion or git, as appropriate, to get the revision number of the given path.

    Returns a pair (revision, is_locally_modified).
    For svn repositories, the revision is just the svn revision number, as a string.
    For git repositories, we return "0xSHA1", i.e. treating the hash for the latest commit
    as a hexadecimal number.
    """
    if os.path.exists(os.path.join(path, '.git')):
        # Git repo
        commit_sha = subprocess.check_output(['git', '-C', path, 'rev-parse', '--short=7','HEAD']).strip()
        retcode = subprocess.call(['git', '-C', path, 'diff-index', '--quiet', 'HEAD', '--'])
        revision = '0x' + commit_sha
        modified = (retcode != 0)
    elif os.path.exists(os.path.join(path, '.svn')):
        # Subversion repo
        revision = subprocess.check_output(['svnversion', path]).strip()
        modified = revision[-1] in 'MSP'
        try:
            # Extract upper end of range, and store modified flag
            while revision[-1] in 'MSP':
                revision = revision[:-1]
            revision = int(revision[1+revision.rfind(':'):])
        except:
            revision = 'UINT_MAX'
    else:
        # In abscence of git, subversion (or ReleaseVersion.txt) we report revision=0
        # This allows a zip download from github to be compiled with SCons
        revision, modified = 0, False
    return (revision, modified)

def GetProjectVersions(projectsRoot, default_revision=None):
    """Return C++ code filling in the project versions map."""
    code = ""
    for entry in sorted(os.listdir(projectsRoot)):
        entry_path = os.path.join(projectsRoot, entry)
        if entry[0] != '.' and os.path.isdir(entry_path):
            revision, _ = GetPathRevision(entry_path)
            if revision == 'Unknown' and default_revision:
                revision = default_revision
            if str(revision).startswith('0x'):
                # Take off the confusing hex symbol
                revision = str(revision[2:])
            code += '%sversions["%s"] = "%s";\n' % (' '*8, entry, revision)
    return code

def GetProjectModified(projectsRoot):
    """Return C++ code filling in the map of whether projects have been modified."""
    code = ""
    for entry in sorted(os.listdir(projectsRoot)):
        entry_path = os.path.join(projectsRoot, entry)
        if entry[0] != '.' and os.path.isdir(entry_path):
            _ , modified = GetPathRevision(entry_path)

            code += '%smodified["%s"] = "%s";\n' % (' '*8, entry, modified)
    return code

def GetVersionCpp(templateFilePath, env):
    """Return the contents of the Version.cpp source file."""
    chaste_root = Dir('#').abspath
    version_file = os.path.join(chaste_root, 'ReleaseVersion.txt')
    wc_modified = False
    if os.path.exists(version_file):
        # Extract just the revision number from the file.
        full_version = open(version_file).read().strip()
        last_component = full_version[1+full_version.rfind('.'):]
        if last_component[0] is not 'r':
            # It's a git SHA
            chaste_revision = '0x' + last_component[1:]
        else:
            # It's an svn revision number
            chaste_revision = last_component
        default_revision_for_projects = chaste_revision
    else:
        chaste_revision, wc_modified = GetPathRevision(chaste_root)
        default_revision_for_projects = None

    time_format = "%a, %d %b %Y %H:%M:%S +0000"
    build_time = time.strftime(time_format, time.gmtime())
    compiler_type = env['build'].CompilerType()
    if compiler_type == "gcc":
        compiler_version = subprocess.Popen( [ "gcc", "-dumpversion" ], stdout=subprocess.PIPE ).communicate()[0].strip()
    elif compiler_type == "intel":
        compiler_version = subprocess.Popen( [ "icpc", "-dumpversion" ], stdout=subprocess.PIPE ).communicate()[0].strip()
    else:
        compiler_version = "unknown"

    xsd_path = os.path.join(chaste_root, env['build'].tools['xsd']) # Account for possible relative path being given
    command = xsd_path + ' version 2>&1'
    xsd_version_string = os.popen(command).readline().strip()
    xsd_version = xsd_version_string[-5:]
    if xsd_version == "ctory": # "No such file or directory" !!
        xsd_version = "undetermined"

    from CheckForCopyrights import current_notice
    licence = current_notice.replace('\nThis file is part of Chaste.\n', '')
    licence = licence.replace('"', '\\"').replace('\n', '\\n')

    subst = {'example': '%(example)s',
             'chaste_root': chaste_root,
             'revision': chaste_revision,
             'wc_modified': str(wc_modified).lower(),
             'project_versions': GetProjectVersions(os.path.join(chaste_root, 'projects'), default_revision_for_projects),
             'project_modified': GetProjectModified(os.path.join(chaste_root, 'projects')),
             'licence': licence,
             'time_format': time_format,
             'time_size': len(build_time)+1,
             'build_time': build_time,
             'uname': ' '.join(os.uname()),
             'build_type': env['build'].build_type,
             'build_dir': env['build'].build_dir,
             'build_info': env['CHASTE_BUILD_INFO'],
             'compiler_type': compiler_type,
             'compiler_version': compiler_version,
             'cc_flags': env['build'].CcFlags(),
             'xsd_version': xsd_version}
    return open(templateFilePath).read() % subst

def GetChasteBuildRootCpp(env):
    """Return the contents of the ChasteBuildRoot.cpp source file."""
    subst = {'chaste_root': Dir('#').abspath,
             'build_type': env['build'].build_type,
             'build_dir': env['build'].build_dir}
    return """
#include "ChasteBuildRoot.hpp"

const char* ChasteBuildRootDir()
{
    return "%(chaste_root)s/";
}

const char* ChasteSourceRootDir()
{
    return "%(chaste_root)s/";
}

std::string ChasteComponentBuildDir(const std::string& rComponent)
{
    return std::string(ChasteBuildRootDir()) + rComponent + "/build/%(build_dir)s/";
}

std::string ChasteBuildDirName()
{
    return "%(build_dir)s";
}

std::string ChasteBuildType()
{
    return "%(build_type)s";
}

""" % subst

def _GenerateCppFromValue(env, target, source):
    """An Action to generate a source file from a value node.

    Use like:
    env.Command('global/src/Version.cpp', [Value(GetVersionCpp(templateFilePath, env))], GenerateCppFromValue)
    or:
    env.Command('global/src/ChasteBuildRoot.cpp', [Value(GetChasteBuildRootCpp())], GenerateCppFromValue)
    """
    out = open(target[0].path, "w")
    out.write(source[0].get_contents())
    out.close()
GenerateCppFromValue = SCons.Action.Action(_GenerateCppFromValue, "Generating $TARGET from build information.")

def RecordBuildInfo(env, build_type, static_libs, use_chaste_libs):
    """Record key build information for the provenance system."""
    # TODO: Add library versions used, etc?
    build_info = build_type
    if use_chaste_libs:
        libtype = ['shared', 'static'][static_libs]
        build_info += ', ' + libtype + ' libraries'
    else:
        build_info += ', no Chaste libraries'
    env['CHASTE_BUILD_INFO'] = build_info

def CreateXsdBuilder(build, buildenv, fakeIt=False):
    """Add a 'builder' for running xsd to generate parser code from an XML schema.

    Adds the converter as a source action to the C/C++ builders.

    If fakeIt is True, the builder doesn't actually run the converter.  This means
    that inner SCons runs for generating dynamically loadable cell models don't try
    regenerating it unnecessarily.
    """
    # Check if  'xsd' is really CodeSynthesis xsd...
    if not SCons.Script.GetOption('clean'):
        command = build.tools['xsd'] + ' version 2>&1'
        xsd_version_string = os.popen(command).readline().strip()
        if xsd_version_string.startswith('XML Schema Definition Compiler'):
            xsd_version = 2
        elif xsd_version_string.startswith('CodeSynthesis XSD XML Schema to C++ compiler'):
            xsd_version = 3 # Or 4, but we don't need the detailed version here
        else:
            print "Unexpected XSD program found:"
            print xsd_version_string
            sys.exit(1)
        # And to assist transitioning to the new builder...
        for path in glob.glob('heart/src/io/ChasteParameters*.?pp'):
            try:
                os.remove(path)
            except OSError:
                pass

    def RunXsd(target, source, env):
        """Action for running XSD."""
        schema_file = str(source[0])
        output_dir = os.path.dirname(target[0].abspath)
        command = [build.tools['xsd'], 'cxx-tree',
                   '--generate-serialization',
                   '--output-dir', output_dir,
                   '--std', 'c++11',
                   '--hxx-suffix', '.hpp', '--cxx-suffix', '.cpp',
                   '--prologue-file', 'heart/src/io/XsdPrologue.txt',
                   '--epilogue-file', 'heart/src/io/XsdEpilogue.txt',
                   '--namespace-regex', 'X.* $Xchaste::parametersX',
                   '--namespace-regex', 'X.* https://chaste.comlab.ox.ac.uk/nss/parameters/(.+)Xchaste::parameters::v$1X',
                   schema_file]
        rc = subprocess.call(command)
        return rc

    if fakeIt:
        XsdAction = buildenv.Action('@echo -n')
    else:
        XsdAction = buildenv.Action(RunXsd)
    def XsdEmitter(target, source, env):
        hpp = os.path.splitext(str(target[0]))[0] + '.hpp'
        if env['INSTALL_FILES']:
            t = env.Install(os.path.join(env['INSTALL_PREFIX'], 'include'), hpp)
            env.Alias('install', t)
        target = target + [hpp]
        env.Precious(target)
        return (target, source)
    # Add XSD as a source of .cpp files
    c_file, cxx_file = SCons.Tool.createCFileBuilders(buildenv)
    cxx_file.add_action('.xsd', XsdAction)
    cxx_file.add_emitter('.xsd', XsdEmitter)

def CreatePyCmlBuilder(build, buildenv):
    """Create a builder for running PyCml to generate C++ source code from CellML.

    PyCml is run to generate as many types of output as we can.  If a .out file is
    present, giving output from Maple, this will include backward Euler code.  A
    -conf.xml file may be given to tune this process somewhat, by specifying extra
    arguments to be passed to ConvertCellModel.py; it will also be used as the
    configuration file for PyCml itself.

    The SCons environment variable PYCML_EXTRA_ARGS may also be used to pass additional
    command-line arguments to ConvertCellModel.py.  This can be set in a project
    SConscript file, by including lines like the following prior to the DoProjectSConscript
    call:
        env = SConsTools.CloneEnv(env)
        env['PYCML_EXTRA_ARGS'] = ['--expose-annotated-variables']
    """
    def IsDynamicSource(source):
        parts = source[0].srcnode().path.split(os.path.sep)
        if parts[0] == 'projects':
            dyn_i = 2
        else:
            dyn_i = 1
        return (parts[dyn_i] == 'dynamic' or
                (parts[dyn_i] == 'build' and parts[dyn_i+2] == 'dynamic'))
    def HasMapleOutput(source):
        out_file = os.path.splitext(source[0].srcnode().abspath)[0] + '.out'
        return os.path.exists(out_file), out_file
    def HasConfigFile(source):
        conf_file = os.path.splitext(source[0].srcnode().abspath)[0] + '-conf.xml'
        return os.path.exists(conf_file), conf_file
    script = os.path.join(Dir('#').abspath, 'python', 'ConvertCellModel.py')
    def GetArgs(target, source, env):
        args = ['-A', '-p', '--output-dir', os.path.dirname(target[0].abspath)]
        args.extend(env.get('PYCML_EXTRA_ARGS', []))
        if SCons.Script.Main.GetOption('silent'):
            args.append('--quiet')
        if IsDynamicSource(source):
            # If we're creating a dynamic library, do things differently:
            # only create a single output .so.  The helper script will recognise
            # the -y flag.
            args.append('-y')
        else:
            args.extend(['--normal', '--opt', '--cvode'])
# Won't work until SCons' C scanner can understand #ifdef
#            if 'CHASTE_CVODE' not in env['CPPDEFINES']:
#                args.remove('--cvode')
            if HasMapleOutput(source)[0]:
                args.append('--backward-euler')
        has_conf, conf_file = HasConfigFile(source)
        if has_conf:
            args.append('--conf=' + conf_file)
        return args
    def RunPyCml(target, source, env):
        args = GetArgs(target, source, env)
        command = [script] + args + [str(source[0])]
        print "Running", command
        rc = subprocess.call(command)
        return rc
    PyCmlAction = buildenv.Action(RunPyCml)
    def PyCmlEmitter(target, source, env):
        args = GetArgs(target, source, env)
        args.append('--show-outputs')
        command = [script] + args + [str(source[0])]
        process = subprocess.Popen(command, stdout=subprocess.PIPE)
        filelist = process.communicate()[0]
        if process.returncode != 0:
            print filelist
            raise IOError('Failed to run PyCml; return code = ' + str(process.returncode))
        # Adjust targets to match what the script will actually create
        target = map(lambda s: s.strip(), filelist.split())
        # Make sure the targets depend on everything they might need
        has_maple, maple_output = HasMapleOutput(source)
        if has_maple:
            env.Depends(target, maple_output)
        has_conf, conf_file = HasConfigFile(source)
        if has_conf:
            env.Depends(target, conf_file)
        # Add dependency on pycml source code (ignoring .pyc files)
        pycml_code = glob.glob(os.path.join(Dir('#/python/pycml').abspath, '*'))
        pycml_code = filter(lambda path: not (path.endswith('.pyc') or path.endswith('protocol.py')), pycml_code)
        env.Depends(target, pycml_code)
        env.Depends(target, script)
        # Install headers if requested
        if not IsDynamicSource(source) and env['INSTALL_FILES']:
            headers = [t for t in target if t.endswith('.hpp')]
            t = env.Install(os.path.join(env['INSTALL_PREFIX'], 'include'), headers)
            env.Alias('install', t)
        return (target, source)

    # Add PyCml as a source of .cpp files
    c_file, cxx_file = SCons.Tool.createCFileBuilders(buildenv)
    cxx_file.add_action('.cellml', PyCmlAction)
    cxx_file.add_emitter('.cellml', PyCmlEmitter)

class PyScanner(SCons.Scanner.Classic):
    """A scanner for import lines in Python source code.

    This partly supports "import module", "from module import ...", package imports and relative imports.
    It doesn't handle cases like "from .package import module" properly, as it makes the package the dependency.
    """
    base = SCons.Scanner.Classic # old-style class so can't super()
    def __init__(self, *args, **kw):
        name = kw.get('name', 'PyScanner')
        suffixes = ['.py']
        path_variable = 'PYINCPATH'
        subst = {'ws': r'[ \t]+',         # required whitespace
                 'ows': r'[ \t]*',        # optional whitespace
                 'nm': r'[a-zA-Z0-9_.]+'} # a module/package/entity name
        regex = r'^%(ows)s(?:import%(ws)s(%(nm)s)|from%(ws)s(%(nm)s)%(ws)simport%(ws)s%(nm)s)(?:%(ws)sas%(ws)s%(nm)s)?%(ows)s$' % subst
        self.base.__init__(self, name, suffixes, path_variable, regex, *args, **kw)

    def find_include_names(self, node):
        # Our regex has 2 groups, only 1 of which will ever match, so join them together
        return [''.join(n) for n in self.cre.findall(node.get_text_contents())]

    def find_include(self, include, source_dir, path):
        if '.' in include:
            # It's a package or relative include
            if include == ".":
                # Special case for "from . import name"
                return (None, include)
            while include[0] == '.':
                include = include[1:]
                if include[0] == '.':
                    source_dir = source_dir.Dir('..')
                path = ()
            include = include.replace('.', os.sep)
            n, i = self.base.find_include(self, include + '.py', source_dir, path)
        else:
            # Simple case
            n, i = self.base.find_include(self, include + '.py', source_dir, path)
        if n is None:
            # The name might be a package, so look for a nested __init__.py
            n, i = self.base.find_include(self, os.path.join(include, '__init__.py'), source_dir, path)
        return n, i


def CreateTexttestBuilder(build, env, otherVars):
    """Create a builder that will run texttest and parse the results."""
    texttest_result_line = re.compile(r'<H2>.*: \d+ tests:( (?P<s>\d+) succeeded)?( (?P<f>\d+) FAILED)?</H2>')
    def TexttestParser(target, source, env):
        """Parse results from texttest to figure out if acceptance tests passed.
        target is a dummy file, since we don't know what we'll output until we're done.
        source is the texttest results file.
        """
        fp = open(str(source[0]))
        fails, succs = 0, 0
        for line in fp:
            m = texttest_result_line.match(line)
            if m:
                fails = int(m.group('f') or 0)
                succs = int(m.group('s') or 0)
                break
        fp.close()
        if fails == 0 and succs == 0:
            status = 'unknown'
        elif fails == 0:
            status = 'OK'
        else:
            status = '%d_%d' % (fails, fails+succs)
        to_file_name = build.output_dir + '/AcceptanceTests.' + status + '.0'
        # Remove any old copies of results from this test
        oldfiles = glob.glob(os.path.join(build.output_dir, 'AcceptanceTests.*'))
        for oldfile in oldfiles:
            os.remove(oldfile)
        # Copy results and update summary dependencies
        env.Execute(Copy(to_file_name, str(source[0])))
        if otherVars['test_summary']:
            env.Alias('test_summary_dependencies', to_file_name)
        return None
    parse_builder = env.Builder(action=TexttestParser)
    env['BUILDERS']['ParseTexttest'] = parse_builder

def RunAcceptanceTests(build, env, appsPath, testsPath, exes, otherVars):
    """
    If appsPath is specifically requested on the command line,
    schedule its acceptance tests for running.
    """
    checkout_dir = Dir('#').abspath
    assert os.path.isabs(appsPath)
    for targ in otherVars['COMMAND_LINE_TARGETS']:
        if not os.path.isabs(targ):
            targ = os.path.join(checkout_dir, targ)
        if targ.endswith(os.path.sep):
            targ = targ[:-len(os.path.sep)]
        if targ == appsPath:
            break
    else:
        # Not found
        return
#    print "Running acceptance tests", map(str, exes)
    tests_dir = Dir(testsPath).abspath
    output_dir = env['ENV']['CHASTE_TEST_OUTPUT']
    if not os.path.isabs(output_dir):
        output_dir = os.path.join(checkout_dir, output_dir)
    texttest = build.tools['texttest'] + ' -d ' + tests_dir
    texttest_output_dir = os.path.join(output_dir, 'texttest_reports', 'chaste')
    time_eight_hours_ago = time.time() - 8*60*60
    canonical_test_date = time.strftime("%d%b%Y", time.localtime(time_eight_hours_ago))
    todays_file = os.path.join(texttest_output_dir, 'test_default_' + canonical_test_date + '.html')
    # The next 2 lines make sure the acceptance tests will get run, and the right results stored
    env.Execute(Delete(todays_file))
    env.Execute(Delete(os.path.join(output_dir, 'texttest_output')))
    env.Command(todays_file, exes,
                ['-' + texttest + ' -b default -c ' + checkout_dir,
                 texttest + ' -b default -c ' + checkout_dir + ' -coll web'])
    dummy = os.path.join(appsPath, 'dummy.texttest')
    env.ParseTexttest(dummy, todays_file)
    env.Depends(appsPath, [todays_file, dummy])
    if otherVars['test_summary']:
        env.Alias('test_summary_dependencies', dummy)

def BuildExes(build, env, appsPath, components, otherVars):
    """Build 'standalone' executables (i.e. with their own main(), not using cxxtest).

    apps_path should refer to a directory structure containing folders:
      src - no subfolders, contains .cpp file(s) each defining main()
      texttest/chaste - definitions for acceptance tests, optional

    components gives the Chaste libraries that these executables link against.
    """
    env = CloneEnv(env)
    src_path = os.path.join(appsPath, 'src')

    if otherVars['static_libs']:
        libpath = '#lib'
        env.Append(LINKFLAGS=' -static ')
        if sys.platform != 'cygwin':
            env.Append(LINKFLAGS='-pthread ')
    else:
        libpath = '#linklib'
    env.Replace(LIBPATH=[libpath] + otherVars['other_libpaths'])
    env.Prepend(CPPPATH=src_path)
    env.Replace(LIBS=components + otherVars['other_libs'])

    exes = []
    for main_cpp in glob.glob(os.path.join(src_path, '*.cpp')):
        exes.append(env.Program(main_cpp))

    if not otherVars['compile_only']:
        # Run acceptance tests if present
        test_path = os.path.join(appsPath, 'texttest', otherVars['texttest_suite'])
        if os.path.isdir(test_path):
            CreateTexttestBuilder(build, env, otherVars)
            RunAcceptanceTests(build, env, appsPath, test_path, exes, otherVars)


def ScheduleTestBuild(env, overrides, testfile, prefix, use_chaste_libs):
    """Set the compilation of a single test.

    This handles the logic of building with or without chaste_libs, and ensures
    the test is added to the default targets.  The behaviour is indentical for
    projects and core components.

    @param env  the main SCons environment to use
    @param overrides  construction variable overrides for building the test exe
    @param testfile  the path of the test .hpp file, relative to the 'test' folder
    @param prefix  testfile without extension
    @param use_chaste_libs  whether to use chaste_libs
    """
    test_hpp = os.path.join('test', testfile)
    runner_cpp = env.Test(prefix+'Runner.cpp', test_hpp)
    runner_exe = env.File(ExeName(env, prefix+'Runner'))
    if use_chaste_libs:
        runner_dummy = None
        runner_exe = env.Program(runner_exe, runner_cpp, **overrides)
        # Make sure we build the test unless the user says otherwise
        env.Default(runner_exe)
    else:
        runner_obj = env.StaticObject(runner_cpp)
        runner_dummy = env.File(prefix+'.dummy')
        try:
            # Work around for some SCons versions
            runner_exe.get_executor()._morph()
        except AttributeError:
            pass
        env.BuildTest(runner_dummy, runner_obj, RUNNER_EXE=runner_exe)
        env.AlwaysBuild(runner_dummy)
        env.Depends(runner_exe, runner_dummy)
        env.Alias('test_exes', runner_dummy)
        # Make sure we build the test unless the user says otherwise
        env.Default(runner_dummy)
    return runner_exe, runner_dummy


def DoDynamicallyLoadableModules(baseEnv, otherVars):
    """Logic to find and schedule for building any dynamically-loadable modules.

    These are generated from source files found in a 'dynamic' folder in a component
    or project.
    """
    # Dynamically loadable modules need different flags
    dyn_env = CloneEnv(baseEnv)
    dyn_env.Append(LINKFLAGS=' '+otherVars['build'].rdynamic_link_flag)
    # Find any source files that should get made into dynamically loadable modules.
    curdir = os.getcwd()
    os.chdir('../..')
    dyn_source, dyn_cpppath = FindSourceFiles(dyn_env, 'dynamic', includeRoot=True, otherVars=otherVars)
    os.chdir(curdir)
    dyn_env.Prepend(CPPPATH=dyn_cpppath)
    # Try to avoid likely conflicts
    for path in [p for p in dyn_env['CPPPATH'] if 'cellml' in str(p)]:
        dyn_env['CPPPATH'].remove(path)
    # Build any dynamically loadable modules
    dyn_libs = []
    for s in dyn_source:
        # Figure out if this is the single requested source, if there is one
        so_dir = os.path.join(curdir, os.pardir, os.pardir, os.path.dirname(s))
        this_dyn_only = (otherVars['dyn_libs_only'] and
                         os.path.realpath(so_dir).startswith(os.path.realpath(otherVars['dyn_folder'])))
        if this_dyn_only and not otherVars['use_chaste_libs']:
            # Warn the user that this combination is unsupported and likely to break
            import SCons.Warnings
            class DynamicLoadingWarning(SCons.Warnings.Warning):
                pass
            SCons.Warnings.enableWarningClass(DynamicLoadingWarning)
            SCons.Warnings.warn(DynamicLoadingWarning,
                                "Using dynamically loaded cell models MAY NOT WORK with chaste_libs=0.\n"
                                "    If you get unresolved symbol errors, use chaste_libs=1 as an argument to scons.")
        if this_dyn_only or not otherVars['dyn_libs_only']:
            # Note: if building direct from CellML, there will be more than 1 target
            dyn_objs = dyn_env.SharedObject(source=s)
            for o in dyn_objs:
                so_lib = dyn_env.OriginalSharedLibrary(source=o)
                dyn_libs.append(dyn_env.Install(so_dir, so_lib))
                if this_dyn_only:
                    # Force dependency on installed version, even if project is a symlink
                    dyn_env.Depends(otherVars['dyn_folder'], dyn_libs[-1])
    return dyn_libs


def DoProjectSConscript(projectName, chasteLibsUsed, otherVars):
    """Main logic for a project's SConscript file.

    The aim of this method is that a project's SConscript file should be able to be as
    simple as:
        import os
        Import("*")
        project_name = os.path.basename(os.path.dirname(os.path.dirname(os.getcwd())))
        chaste_libs_used = ['core'] # or ['heart'], ['cell_based', 'crypt'], etc.
        result = SConsTools.DoProjectSConscript(project_name, chaste_libs_used, globals())
        Return("result")
    """
    if otherVars['debug']:
        print "Executing SConscript for project", projectName
    # Commonly used variables
    env = CloneEnv(otherVars['env'])
    use_chaste_libs = otherVars['use_chaste_libs']
    project_path = os.path.join('projects', projectName)
    # Store our dependencies
    env['CHASTE_COMP_DEPS'][project_path] = chasteLibsUsed
    # Note that because we are using SCons' variant dir functionality, and the build
    # dir is created before the SConscript files are executed, that the working dir
    # will be set to <project>/build/<something>.
    curdir = os.getcwd()
    # Look for .cpp files within the project's src folder
    os.chdir('../..') # This is so .o files are built in <project>/build/<something>/
    files, cpppath = FindSourceFiles(env, 'src', ignoreDirs=['broken'], includeRoot=True, otherVars=otherVars)
    pydirs = FindPyDirs(env, otherVars)
    # Look for source files that tests depend on under <project>/test/.
    testsource, test_cpppath = FindSourceFiles(env, 'test', ignoreDirs=['data'], otherVars=otherVars)
    # Move back to the build dir
    os.chdir(curdir)

    # Update CPPPATH for folders in this project, with a hook to pull in its dependencies later
    env['CHASTE_CPP_PATHS'][projectName] = CompleteFolderPaths(cpppath)
    all_cpppath = cpppath + test_cpppath
    if all_cpppath:
        env.Prepend(CPPPATH=all_cpppath)
    env.Append(CPPPATH="${CHASTE_CPP_PATH['%s']}" % projectName)
    # Python search path is handled similarly
    AddPyDirs(env, project_path, pydirs, projectName)

    # Build any dynamically loadable modules
    dyn_libs = DoDynamicallyLoadableModules(env, otherVars)
    if otherVars['dyn_libs_only']:
        # Short-circuit most stuff, and just build the requested .so
        return []

    # Look for files containing a test suite
    # A list of test suites to run will be found in a test/<name>TestPack.txt file, one per line.
    # Alternatively, a single test suite may have been specified on the command line.
    testfiles = FindTestsToRun(env, otherVars['build'], otherVars['BUILD_TARGETS'],
                               otherVars,
                               project=projectName)
    if otherVars['debug']:
        print "  Will run tests:", map(str, testfiles)

    # Build and install the library for this project
    project_lib = test_lib = libpath = None
    set_env_chaste_libs = True
    if use_chaste_libs:
        if otherVars['static_libs']:
            if files:
                project_lib = env.StaticLibrary(projectName, files)
                project_lib = env.Install('#lib', project_lib)[0]
            libpath = '#lib'
        else:
            if files:
                project_lib = env.SharedLibrary(projectName, files, libIdPrefix='projects/')[0]
                set_env_chaste_libs = False # Done by FasterSharedLibrary
            libpath = '#linklib'
        # Build the test library for this project
        if testsource:
            test_lib = env.StaticLibrary('test'+projectName, testsource)[0]
    else:
        # Build the object files for this project
        for source_file in files + testsource:
            objs = env.StaticObject(source_file)
            key = os.path.join(project_path, source_file)
            RegisterObjects(env, key, objs)
    if set_env_chaste_libs:
        env['CHASTE_LIBRARIES'][project_path] = project_lib

    # Libraries to link against
    chaste_libs = ["${CHASTE_COMP_DEPS['%s']}" % project_path]
    if project_lib:
        chaste_libs[0:0] = [projectName]
    if test_lib:
        chaste_libs[0:0] = ['test' + projectName]
    all_libs = chaste_libs + otherVars['other_libs']

    # Make test output depend on shared libraries, so if implementation changes then tests are re-run.
    lib_deps = []
    if test_lib:
        lib_deps.append(test_lib) # Source code used by our tests
    #lib_deps.extend(map(lambda lib: '#lib/lib%s.so' % lib, chasteLibsUsed)) # all Chaste libs used
    if project_lib:
        lib_deps.append(project_lib)

    test_log_files = ScheduleTests(env, all_libs, dyn_libs, lib_deps, libpath, testfiles, otherVars)

    # Any executables to build?
    if otherVars['build_exes']:
        apps_path = os.path.join(Dir('#').abspath, 'projects', projectName, 'apps')
        if os.path.isdir(apps_path):
            BuildExes(otherVars['build'], env, apps_path, components=chaste_libs, otherVars=otherVars)

    return test_log_files

def CheckForSpecialFiles(env, component, files, otherVars):
    """Schedule compiles of source files that need special compilation flags."""
    special_files = [('mesh', 'src/3rdparty/tetgen1.4.2/predicates.cpp')]
    special_objs = []
    overrides = {'CCFLAGS': '-O0'}
    for special_comp, special_file in special_files:
        if special_comp == component and special_file in files:
            # Remove from files
            files.remove(special_file)
            # Handle specially
            if otherVars['use_chaste_libs']:
                if otherVars['static_libs']:
                    special_objs.extend(env.StaticObject(special_file, **overrides))
                else:
                    special_objs.extend(env.SharedObject(special_file, **overrides))
            else:
                objs = env.StaticObject(special_file, **overrides)
                key = os.path.join(component, str(special_file))
                RegisterObjects(env, key, objs)
                special_objs.extend(objs)
    return special_objs

def FindPyDirs(env, otherVars):
    """Find folders containing Python sources in this component/project."""
    pydirs = FindSourceFiles(env, 'src', sourceExts=['.py'], includeRoot=True, checkSourcesExist=True, dirsOnly=True, otherVars=otherVars)
    # Clear out any stale .pyc files (i.e. those without associated .py sources)
    for pydir in pydirs:
        for pyc in glob.glob(os.path.join(pydir, '*.pyc')):
            if not os.path.exists(pyc[:-1]):
                env.Execute(Delete(pyc))
    return pydirs

def AddPyDirs(env, basePath, pydirs, componentName):
    """Add folders containing Python sources, if any, to the module search path.

    Also sets up build data structures so other components can use them.
    """
    pydir_paths = map(lambda d: os.path.join('#', basePath, d), pydirs)
    if pydirs:
        env.Append(PYINCPATH=pydir_paths)
    # Set the paths so other components/projects can use our Python code
    env['CHASTE_PYINCPATHS'][componentName] = pydir_paths
    env.Append(PYINCPATH='${CHASTE_PYINCPATH["%s"]}' % componentName)


def DoComponentSConscript(component, otherVars):
    """Main logic for a Chaste component's SConscript file.

    The aim of this method is that a component's SConscript file should be able to be as
    simple as:
        import os
        Import("*")
        toplevel_dir = os.path.basename(os.path.dirname(os.path.dirname(os.getcwd())))
        result = SConsTools.DoComponentSConscript(toplevel_dir, globals())
        Return("result")
    """
    if otherVars['debug']:
        print "Executing SConscript for component", component
    # Commonly used variables
    env = CloneEnv(otherVars['env'])
    use_chaste_libs = otherVars['use_chaste_libs']

    # Note that because we are using SCons' variant dir functionality, and the build
    # dir is created before the SConscript files are executed, that the working dir
    # will be set to <component>/build/<something>.
    curdir = os.getcwd()
    # Look for source files within the <component>/src folder
    os.chdir('../..') # This is so .o files are built in <component>/build/<something>/
    files, cpppath = FindSourceFiles(env, 'src', ignoreDirs=['broken'], includeRoot=True, otherVars=otherVars)
    pydirs = FindPyDirs(env, otherVars)
    # Look for source files that tests depend on under test/.
    # We also need to add any subfolders to the CPPPATH, so they are searched
    # for #includes.
    testsource, test_cpppath = FindSourceFiles(env, 'test', ignoreDirs=['data'], includeRoot=True, otherVars=otherVars)
    # Install headers?
    if env['INSTALL_FILES']:
        headers, _ = FindSourceFiles(env, 'src', sourceExts=['.hpp'], otherVars=otherVars)
        t = env.Install(os.path.join(otherVars['install_prefix'], 'include'), headers)
        env.Alias('install', t)
    # Move back to the buid dir
    os.chdir(curdir)

    # Update CPPPATH for folders in this component, with a hook to pull in its dependencies later
    env['CHASTE_CPP_PATHS'][component] = CompleteFolderPaths(cpppath)
    all_cpppath = cpppath + test_cpppath
    if all_cpppath:
        env.Prepend(CPPPATH=all_cpppath)
    env.Append(CPPPATH='${CHASTE_CPP_PATH["%s"]}' % component)
    # Python search path is handled similarly
    AddPyDirs(env, component, pydirs, component)

    # Build any dynamically loadable modules
    dyn_libs = DoDynamicallyLoadableModules(env, otherVars)
    if otherVars['dyn_libs_only']:
        # Short-circuit most stuff, and just build the requested .so
        return []

    # Look for files containing a test suite
    # A list of test suites to run will be found in a test/<name>TestPack.txt
    # file, one per line.
    # Alternatively, a single test suite may have been specified on the command line.
    testfiles = FindTestsToRun(env, otherVars['build'], otherVars['BUILD_TARGETS'],
                               otherVars,
                               component=component)
    if otherVars['debug']:
        print "  Will run tests:", map(str, testfiles)

    special_objects = CheckForSpecialFiles(env, component, files, otherVars)

    # Build and install the library for this component
    set_env_chaste_libs = True
    lib = test_lib = libpath = None
    if use_chaste_libs:
        if otherVars['static_libs']:
            if files:
                lib = env.StaticLibrary(component, files + special_objects)
                # Add explicit dependencies on the Chaste libraries we use too,
                # so this library and its tests will be re-linked if they change
                env.Depends(lib, map(lambda lib: '#lib/lib%s.a' % lib,
                                     env['CHASTE_COMP_DEPS'][component]))
                lib = env.Install('#lib', lib)[0]
            libpath = '#lib'
        else:
            if files:
                lib = env.SharedLibrary(component, files + special_objects)
                set_env_chaste_libs = False # Done by FasterSharedLibrary
            libpath = '#linklib'
        # Build the test library for this component
        if testsource:
            test_lib = env.StaticLibrary('test'+component, testsource)
        # Install libraries?
        if lib and env['INSTALL_FILES']:
            t = env.Install(os.path.join(otherVars['install_prefix'], 'lib'), lib)
            env.Alias('install', t)
    else:
        # Don't build libraries - tests will link against object files directly
        for source_file in files + testsource:
            objs = env.StaticObject(source_file)
            key = os.path.join(component, str(source_file))
            RegisterObjects(env, key, objs)
    if set_env_chaste_libs:
        env['CHASTE_LIBRARIES'][component] = lib

    # Determine libraries to link against.
    # Note that order does matter!
    chaste_libs = ["${CHASTE_COMP_DEPS['%s']}" % component]
    if lib:
        chaste_libs = [component] + chaste_libs
    if test_lib:
        chaste_libs = ['test' + component] + chaste_libs
    all_libs = chaste_libs + otherVars['other_libs']

    # Make test output depend on shared libraries, so if implementation changes then tests are re-run.
    lib_deps = []
    if test_lib:
        lib_deps.append(test_lib) # Source code used by our tests
    #lib_deps.extend(map(lambda lib: '#lib/lib%s.so' % lib, chaste_libs)) # all Chaste libs used
    if lib:
        lib_deps.append(lib)

    test_log_files = ScheduleTests(env, all_libs, dyn_libs, lib_deps, libpath, testfiles, otherVars)

    if otherVars['build_exes']:
        apps_path = os.path.join(Dir('#').abspath, component, 'apps')
        if os.path.isdir(apps_path):
            BuildExes(otherVars['build'], env, apps_path, components=chaste_libs, otherVars=otherVars)
        if component == 'heart':
            # Legacy special case!
            BuildExes(otherVars['build'], env, Dir('#/apps').abspath, components=chaste_libs, otherVars=otherVars)

    return test_log_files


def ScheduleTests(env, allLibs, dynLibs, libDeps, libPath, testFiles, otherVars):
    """Schedule the test builds & runs for this component/project."""
    use_chaste_libs = otherVars['use_chaste_libs']
    # Collect a list of test log files to use as dependencies for the test summary generation
    test_log_files = []

    # Build and run tests of this component/project
    if testFiles:
        if not use_chaste_libs:
            overrides = {'LIBS': otherVars['other_libs'],
                         'LIBPATH': otherVars['other_libpaths']}
            env['TestBuilder'] = lambda target, source: env.Program(target, source, **overrides)
        else:
            overrides = {'LIBS': allLibs,
                         'LIBPATH': [libPath, '.'] + otherVars['other_libpaths']}

    for testfile in testFiles:
        prefix, ext = os.path.splitext(testfile)
        is_py_test = ext == '.py'
        if not is_py_test:
            (runner_exe, runner_dummy) = ScheduleTestBuild(env, overrides, testfile, prefix, use_chaste_libs)
            if otherVars['offline_mode']:
                otherVars['offline_test_runners'].append(runner_exe)
        if not otherVars['compile_only']:
            log_file = env.File(prefix+'.log')
            test_log_files.append(log_file)
            if is_py_test:
                env.PyRunTest(log_file, os.path.join('test', testfile))
            else:
                if use_chaste_libs:
                    env.Depends(log_file, libDeps)
                else:
                    env.Depends(log_file, runner_dummy)
                if dynLibs:
                    # All tests should depend on dynamically loadable modules, just in case
                    env.Depends(log_file, dynLibs)
                env.RunTest(log_file, runner_exe)
            if otherVars['force_test_runs']:
                env.AlwaysBuild(log_file)

    return test_log_files
