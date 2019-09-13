
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

"""Host-specific configuration.

This module contains the logic to set various build parameters based on a
host configuration file.  Configuration files for developer machines may be
stored in a 'machines' package.  If a machine doesn't have a specific file,
default.py is used.

Each configuration file should be a Python program providing certain variables
specifying where to find libraries and tools.  In each case, if the library or
tool is not present, the variable should be set to None.  These variables are:
 * petsc_2_2_path - path to PETSc 2.2 install
 * petsc_2_3_path - path to PETSc 2.3 install
 * petsc_3_0_path - path to PETSc 3.0 install
 * petsc_path     - path to PETSc install for versions 3.1 and newer
 * dealii_path    - path to Deal.II install
 * metis_path     - path to METIS install
 * intel_path     - path to Intel compiler installation

These flags indicate whether to use certain optional external libraries.
Eventually we will use SCons' Configure functionality to make this unnecessary.
If the variable is not present, it will default to False.  If the variable is
set to True, then the appropriate paths (see below) should be specified, too.
 * use_cvode  - whether to use CVODE
 * use_vtk - whether to use VTK development libraries
 * use_adaptivity - whether to use adaptivity library
 
 * other_includepaths - list of paths containing other header files
 * other_libpaths     - list of paths containing other libraries, including
                        metis, xsd, and boost
 * other_libraries    - list of other libraries to link against
 * blas_lapack        - the names of the blas and lapack libraries on this system
                        (a 2 element list)

 * ccflags - any extra compiler flags needed, as a string.

 * tools - a dictionary mapping executable names to their absolute paths, for tools
           not found on $PATH.
 * icpc - a special case tool: allows you to change the Intel compiler command, e.g.
          to add flags specific to 64-bit machines.

Any non-absolute paths will be considered relative to the root of the Chaste install.
"""

import glob
import imp
import os
import re
import sys
import types

# A local.py file in this folder overrides everything
try:
    import local as conf
except ImportError:
    # Do we have any machine-specific config for our host name?
    try:
        import machines
        conf = machines.config_module()
    except ImportError:
        # How about distro-specific config?
        try:
            if sys.platform == 'darwin':
                distro = 'macosx'
            else:
                # Get distro via Python
                import platform
                distro = platform.linux_distribution()[0].lower()
                # Get distro manually
                if distro == '':
                    fp = open('/etc/issue')
                    distro = fp.read().split()[0].lower()
                    fp.close()
            (file, pathname, desc) = imp.find_module(distro, ['python/hostconfig'])
            try:
                conf = imp.load_module(distro, file, pathname, desc)
            finally:
                file.close()
        except (ImportError, IOError):
            # Fall-back to default settings
            import default as conf


# This is a bit ugly at present: SConstruct calls configure() to fill
# these global variables in, then reads them directly.
# Note that the order in which things are added to these lists often matters!
libpaths = []
incpaths = []
libraries = []

# The functions below are supplied to machine config files for their use

def ConfigError(message):
    """Print an error message and quit."""
    raise ValueError(message)
    
def CheckPathExists(path, component):
    """Check that a path exists, dieing if not."""
    if not os.path.exists(path):
        ConfigError("Cannot use %s from %s: path does not exist." % (component, path))

def RemoveFromPath(pathList, searchString):
    """Remove path entries from pathList that contain searchString."""
    for path in pathList[:]:
        if path.find(searchString) != -1:
            pathList.remove(path)
    return

def AddBoost(basePath, version, forceUseSystem=False):
    """Use Boost installed in a non-standard location.

    Expects basePath to point to a folder containing include and lib folders,
    and version to be of the form '1.36' or '1.33.1'.

    Will automatically account for extended Boost library naming schemes.
    Can also handle boost libraries already appearing in other_libpaths etc.
    """
    CheckPathExists(basePath, 'Boost')
    # Remove existing Boost libs
    for lib in conf.other_libraries[:]:
        if lib.startswith('boost_'):
            conf.other_libraries.remove(lib)
    RemoveFromPath(conf.other_includepaths, 'boost')
    RemoveFromPath(conf.other_libpaths, 'boost')
    # Add libs from new location
    if float(version[:4]) >= 1.40:
        # The [:4] is to cope with versions like '1.33.1'
        inc = ''
    else:
        inc = 'boost-' + version.replace('.', '_')
    conf.other_includepaths.append(os.path.join(basePath, 'include', inc))
    libpath = os.path.join(basePath, 'lib')
    conf.other_libpaths.append(libpath)
    boost_libs = ['boost_serialization', 'boost_filesystem']
    if float(version[:4]) >= 1.42 or forceUseSystem:
        boost_libs.append('boost_system')
    testlib = boost_libs[0]
    base = os.path.join(libpath, 'lib' + testlib)
    matches = glob.glob(base + '*.so')
    if not matches:
        # We might only have static libraries
        matches = glob.glob(base + '*.a')
        if not matches:
            raise ValueError('Boost library ' + testlib + ' not found in ' + basePath)
    suffix = os.path.splitext(matches[0][len(base):])[0]
    for lib in boost_libs:
        conf.other_libraries.append(lib + suffix)
    return

def AddHdf5(basePath):
    """Use HDF5 from a non-standard location."""
    CheckPathExists(basePath, 'HDF5')
    # Remove existing paths
    RemoveFromPath(conf.other_includepaths, 'hdf5')
    RemoveFromPath(conf.other_libpaths, 'hdf5')
    # Add new location
    conf.other_includepaths.append(os.path.join(basePath, 'include'))
    conf.other_libpaths.append(os.path.join(basePath, 'lib'))
    if not 'hdf5' in conf.other_libraries:
        conf.other_libraries.extend(['hdf5', 'z'])
    return

def AddXerces(basePath):
    """Use Xerces from a non-standard location."""
    CheckPathExists(basePath, 'Xerces')
    # Remove existing paths
    RemoveFromPath(conf.other_includepaths, 'xerces')
    RemoveFromPath(conf.other_libpaths, 'xerces')
    RemoveFromPath(conf.other_libraries, 'xerces') # The default library is called xerces-c-3.1 on some systems
    # Add new location
    conf.other_includepaths.append(os.path.join(basePath, 'include'))
    conf.other_libpaths.append(os.path.join(basePath, 'lib'))
    if not 'xerces-c' in conf.other_libraries:
        conf.other_libraries.append('xerces-c')
    return

def AddParMetis(basePath):
    """Use ParMetis from a non-standard location."""
    CheckPathExists(basePath, 'ParMetis')
    # Remove existing paths
    RemoveFromPath(conf.other_includepaths, 'ParMetis')
    RemoveFromPath(conf.other_libpaths, 'ParMetis')
    # Add new location
    conf.other_includepaths.append(basePath)
    conf.other_libpaths.append(basePath)
    if not 'parmetis' in conf.other_libraries:
        conf.other_libraries.extend(['parmetis', 'metis'])
    return

def AddXsd(basePath):
    """Use CodeSynthesis XSD from a non-standard location."""
    CheckPathExists(basePath, 'XSD')
    # Remove existing include path
    RemoveFromPath(conf.other_includepaths, 'libxsd')
    # Add new location
    conf.other_includepaths.append(os.path.join(basePath, 'libxsd'))
    conf.tools['xsd'] = os.path.join(basePath, 'bin', 'xsd')
    return

def DetermineCvodeVersion(includePath):
    """Figure out which version of CVODE is installed by grepping its version header."""
    CheckPathExists(includePath, 'CVODE')
    config_h = os.path.join(includePath, 'sundials', 'sundials_config.h')
    if os.path.exists(config_h):
        version_re = re.compile(r'#define SUNDIALS_(?:PACKAGE_)?VERSION "?(\d+).(\d+).(\d+)"?')
        for line in open(config_h):
            m = version_re.match(line)
            if m:
                major = int(m.group(1))
                minor = int(m.group(2))
                patch = int(m.group(3))
                conf.cvode_version = str(major*10000 + minor*100 + patch)
                break
        else: 
            ConfigError('Unable to find version string within sundials_config.h')
    else:
        ConfigError('Unable to find sundials_config.h to read CVODE version')

def TryRemove(pathGlob):
    """Try to remove files matching the given glob pattern, ignoring errors."""
    for path in glob.glob(pathGlob):
        try:
            os.remove(path)
        except OSError:
            pass

# Supply the above functions to the config module
for name in dir():
    if name[0] != '_':
        item = globals()[name]
        if type(item) == types.FunctionType:
            setattr(conf, name, item)

if sys.version_info[:2] < (2,7):
    # Add subprocess.check_output to earlier Python versions
    import subprocess
    def check_output(*popenargs, **kwargs):
        """Run command with arguments and return its output as a string. Copied from Python 2.7"""
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            print >>sys.stderr, "Called process failed; output was:"
            print >>sys.stderr, output
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = check_output


def EnsureVariablesDefined():
    """Ensure that the other_* variables are defined; default to empty lists if not."""
    conf.other_libraries = getattr(conf, 'other_libraries', [])
    conf.other_libpaths = getattr(conf, 'other_libpaths', [])
    conf.other_includepaths = getattr(conf, 'other_includepaths', [])

def DoPetsc(version, optimised, profile=False, production=False, includesOnly=False):
    """Determine PETSc include and library paths.

    The locations vary depending on the version of PETSc, and possibly
    whether optimised libraries are to be used.

    The version can be given as 2.2, 2.3, or 3.0 etc. to choose PETSc version.
    If a host doesn't support the requested version we try older ones.
    A ValueError is raised if no suitable version is present.

    Set optimised to True to use optimised builds of the libraries rather than debug builds.
    Set profile to True to use profile builds of PETSc.
    """
    if os.environ.get('CHASTE_LOAD_ENV', ''):
        return

    conf.petsc_2_2_path = getattr(conf, 'petsc_2_2_path', None)
    conf.petsc_2_3_path = getattr(conf, 'petsc_2_3_path', None)
    conf.petsc_3_0_path = getattr(conf, 'petsc_3_0_path', None)
    conf.petsc_path = getattr(conf, 'petsc_path', None)
    requested_version = version
    version_number = map(int, version.split('.'))
    if version_number > [3,0] and (conf.petsc_path is None or not os.path.isdir(conf.petsc_path)):
        # Use 3.0 instead
        version = '3.0'
    if version == '3.0' and (conf.petsc_3_0_path is None or 
                             not os.path.isdir(conf.petsc_3_0_path)):
        # Use 2.3 instead
        version = '2.3'
    if version == '2.3' and (conf.petsc_2_3_path is None or 
                             not os.path.isdir(conf.petsc_2_3_path)):
        # Use 2.2 instead
        version = '2.2'
    if version == '2.2' and (conf.petsc_2_2_path is None or 
                             not os.path.isdir(conf.petsc_2_2_path)):
        # Raise a friendly error
        ConfigError('PETSc %s requested, but no path for this or an earlier version given in the host config.' % requested_version)
    version_number = map(int, version.split('.')) # The actual version being used, now
    
    def GetBuildNameList():
        build_names = []
        if production:
            build_names.append(getattr(conf, 'petsc_build_name_production', conf.petsc_build_name))
        elif profile:
            build_names.append(getattr(conf, 'petsc_build_name_profile', conf.petsc_build_name))
        elif optimised:
            build_names.append(getattr(conf, 'petsc_build_name_optimized', conf.petsc_build_name))
        # Always have a fallback option to the default build name
        build_names.append(conf.petsc_build_name)
        # Some systems may only have optimised PETSc, so have that as final fallback
        build_names.append(getattr(conf, 'petsc_build_name_optimized', conf.petsc_build_name))
        return build_names
    
    if version == '2.2':
        petsc_base = os.path.abspath(conf.petsc_2_2_path)
        # Gracefully fall back to optimised/non-opt if the requested one isn't there
        if optimised:
            dirs = ['libO_c++', 'libg_c++']
        else:
            dirs = ['libg_c++', 'libO_c++']
        for d in dirs:
            libpath = os.path.join(petsc_base, 'lib', d, conf.petsc_build_name)
            if os.path.exists(libpath): break
        else:
            ConfigError('No PETSc 2.2 libraries found.')
        incpaths.append(os.path.join(petsc_base, 'bmake', conf.petsc_build_name))
    elif version == '2.3':
        petsc_base = os.path.abspath(conf.petsc_2_3_path)
        CheckPathExists(petsc_base, 'PETSc')
        for build_name in GetBuildNameList():
            libpath = os.path.join(petsc_base, 'lib', build_name)
            if os.path.isdir(libpath): break
        incpaths.append(os.path.join(petsc_base, 'bmake', build_name))
    elif version_number >= [3,0]:
        if version == '3.0':
            petsc_base = os.path.abspath(conf.petsc_3_0_path)
        else:
            petsc_base = os.path.abspath(conf.petsc_path)
        CheckPathExists(petsc_base, 'PETSc')
        for build_name in GetBuildNameList():
            libpath = os.path.join(petsc_base, build_name, 'lib')
            if os.path.isdir(libpath): break
        incpaths.append(os.path.join(petsc_base, build_name, 'include'))
        # PETSc 3 allows us to automatically download openmpi.
        # If we do, make sure to use the correct mpicxx/mpirun.
        binpath = os.path.join(petsc_base, build_name, 'bin')
        if not hasattr(conf, 'tools'):
            conf.tools = {}
        if os.path.exists(os.path.join(binpath, 'mpicxx')):
            conf.tools['mpicxx'] = os.path.abspath(os.path.join(binpath, 'mpicxx'))
        if os.path.exists(os.path.join(binpath, 'mpirun')):
            conf.tools['mpirun'] = os.path.abspath(os.path.join(binpath, 'mpirun'))
    else:
        ConfigError('Unrecognised PETSc version requested: ' + version)
    incpaths.append(os.path.join(petsc_base, 'include'))
    if not os.path.isdir(incpaths[-2]):
        ConfigError('PETSc headers directory %s not found.' % incpaths[-2])
    if not includesOnly:
        if not os.path.isdir(libpath):
            ConfigError('PETSc libraries directory %s not found.' % libpath)
        libpaths.append(libpath)
        if version_number >= [3,1]:
            if os.path.exists(os.path.join(libpath, 'libpetsc_real.so')):
                libraries.append('petsc_real')
            else:
                libraries.append('petsc')
        else:
            libraries.extend(['petscts', 'petscsnes', 'petscksp', 'petscdm', 
                              'petscmat', 'petscvec', 'petsc'])
            if sys.platform == 'cygwin':
                libraries.extend(['gdi32', 'user32', 'advapi32', 'kernel32', 'dl'])

def DoDealii(build):
    """Add Deal.II include & library paths, and libraries.

    Deal.II uses different library *names* to distinguish optimised versions.
    """
    conf.dealii_path = getattr(conf, 'dealii_path', None)
    if conf.dealii_path is None:
        ConfigError('Deal.II required, but no path given in the host config.')
    base = os.path.abspath(conf.dealii_path)
    CheckPathExists(base, 'Deal.II')
    # Check Deal.II version
    version = open(os.path.join(base, 'Version')).read()
    if not version.startswith('6.'):
        # Older Deal.II requires PETSc 2.2
        do_petsc('2.2', build.is_optimised)
    else:
        # Just pick up the header files
        do_petsc('2.3', build.is_optimised, includes_only=True)
    # Add Deal.II libraries
    libpaths.append(os.path.join(base, 'lib'))
    relative_incpaths = ['base/include', 'lac/include', 'deal.II/include']
    incpaths.extend(map(lambda relpath: os.path.join(base, relpath),
                        relative_incpaths))
    libs = ['deal_II_1d', 'deal_II_2d', 'deal_II_3d', 'lac', 'base']
    if version.startswith('6.'):
        libs.append('petscall')
    if build.dealii_debugging:
        libs = map(lambda s: s + '.g', libs)
    libraries.extend(libs)

def Flatten(iterable):
    """Flatten a list which may have nested lists.
    A utility method for OptionalLibraryDefines."""
    it = iter(iterable)
    for e in it:
        if isinstance(e, (list, tuple)):
            for f in Flatten(e):
                yield f
        else:
            yield e

def OptionalLibraryDefines():
    """
    Work out what optional libraries have been asked for,
    and return the appropriate #define names, as a list.
    """
    # Todo #2367 take out adaptivity, and replace with warning/error?
    possible_flags = {'cvode': 'CHASTE_CVODE', 'vtk': 'CHASTE_VTK', 'adaptivity': 'CHASTE_ADAPTIVITY', 'gmsh': 'CHASTE_GMSH', 'complex_bessel': 'CHASTE_COMPLEX_BESSEL'}
    actual_flags = []
    if getattr(conf, 'use_cvode', False):
        # Need to set a define for CVODE version.  Assume 2.3.0 if not specified.
        actual_flags.append('CHASTE_SUNDIALS_VERSION=' + getattr(conf, 'cvode_version', '20300'))
    for libname, symbol in possible_flags.iteritems():
        if getattr(conf, 'use_' + libname, False):
            actual_flags.append(symbol)
    for lib in Flatten(conf.other_libraries):
        if lib.startswith('xerces-c'):
            actual_flags.append('CHASTE_XERCES')
            break
    return actual_flags

def Configure(build):
    """Given a build object (BuildTypes.BuildType instance), configure the build."""
    if not build.quiet:
        print "Using hostconfig settings from", os.path.splitext(conf.__file__)[0] + ".py"
    prefs = build.GetPreferedVersions()
    if hasattr(conf, 'Configure') and callable(conf.Configure):
        # The machine config has a method to do its configuration, so call that first.
        conf.Configure(prefs, build)
    if build.using_dealii:
        EnsureVariablesDefined()
        DoDealii(build)
        libraries.extend(conf.other_libraries) # Some of "other_libraries" may depend on BLAS/LAPACK, make sure they are included before them.
        libraries.extend(['blas', 'lapack']) # Use versions provided with Deal.II
    else:
        if prefs:
            if hasattr(conf, 'SetPreferedVersions') and callable(conf.SetPreferedVersions):
                conf.SetPreferedVersions(prefs)
            elif not (hasattr(conf, 'Configure') and callable(conf.Configure)):
                raise ValueError('Machine configuration has no support for setting prefered library versions.')
            petsc_version = prefs.get('petsc', '3.1')[:3]
        else:
            petsc_version = '3.1'
        EnsureVariablesDefined()
        DoPetsc(petsc_version, build.is_optimised, build.is_profile, build.is_production) # PETSc links against some objects defined in "other_libraries"
        libraries.extend(conf.other_libraries) # Some of "other_libraries" may depend on BLAS/LAPACK, make sure they are included before them.
        if build.is_production:
            libraries.extend(conf.blas_lapack_production)
        else:
            libraries.extend(conf.blas_lapack)
    if build.CompilerType() == 'intel':
        intel_path = os.path.abspath(conf.intel_path)
        CheckPathExists(intel_path, 'Intel compilers')
        intel_lib_path = os.path.join(intel_path, 'lib')
        CheckPathExists(intel_lib_path, 'Intel compilers')
        libpaths.append(intel_lib_path)
        # This is a bit of a hack to get version 12 working for #2059
        intel_lib_path = os.path.join(intel_lib_path, 'intel64')
        if os.path.exists(intel_lib_path):
            libpaths.append(intel_lib_path)
    incpaths.extend(conf.other_includepaths)
    libpaths.extend(map(os.path.abspath, conf.other_libpaths))
    # Needed for dynamically loaded cell models
    libraries.append('dl')

    build.tools.update(conf.tools)

    if build.CompilerType() == 'intel':
        # Switch to use Intel toolchain
        if hasattr(conf, 'icpc'):
            build.tools['mpicxx'] += ' -CC="'+conf.icpc+'"'
        intel_bin_path = os.path.join(intel_path, 'bin')
        # This is a bit of a hack to get version 12 working for #2059
        if os.path.exists(os.path.join(intel_bin_path, 'intel64')):
            intel_bin_path = os.path.join(intel_bin_path, 'intel64')
        CheckPathExists(intel_bin_path, 'Intel compilers')
        build.tools['cxx'] = os.path.join(intel_bin_path, 'icpc')
        build.tools['ar'] = os.path.join(intel_bin_path, 'xiar')

    if hasattr(conf, 'ModifyBuild') and callable(conf.ModifyBuild):
        conf.ModifyBuild(build)
    
    if build.debug:
        for name in dir(conf):
            item = getattr(conf, name)
            if name[0] != '_' and type(item) != types.FunctionType:
                print name, '=', item
        print "Libraries:", libraries
        print "Library paths:", libpaths
        print "Include paths:", incpaths
        build.DumpDebugInfo()

def CppDefines():
    """Return a list of extra C preprocessor defines."""
    return OptionalLibraryDefines()

def CcFlags():
    """Return a string containing extra flags for the C++ compiler."""
    return getattr(conf, 'ccflags', '')

def LdFlags():
    """Return a string containing extra flags for the linker."""
    flags = getattr(conf, 'ldflags', '')
    if sys.platform == 'darwin':
        flags += ' -Wl,-undefined,dynamic_lookup'
    return flags
