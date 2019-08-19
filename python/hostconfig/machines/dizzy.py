# Configuration

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

import glob
import os
import subprocess

# Check which version of Ubuntu this is
fp = open('/etc/issue')
ubuntu_ver = fp.read().split()[1]
fp.close()
ubuntu_ver = map(int, ubuntu_ver.split('.')[0:2]) 

# We've installed PETSc ourselves, and got it to install OpenMPI, ParMETIS, and HDF5
if os.path.isdir('/home/bob/petsc-3.6.4'):
    petsc_ver = [3,6]
    petsc_path = '/home/bob/petsc-3.6.4'
elif os.path.isdir('/home/bob/petsc-3.7.7'):
    petsc_ver = [3,7]
    petsc_path = '/home/bob/petsc-3.7.7'
else:
    petsc_ver = [3,1]
    petsc_path = '/home/bob/petsc-3.1-p8/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu-profile'
petsc_build_name_optimized = 'linux-gnu-opt'

# Intel paths etc.
petsc_build_name_production = 'linux-intel-opt-mkl'
intel_path = '/opt/intel'
blas_lapack_production = ['mkl_intel_lp64', 'mkl_core', 'mkl_sequential']
other_libpaths = [os.path.join(petsc_path, 'external_packages', 'f2cblaslapack', petsc_build_name), '/opt/intel/composerxe/mkl/lib/intel64']

# Includes should all be under /usr/include or PETSc
other_includepaths = []

# Default extra compiler & linker flags
ccflags = ''
ldflags = ''

# Basic boost libraries we use
boost_libs = ['boost_serialization', 'boost_filesystem']

# Add this to enable compression of boost archives
#boost_libs.append('boost_iostreams')
# Add an extra flag to tell project that boost_iostreams is available.
# ccflags+=' -DCHASTE_BOOST_IOSTREAMS '

if ubuntu_ver >= [10,10]:
    boost_libs.append('boost_system')
if ubuntu_ver >= [9,10] and ubuntu_ver <= [12,10]:
    boost_libs = map(lambda l: l+'-mt', boost_libs)

# Determine installed version of Xerces
try:
    xerces3 = bool(subprocess.check_output(['dpkg-query', '-W', '-f', '${version}', 'libxerces-c-dev'], stderr=subprocess.STDOUT))
except:
    xerces3 = False
if xerces3:
    if ubuntu_ver >= [18,04]:
        # Xerces 3.2
        xerces_lib = 'xerces-c-3.2'
    else:
        # Xerces 3.1
        xerces_lib = 'xerces-c-3.1'
else:
    # Xerces 2.8
    xerces_lib = 'xerces-c'

other_libraries = boost_libs + [xerces_lib,
                                'hdf5', 'z',
                                'parmetis', 'metis']

def MaybeAdd(libname, libpath, incpath):
    found = False
    if os.path.isdir(libpath) and os.path.isdir(incpath):
        other_libraries.append(libname)
        other_libpaths.append(libpath)
        other_includepaths.append(incpath)
        found = True
    return found

# Use R if it's available
if (MaybeAdd('R', '/usr/lib/R', '/usr/share/R/include/') and
    MaybeAdd('Rcpp', '/usr/local/lib/R/site-library/Rcpp/lib', '/usr/local/lib/R/site-library/Rcpp/include/') and
    MaybeAdd('RInside', '/usr/local/lib/R/site-library/RInside/lib', '/usr/local/lib/R/site-library/RInside/include/')):
    # Add an extra flag to tell Gary's project that R is available.
    ccflags+=' -DCHASTE_R '

# Figure out which lapack/blas packages are actually installed!
if os.path.exists('/usr/lib/liblapack-3.so'):
    blas_lapack = ['lapack-3', 'blas-3']
else:
    blas_lapack = ['lapack', 'blas']

tools = {'texttest': '../../../texttest-3.19/source/bin/texttest.py',
         'xsd': '/usr/bin/xsdcxx'}

if os.path.exists('/usr/bin/google-pprof'):
    tools['pprof'] = '/usr/bin/google-pprof'

do_inf_tests = 1

def Configure(prefs, build):
    """Allow the user to specify which version of various libraries we're using."""
    import socket
    if socket.getfqdn().startswith('scoop'):
        BASE = '/home/lofty/'
    else:
        BASE = '/home/robert/'

    if build.CompilerType() == 'intel':
        # Might need a different PETSc build name!
        global petsc_build_name, petsc_build_name_optimized
        petsc_build_name = 'linux-intel'
        petsc_build_name_optimized = 'linux-intel'

    # PETSc
    global petsc_ver
    global petsc_path
    global petsc_3_0_path
    global petsc_2_3_path
    external_packages = 'externalpackages'
    # We default to the PETSc set above if present for this user
    petsc_set = 'petsc' in prefs
    if petsc_set:
        petsc_ver = map(int, prefs['petsc'].split('.')[:2])
    else:
        prefs['petsc'] = '%d.%d' % tuple(petsc_ver)

    # Note: these paths assume that PETSc 3.1 and above have their own MPI or use the system version, and earlier versions used BASE/mpi
    if not petsc_set:
        pass
    elif prefs['petsc'] == '3.7':
        petsc_path = BASE+'petsc-3.7.2'
    elif prefs['petsc'] == '3.6':
        petsc_path = BASE+'petsc-3.6.1'
    elif prefs['petsc'] == '3.5':
        petsc_path = BASE+'petsc-3.5.2'
    elif prefs['petsc'] == '3.4':
        petsc_path = BASE+'petsc-3.4.2'
    elif prefs['petsc'] == '3.3':
        petsc_path = BASE+'petsc-3.3-p5'
    elif prefs['petsc'] == '3.2':
        petsc_path = BASE+'petsc-3.2-p7'
    elif prefs['petsc'] == '3.1':
        petsc_path == '/home/bob/petsc-3.1-p8/'
    elif prefs['petsc'] == '3.0':
        petsc_3_0_path = BASE+'petsc-3.0.0-p12'
        petsc_path = petsc_3_0_path
        tools['mpicxx'] = BASE+'mpi/bin/mpicxx'
        tools['mpirun'] = BASE+'mpi/bin/mpirun'
    elif prefs['petsc'] == '2.3.3':
        petsc_2_3_path = BASE+'petsc-2.3.3-p15' #For testing SUBMINOR=3
        petsc_path = petsc_2_3_path
        tools['mpicxx'] = BASE+'mpi/bin/mpicxx'
        tools['mpirun'] = BASE+'mpi/bin/mpirun'
    elif prefs['petsc'] == '2.3':
        petsc_2_3_path = BASE+'petsc-2.3.2-p4'
        petsc_path = petsc_2_3_path
        tools['mpicxx'] = BASE+'mpi/bin/mpicxx'
        tools['mpirun'] = BASE+'mpi/bin/mpirun'
    else:
        raise ValueError('Unsupported PETSc version "%s" requested' % prefs['petsc'])

    if petsc_ver >= [3,4]:
        # Add a path to the local version of h5dump for acceptance testing
        # h5dump should have the correct indentation (currently 1.8.10 or above) and the binary should match with any LD_LIBRARY_PATH .so files
        binpath = os.path.join(petsc_path, petsc_build_name, 'bin')
        assert os.path.isdir(binpath)
        os.environ['PATH'] = binpath + ':' + os.environ['PATH']

    if build.is_optimised:
        build_name = petsc_build_name_optimized
    else:
        build_name = petsc_build_name
    # This assumes that other_libpaths starts with a location (might be fake) to a PETSc external packages location
    assert 'f2cblaslapack' in other_libpaths[0]
    other_libpaths[0] = os.path.join(petsc_path, external_packages, 'f2cblaslapack', build_name)

    # Boost library version
    installed_boost = subprocess.check_output(['dpkg-query', '-W', '-f', '${version}', 'libboost-serialization-dev'], stderr=subprocess.STDOUT)[:4]
    prefs['boost'] = prefs.get('boost', installed_boost)
    if prefs['boost'] != installed_boost:
        boost_dir = BASE+'boost_' + prefs['boost'].replace('.', '_')
        if os.path.isdir(boost_dir):
            AddBoost(boost_dir, prefs['boost'], forceUseSystem=True)
        else:
            raise ValueError('Unsupported Boost version "%s" requested' % prefs['boost'])

    # HDF5.  Note that sometimes this is installed by PETSc, and so this setting will be ignored!
    if 'hdf5' in prefs:
        if prefs['hdf5'] == '1.8':
            AddHdf5(BASE+'hdf5_1_8_9/')
        elif prefs['hdf5'] == '1.8.13':
            assert prefs['petsc'] == '3.2', 'HDF5 1.8.13 is only available with PETSc 3.2'
        elif prefs['hdf5'] == '1.6':
            AddHdf5(BASE+'hdf5/')
        else:
            raise ValueError('Unsupported HDF5 version "%s" requested' % prefs['hdf5'])

    # Xerces
    if 'xerces' in prefs:
        if prefs['xerces'] == '2.8' and xerces3:
            AddXerces(BASE+'xercesc_2_8')
        elif prefs['xerces'] == '3.1' and not xerces3:
            AddXerces(BASE+'xercesc_3_1')

    # XSD
    #installed_xsd = subprocess.check_output(['dpkg-query', '-W', '-f', '${version}', 'xsdcxx'], stderr=subprocess.STDOUT)[:3]
    if 'xsd' in prefs:
        if prefs['xsd'] == '4.0':
            AddXsd(BASE+'xsd-4.0.0-x86_64-linux-gnu')
        elif prefs['xsd'] == '3.3':
            AddXsd(BASE+'xsd-3.3.0-x86_64-linux-gnu')
        elif prefs['xsd'] == '3.2':
            AddXsd(BASE+'xsd-3.2.0-x86_64-linux-gnu')
        elif prefs['xsd'] == '3.1':
            AddXsd(BASE+'xsd-3.1.0-x86_64-linux-gnu')
        elif prefs['xsd'] == '2.3':
            AddXsd(BASE+'xsd-2.3.1-x86_64-linux-gnu')
        else:
            raise ValueError('Unsupported XSD version "%s" requested' % prefs['xsd'])

    # Parmetis.  We often don't switch the libraries explicitly for this, since Parmetis 4 is obtained via PETSc 3.3, and Parmetis 3 by PETSc 3.1 & 3.2.
    # Instead we check that the requested PETSc version is suitable.
    if petsc_ver <= [3, 0]:
        prefs['parmetis'] = prefs.get('parmetis', '3')
    if 'parmetis' in prefs:
        global ccflags
        parmetis = prefs['parmetis']
        ccflags += ' -DCHASTE_PARMETIS_REQUIRED=' + str(parmetis)
        if parmetis == '3' and petsc_ver >= [3, 3]:
            raise ValueError('To use ParMETIS 2 you must install it using an older PETSc (<= 3.2)')
        if parmetis == '4' and petsc_ver < [3, 3]:
            raise ValueError('To use ParMETIS 4 you must install it using a recent PETSc (>= 3.3)')
        if parmetis == '3' and petsc_ver <= [3, 0]:
            other_includepaths.append(BASE+'ParMetis-3.1')
            other_libpaths.append(BASE+'ParMetis-3.1')

    # Determine the version of CVODE to be used
    global use_cvode
    use_cvode = int(prefs.get('use-cvode', 1))
    if use_cvode:
        if 'cvode' not in prefs:
            cvode_path = '/usr/' # Use the system installation
        elif prefs['cvode'] == '2.5':
            cvode_path = BASE+'cvode'
        elif prefs['cvode'] == '2.6':
            cvode_path = BASE+'cvode_2_6'
        elif prefs['cvode'] == '2.7':
            cvode_path = BASE+'cvode_2_7'
        elif prefs['cvode'] == '2.8':
            if os.path.exists('/opt/cvode-libs') and build.CompilerType() == 'intel':
                # CVODE 2.8.1 (Sundials 2.6.1) specially downloaded and compiled with intel optimisations on this machine,
                # these shared libraries link to intel libraries, so need to have intel paths set up to work...
                # N.B. we needed to hack the CMakeLists.txt slightly to get the version reported correctly - see #2479 for a copy.
                cvode_path = '/opt/cvode-libs'
            else:
                cvode_path = BASE+'cvode_2_8'
        elif prefs['cvode'] == '2.9':
            cvode_path = BASE+'cvode_2_9'
        elif prefs['cvode'] == '3.1':
            cvode_path = BASE+'cvode_3_1'
        elif prefs['cvode'] == '4.1':
            cvode_path = BASE+'cvode_4_1'
        else:
            raise ValueError('Unsupported CVODE version "%s" requested' % prefs['cvode'])
        cvode_inc = os.path.join(cvode_path, 'include')
        DetermineCvodeVersion(cvode_inc)
        if cvode_path != '/usr/':
            other_includepaths.append(cvode_inc)
            other_libpaths.append(os.path.join(cvode_path, 'lib'))
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])

    # Extra libraries for VTK output
    global use_vtk
    use_vtk = int(prefs.get('use-vtk', 1))
    if use_vtk:
        # Pick version 1) if the user has asked for a version, 2) if there is a system 6.2 version on the machine, 3) fall back to system version 5
        default_vtk_version = '5'                             # (3) Fall back to version 5.x (Ubuntu 14.04)
        vtk_62_include_path = '/usr/include/vtk-6.2'
        vtk_63_include_path = '/usr/include/vtk-6.3'
        if (os.path.isdir(vtk_62_include_path)):
            default_vtk_version = '6.2'                       # (2) System-wide 6.2 (Ubuntu 16.04) 
        vtk_71_include_path = '/usr/include/vtk-7.1'
        if (os.path.isdir(vtk_71_include_path)):
            default_vtk_version = '7.1'                       # Ubuntu 18,04
        if (os.path.isdir(vtk_63_include_path)):              # (2) System-wide 6.3 (Ubuntu 18.04 - optional dependency)
            default_vtk_version = '6.3'
        prefs['vtk'] = prefs.get('vtk', default_vtk_version)  # (1) User's choice
        # Assume that the system wide version of VTK is 6.2 (on 16.04) and use this in preference
        if prefs['vtk'] == '7.1':
            if (not os.path.isdir(vtk_71_include_path)):
                raise ValueError("No system headers for VTK 7 found at "+vtk_71_include_path)
            other_includepaths.append(vtk_71_include_path)
            vtk_libs = ['CommonCore','CommonDataModel', 'IOParallelXML', 'IOXML','IOGeometry','CommonExecutionModel','FiltersCore','FiltersGeometry','FiltersModeling','FiltersSources']
            vtk_libs = map(lambda l: 'vtk' + l + '-7.1', vtk_libs)
            other_libraries.extend(vtk_libs)
        # Assume that the system wide version of VTK is 6.2 (on 16.04) and use this in preference
        elif prefs['vtk'] == '6.2':
            if (not os.path.isdir(vtk_62_include_path)):
                raise ValueError("No system headers for VTK 6 found at "+vtk_62_include_path)
            other_includepaths.append(vtk_62_include_path)
            vtk_libs = ['CommonCore','CommonDataModel', 'IOParallelXML', 'IOXML','IOGeometry','CommonExecutionModel','FiltersCore','FiltersGeometry','FiltersModeling','FiltersSources']
            vtk_libs = map(lambda l: 'vtk' + l + '-6.2', vtk_libs)
            other_libraries.extend(vtk_libs)
        # The rest of this *may* switch versions to those in old /home/lofty folders
        elif prefs['vtk'][0] == '6':
            vtk_ver = map(int, prefs['vtk'].split('.')[:2])
            subst = {'v': prefs['vtk']}
            other_includepaths.append(BASE+'vtk-%(v)s/include/vtk-%(v)s' % subst)
            other_libpaths.append(BASE+'vtk-%(v)s/lib' % subst)
            vtk_libs = ['CommonCore','CommonDataModel','IOXML','IOGeometry','CommonExecutionModel','FiltersCore','FiltersGeometry','FiltersModeling','FiltersSources']
            if vtk_ver >= [6,2]:
                vtk_libs[2:2] = ['IOParallelXML']
            vtk_libs = map(lambda l: 'vtk' + l + '-' + prefs['vtk'], vtk_libs)
            other_libraries.extend(vtk_libs)
        elif prefs['vtk'] == '5.10':
            other_includepaths.append(BASE+'vtk-5.10/include/vtk-5.10')
            other_libpaths.append(BASE+'vtk-5.10/lib/vtk-5.10')
            other_libraries.extend(['vtkGraphics','vtkFiltering','vtkIO','vtkCommon', 'vtksys', 'vtkexpat', 'vtkzlib'])
        else:
            # Included for completeness: assume an older version of Ubuntu (and that the vtk version has been set)
            vtk_include_path = filter(os.path.isdir, glob.glob('/usr/include/vtk-5*'))
            if len(vtk_include_path) != 1:
                raise ValueError("No or multiple system headers for VTK 5 found")
            vtk_version = float(vtk_include_path[0].split('-')[1])
            other_includepaths.append(vtk_include_path[0])
            other_libraries.extend(['vtkIO','vtkGraphics','vtkCommon', 'z'])
            if vtk_version >= 5.6:
                other_libraries.append('vtkFiltering')

    # Check that adaptivity is NOT requested, as it is no longer supported (#2367)
    if int(prefs.get('use-adaptivity', 0)):
        raise ValueError('The adaptivity library is no longer supported (see ticket:2367)')

    # Support for GoogleProfile build
    if build.build_dir == 'google_profile':
        other_libraries.extend(['profiler'])
