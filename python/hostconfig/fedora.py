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

# If you are an active developer, committing back to the trunk, 
# please uncomment the following line to run tests on duplicate file names, 
# orphaned tests and copyright notices.
# do_inf_tests = 1

# Find PETSC, assuming the installation guide was followed so that CHASTE_LIBS
# is set and 3.6.2 is installed
chaste_libs = os.environ.get('CHASTE_LIBS')
petsc_path = None
if chaste_libs is not None:
    petsc_ver = 3.6
    petsc_path = os.path.join(chaste_libs, 'petsc-3.6.2/')

petsc_build_name = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu'
petsc_build_name_profile = petsc_build_name

dealii_path = None

other_includepaths = ['/usr/include/metis/']
other_libpaths = []

# Figure out if Intel compilers are available, here we assume they are installed under the 
# default location of /opt/intel, if not you will need to change some hardcoded paths below.
if os.path.isdir('/opt/intel'):
    blas_lapack_production = ['mkl_intel_lp64', 'mkl_core', 'mkl_sequential']
    if os.path.isdir('/opt/intel/composerxe'):
        intel_path = '/opt/intel/composerxe'
        if os.path.isdir(intel_path + '/mkl/lib/intel64'):
            other_libpaths.append(intel_path + '/mkl/lib/intel64')
    else:
        options = glob.glob('/opt/intel/compilers_and_libraries_*/linux')
        options.sort()
        if options:
            intel_path = options[-1]
        else:
            intel_path = None
else:
    intel_path = None

# 2017-06-27: Removed openmpi to make scons run
#hdf5_lib = 'hdf5_openmpi'
#other_includepaths.append('/usr/include/hdf5/openmpi/')
hdf5_lib = 'hdf5'

libs_for_petsc = ['petsccontrib', 'X11',
# 2017-06-27 Remove HYPRE and spooles to make scons run
#                  'HYPRE', 'spooles',
                  'superlu',
                  'umfpack', 'amd' # Both for Umfpack
                  ]

# 2017-06-27 Remove HYPRE to make scons run
#libs_for_petsc.append(['HYPRE_utilities', 
#                   'HYPRE_struct_mv', 'HYPRE_struct_ls',  
#                   'HYPRE_sstruct_mv', 'HYPRE_sstruct_ls', 
#                   'HYPRE_IJ_mv', 'HYPRE_parcsr_ls', 'dmumps'])

# 2017-06-27 Removed scotch
#libs_for_petsc.append('scotch')

if petsc_ver >= 3.1:
    libs_for_petsc.remove('petsccontrib')

boost_libs = ['boost_serialization', 'boost_filesystem']
boost_libs.append('boost_system')

# Determine installed version of Xerces
try:
    xerces3 = bool(subprocess.check_output(['dpkg-query', '-W', '-f', '${version}', 'libxerces-c-dev'], stderr=subprocess.STDOUT))
except:
    xerces3 = False
if xerces3:
    # Xerces 3.1
    xerces_lib = 'xerces-c-3.1'
else:
    # Xerces 2.8
    xerces_lib = 'xerces-c'

other_libraries = libs_for_petsc + boost_libs + \
                  [xerces_lib,
                   hdf5_lib, 'z',
                   'parmetis', 'metis']

# Figure out which lapack/blas packages are actually installed!
if os.path.exists('/usr/lib/liblapack-3.so'):
    blas_lapack = ['lapack-3', 'blas-3']
else:
    blas_lapack = ['lapack', 'blas']

tools = {'xsd': '/usr/bin/xsdcxx',
         'mpirun': '/usr/bin/mpirun.openmpi',
         'mpicxx': '/usr/bin/mpic++.openmpi'}

def Configure(prefs, build):
    """Set up the build configuring.
    
    prefs can specify which version of various libraries we should use, and which optional libraries.
    VTK and CVODE support default on if they are installed.
    
    build is an instance of BuildTypes.BuildType.
    """
    global use_cvode
    global use_vtk

    # Whether to use the Intel compiler
    if build.CompilerType() == 'intel':
        tools['mpicxx'] = 'OMPI_CXX=icpc ' + tools['mpicxx']

    # Extra libraries for VTK output
    vtk_base = '/usr/include/vtk-'
    vtk5_include_path = filter(os.path.isdir, glob.glob(vtk_base + '5*'))
    vtk6_include_path = filter(os.path.isdir, glob.glob(vtk_base + '6*'))
    if vtk5_include_path:
        vtk_include_path = vtk5_include_path[0]
    elif vtk6_include_path:
        vtk_include_path = vtk6_include_path[0]
    else:
        vtk_include_path = ''
    use_vtk = int(prefs.get('use-vtk', True))
    use_vtk = use_vtk and bool(vtk_include_path)
    if use_vtk:
        # Note: 10.10 uses VTK 5.4, 10.04 uses 5.2, and early use 5.0.
        # Some systems may have VTK6 but not VTK5.
        vtk_version = vtk_include_path[len(vtk_base):]
        other_includepaths.append(vtk_include_path)
        vtk_libs = ['CommonCore','CommonDataModel','IOXML','IOGeometry','CommonExecutionModel','FiltersCore','FiltersGeometry','FiltersModeling','FiltersSources']
        vtk_ver = map(int, vtk_version.split('.')[:2])
        if vtk_ver >= [6,2]:
            vtk_libs[2:2] = ['IOParallelXML']
        vtk_libs = map(lambda l: 'vtk' + l + '-' + vtk_version, vtk_libs)
        other_libraries.extend(vtk_libs)

    # Is CVODE installed?
    use_cvode = int(prefs.get('use-cvode', True))
    use_cvode = use_cvode and os.path.exists('/usr/lib64/libsundials_cvode.so')
    if use_cvode:
        DetermineCvodeVersion('/usr/include')
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
