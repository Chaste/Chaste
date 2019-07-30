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

# If you are an active developer, committing back to the trunk, 
# please uncomment the following line to run tests on duplicate file names, 
# orphaned tests and copyright notices.
do_inf_tests = 1

# Check which version of Ubuntu this is
fp = open('/etc/issue')
ubuntu_ver = fp.read().split()[1]
fp.close()

# First deal with special cases for beta releases etc.
if ubuntu_ver == 'Trusty':
    ubuntu_ver = [14,04]
else:
    ubuntu_ver = map(int, ubuntu_ver.split('.')[0:2])

if ubuntu_ver >= [14,04]:
    petsc_ver = 3.4
    petsc_path = '/usr/lib/petscdir/3.4.2/'
elif ubuntu_ver >= [12,10]:
    petsc_ver = 3.2
    petsc_path = '/usr/lib/petscdir/3.2/'
elif ubuntu_ver >= [10,10]:
    petsc_ver = 3.1
    petsc_path = '/usr/lib/petscdir/3.1/'
elif ubuntu_ver >= [9,10]:
    petsc_ver = 3
    petsc_3_0_path = '/usr/lib/petscdir/3.0.0/'
else:
    petsc_ver = 2
    petsc_2_3_path = '/usr/lib/petscdir/2.3.3/'

petsc_2_2_path = ''
petsc_build_name = 'linux-gnu-c-debug'
petsc_build_name_profile = petsc_build_name
petsc_build_name_optimized = 'linux-gnu-c-opt'

dealii_path = None
intel_path = None
icpc = 'icpc'

# The following c++ compiler flag switches on POSIX multithreading.
# UPDATE - this doesn't seem to be required - is it default, or already pulled in by something else?
#ccflags = '-pthread' 

other_includepaths = ['/usr/include/metis/']
other_libpaths = []
libs_for_petsc = ['petsccontrib', 'X11',
                  'HYPRE', 'spooles', 'superlu',
                  'umfpack', 'amd' # Both for Umfpack
                  ]

#Fixes (possibly temporary) for Natty
if ubuntu_ver >= [11,04]:
    libs_for_petsc.append(['HYPRE_utilities', 
                           'HYPRE_struct_mv', 'HYPRE_struct_ls',  
                           'HYPRE_sstruct_mv', 'HYPRE_sstruct_ls', 
                           'HYPRE_IJ_mv', 'HYPRE_parcsr_ls', 'dmumps'])
if petsc_ver >= 3:
    libs_for_petsc.append('scotch')
else:
    libs_for_petsc.append('sidl')
if petsc_ver >= 3.1:
    libs_for_petsc.remove('petsccontrib')

ccflags=''

boost_libs = ['boost_serialization', 'boost_filesystem','boost_iostreams','boost_program_options']
# New library dependency
ccflags += '-DCHASTE_BOOST_IOSTREAMS -DCHASTE_BOOST_PROGRAM_OPTIONS '

if ubuntu_ver >= [10,10]:
    boost_libs.append('boost_system')
if ubuntu_ver >= [9,10] and ubuntu_ver <= [12,10]:
    boost_libs = map(lambda l: l+'-mt', boost_libs)

other_libraries = libs_for_petsc + boost_libs + \
                  ['xerces-c',
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
    ccflags += '-DCHASTE_R '

# Use QHull if it is available too, requires a couple of fiddly options...
if (MaybeAdd('qhullcpp', '/home/garmir/git/qhull/lib' , '/home/garmir/git/qhull/src/libqhullcpp/') and
    MaybeAdd('qhull', '/home/garmir/git/qhull/lib' , '/home/garmir/git/qhull/src/') ):    
    other_includepaths.append('/home/garmir/git/qhull/src/libqhull/')
    ccflags += '-DCHASTE_QHULL -D qh_QHpointer '   

if (MaybeAdd('sundials_ida','/usr/lib','/usr/include')):
    ccflags += '-DCHASTE_SUNDIALS_IDA '

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
    
    # Extra libraries for VTK output
    vtk_include_path = filter(os.path.isdir, glob.glob('/usr/include/vtk-5*'))
    use_vtk = int(prefs.get('use-vtk', True))
    use_vtk = use_vtk and bool(vtk_include_path)
    if use_vtk:
        # Note: 10.10 uses VTK 5.4, 10.04 uses 5.2, and early use 5.0
        other_includepaths.extend(vtk_include_path)
        other_libraries.extend(['vtkIO', 'vtkCommon', 'vtkGraphics', 'z'])
        if ubuntu_ver >= [11,10]: # 11.10 uses VTK 5.6
            other_libraries.extend(['vtkFiltering'])

    # Is CVODE installed?
    use_cvode = int(prefs.get('use-cvode', True))
    use_cvode = use_cvode and os.path.exists('/usr/lib/libsundials_cvode.so')
    if ubuntu_ver <= [9,04]:
        # We don't support CVODE 2.4
        use_cvode = False
    if use_cvode:
        DetermineCvodeVersion('/usr/include')
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])

