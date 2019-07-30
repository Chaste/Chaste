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

import os
import sys

############################################################
# TO CONFIGURE YOUR MACHINE: if you have followed the manual
# installation instructions, edit the definition of
# chaste_libs_path below to point to where you installed the
# dependencies.
############################################################

# If you are an active developer, committing back to the trunk, 
# please uncomment the following line to run tests on duplicate file names, 
# orphaned tests and copyright notices.
# do_inf_tests = 1

#EDIT HERE
#For a simple installation all paths will be below this directory
chaste_libs_path =  '/Users/administrator1/Documents/chaste-libs'
#EDIT HERE

if not os.path.exists(chaste_libs_path) or not os.path.isdir(chaste_libs_path):
    print >>sys.stderr, "Chaste dependencies folder", chaste_libs_path, \
        "not found; please edit python/hostconfig/default.py"
    sys.exit(1)

petsc_3_0_path = chaste_libs_path+'petsc-3.0.0-p8/'
petsc_2_3_path = chaste_libs_path+'petsc-2.3.3-p15/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
dealii_path = ''
parmetis_path = chaste_libs_path+'/ParMetis-3.1/'
# If you have the Intel compiler installed, set this to the folder where it lives
intel_path = '/opt/intel/cc/9.1.039/'
# You may need to edit this to ensure that the intel compiler finds the right gcc libraries, e.g.
#icpc = 'icpc -gcc-version=410 -I /usr/include/c++/4.1.3/x86_64-linux-gnu/ -I/usr/include/c++/4.1.3/    -I/usr/include/c++/4.1.3/backward'
icpc = 'icpc'

if os.uname()[4] == 'x86_64':
    xsd_path = chaste_libs_path + 'xsd-3.2.0-x86_64-linux-gnu/'
else:
    xsd_path = chaste_libs_path+'xsd-3.2.0-i686-linux-gnu/'

other_includepaths = [chaste_libs_path+'hdf5/include',
                      chaste_libs_path+'xerces/include',
                      chaste_libs_path+'boost/include/boost-1_34_1',
                      xsd_path + 'libxsd',
                      parmetis_path]

other_libpaths = [chaste_libs_path+'lib',
                  chaste_libs_path+'boost/lib', 
                  chaste_libs_path+'xerces/lib',
                  chaste_libs_path+'hdf5/lib',
                  chaste_libs_path+'rdf/lib',
                  os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu'),
                  parmetis_path]

# The order of libraries in these lists matters!
# Note that boost serialization sometimes has a different name: eg boost_serialization-gcc41
other_libraries = ['boost_serialization', 'boost_filesystem', 'xerces-c', 'hdf5', 'z', 'parmetis', 'metis']

# If using Boost >= 1.42, uncomment the following line
#other_libraries.append('boost_system')
if os.path.exists(petsc_3_0_path):
    # Assume user has followed INSTALLATION.txt and has Fortran blas & lapack
    blas_lapack = ['flapack', 'fblas', 'gfortran', 'gfortranbegin']
    other_libraries.append('HYPRE')
else:
    blas_lapack = ['f2clapack', 'f2cblas'] # Note: Lapack before BLAS
  # Note: parmetis before metis, hdf5 before z.

# Note that boost serialization sometimes has a different name:
# other_libraries = ['boost_serialization-gcc41', 'xerces-c', 'hdf5', 'z', 'parmetis', 'metis']

tools = {'mpirun': chaste_libs_path+'mpi/bin/mpirun',
         'mpicxx': chaste_libs_path+'mpi/bin/mpicxx',
         'xsd': xsd_path+'bin/xsd'}


def Configure(prefs, build):
    """Set up the build configuring.
    
    prefs can specify which version of various libraries we should use, and which optional libraries.
    
    build is an instance of BuildTypes.BuildType.
    """
    global use_cvode
    global use_vtk
    
    # use_vtk defaults to True. Change to False if VTK development libraries are not available.
    use_vtk = int(prefs.get('use-vtk', True))
    
    # VTK is required for adaptivity to work, so if vtk is turned off, turn off adaptivity too.
    # See also https://chaste.cs.ox.ac.uk/trac/wiki/InstallAdaptivityLibrary
    use_adaptivity = int(prefs.get('use-adaptivity', False)) and use_vtk
    if use_adaptivity:
        other_includepaths.append(chaste_libs_path+'libadaptivity/include')
        other_libpaths.append(chaste_libs_path+'libadaptivity/lib')
        other_libraries.extend(['adaptivity', 'gfortran', 'gfortranbegin'])

    # Extra libraries for VTK output
    # This has to come after the 'if use_adaptivity' block, because the libraries there depend on these
    if use_vtk:
        other_includepaths.append(chaste_libs_path+'Vtk5/include/vtk-5.2')
        other_libpaths.append(chaste_libs_path+'Vtk5/lib/vtk-5.2')
        other_libraries.extend(['vtkFiltering', 'vtkIO', 'vtkCommon', 'vtksys', 'vtkzlib', 'vtkexpat', 'vtkGraphics'])
    
    # Chaste may also optionally link against CVODE.
    use_cvode = int(prefs.get('use-cvode', False))
    if use_cvode:
        other_includepaths.append(chaste_libs_path+'cvode/include')
        DetermineCvodeVersion(other_includepaths[-1])
        other_libpaths.append(chaste_libs_path+'cvode/lib')
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
