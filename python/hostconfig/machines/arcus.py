# Configuration for arcus.arc.ox.ac.uk

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

chaste_libs_path = '/system/software/recipes/chaste/3.2/lib'

# PETSc (assumes you have done module load)
petsc_path = os.environ.get("PETSC_DIR")
if petsc_path == None:
  print "*** WARNING: PETSC_DIR not set.  Please do \"module add PETSc\" for a compatible PETSc." 
petsc_build_name = ''
petsc_build_name_profile = ''
petsc_build_name_optimized = 'linux-intel-opt'
dealii_path = ''

# Intel compilers
intel_path = '/system/software/linux-x86_64/compilers/intel/intelCS-2013'
intel_bin_path = '/system/software/linux-x86_64/compilers/intel/intelCS-2013/bin'
#icpc = 'icpc'
#ccflags = '-gcc-version=440'

# other libraries
hdf5_path      = '/system/software/arcus/lib/hdf5/1.8.9_openmpi'
xerces_path    = '/system/software/linux-x86_64/xerces-c/3.3.1'
boost_path     = '/system/software/linux-x86_64/lib/boost/1_56_0'
# Suspicious
parmetis_path  = '/system/software/arcus/lib/parmetis/4.0.3/ompi-1.8.3__intel-2013'
xsd_path       = '/system/software/linux-x86_64/lib/xsd/3.3.0-1'
szip_path      = '/system/software/linux-x86_64/lib/szip/2.1'
vtk_path       = '/system/software/linux-x86_64/lib/vtk/5.10.1__python2.7'
# Suspicious
sundials_path = '/system/software/arcus/lib/sundials/openmpi-1.8.3/2.5.0/double'


other_includepaths = [intel_path + '/include',
                      petsc_path + '/include',
                      hdf5_path + '/include',
                      xerces_path + '/include',
                      boost_path + '/include',
                      xsd_path + '/include',
                      parmetis_path + '/include',
                      szip_path + '/include',
                      '/usr/include']

other_includepaths.extend(['/usr/include/c++/4.4.4/x86_64-redhat-linux', '/usr/include/c++/4.4.4',
                           '/usr/include/c++/4.4.4/backward', '/usr/lib/gcc/x86_64-redhat-linux/4.4.4/include'])

other_libpaths = [intel_path + '/lib/intel64',
                  intel_path + '/mkl/lib/intel64',
                  petsc_path + '/lib',
                  hdf5_path + '/lib',
                  xerces_path + '/lib',
                  boost_path + '/lib',
                  xsd_path + '/lib',
                  parmetis_path + '/lib',
                  szip_path + '/lib',
                  '/system/software/linux-x86_64/python/2.7.8/lib']

other_libraries = [ 'dmumps', 'mumps_common', 'mkl_scalapack_lp64', 'mkl_blacs_openmpi_lp64', 'pord','boost_serialization', 'boost_filesystem', 'boost_program_options', 'xerces-c', 'hdf5', 'sz', 'z', 'parmetis', 'metis', 'HYPRE']
blas_lapack = ['mkl_intel_lp64', 'mkl_core', 'mkl_sequential', 'pthread', 'ifcore']

# If using Boost >= 1.42, uncomment the following line
other_libraries.append('boost_system')

tools = {'mpirun': 'mpirun',
         'mpicxx': 'mpicxx',
         'xsd': xsd_path + '/bin/xsd'}

do_inf_tests=0

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
        other_includepaths.append(chaste_libs_path + '/libadaptivity/include')
        other_libpaths.append(chaste_libs_path + '/libadaptivity/lib')
        other_libraries.extend(['adaptivity', 'gfortran', 'gfortranbegin'])

    # Extra libraries for VTK output
    # This has to come after the 'if use_adaptivity' block, because the libraries there depend on these
    if use_vtk:
        other_includepaths.append(vtk_path + '/include/vtk-5.10')
        other_libpaths.append(vtk_path + '/lib/vtk-5.10')
        other_libraries.extend(['vtkGraphics','vtkFiltering','vtkIO','vtkCommon', 'vtksys', 'vtkexpat', 'vtkzlib'])
        # 5.8 was other_libraries.extend(['vtkFiltering', 'vtkIO', 'vtkCommon', 'vtksys', 'vtkzlib', 'vtkexpat', 'vtkGraphics'])
    
    # Chaste may also optionally link against CVODE.
    use_cvode = int(prefs.get('use-cvode', True))
    if use_cvode:
        other_includepaths.append(sundials_path + '/include/')
        DetermineCvodeVersion(other_includepaths[-1])
        other_libpaths.append(sundials_path + '/lib')
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
