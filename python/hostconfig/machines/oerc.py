# Configuration for HAL or SAL

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


petsc_3_0_path = '/home/system/software/hal/lib/PETSc/petsc-3.0.0-p12/icc-2011/'
petsc_2_3_path = ''
petsc_build_name = ''
petsc_build_name_profile = ''
petsc_build_name_optimized = 'linux-intel-opt'
dealii_path = ''
parmetis_path = '/home/system/software/hal/lib/parmetis/3.1.1/sgimpt__intel-11.1/'
# If you have the Intel compiler installed, set this to the folder where it lives
intel_path = '/system/software/linux-x86_64/compilers/intel/intelCS-2013/'
intel_bin_path= '/system/software/linux-x86_64/compilers/intel/intelCS-2013/bin'

xsd_path = '/system/software/hal/lib/xsd/3.3.0-1/'

other_includepaths = ['/system/software/linux-x86_64/compilers/intel/intelCS-2013/mkl/include/intel64/',
                      '/system/software/linux-x86_64/compilers/intel/intelCS-2013/include/intel64/',
                      '/system/software/hal/lib/hdf5/hdf5-1.6.9_parallel/include',
                      '/home/system/software/redqueen/libs/szip-2.1/include',
                      '/system/software/linux-x86_64/xerces-c/3.3.1/include',
                      '/system/software/redqueen/libs/boost-1_45_0/include',
                      xsd_path + '/include',
                      parmetis_path + 'include',
                      petsc_3_0_path + 'include'] 

other_libpaths = ['/usr/lib64',
                  '/system/software/linux-x86_64/compilers/intel/intelCS-2013/mkl/lib/intel64/',
                  '/system/software/linux-x86_64/compilers/intel/intelCS-2013/lib/intel64/',
                  '/system/software/redqueen/libs/boost-1_45_0/lib',
                  '/system/software/linux-x86_64/xerces-c/3.3.1/lib',
                  '/system/software/hal/lib/hdf5/hdf5-1.6.9_parallel/lib',
                  '/home/system/software/redqueen/libs/szip-2.1/lib',
                  '/system/software/hal/lib/python2.6/site-packages/rdflib',
                  os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu'),
                  parmetis_path + 'lib',
                  petsc_3_0_path + 'lib'] # own inclusion


# The order of libraries in these lists matters!
other_libraries = [ 'dmumps', 'mumps_common', 'scalapack', 'blacs', 'pord', 'mpi++','boost_serialization', 'boost_filesystem', 'boost_program_options', 'xerces-c', 'hdf5', 'sz', 'z', 'parmetis', 'metis']
if os.path.exists(petsc_3_0_path):
    # Assume user has followed INSTALLATION.txt and has Fortran blas & lapack
    blas_lapack = ['mkl_intel_lp64',
                   'mkl_intel_thread','mkl_core','iomp5','pthread', 'ifcore']
    other_libraries.append('HYPRE')
else:
    blas_lapack = ['f2clapack', 'f2cblas'] # Note: Lapack before BLAS
  # Note: parmetis before metis, hdf5 before z.
# Note that boost serialization sometimes has a different name:
# other_libraries = ['boost_serialization-gcc41', 'xerces-c', 'hdf5', 'z', 'parmetis', 'metis']

blas_lapack_production = ['mkl_solver_lp64_sequential', 'mkl_intel_lp64', 'mkl_sequential', 'mkl_core']

tools = {'mpirun': '/usr/bin/mpirun',
         'mpicxx': 'icpc' + ' -static-intel -wr191 -wr304 -wr981 -wr383 -wr1419  -wr82 -wr2304' + ' -wr2026' + ' -lmpi -lmpi++',
         'xsd': xsd_path+'bin/xsd'}


def Configure(prefs, build):
    """Set up the build configuring.
    
    prefs can specify which version of various libraries we should use, and which optional libraries.
    
    build is an instance of BuildTypes.BuildType.
    """
    global use_cvode
    global use_vtk
    
    # use_vtk defaults to false. Change to True if VTK development libraries are available.
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
        #Todo VTK is not installed yet
        other_includepaths.append('/system/software/linux-x86_64/lib/vtk/5.10.1/include/vtk-5.10')
        other_libpaths.append('/system/software/linux-x86_64/lib/vtk/5.10.1/lib/vtk-5.10')
        other_libraries.extend(['vtkFiltering', 'vtkIO', 'vtkCommon', 'vtksys', 'vtkzlib', 'vtkexpat', 'vtkGraphics'])
    
    # Chaste may also optionally link against CVODE.
    use_cvode = int(prefs.get('use-cvode', False))
    if use_cvode:
        #Todo CVODE is not installed yet
        other_includepaths.append(chaste_libs_path+'cvode/include')
        DetermineCvodeVersion(other_includepaths[-1])
        other_libpaths.append(chaste_libs_path+'cvode/lib')
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
