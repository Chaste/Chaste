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

# If you are an active developer, committing back to the trunk, 
# please uncomment the following line to run tests on duplicate file names, 
# orphaned tests and copyright notices.
do_inf_tests = 1

chaste_libs_path = '/usr/chaste_deps/'

petsc_2_2_path = None
petsc_2_3_path = None
petsc_3_0_path = chaste_libs_path+'petsc-3.0.0-p8/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_production = 'linux-intel-opt'
intel_path = '/usr/intel'

other_includepaths = [chaste_libs_path+'xerces/include',
                      '/users/nejzem/boost_1_37/include/boost-1_37',
                      chaste_libs_path+'xsd-3.2.0-i686-linux-gnu/libxsd',
                      chaste_libs_path+'petsc-3.0.0-p8/include']

other_libpaths = [chaste_libs_path+'lib',
                  chaste_libs_path+'rdf/lib',
                  '/users/nejzem/boost_1_37/lib', 
                  chaste_libs_path+'xerces/lib',
                 '/usr/intel/mkl/lib/em64t',
                 '/usr/intel/lib/intel64']

blas_lapack = ['flapack', 'fblas']
blas_lapack_production = ['mkl_solver_lp64_sequential', 'mkl_intel_lp64', 'mkl_sequential', 'mkl_core']
other_libraries = ['boost_serialization-gcc43-mt', 'boost_filesystem-gcc43-mt', 'xerces-c', 'hdf5', 'z', 'libgfortran', 'parmetis', 'metis', 'HYPRE']

tools = {'mpirun': chaste_libs_path+'mpi-intel/bin/mpirun', # this will be only used for IntelProduction, since other PETSc builds have OpenMPI built in.
         'mpicxx': chaste_libs_path+'mpi-intel/bin/mpicxx',# this will be only used for IntelProduction, since other PETSc builds have OpenMPI built in.
         'xsd': chaste_libs_path+'xsd-3.2.0-i686-linux-gnu/bin/xsd',
         'texttest': ''}

use_vtk = True
if use_vtk:
    other_includepaths.append('/usr/include/vtk')
    other_libraries.extend(['vtkIO','vtkGraphics','vtkCommon', 'z'])
