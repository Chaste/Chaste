# Configuration for core Chaste machines

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

chaste_libs_dir = '/home/scratch/chaste-libs'

petsc_3_0_path = os.path.join(chaste_libs_dir, 'petsc-3.0.0-p9')
petsc_path = os.path.join(chaste_libs_dir, 'petsc-3.2-p7')
petsc_build_name = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_profile = 'linux-gnu-profile'

#intel_path = '/opt/intel/cc/9.1.039/lib'
#icpc='icpc'

other_includepaths = [os.path.join(chaste_libs_dir, 'include')]
other_libpaths = [os.path.join(chaste_libs_dir, 'lib'),
                  '/usr/lib64/openmpi/lib']
blas_lapack = ['flapack', 'fblas', 'gfortran', 'gfortranbegin']
other_libraries = ['boost_serialization', 'boost_filesystem', 'xerces-c', 'hdf5', 'z', 'parmetis', 'metis', 'HYPRE']

use_cvode = True
if use_cvode:
    other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])

tools = {'mpicxx': '/usr/lib64/openmpi/bin/mpicxx',
         'mpirun': '/usr/lib64/openmpi/bin/mpirun',
         'xsd': 'xsdcxx',
         'texttest': os.path.join(chaste_libs_dir, 'texttest-3.19/source/bin/texttest.py')}

do_inf_tests = 1

