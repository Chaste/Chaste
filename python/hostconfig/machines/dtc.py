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

petsc_2_3_path = '/usr/apps/chaste-plugins/petsc-2.3.2-p4/'
petsc_build_name = 'linux-gnu'
icpc='icpc'

#other_includepaths = ['../../../xsd-2.3.1-i686-linux-gnu/libxsd', '/usr/apps/chaste-plugins/hdf5/include']
other_includepaths = ['/usr/apps/chaste-plugins/hdf5/include']
other_libpaths = [os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu/'),
                  '/usr/apps/chaste-plugins/hdf5/lib']
blas_lapack = ['f2clapack', 'f2cblas']
other_libraries = ['boost_serialization', 'boost_filesystem', 'xerces-c', 'z', 'hdf5']

tools = {'mpicxx': '/usr/apps/chaste-plugins/mpi/bin/mpicxx',
         'mpirun': '/usr/apps/chaste-plugins/mpi/bin/mpirun'}

#These are unused - possibly set erroneously
petsc_2_2_path = '/usr/petsc-2.2.1/'
petsc_build_name_profile = 'linux-gnu-profile'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_production = 'linux-intel-opt-mkl'
dealii_path = '../../../deal.II/'
metis_path = '../../../metis-4.0/'
intel_path = '/opt/intel/cc/9.1.039'
blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']
#tools = {'texttest': '/home/chaste/texttest-3.10/source/bin/texttest.py'}
