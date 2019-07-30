# Configuration for Joe's machines

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


petsc_2_2_path = None
petsc_2_3_path = '../../petsc-2.3.2-p10/'
petsc_3_0_path = None
petsc_build_name = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_production = 'linux-intel-opt-mkl'
dealii_path = None
parmetis_path = '../../ParMetis-3.1'
intel_path = ''
icpc = 'icpc'


other_includepaths = ['../../hdf5/include', 
'../../xsd-2.3.1-i686-linux-gnu/libxsd', parmetis_path]
other_libpaths = [os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu-opt/'),  
                    '../../hdf5/lib', parmetis_path]
blas_lapack = ['f2clapack', 'f2cblas']
blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']
other_libraries = ['boost_serialization', 'boost_filesystem', 'xerces-c', 'hdf5', 'z', 'parmetis', 'metis']

tools = {'mpirun': '../../mpi/bin/mpirun',
         'mpicxx': '../../mpi/bin/mpicxx',
            'xsd': '../../xsd-2.3.1-i686-linux-gnu/bin/xsd'}


use_vtk = True
if use_vtk:
    other_libraries.extend(['vtkGraphics', 'vtkFiltering', 'vtkIO', 'vtkCommon',  'z'])
    other_includepaths.extend(['/usr/include/vtk-5.0/'])

do_inf_tests = 1

    
