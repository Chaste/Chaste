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

"""
Add to your .bashrc:
#For the Intel compilers icc (C) and icpc (C++)
. /opt/intel/bin/compilervars.sh intel64

#For OpenMPI binary and library
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib/
export PATH=${PATH}:/usr/lib64/openmpi/bin/
"""

petsc_3_0_path = '/usr/local/petsc-3.0.0-p8/'
petsc_build_name = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
#petsc_build_name_production = 'linux-intel-opt-mkl'
#intel_path = '/opt/intel/cce/10.0.025/'
#icpc = 'icpc -gcc-version=413 -I /usr/include/c++/4.1.3/x86_64-linux-gnu/ -I/usr/include/c++/4.1.3/'

other_includepaths = ['/usr/local/boost/include']
other_libpaths = ['/usr/local/boost/lib', '/usr/lib64/openmpi/lib']

blas_lapack = ['flapack', 'fblas']
blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']
other_libraries = ['boost_serialization', 'boost_filesystem', 'xerces-c', 'z', 'hdf5', 'parmetis', 'metis', 'gfortran', 'HYPRE']

do_inf_tests = 1

use_cvode = True
if use_cvode:
  other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])

use_vtk = True
if use_vtk:
  other_includepaths.append('/usr/include/vtk')
  other_libpaths.append('/usr/lib64/vtk')
  other_libraries.extend(['vtkFiltering','vtkIO',  'vtkCommon', 'vtksys', 'vtkGraphics'])
