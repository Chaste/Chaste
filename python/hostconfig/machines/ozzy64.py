# Configuration for Ozzy's Comlab machine

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

petsc_2_2_path = None
petsc_2_3_path = '/home/ozzy/petsc-2.3.3-p15/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'
dealii_path = None
metis_path = None
intel_path = None #'/opt/intel/cce/10.0.025/'
#icpc = 'icpc -gcc-version=413 -I /usr/include/c++/4.1.3/x86_64-linux-gnu/ -I/usr/include/c++/4.1.3/'

ldflags = ' /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf '

other_includepaths = [] 
other_libpaths = ['/usr/lib', '/home/ozzy/lib']
blas_lapack = []
blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']
other_libraries = ['boost_serialization', 'boost_filesystem', 'xerces-c', 'hdf5']

use_cvode = False
if use_cvode:
    other_includepaths.append('/home/ozzy/cvode/include')
    other_libpaths.append('home/ozzy/cvode/lib')
    other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])

tools = {'texttest': '/home/chaste/texttest-3.10/source/bin/texttest.py',
         'mpirun': '/home/ozzy/mpi/bin/mpirun',
         'mpicxx': '/home/ozzy/mpi/bin/mpicxx',
         'tetgen': '/home/ozzy/bin'}

do_inf_tests = 1


