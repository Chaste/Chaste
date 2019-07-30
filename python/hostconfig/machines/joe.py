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
import hostconfig

petsc_3_0_path = '/home/jmpf/petsc-3.0.0-p8/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu-profile'
petsc_build_name_optimized = 'linux-gnu-opt'
petsc_build_name_production = 'linux-intel-opt-mkl'
intel_path = '/opt/intel/composerxe/'
#icpc = 'icpc -gcc-version=413 -I /usr/include/c++/4.1.3/x86_64-linux-gnu/ -I/usr/include/c++/4.1.3/'

other_includepaths = ['/home/jmpf/boost/include']
other_libpaths = ['/home/jmpf/boost/lib']
#other_includepaths = ['../../../xsd-2.3.1-i686-linux-gnu/libxsd', 
#                     '../../../hdf5/include', 
#                     os.path.join(petsc_2_3_path, 'externalpackages/hypre-1.11.1b/linux-gnu/include/'),
#                  parmetis_path]
#other_libpaths = [os.path.join(petsc_2_3_path, 'externalpackages/f2cblaslapack/linux-gnu/'),  
#               '/opt/intel/mkl/9.1.023/lib/em64t',
#               '../../../hdf5/lib',
#               os.path.join(petsc_2_3_path, 'externalpackages/hypre-1.11.1b/linux-gnu/lib/'),
#               parmetis_path]


#blas_lapack_production = ['mkl_lapack', 'mkl', 'svml']
blas_lapack_production = ['mkl_solver_lp64_sequential', 'mkl_intel_lp64', 'mkl_sequential', 'mkl_core']
other_libraries = ['boost_serialization', 'boost_filesystem', 'boost_system', 'xerces-c', 'z', 'hdf5', 'parmetis', 'metis', 'HYPRE']
blas_lapack = ['flapack', 'fblas', 'gfortran']

tools = {'xsd': 'xsdcxx'}

do_inf_tests = 1

def Configure(prefs, build):
    """Set up the build configuring.
       
    prefs can specify which version of various libraries we should use, and which optional libraries.
        
    build is an instance of BuildTypes.BuildType.
    """
    global use_cvode
    global use_vtk
                                
    use_cvode = True
    if use_cvode:
      DetermineCvodeVersion('/usr/include')
      other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
  
    use_vtk = True
    if use_vtk:
      other_includepaths.append('/usr/include/vtk-5.8')
      other_libpaths.append('/usr/lib64/vtk-5.8')
      other_libraries.extend(['vtkFiltering','vtkIO',  'vtkCommon', 'vtksys', 'vtkGraphics'])
        



""" Some instructions on what produced my configuration:

#Boost 1.42 (Boost 1.46 is installed under Chaste dependencies but
#1.There's the segfault issue
#2.There's a problem with the Intel (12) compiler and 1.46
wget http://sourceforge.net/projects/boost/files/boost/1.42.0/boost_1_42_0.tar.gz
tar xvfz boost_1_42_0.tar.gz 
cd boost_1_42_0/
bjam "-sTOOLS=gcc" --prefix=$HOME/boost install

#PETSc etc.
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.0.0-p8.tar.gz
tar xzvf petsc-3.0.0-p8.tar.gz 
rm petsc-3.0.0-p8.tar.gz
cd petsc-3.0.0-p8/
export PETSC_DIR=`pwd`

export PETSC_ARCH=linux-gnu
./config/configure.py  --download-f-blas-lapack=yes --with-mpi-dir=/usr/lib64/openmpi \
--download-parmetis=yes --download-hypre=yes -with-x=false --with-clanguage=cxx
make all
export PETSC_ARCH=linux-gnu-opt
./config/configure.py  --download-f-blas-lapack=yes --with-mpi-dir=/usr/lib64/openmpi \
--download-parmetis=yes --download-hypre=yes -with-x=false --with-clanguage=cxx --with-debugging=0
make all

export PETSC_ARCH=linux-intel
./config/configure.py --with-cxx=icpc --with-cc=icc --with-vendor-compiler=intel  \
--download-c-blas-lapack=yes \
--with-fortran=0 --with-x=false --with-clanguage=cxx --download-openmpi=yes \
--CXXFLAGS=-fPIC --CFLAGS=-fPIC --download-hdf5=yes --download-parmetis=yes  --with-shared=0
make all

#The following is a LIE --with-blas-lapack-dir cannot cope with the newer Intel MKL names.
#We have to have a system blas installed and then hope for the best!
export PETSC_ARCH=linux-intel-opt-mkl
./config/configure.py --with-cxx=icpc --with-cc=icc --with-vendor-compiler=intel \
--with-blas-lapack-dir=/opt/intel/composerxe/mkl/lib/intel64 \
--with-fortran=0 --with-x=false --with-clanguage=cxx --download-openmpi=yes \
--CXXFLAGS=-fPIC --CFLAGS=-fPIC --with-debugging=0 --download-hdf5=yes --download-parmetis=yes  --with-shared=0
make all
#Cannot do Fortran/HYPRE in the prescence of the Intel compiler yet...
#--with-fortran=1 --download-hypre=yes 

"""

