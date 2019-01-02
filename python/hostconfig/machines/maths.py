# Configuration for Oxford's Maths Institute

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

#EDIT HERE
#For a simple installation all paths will be below this directory
chaste_libs_path = '/scratch/chaste/'
#EDIT HERE

if not os.path.exists(chaste_libs_path) or not os.path.isdir(chaste_libs_path):
    print >>sys.stderr, "Chaste dependencies folder", chaste_libs_path, \
        "not found; please edit python/hostconfig/default.py"
    sys.exit(1)

petsc_2_2_path = None
petsc_2_3_path = None
petsc_3_0_path = None
petsc_path = '/usr/lib/petscdir/3.4.2/'
petsc_build_name = 'linux-gnu-c-opt'
petsc_build_name_profile = 'linux-gnu-c-opt'
petsc_build_name_optimized = 'linux-gnu-c-opt'

parmetis_path = chaste_libs_path+'ParMetis-3.1/'
intel_path = chaste_libs_path+'intel/cc/9.1.039/'
icpc = 'icpc'

other_includepaths = [chaste_libs_path+'hdf5/include',
                      chaste_libs_path+'include',
                      chaste_libs_path+'boost/include/boost-1_34_1',
                      chaste_libs_path+'xsd-2.3.1-i686-linux-gnu/libxsd',
                      parmetis_path,
                      chaste_libs_path+'metis/include']

other_libpaths = [chaste_libs_path+'lib',
                  chaste_libs_path+'xerces/lib',
                  chaste_libs_path+'hdf5/lib',
                  parmetis_path,
                  '/usr/lib',
                  chaste_libs_path+'metis/lib']

blas_lapack = ['lapack', 'blas']
other_libraries = ['boost_serialization','boost_system', 'boost_filesystem', 'xerces-c', 'z', 'hdf5', 'parmetis', 'metis']

# use_vtk set to false initially. Change to True if VTK development libraries are available.
use_vtk = True

tools = {'xsd': chaste_libs_path+'xsd-2.3.1-i686-linux-gnu/bin/xsd'}

#Extra libraries for VTK output
if use_vtk:
    other_libraries.extend(['vtkGraphics', 'vtkFiltering', 'vtkIO', 'vtkCommon', 'z'])
    other_includepaths.extend(['/usr/include/vtk-5.8'])

use_cvode = True
if use_cvode:
    other_includepaths.append('/scratch/chaste/cvode/include')
    other_libpaths.append('/scratch/chaste/cvode/lib')
    other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
ccflags = ''
