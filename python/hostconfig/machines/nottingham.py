# Configuration for James Preston in Nottingham

"""Copyright (c) 2005-2018, University of Oxford.
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

# If you are an active developer, committing back to the trunk, 
# please uncomment the following line to run tests on duplicate file names, 
# orphaned tests and copyright notices.
do_inf_tests = 1

#EDIT HERE
#For a simple installation all paths will be below this directory
chaste_libs_path = '/local/pmxjp6/chaste_libraries/'
#EDIT HERE

if not os.path.exists(chaste_libs_path) or not os.path.isdir(chaste_libs_path):
    print >>sys.stderr, "Chaste dependencies folder", chaste_libs_path, \
        "not found; please edit python/hostconfig/local.py"
    sys.exit(1)

petsc_ver = 3.4
petsc_path = chaste_libs_path+'petsc-3.7.4/'
petsc_build_name = 'linux-gnu'
petsc_build_name_profile = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'

# If you have the Intel compiler installed, set this to the folder where it lives
#intel_path = '/opt/intel/cc/9.1.039/'
# You may need to edit this to ensure that the intel compiler finds the right gcc libraries, e.g.
#icpc = 'icpc -gcc-version=410 -I /usr/include/c++/4.1.3/x86_64-linux-gnu/ -I/usr/include/c++/4.1.3/    -I/usr/include/c++/4.1.3/backward'
#icpc = 'icpc'

if os.uname()[4] == 'x86_64':
    xsd_path = chaste_libs_path + 'xsd-3.3.0-x86_64-linux-gnu/'
else:
    xsd_path = chaste_libs_path + 'xsd-3.3.0-i686-linux-gnu/'

other_includepaths = [chaste_libs_path+'include',
                      xsd_path + 'libxsd']

other_libpaths = [chaste_libs_path+'lib']

# The order of libraries in these lists matters!
# Note that boost serialization sometimes has a different name: eg boost_serialization-gcc41
other_libraries = ['boost_serialization', 'boost_filesystem', 'boost_system', 'xerces-c', 'hdf5', 'parmetis', 'metis']
# Note: parmetis before metis, hdf5 before z.

# Assuming you installed PETSc with HYPRE support, you'll want this line:
other_libraries.append('HYPRE')

blas_lapack = ['f2clapack', 'f2cblas']

tools = {'xsd': xsd_path+'bin/xsd'}


def Configure(prefs, build):
    """Set up the build configuring.
    
    prefs can specify which version of various libraries we should use, and which optional libraries.
    
    build is an instance of BuildTypes.BuildType.
    """
    global use_cvode
    global use_vtk
    
    # use_vtk defaults to True. Change to False if VTK development libraries are not available.
    use_vtk = int(prefs.get('use-vtk', True))
    
    # Extra libraries for VTK output
    if use_vtk:
        other_includepaths.append(chaste_libs_path + '/include/vtk-5.10')
        other_libpaths.append(chaste_libs_path + '/lib/vtk-5.10')

        # VTK 5:
        #other_libraries.extend(['vtkFiltering', 'vtkIO', 'vtkCommon', 'vtksys', 'vtkzlib', 'vtkexpat', 'vtkGraphics'])
        
        # VTK 6:
        other_libraries.extend(['vtkGraphics','vtkFiltering','vtkIO','vtkCommon', 'vtksys', 'vtkexpat', 'vtkzlib'])
    
    # Chaste may also optionally link against CVODE.
    use_cvode = int(prefs.get('use-cvode', True))
    if use_cvode:
        other_includepaths.append(petsc_path + 'linux-gnu/include/')
        DetermineCvodeVersion(other_includepaths[-1])
        other_libpaths.append(chaste_libs_path+'cvode/lib')
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
