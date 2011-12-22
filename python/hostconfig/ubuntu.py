# Configuration

"""Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.
"""

import glob
import os

# Check which version of Ubuntu this is
fp = open('/etc/issue')
ubuntu_ver = fp.read().split()[1]
fp.close()
#First deal with special cases for beta releases etc.
if ubuntu_ver == 'natty':
    ubuntu_ver = [11,04]
else:
    ubuntu_ver = map(int, ubuntu_ver.split('.')[0:2]) 

if ubuntu_ver >= [10,10]:
    petsc_ver = 3.1
    petsc_path = '/usr/lib/petscdir/3.1/'
elif ubuntu_ver >= [9,10]:
    petsc_ver = 3
    petsc_3_0_path = '/usr/lib/petscdir/3.0.0/'
else:
    petsc_ver = 2
    petsc_2_3_path = '/usr/lib/petscdir/2.3.3/'

petsc_2_2_path = ''
petsc_build_name = 'linux-gnu-c-debug'
petsc_build_name_profile = petsc_build_name
petsc_build_name_optimized = 'linux-gnu-c-opt'

dealii_path = None
intel_path = None
icpc = 'icpc'

other_includepaths = ['/usr/include/metis/']
other_libpaths = []
libs_for_petsc = ['petsccontrib', 'X11',
                  'HYPRE', 'spooles', 'superlu',
                  'umfpack', 'amd' # Both for Umfpack
                  ]
#Fixes (possibly temporary) for Natty
if ubuntu_ver >= [11,04]:
    libs_for_petsc.append(['HYPRE_utilities', 
		'HYPRE_struct_mv', 'HYPRE_struct_ls',  
		'HYPRE_sstruct_mv', 'HYPRE_sstruct_ls', 
		'HYPRE_IJ_mv', 'HYPRE_parcsr_ls', 'dmumps'])
if petsc_ver >= 3:
    libs_for_petsc.append('scotch')
else:
    libs_for_petsc.append('sidl')
if petsc_ver >= 3.1:
    libs_for_petsc.remove('petsccontrib')
if ubuntu_ver >= [9,10]:
    boost_suffix = '-mt'
else:
    boost_suffix = ''

other_libraries = libs_for_petsc + \
                  ['boost_serialization'+boost_suffix, 'xerces-c',
                   'hdf5', 'z',
                   'parmetis', 'metis']

# Figure out which lapack/blas packages are actually installed!
if os.path.exists('/usr/lib/liblapack-3.so'):
    blas_lapack = ['lapack-3', 'blas-3']
else:
    blas_lapack = ['lapack', 'blas']

tools = {'xsd': '/usr/bin/xsdcxx',
         'mpirun': '/usr/bin/mpirun.openmpi',
         'mpicxx': '/usr/bin/mpic++.openmpi'}

def Configure(prefs, build):
    """Set up the build configuring.
    
    prefs can specify which version of various libraries we should use, and which optional libraries.
    VTK and CVODE support default on if they are installed.
    
    build is an instance of BuildTypes.BuildType.
    """
    global use_cvode
    global use_vtk
    
    # Extra libraries for VTK output
    vtk_include_path = filter(os.path.isdir, glob.glob('/usr/include/vtk-5*'))
    use_vtk = int(prefs.get('use-vtk', True))
    use_vtk = use_vtk and bool(vtk_include_path)
    if use_vtk:
        # Note: 10.10 uses VTK 5.4, 10.04 uses 5.2, and early use 5.0
        other_includepaths.extend(vtk_include_path)
        other_libraries.extend(['vtkIO', 'vtkCommon', 'vtkGraphics', 'z'])
	if ubuntu_ver >= [11,10]: # 11.10 uses VTK 5.6
		other_libraries.extend(['vtkFiltering'])

    # Is CVODE installed?
    use_cvode = int(prefs.get('use-cvode', True))
    use_cvode = use_cvode and os.path.exists('/usr/lib/libsundials_cvode.so')
    if ubuntu_ver <= [9,04]:
        # We don't support CVODE 2.4
        use_cvode = False
    if use_cvode:
        DetermineCvodeVersion('/usr/include')
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
