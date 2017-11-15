"""Copyright (C) University of Oxford, 2005-2013

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

import os

petsc_ver = 3.7
petsc_path='/usr/local/Cellar/petsc/3.7.4_1'
petsc_build_name = ''
petsc_build_name_profile = ''
petsc_build_name_optimized = ''

noccache = "true"

other_includepaths = ['/usr/local/opt/libxsd',
                      '/usr/local/include/']

other_libpaths = [ '/usr/X11/lib',
                  '/usr/local/lib/']

blas_lapack = []
other_libraries = ['X11', 'boost_serialization-mt', 'boost_filesystem-mt', 'boost_system-mt','xerces-c', 'z', 'hdf5', 'parmetis','metis']


ldflags='-framework Accelerate'


def Configure(prefs, build):
    """Set up the build configuring.
    
    prefs can specify which version of various libraries we should use, and which optional libraries.
    
    build is an instance of BuildTypes.BuildType.
    """

    # VTK setup
    global use_vtk
    use_vtk = True
    if use_vtk:
        #homebrew installs vtk5 without simlinks
        other_includepaths.append('/usr/local/opt/vtk5/include/vtk-5.10/')
        other_libpaths.append('/usr/local/opt/vtk5/lib/vtk-5.10')
        other_libraries.extend(['vtkFiltering', 'vtkIO', 'vtkCommon', 'vtksys', 'vtkGraphics'])

    # CVODE setup
    global use_cvode
    use_cvode = True
    if use_cvode:
        DetermineCvodeVersion('/usr/local/include')
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
