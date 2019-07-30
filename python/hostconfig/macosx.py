# This configuration is for Mac OSX Mountain Lion
# with package installed via Homebrew
# See https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/ChasteInstallationOnMountainLion

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

petsc_ver = 3.3
petsc_path='/usr/local/Cellar/petsc/3.3-p5/'
petsc_build_name = ''
petsc_build_name_profile = ''
petsc_build_name_optimized = ''

noccache = "true"

other_includepaths = ['/usr/local/opt/libxsd',
                      '/usr/local/include/']

other_libpaths = [ '/usr/X11/lib',
                  '/usr/local/lib/']
blas_lapack = []
other_libraries = ['X11', 'boost_serialization-mt', 'boost_filesystem-mt', 'boost_system-mt','xerces-c', 'z', 'hdf5', 'parmetis', 'metis']


# Location of Apple's versions of BLAS and LAPACK: -framework vecLib
ldflags='-framework vecLib'


def Configure(prefs, build):
    """Set up the build configuring.
    
    prefs can specify which version of various libraries we should use, and which optional libraries.
    
    build is an instance of BuildTypes.BuildType.
    """

    # VTK setup
    global use_vtk
    use_vtk = True
    if use_vtk:
        other_includepaths.append('/usr/local/include/vtk-5.10')
        other_libpaths.append('/usr/local/lib/vtk-5.10')
        other_libraries.extend(['vtkFiltering', 'vtkIO', 'vtkCommon', 'vtksys', 'vtkzlib', 'vtkexpat', 'vtkGraphics'])

    # CVODE setup
    global use_cvode
    use_cvode = True
    if use_cvode:
        DetermineCvodeVersion('/usr/local/include')
        other_libraries.extend(['sundials_cvode', 'sundials_nvecserial'])
