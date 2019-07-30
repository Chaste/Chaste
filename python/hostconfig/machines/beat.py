# Configuration for beat.cs.ox.ac.uk, see
# https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/Beat

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

petsc_path = '/opt/petsc-3.5.2'
petsc_build_name = 'linux-gnu'
petsc_build_name_optimized = 'linux-gnu-opt'

use_cvode = True
use_vtk = True
other_includepaths = ['/usr/include/vtk/',
                      '/opt/xsd-3.3.0-x86_64-linux-gnu/libxsd']
other_libpaths = ['/usr/lib64/vtk/']

blas_lapack = ['flapack', 'fblas']

other_libraries = ['boost_system', 'boost_serialization', 'boost_filesystem', 'xerces-c', 'hdf5', 'z', 'parmetis', 'metis']
other_libraries.append(['vtkCommonCore','vtkCommonDataModel','vtkIOXML','vtkCommonExecutionModel','vtkFiltersCore','vtkFiltersGeometry','vtkFiltersModeling','vtkFiltersSources'])
other_libraries.append(['sundials_cvode', 'sundials_nvecserial'])

tools = {'xsd': '/opt/xsd-3.3.0-x86_64-linux-gnu/bin/xsd'}

def Configure(prefs, build):
    DetermineCvodeVersion(os.path.join(petsc_path,petsc_build_name,'include'))

