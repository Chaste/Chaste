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

chaste_libs = os.environ.get("CHASTE_LIBS")
compiler = os.environ.get("PE_ENV")
if compiler != 'INTEL':
    print "Only the Intel compiler is supported. Load PrgEnv-intel module."
    sys.exit(1)

# The modules system deals with most include/lib paths transparently.
other_includepaths = [os.path.join(chaste_libs, 'xsd/libxsd')]
other_libpaths = []

# Provided by cray-libsci
blas_lapack = []
blas_lapack_production = []

other_libraries = ['boost_serialization', 'boost_system', 'boost_filesystem', 'hdf5', 'xerces-c', 'z']

use_vtk = True
if use_vtk:
    other_libraries.extend(['vtkGraphics', 'vtkFiltering', 'vtkIO', 'vtkCommon', 'vtksys', 'vtkexpat', 'vtkzlib'])

tools = {'mpicxx': 'CC',
         'xsd': os.path.join(chaste_libs, 'bin/xsd')}

def Configure(prefs, build):
    global use_cvode
    use_cvode = int(prefs.get('use-cvode', 1))
    if use_cvode:
        DetermineCvodeVersion(os.path.join(os.environ.get('CRAY_TPSL_PREFIX_DIR'), 'include'))
