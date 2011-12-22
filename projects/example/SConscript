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

# Example SCons build script for user projects.

import os

# Get variables exported from the SConstruct file
Import("*")

# Get our project directory name.
# This assumes that this project is located at <chasteRoot>/projects/<project>.
# Other things will probably go wrong if this is not the case...
project_name = os.path.basename(os.path.dirname(os.path.dirname(os.getcwd())))

# Determine Chaste libraries to link against.
# Note that order does matter!
# Select which line to uncomment based on what your project needs.
chaste_libs_used = comp_deps['core']
#chaste_libs_used = ['cell_based'] + comp_deps['cell_based']
#chaste_libs_used = ['heart'] + comp_deps['heart']
#chaste_libs_used = ['cell_based', 'heart'] + comp_deps['heart']

# Do the build magic
result = SConsTools.DoProjectSConscript(project_name, chaste_libs_used, globals())
Return("result")
