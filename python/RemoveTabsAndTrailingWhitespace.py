#!/usr/bin/env python

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

"""Script to 
i) replace tabs with 4 spaces - calls "sed -i 's/\t/    /g'" on all source files.
ii) remove any trailing whitespace - calls " sed -i 's/[ \\t]*$//' "
"""


import os, sys

exts = ['.cpp', '.hpp']
dir_ignores = ['build', 'cxxtest', 'testoutput', 'docs', 'doxygen', 'projects', 'data']
tab_spaces = ' ' * 4
chaste_dir = '.'
    
for root, dirs, files in os.walk(chaste_dir):
    for dir in dir_ignores:
        if dir in dirs:
            dirs.remove(dir)
    # Check for source files
    for file in files:
        name, ext = os.path.splitext(file)
        if ext in exts:
            file_name = os.path.join(root, file)
            command = "sed -i 's/\\t/%s/g' %s" % (tab_spaces, file_name)
            print "Checking " + file_name
            os.system(command)
            ### for removing trailing whitespace
            command = " sed -i 's/[ \\t]*$//' " + file_name
            os.system(command)
