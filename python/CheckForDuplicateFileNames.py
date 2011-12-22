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


# Check through all the Chaste source directories, noting any source files
# that have the same name.

exts = ['.cpp', '.hpp']
dir_ignores = ['build', 'dynamic', 'data']
startchar_ignores = ['_', '.']
chaste_dir = '.'

import glob
import os

# Dictionary mapping file names to locations
source_files = {}

def DoWalk(root_dir):
    for root, dirs, files in os.walk(root_dir):
        # Check for ignored dirs
        for dirname in dirs[:]:
            if dirname in dir_ignores or dirname[0] in startchar_ignores:
                dirs.remove(dirname)
        # Check for source files
        for file in files:
            name, ext = os.path.splitext(file)
            if ext in exts:
                if source_files.has_key(file):
                    # We've already found a file with this name
                    source_files[file].append(os.path.join(root, file))
                else:
                    # This is the first occurence of this name
                    source_files[file] = [os.path.join(root, file)]

components = os.path.join(chaste_dir, '*', '')
projects = os.path.join(chaste_dir, 'projects', '*', '')
for root_dir in (glob.glob(components + 'src') + glob.glob(components + 'test')
                 + glob.glob(projects + 'src') + glob.glob(projects + 'test')):
    DoWalk(root_dir)

# Now check dictionary for duplicates
num_found_dups = 0
for file in source_files:
    if len(source_files[file]) > 1:
        print "Duplicate occurrences of", file, ":"
        for loc in source_files[file]:
            print "   ", loc
        num_found_dups += 1

# Let the test summary script know
if num_found_dups > 0:
    print
    print "The next line is for the benefit of the test summary scripts."
    print "Failed", num_found_dups, "of", len(source_files), "tests"

    # Return a non-zero exit code if orphans were found
    import sys
    sys.exit(num_found_dups)
else:
    print "Infrastructure test passed ok."
