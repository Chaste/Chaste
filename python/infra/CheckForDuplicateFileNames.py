#!/usr/bin/env python


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
