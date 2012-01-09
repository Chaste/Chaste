#!/usr/bin/env python


"""Copyright (C) University of Oxford, 2005-2012

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


# This script checks all the test folders to see if there are any
# Test*.hpp files that aren't listed in a test pack.
# It expects to be run from the trunk of the Chaste distribution.

import glob
import os
import re
import sys

sys.path[0:0] = ['python']
import BuildTools
set = BuildTools.set

chaste_dir = '.'

suite_res = {'.hpp': re.compile(r'class\s+(\w+)\s*:\s*public\s+((::)?\s*CxxTest\s*::\s*)?\w*TestSuite\s+$'),
             '.py': re.compile(r'class\s+\w+\(unittest\.TestCase\):\s+$')}

def IsTestFile(test_dir, test_file_path):
    """Does the given file define a test suite?"""
    is_test = False
    test_file = os.path.basename(test_file_path)
    test_ext = os.path.splitext(test_file)[1]
    if test_file[:4] == 'Test' and test_ext in suite_res.keys():
        fp = open(os.path.join(test_dir, test_file_path))
        for line in fp:
            m = suite_res[test_ext].match(line)
            if m:
                is_test = True
                break
        fp.close()
    #print test_dir, test_file, test_ext, is_test
    return is_test

test_packs  = set()  # Names of test packs found
orphans     = set()  # Names of any orphaned test files
found_tests = set()  # Names of tests found in test packs
test_dirs = glob.glob('*/test') + glob.glob('projects/*/test')

# First get a list of all tests in all test packs
for test_dir in test_dirs:
    tf, pn = BuildTools.GetTestsInTestPacks(test_dir, returnFoundPacks=True)
    found_tests.update(tf)
    test_packs.update(pn)

# Now check for orphaned tests in each top-level dir
local_found_tests = {} # Names of tests found in test packs in each folder

for test_dir in test_dirs:
    local_found_tests[test_dir] = BuildTools.GetTestsInTestPacks(test_dir)
    #print test_dir, local_found_tests[test_dir]
    # Check for orphans in this folder
    for dirpath, dirnames, filenames in os.walk(test_dir):
        for dirname in dirnames[:]:
            if dirname in ['.svn', 'data']:
                dirnames.remove(dirname)
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            filepath = BuildTools.relpath(filepath, test_dir)
            if IsTestFile(test_dir, filepath):
                if not filepath in local_found_tests[test_dir]:
                    orphans.add(os.path.join(test_dir, filepath))
                else:
                    local_found_tests[test_dir].remove(filepath)

# Output the names of test packs found
if test_packs:
    print "Test packs found:"
    for test_pack in test_packs:
        print "   ", test_pack
    print

# Compute a list of tests listed in test packs without .hpp files
not_found = []
for test_dir in local_found_tests.keys():
    for test_file in local_found_tests[test_dir]:
        not_found.append(test_dir + test_file)

# Display results
if orphans or not_found:
    if orphans:
        print "Orphaned tests found:"
        for orphan in sorted(orphans):
            print "   ", orphan
        print
    if not_found:
        print "Tests that don't exist:"
        for test in sorted(not_found):
            print "   ", test
        print
    print "The next line is for the benefit of the test summary scripts."
    n_orphans, n_found = len(orphans), len(found_tests)
    print "Failed", n_orphans, "of", n_orphans+n_found, "tests"

    # Return a non-zero exit code if problems were found
    sys.exit(n_orphans + len(not_found))
else:
    print "Infrastructure test passed ok."
  
