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


# This script checks all the test folders to see if there are any
# Test*.hpp files that aren't listed in a test pack.
# It expects to be run from the trunk of the Chaste distribution.

import glob
import os
import re
import sys

sys.path[0:0] = ['python/infra']
import BuildTools
set = BuildTools.set

chaste_dir = '.'

suite_res = {'.hpp': re.compile(r'class\s+(\w+)\s*:\s*public\s+((::)?\s*CxxTest\s*::\s*)?\w*TestSuite\s+$'),
             '.py': re.compile(r'class\s+\w+\(unittest\.TestCase\):\s+$')}

# Should we check for orphans in projects too?
projects_to_check = sys.argv[1:]

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
parallel_not_tested_continuous_too = [];
test_dirs = glob.glob('*/test')
test_dirs.extend(map(lambda p: os.path.join(p, 'test'), projects_to_check))

# First get a list of all tests in all test packs
print "Test folders checked:"
for test_dir in sorted(test_dirs):
    print "   ", test_dir
    tf, pn = BuildTools.GetTestsInTestPacks(test_dir, returnFoundPacks=True)
    found_tests.update(tf)
    test_packs.update(pn)
    
    # Make a list of parallel and continuous tests
    parallel = BuildTools.GetTestsInTestPacks(test_dir,['Parallel'])
    continuous = BuildTools.GetTestsInTestPacks(test_dir,['Continuous'])
    
    # Check that parallel tests are a subset of continuous (coverage just runs continuous in parallel!)
    if not parallel.issubset(continuous):
        not_tested_continuous_too = parallel.difference(continuous)
        for particular_test in not_tested_continuous_too:
            parallel_not_tested_continuous_too.append(particular_test)

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
    for test_pack in sorted(test_packs):
        print "   ", test_pack
    print
    
# Compute a list of tests listed in test packs without .hpp files
not_found = []
for test_dir in local_found_tests.keys():
    for test_file in local_found_tests[test_dir]:
        not_found.append(os.path.join(test_dir,test_file))

# Display results
if orphans or not_found or parallel_not_tested_continuous_too:
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
    if parallel_not_tested_continuous_too:
        print "Parallel tests that need to be in a continuous test pack too:"
        for test in sorted(parallel_not_tested_continuous_too):
            print "   ", test
        print
    print "The next line is for the benefit of the test summary scripts."
    n_orphans, n_found = len(orphans)+len(parallel_not_tested_continuous_too), len(found_tests),
    print "Failed", n_orphans, "of", n_orphans+n_found, "tests"

    # Return a non-zero exit code if problems were found
    sys.exit(n_orphans + len(not_found) + len(parallel_not_tested_continuous_too))
else:
    print "Infrastructure test passed ok."
  
