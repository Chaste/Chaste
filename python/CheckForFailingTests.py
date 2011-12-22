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

"""
Check whether there are any tests listed in the Failing test pack.

This is intended to be run as a script from the root of the Chaste
distribution.
"""

import glob
import sys

sys.path[0:0] = ['python']
import BuildTools

failing_tests = BuildTools.set()
for test_dir in glob.glob('*/test') + glob.glob('projects/*/test'):
    failing_tests.update(BuildTools.GetTestsInTestPacks(test_dir, ['Failing']))

# Display results
if failing_tests:
    print "Failing tests found:"
    for test in sorted(failing_tests):
        print "   ", test
    print "The next line is for the benefit of the test summary scripts."
    print "Failed", len(failing_tests), "of", len(failing_tests), "tests"
else:
    print "Infrastructure test passed ok."
