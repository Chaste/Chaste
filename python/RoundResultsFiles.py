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

"""
This script is used as a filter on acceptance test results files, to
round any floating point numbers to a given number of (decimal) digits.
This allows us to ignore any variations below a tolerance level.

The script takes a single argument (the tolerance), reads from stdin,
and writes to stdout.
"""

import re
import sys

if len(sys.argv) != 2:
    print >> sys.stderr, "Usage:", sys.argv[0], " <tolerance>"
    sys.exit(1)
tolerance = int(sys.argv[1])
print 'Warning: These have been filtered by ',sys.argv[0],'.'
print 'This means that all floating point numbers have been rounded to ', tolerance, ' decimal places.'
def Replace(matchobj):
    """Given a match, round the number to the tolerance."""
    return str(round(float(matchobj.group(0)), tolerance))

number = re.compile(r'((\+|-)?[0-9]+\.[-+0-9e]+)')

for line in sys.stdin:
    print re.sub(number, Replace, line),
