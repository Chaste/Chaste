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
print 'This means that all floating point numbers have been rounded to', tolerance, 'decimal places.'
def Replace(matchobj):
    """Given a match, round the number to the tolerance."""
    number = float(matchobj.group(0))
    rounded = round(number, tolerance)
    discriminant = (rounded - number)*10**(tolerance+1) # How much we have moved the truncated digit
    discriminant = round(discriminant, 10)
    if (number<0 and discriminant == 5.0):
        # See: http://docs.python.org/2/library/functions.html#round
        # Negative number has been arbitrarily rounded towards zero - round it away
        rounded -= 10**(-tolerance)
    elif (number>0 and discriminant == -5.0):
        # Positive number has been arbitrarily rounded towards zero - round it away
        rounded += 10**(-tolerance)
    return str(rounded)

#Number must either have a decimal point (and no "e") or match scientific notation
decimal = '(\+|-)?\d+\.\d+'
#scientific notation usually has a decimal point
scientific = '(\+|-)?\d*\.\d*e(\+|-)?\d+'
#scientific notation might be of the form -8e-10 (a single digit with no decimal point)
scientific_no_point = '(\+|-)?\de(\+|-)?\d+'

#Match scientific first
number = re.compile('(' + scientific +'|' + decimal +'|' + scientific_no_point + ')')

for line in sys.stdin:
    print re.sub(number, Replace, line),
