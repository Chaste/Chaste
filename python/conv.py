
"""Copyright (c) 2005-2012, University of Oxford.
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

in_name = sys.argv[1]
out_name = sys.argv[2]

max_row, max_col, nnz = 0,0,0

fp = open(in_name)

i = open('/tmp/i.m', 'w')
j = open('/tmp/j.m', 'w')
sv = open('/tmp/sv.m', 'w')
i.write('i = [\n')
j.write('j = [\n')
sv.write('sv = [\n')

for line in fp:
    items = line.split('(')
    row = int(items[0].split()[1][:-1]) + 1
    max_row = max(row, max_row)
    items = items[1:]
    nnz += len(items)
    for item in items:
        col, rest = item.split(',')
        col = int(col) + 1
        max_col = max(col, max_col)
        val = rest.split(')')[0].strip()
        i.write('%d\n' % row)
        j.write('%d\n' % col)
        sv.write('%s\n' % val)
        #out.write('A(%d,%d) = %s;\n' % (row,col,val))

i.write('];\n')
j.write('];\n')
sv.write('];\n')

i.close()
j.close()
sv.close()

os.system('cat /tmp/i.m /tmp/j.m /tmp/sv.m > %s' % out_name)
out = open(out_name, 'a+')
out.write("A = sparse(i, j, sv, %d,%d);\n" % (max_row,max_col))
out.close()
fp.close()

print max_row, max_col, nnz
