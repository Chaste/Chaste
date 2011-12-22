
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
