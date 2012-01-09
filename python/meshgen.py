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
import os,sys

if len(sys.argv) < 2:
  sys.exit(1)

pref = sys.argv[1]
eles = 20
nodes = eles+1

f=file(pref+'.node', 'w')
f.write("%d\t1\t0\t1\n" % (nodes))
for i in range(0,nodes):
  b = 0
  if i == 0 or i == nodes: b=1
  f.write("%d %f %d\n" % (i, 1.0*i/eles, b))
f.close()

f=file(pref+'.ele', 'w')
f.write("%d\t2\t0\n" % (eles))
for i in range(0, eles):
  f.write("%d %d %d\n" % (i, i, i+1))
f.close()
