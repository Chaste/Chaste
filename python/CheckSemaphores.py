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


# This test is specific to a semaphore-based implementation of MPI.
# Its purpose it to provide a early warning when system semaphore are being
# used up, due to parallel tests aborting.

# CLEANING UP
# Look at  testoutput/$HOSTNAME.BUILD_TYPE/Semaphores...
# 
# There are 16 semaphores open, owned by 1 users:
#	chaste (16)
# etc.
# On the relevant machine and with the relevant user:
# ~/mpi/sbin/cleanipcs
# e.g. If you want to clean "spud" on the build server
# sudo su spud ~/mpi/sbin/cleanipcs

import os
from pwd import getpwuid


max_semaphore_arrays=int(os.popen('ipcs -l | awk \'/max number of arrays/\'').read().split()[-1])
#On the build server (and other Ubuntu machines) the value is max_semaphore_arrays==128
#MPI starts to fail when the number of semaphores reaches 128
semaphore_limit = max_semaphore_arrays/2

semaphore_data=os.popen('tail -n +2 /proc/sysvipc/sem  |awk \'{print $5}\'| sort | uniq -c').readlines()

total_open = 0
names = []
for entry in semaphore_data:
  total_open += int(entry.split()[0])
  names.append( getpwuid(int(entry.split()[1]))[0] + ' (' + entry.split()[0] +')' )
  
# Let the test summary script know
print "There are", total_open,"semaphores open, owned by",len(semaphore_data), "users:"
for name in names:
    print "\t", name
print "The next line is for the benefit of the test summary scripts."
if total_open > semaphore_limit:
    print "Failed",total_open,"of",total_open,"tests"
else:
    print "Infrastructure test passed ok."
