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


# This test is specific to a semaphore-based implementation of MPI.
# Its purpose it to provide a early warning when system semaphore are being
# used up, due to parallel tests aborting.

# CLEANING UP
# Look at  testoutput/$HOSTNAME.BUILD_TYPE/Semaphores...
# 
# There are 16 semaphores open, owned by 1 users:
# chaste (16)
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
