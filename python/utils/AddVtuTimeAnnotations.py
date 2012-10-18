#!/usr/bin/env python
"""Copyright (C) Victor Chang Cardiac Research Institute, 2012. (Original author A.Sadrieh)."""

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
 
"""Convert a .vtu file from a heart simulation to have time annotations.
   This changes lines like
   <DataArray type="Float64" Name="V_000019" ...
   <DataArray type="Float64" Name="Phi_e_000020" ...
   to
   <DataArray type="Float64" Name="V" TimeStep="19" ...
   <DataArray type="Float64" Name="Phi_e" TimeStep="20" ...
   
   It also needs to add a TimeValues annotation to the <UnstructuredGrid> component:
   <UnstructuredGrid TimeValues="0 1 ..." > 
"""
import fileinput
import sys
import re

"""
    Check that the time steps start at zero and are monotonically increasing by 1 each time
    Note that there are likely to be multiple copies of each time step.  
    We could count the number of copies and check that they're all the same...
    
    Returns the value of the final time step
"""
def QueryVtuFile(inputFileName, nameTimePair):
    last_time = 0
    for line in fileinput.input(inputFileName):
       match = nameTimePair.search(line)
       if (match):
         time_step = int(match.group(2))
         if ( time_step != last_time ):
             assert( time_step == last_time + 1)
             last_time = time_step
    return time_step

"""
Convert from one .vtu file to another adding time step annotations
nameTimePair is a regular expression
lastTimeStep is calculated during an earlier pass (QueryVtuFile)
"""
def AnnotateVtuFile(inputFileName, nameTimePair, outputFileName, lastTimeStep):
    timevalues = range(0, lastTimeStep+1, 1)
    header_mode = True
    out_fp = file(outputFileName, 'w')
    for line in fileinput.input(inputFileName):
        if (header_mode):
           if (line.find('<UnstructuredGrid>')>0):
              line = '<UnstructuredGrid TimeValues="' 
              for time in timevalues:
                  line += str(time)+ ' '
              line += '">\n'
           if (line.find('</UnstructuredGrid>')>0):
              header_mode = False
           match = nameTimePair.search(line)
           if (match):
             var_name = match.group(1)
             time_step_name = match.group(2)
             #Strip the prefix of zeros from the time_step
             time_step = str(int(time_step_name))
             line = line.replace(  '_'+time_step_name+'\"', '\" TimeStep=\"'+time_step+'\" ')
        # Line (possibly amended) gets written to the new file
        out_fp.write(line)

if __name__ == "__main__":
    #Checking command line arguments
    if len(sys.argv) != 3:
        print >> sys.stderr, "Usage:", sys.argv[0], "<input_vtu_file> <output_vtu_file>"
        sys.exit(1)
    #Reading in command line arguments
    input_name = sys.argv[1]
    output_name = sys.argv[2]
    #If input name doesn't end in .vtu then add it
    if (input_name[-4:] != '.vtu'):
        input_name = input_name+'.vtu'
    #If output_name doesn't end in .vtu then add it
    if (output_name[-4:] != '.vtu'):
        output_name=output_name+'.vtu'
        
    print'Reading from', input_name, 'and writing to', output_name
    
    #Regular expression for  ... Name="some_var_name_000020" ... as Chaste writes it.
    name_time_pair = re.compile('Name=\"(.*)_([0-9]{6})\"')
    
    # Check the time steps in the vtu file and query for the final time step 
    last_time_step = QueryVtuFile(input_name, name_time_pair)

    #Write a copy of the file but with annotations
    AnnotateVtuFile(input_name, name_time_pair, output_name, last_time_step)

