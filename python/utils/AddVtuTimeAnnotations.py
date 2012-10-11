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
   <DataArray type="Float64" Name="V_000020" ...
   to
   <DataArray type="Float64" Name="V" TimeStep="19" ...
   <DataArray type="Float64" Name="V" TimeStep="20" ...
   
   It also needs to add a TimeValues annotation to the <UnstructuredGrid> component:
   <UnstructuredGrid TimeValues="0 1 ..." > 
"""
import fileinput
import sys
import re
import shutil

def getTimeSteps(fname):
   fname += '_times.info'
   filehandle = open (fname, 'r')
   result =  {'time_step_cnt': 0, 'time_step_cnt': 0, 'first': 0, 'last': 0}
   for line in filehandle:
       temp = line.split()
       if temp[0] == 'Number':
          result['time_step_cnt']=(int(temp[3]))
#       elif temp[0] == 'timestep':
#          result['time_step_len']=(int(temp[1]))
#       elif temp[0] == 'First':
#          result['first']=(int(temp[2]))
#       elif temp[0] == 'Last':
#          result['last']=(int(temp[2]))

# \todo #2245 This is a hack because the original assumption was that the actual times were going into the
# .vtu rather than the timestep numbers.
   result['time_step_len']=1
   result['first']=0
   result['last']=result['time_step_cnt']-1
   #print result
   return result


"""\todo #2245 Make this code generic to look for 
                 Name="some_name_with_possible_underscores_and_numbers_000000"
             and make it have the TimeStep annotation
 \todo #2245 We need to investigate RangeMin and RangeMax to get Paraview to automatically set the range for all time
"""
def AnnotateVtuFile(fname,timestep_info):
    timevalues = range(timestep_info['first'],timestep_info['last']+1,timestep_info['time_step_len'])
    header_mode = True
    for line in fileinput.input(fname, inplace=1):
        if (header_mode):
           if (line.find('<UnstructuredGrid>')>0):
              line = '<UnstructuredGrid TimeValues="' 
              for time in timevalues:
                  line += str(time)+ ' '
              line += '">\n'
           if (line.find('</UnstructuredGrid>')>0):
              header_mode = False
           if (line.find('DataArray type="Float64" Name="V_')>0):
              search_result = re.search('V_[0-9]{6}',line)
              v_str = search_result.group(0)
#              match = re.search('_[0-9]{6}', line)
#              print match.group(0), '--', match.group(0)  
              line = line.replace(v_str,'V')
              v_ind = int(v_str[2:])
              line = line.replace('Name="V"','Name="V" TimeStep="'+str(v_ind)+'" '  )
        # possibly a debug line from A.Sadrieh
        sys.stdout.write(line)

if __name__ == "__main__":
    #Checking command line arguments
    if len(sys.argv) != 3:
        print >> sys.stderr, "Usage:", sys.argv[0], "<path_to_original_vtu> <output_file>"
        sys.exit(1)
    #Reading in command line arguments
    fnameprefix = sys.argv[1]
    output_name = sys.argv[2]
    #If input name ends in .vtu then strip back to the prefix
    if (fnameprefix[-4:] == '.vtu'):
        fnameprefix=fnameprefix[:-4]
    #If output_name doesn't end in .vtu then add it
    if (output_name[-4:] != '.vtu'):
        output_name=output_name+'.vtu'
        
    print'Reading from', fnameprefix, 'and writing to', output_name
    #Read _times.info file    
    tstepinfo  = getTimeSteps(fnameprefix)
    #Take copy of original file
    shutil.copy(fnameprefix+'.vtu',output_name)
    #Annotate the copy
    AnnotateVtuFile(output_name,tstepinfo)

