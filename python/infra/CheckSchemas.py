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

import os
import sys
import difflib

# Check XSD files for consistent schemas

#dir_ignores = ['build', 'cxxtest', 'testoutput', 'doc', 'projects', 'global', 'linalg', 'ode', 'pde', 'docs', 'anim', 'link', 'linklib', 'cell_based', 'notforrelease_cell_based', 'notforrelease_lung', 'python']
# Important directories are heart, apps and notforrelease
dir_ignores = ['build', 'cxxtest', 'testoutput', 'python']
startchar_ignores = ['_', '.']
exception_tests = ['BrokenSchema.xsd',              #Used in TestHeartConfig.hpp
                   'EmptyRoot.xsd',                 #Used in TestHeartConfig.hpp
                   'schema with spaces.xsd',        #Used in TestHeartConfig.hpp
                   'ChasteParametersRelease1.xsd',  #Used specifically in TestHeartConfig.hpp as it's before name-spaces
                   'ChasteParametersRelease1_1.xsd' #Used specifically in TestHeartConfig.hpp for testing without schema location
                  ]

chaste_dir = '.'
num_bad = 0
num_good = 0
reference_path='./heart/src/io'

def FilesAreDifferent(pathToFile1, pathToFile2):
    """Find differences between to files"""
    file_string_1 = open(pathToFile1).readlines()
    file_string_2 = open(pathToFile2).readlines()
    #Produce generator for unified diff
    udiff = difflib.unified_diff(file_string_1, file_string_2, pathToFile1, pathToFile2)
    try:
        n = 10
        #Print the first n lines of the unified diff
        for _ in range(0, n):
            print udiff.next()     
        return True     
    except StopIteration:
        #The generator has little content and so was probably empty
        return (False)

#Build a list of reference files (in './heart/src/io') and a list of files which need to be verified
reference_files=[]
files_to_check=[]
for root, dirs, files in os.walk(chaste_dir):
    # Check for ignored dirs
    for dirname in dirs[:]:
        if dirname in dir_ignores or dirname[0] in startchar_ignores:
            dirs.remove(dirname)
    # Check for source files
    for file in files:
        name, ext = os.path.splitext(file)
        #Find schemas
        if (ext=='.xsd'):
            if root == reference_path:
                reference_files.append( file )
            else:
                files_to_check.append( (root,file) )

for (path, file) in files_to_check:
    file_name = os.path.join(path, file)
    if file in reference_files:
        if (FilesAreDifferent(os.path.join(reference_path, file), file_name)):
          num_bad+=1
        else:
          num_good+=1
    else:
        if file in exception_tests:
          num_good+=1
        else:  
          print file_name,' is an orphan.  There is no reference schema with that name.'
          num_bad+=1

print "Schema test run over (",num_bad+num_good,") files"
if num_bad > 0:
    print
    print "The next line is for the benefit of the test summary scripts."
    print "Failed",num_bad,"of",num_bad+num_good,"tests"

    # Return a non-zero exit code if orphans or bad schema were found
    sys.exit(num_bad)
else:
    print "Infrastructure test passed ok."

