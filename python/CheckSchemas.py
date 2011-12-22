#!/usr/bin/env python

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
import difflib

# Check XSD files for consistent schemas

#dir_ignores = ['build', 'cxxtest', 'testoutput', 'doc', 'projects', 'global', 'linalg', 'ode', 'pde', 'docs', 'anim', 'link', 'linklib', 'cell_based', 'notforrelease_cell_based', 'python']
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

