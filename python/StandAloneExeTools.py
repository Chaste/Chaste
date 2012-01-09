
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

"""
Support methods for generating stand-alone executable packages using
shared libraries.
"""

import glob
import os
import re
import shutil
import subprocess
import sys

def CopyStrip(inFileName, outFileName, stripText):
    """Strip stripText from inFileName and write the result to outFileName.
    
    This can be used to, for example, strip the XML Schema location path from
    XML configuration files.
    """
    out_file = open(outFileName, 'w')
    for line in open(inFileName):
        out_file.write(line.replace(stripText, ''))
    out_file.close()

def GetArchBits():
    """Get the architecture type - returns '32' or '64' depending on the word size."""
    arch = os.uname()[-1] #Architecture is the last in the 5-tuple produced by uname
    if arch == 'x86_64':
        arch_bits = '64'
    else:
        arch_bits = '32'
    return arch_bits

def CompileChaste(target, build='GccOpt'):
    """Compile the code and bail out if necessary.
    
    target gives the target to build.
    """
    print 'Compiling dynamically-linked executable'
    if os.system('scons build=' + build + ' static=0 chaste_libs=1 exe=1 compile_only=1 ' + target):
        print "General build failure.  You aren't ready to release!"
        sys.exit(1)

def CopySharedLibraries(executablePath, librariesPath):
    """Copy shared libraries needed by executablePath into librariesPath, which will be created."""
    print 'Making a subdirectory of library dependencies'
    found_chaste_libraries = False
    ldd = subprocess.Popen(['ldd', executablePath], stdout=subprocess.PIPE ).communicate()[0]
    os.mkdir(librariesPath)
    #Look for "LIB WHITESPACE => WHITESPACE PATH WHITESPACE"
    lib_location = re.compile(r'(\S*)\s*=>\s*(\S*)\s*')
    for lib_pair in lib_location.findall(ldd):
        if lib_pair[1][0] == '(' or lib_pair[1] == 'not':
            print 'No library found for ', lib_pair[0]
        elif lib_pair[0].startswith('libc.so') or lib_pair[0].startswith('libpthread.so'):
            print 'Ignoring library ', lib_pair[0], 'for compatibility'
        else:
            if lib_pair[0] == 'libglobal.so':
                found_chaste_libraries = True
            shutil.copy(lib_pair[1], librariesPath)
    
    if not found_chaste_libraries:
      print 'Could not find Chaste libraries (e.g. libglobal.so).  Please set LD_LIBRARY_PATH.'
      sys.exit(1)

def CopyXmlFiles(destPath):
    """Copy the XML parameter files and schemas to the destination folder."""
    # Copy the xml files (with the hard coded location erased)
    CopyStrip('ChasteParameters.xml', os.path.join(destPath, 'ChasteParameters.xml'), 'heart/src/io/')
    CopyStrip('heart/test/data/xml/ChasteParametersFullFormat.xml', os.path.join(destPath, 'ChasteParametersFullFormat.xml'), '../../../src/io/')
    CopyStrip('heart/test/data/xml/ChasteParametersResumeSimulationFullFormat.xml', os.path.join(destPath, 'ChasteParametersResumeSimulationFullFormat.xml'), '../../../src/io/')
    # Copy all xsd files
    for xsd in glob.glob('heart/src/io/ChasteParameters*.xsd'):
        shutil.copy(xsd, destPath)

def CopyTree(sourcePath, destPath):
    """Copy a tree, ignoring .svn folders and editor roll-backs."""
    shutil.copytree(sourcePath, destPath)
    for dirpath, dirnames, filenames in os.walk(destPath):
        if '.svn' in dirnames:
            dirnames.remove('.svn')
            shutil.rmtree(os.path.join(dirpath, '.svn'))
        for filename in filenames:
            if filename[-1] == '~':
                os.remove(os.path.join(dirpath, filename))

def CreateWrapperScript(exeName, folderPath):
    """Create a wrapper script to run the given executable."""
    script_path = os.path.join(folderPath, exeName + '.sh')
    script = open(script_path, 'w')
    script.write(r"""#!/bin/bash

# A wrapper script for Chaste that can figure out where it really lives on the
# filesystem, set the LD_LIBRARY_PATH to the correct subfolder and run Chaste.

#IFS=" \t\n"
#declare -x PATH=/bin:/usr/bin

script="$0"

# Have we been called via a symlink?
while [[ -L "$script" ]]; do script=$(readlink -n "$script"); done

# Figure out the folder $script is in
script_path=$(2>/dev/null cd "${script%%/*}" >&2; echo "`pwd -P`/${script##*/}")
script_dir=$(dirname "$script_path")

#Set the LD_LIBRARY_PATH and run Chaste
export LD_LIBRARY_PATH="$script_dir/libs"

#Inform the user where the output will appear
if [ -z "$CHASTE_TEST_OUTPUT" ]; then
  echo "\$CHASTE_TEST_OUTPUT is currently unset.  Your output will appear in ./testoutput"
else
  echo "\$CHASTE_TEST_OUTPUT is currently set to " $CHASTE_TEST_OUTPUT. 
fi

# This line actually run Chaste with the given arguments
"$script_dir/%s" "$@"
""" % exeName)
    script.close()
    os.chmod(script_path, 0755)
