
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
