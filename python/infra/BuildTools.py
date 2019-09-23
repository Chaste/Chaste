
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
Useful functions for use by the build system, and utility scripts.

These have been extracted from SConsTools.py so they can be used by scripts
run externally to an SConscript.
"""

import os

# Compatability with Python 2.3
try:
    set = set
except NameError:
    import sets
    set = sets.Set

# Pre-2.6 compatibility (on posix)
try:
    relpath = os.path.relpath
except AttributeError:
    def relpath(path, start=os.path.curdir):
        """Return a relative version of a path"""
    
        if not path:
            raise ValueError("no path specified")
        
        start_list = os.path.abspath(start).split(os.path.sep)
        path_list = os.path.abspath(path).split(os.path.sep)
        
        # Work out how much of the filepath is shared by start and path.
        i = len(os.path.commonprefix([start_list, path_list]))
    
        rel_list = [os.path.pardir] * (len(start_list)-i) + path_list[i:]
        return os.path.join(*rel_list)



def GetTestsInTestPacks(testRootDir, packNames=[], returnFoundPacks=False, subfolder=''):
    """Generate a set of all test files listed in test pack files under the given folder.
    
    If packNames is non-empty, only test packs with matching names will be considered.
    If subfolder is given, only test packs under that folder will be considered.
    """
    pack_suffix = 'TestPack.txt'
    suffix_len = len(pack_suffix)
    testfiles = set()
    found_packs = set()
    for dirpath, dirnames, filenames in os.walk(testRootDir):
        for dirname in dirnames[:]:
            if dirname in ['.svn', 'data']:
                dirnames.remove(dirname)
        for filename in filenames:
            if filename.endswith(pack_suffix):
                pack_name = filename[:-suffix_len]
                found_packs.add(pack_name)
                if not packNames or pack_name in packNames:
                    # Process this test pack file
                    try:
                        pack_file = file(os.path.join(dirpath, filename), 'r')
                        for rel_testfile in pack_file:
                            # Ignore empty lines and duplicates
                            rel_testfile = rel_testfile.strip()
                            if rel_testfile:
                                rel_test_path = relpath(os.path.join(dirpath, rel_testfile.strip()),
                                                        testRootDir)
                                if rel_test_path.startswith(subfolder):
                                    testfiles.add(rel_test_path)
                        pack_file.close()
                    except IOError:
                        pass
    if returnFoundPacks:
        return testfiles, found_packs
    else:
        return testfiles
