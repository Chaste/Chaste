
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



def GetTestsInTestPacks(testRootDir, packNames=[], returnFoundPacks=False):
    """Generate a set of all test files listed in test pack files under the given folder.
    
    If packNames is non-empty, only test packs with matching names will be considered.
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
                                testfiles.add(relpath(os.path.join(dirpath, rel_testfile.strip()),
                                                      testRootDir))
                        pack_file.close()
                    except IOError:
                        pass
    if returnFoundPacks:
        return testfiles, found_packs
    else:
        return testfiles
