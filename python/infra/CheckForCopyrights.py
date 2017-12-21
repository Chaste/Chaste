#!/usr/bin/env python

"""Copyright (c) 2005-2017, University of Oxford.
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
import re
import sys


deprecated_notice = re.compile(r"""Copyright \(c\) 2005-\d{4}, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.
((
This file is part of Chaste.
)?)
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 \* Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 \* Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 \* Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES \(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION\)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT \(INCLUDING NEGLIGENCE OR OTHERWISE\) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
""", re.MULTILINE)

deprecated_notice_GPL = re.compile(r"""Copyright \(C\) University of Oxford, 2005-\d{4}

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.
((
This file is part of Chaste.
)?)
Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
\(at your option\) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.
""", re.MULTILINE)

current_notice="""Copyright (c) 2005-2017, University of Oxford.
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

py_current_notice='"""'+current_notice+'"""\n'
cpp_current_notice='/*\n\n'+current_notice+'\n*/'
cpp_notice_to_add = cpp_current_notice + "\n\n"

# This is used when replacing a deprecated notice with the latest version,
# to account for the optional text.
replacement_notice = current_notice.replace("\nThis file is part of Chaste.\n", r"\1")

output_notice=current_notice.replace("\nThis file is part of Chaste.\n", "")
boost_random_distribution_notice = """
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 """
pycml_notice=" Processed by pycml - CellML Tools in Python"
xsd2_notice="// Copyright (C) 2005-2007 Code Synthesis Tools CC"
xsd3_notice="// Copyright (C) 2005-2008 Code Synthesis Tools CC"
triangle_notice="""/*  Copyright 1993, 1995, 1997, 1998, 2002, 2005                             */
/*  Jonathan Richard Shewchuk                                                */"""
tetgen_notice="""///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TetGen                                                                    //
//                                                                           //
// A Quality Tetrahedral Mesh Generator and 3D Delaunay Triangulator         //
//                                                                           //
// Version 1.4                                                               //
// April 16, 2007                                                            //
//                                                                           //
// Copyright (C) 2002--2007                                                  //
// Hang Si                                                                   //
// Research Group Numerical Mathematics and Scientific Computing             //
// Weierstrass Institute for Applied Analysis and Stochastics                //
// Mohrenstr. 39, 10117 Berlin, Germany                                      //
// si@wias-berlin.de                                                         //
//                                                                           //
// TetGen is freely available through the website: http://tetgen.berlios.de. //
//   It may be copied, modified, and redistributed for non-commercial use.   //
//   Please consult the file LICENSE for the detailed copyright notices.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
"""
tetgen_predicates_notice="""/*****************************************************************************/
/*                                                                           */
/*  Routines for Arbitrary Precision Floating-point Arithmetic               */
/*  and Fast Robust Geometric Predicates                                     */
/*  (predicates.c)                                                           */
/*                                                                           */
/*  May 18, 1996                                                             */
/*                                                                           */
/*  Placed in the public domain by                                           */
/*  Jonathan Richard Shewchuk                                                */
/*  School of Computer Science                                               */
/*  Carnegie Mellon University                                               */
/*  5000 Forbes Avenue                                                       */
/*  Pittsburgh, Pennsylvania  15213-3891                                     */
/*  jrs@cs.cmu.edu                                                           */
"""
py_lgpl_notice = """# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details."""



def CheckForCopyrightNotice(findStrOrRe, fileIn):
    """Test if the (possibly multi-line) string/regexp findStr is contained anywhere in fileIn."""
    fileIn.seek(0)
    file_text = fileIn.read()
    if isinstance(findStrOrRe, type('')):
        found = file_text.find(findStrOrRe) >= 0
    else:
        found = findStrOrRe.search(file_text) is not None
    return found

def UpdateFile(oldFilePath, newFilePath):
    """Replace the contents of oldFilePath with newFilePath.

    This removes the old file and renames the new to match, but also
    transfers permissions etc.
    """
    perm = os.stat(oldFilePath).st_mode
    os.rename(newFilePath, oldFilePath)
    os.chmod(oldFilePath, perm)

def ReplaceStringInFile(findRe, repStr, filePath):
    """Replaces all strings matching findRe by repStr in file filePath."""
    tempName = filePath+'~'
    input = open(filePath)
    output = open(tempName, 'w')
    s = input.read()
    output.write(findRe.sub(repStr, s))
    output.close()
    input.close()
    UpdateFile(filePath, tempName)
    print 'Notice: replaced deprecated copyright notice in', filePath

def HeadAppendStringInFile(appendString, filePath):
    """Adds appendStr to the top of file filePath"""
    tempName = filePath+'~'
    input = open(filePath)
    output = open(tempName, 'w')
    s = input.read()
    output.write(appendString)
    output.write(s)
    output.close()
    input.close()
    UpdateFile(filePath, tempName)
    print 'Notice: applied copyright notice in ', filePath


def InspectFile(fileName):
    file_in = open(fileName)
    if fileName[-21:] == 'CheckForCopyrights.py':
        #Can't really check this one, since it knows all the licences
        return True
    valid_notice = False
    if (CheckForCopyrightNotice(cpp_current_notice, file_in) or
        CheckForCopyrightNotice(py_current_notice, file_in) or
        CheckForCopyrightNotice(output_notice, file_in)):
        #print 'Found current notice in '+file_name
        valid_notice=True
    if (CheckForCopyrightNotice(pycml_notice, file_in) or
        CheckForCopyrightNotice(boost_random_distribution_notice, file_in) or
        CheckForCopyrightNotice(xsd2_notice, file_in) or
        CheckForCopyrightNotice(xsd3_notice, file_in) or
        CheckForCopyrightNotice(triangle_notice, file_in) or
        CheckForCopyrightNotice(tetgen_predicates_notice, file_in) or
        CheckForCopyrightNotice(tetgen_notice, file_in) or
        CheckForCopyrightNotice(py_lgpl_notice, file_in)):
        #print 'Found 3rd party notice in '+file_name
        if valid_notice:
            print "Multiple notices on", file_name
            return False
        else:
            return True

    if valid_notice:
        return True

    if CheckForCopyrightNotice(deprecated_notice, file_in):
        print 'Found deprecated copyright notice for', fileName
        if apply_update:
            ReplaceStringInFile(deprecated_notice, replacement_notice, fileName)
            return True
        else:
            print 'Fix this by doing:',sys.argv[0],'-update'
            return False

    if CheckForCopyrightNotice(deprecated_notice_GPL, file_in):
        print 'Found deprecated GPL copyright notice for', fileName
        if apply_update:
            ReplaceStringInFile(deprecated_notice_GPL, replacement_notice, fileName)
            return True
        else:
            print 'Fix this by doing:',sys.argv[0],'-update'
            return False

    print 'Found no copyright notice for', fileName
    if apply_new:
        if fileName[-3:] == '.py':
            print 'Not implemented for .py files'
            return False
        else:
            HeadAppendStringInFile(cpp_notice_to_add, fileName)
        return True
    else:
        print 'Fix this by doing:',sys.argv[0],'-new'
        return False


if __name__ == '__main__':
    # Check, apply or modify the copyright notices.
    # .cpp, .hpp., .py, .java are C++, Python and Java code.
    exts = ['.cpp', '.hpp', '.py', '.java']

    # SCons files
    # output.chaste files in acceptance tests (all Chaste executables should output the valid copyright notice)
    # Version.cpp.in is the provenance file
    named_files = ['SConscript', 'SConstruct', 'output.chaste', 'Version.cpp.in']

    dir_ignores = ['Debug', 'Release', 'build', 'cxxtest', 'testoutput', 'doc', 'projects', 'hierwikiplugin']
    startchar_ignores = ['_', '.']
    exclusions = ['python/pycml/enum.py', 'python/pycml/pyparsing.py', 'python/pycml/schematron.py']

    apply_update = '-update' in sys.argv
    apply_new = '-new' in sys.argv

    chaste_dir = '.'
    if '-dir' in sys.argv:
        i = sys.argv.index('-dir')
        chaste_dir = os.path.realpath(sys.argv[i+1])

    num_no_copyrights = 0
    num_copyrights = 0
    chaste_dir_len = len(os.path.join(chaste_dir, ''))
    for root, dirs, files in os.walk(chaste_dir):
        relative_root = root[chaste_dir_len:]
        # Check for ignored dirs
        for dirname in dirs[:]:
            if dirname in dir_ignores or dirname[0] in startchar_ignores:
                dirs.remove(dirname)
        # Check for source files
        for file in files:
            relative_path = os.path.join(relative_root, file)
            name, ext = os.path.splitext(file)
            if ((ext in exts or file in named_files) and
                relative_path not in exclusions):
                file_name = os.path.join(root, file)
                if InspectFile(file_name) == False:
                    num_no_copyrights += 1
                else:
                    num_copyrights += 1

    # Let the test summary script know
    if chaste_dir == ".":
        dir = os.getcwd()
    else:
        dir = chaste_dir

    print "Copyright test run over ",dir," (",num_no_copyrights+num_copyrights,") files"
    if num_no_copyrights > 0:
        print
        print "The next line is for the benefit of the test summary scripts."
        print "Failed",num_no_copyrights,"of",num_no_copyrights+num_copyrights,"tests"

        # Return a non-zero exit code if orphans were found
        sys.exit(num_no_copyrights)
    else:
        print "Infrastructure test passed ok."
