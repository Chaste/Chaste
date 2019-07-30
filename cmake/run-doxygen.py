#!/usr/bin/python
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
import shutil
import sys

error_log = 'doxygen-error.log'
output_log = 'doxygen-output.log'


def run_doxygen(revision, hide_undoc_classes=True):
    """Run Doxygen on the given checkout."""
    # Run Doxygen
    cmd = '( cat Doxyfile ; echo "PROJECT_NUMBER=Build:: ' + str(revision) + '"'
    if hide_undoc_classes:
        cmd += '; echo "HIDE_UNDOC_CLASSES = YES" '
    cmd += ' ) | doxygen - 2>'+error_log+' 1>'+output_log

    print(cmd)

    exitcode = os.system(cmd)
    assert exitcode == 0, "Doxygen returned non-zero exit code"
    assert os.path.exists('doxygen/html'), "What happened to our output???"


if __name__ == "__main__":
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        print >>sys.stderr, "Usage:", sys.argv[0], "<source_root>", "<docs_folder>", "<revision>", "[check_coverage]"
        sys.exit(1)
    SOURCE_DIR = sys.argv[1]
    DOCS_DIR = sys.argv[2]
    build_version = sys.argv[3]
    if len(sys.argv) == 5:
        check_coverage = bool(sys.argv[4])
    else:
        check_coverage = False

    CWD = os.getcwd()
    os.chdir(SOURCE_DIR)

    # Run Doxygen
    run_doxygen(build_version, hide_undoc_classes=check_coverage)

    # Copy the files to DOCS_DIR
    if os.path.exists(DOCS_DIR):
        shutil.rmtree(DOCS_DIR)
    shutil.copytree('doxygen/html', DOCS_DIR)

    if check_coverage:
        sys.path.insert(0, os.path.join(SOURCE_DIR, 'python', 'infra'))
        from ParseDoxygen import parse_doxygen
        parse_doxygen(os.path.join(SOURCE_DIR, output_log), os.path.join(SOURCE_DIR, error_log), DOCS_DIR)

    os.chdir(CWD)
