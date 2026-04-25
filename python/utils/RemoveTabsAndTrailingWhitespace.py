#!/usr/bin/env python3

"""Copyright (c) 2005-2026, University of Oxford.
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

"""Script to
i) replace tabs with 4 spaces - calls "sed -i 's/\t/    /g'" on all source files.
ii) remove any trailing whitespace - calls " sed -i 's/[ \\t]*$//' "
"""


import subprocess
from pathlib import Path

# List of files to include, and the number of spaces to replace tabs with (None means ignore tabs)
include_filenames: dict[str, str] = {
    'CMakeLists.txt': ' ' * 4,
    'Makefile': None,
    '.gitignore': None,
    'Dockerfile': None,
    'Doxyfile': None,
    'ChasteConfig.cmake.in': ' ' * 4,
    'ChasteBuildInfo_cmake.cpp.in': ' ' * 4,
    'Version.cpp.in': ' ' * 4,
    'CITATION.cff': ' ' * 2,
    'README': ' ' * 4,
    '.clang-tidy': None,
    '.clang-format': None,
}

# List of file extensions to include, and the number of spaces to replace tabs with (None means ignore tabs)
include_extensions: dict[str, str] = {
    '.hpp': ' ' * 4,
    '.cpp': ' ' * 4,
    '.py': ' ' * 4,
    '.ipynb': ' ' * 4,
    '.md': ' ' * 4,
    '.cmake': ' ' * 4,
    '.xml': ' ' * 4,
    '.json': ' ' * 4,
    '.yml': ' ' * 2,
    '.yaml': ' ' * 2,
    '.rst': ' ' * 4,
    '.rst': ' ' * 4,
    '.sh': ' ' * 2,
    '.xsd': ' ' * 2,
    '.m': ' ' * 3,
    '.java': ' ' * 4,
    '.txt': None,
    '.cellml': 3,
    '.cellml': 3,
    '.toml': None,
    '.cfg': None,
}



dir_ignores = ['cxxtest', 'docs', 'data', 'third_party_libs', '3rdparty', 'external', 'texttest']

tab_spaces = ' ' * 4
tracked_files = subprocess.check_output(['git', 'ls-files'], text=True).splitlines()

exts = set()

for file_name in tracked_files:
    file_path = Path(file_name)

    if any(dir_ignore in file_path.parts for dir_ignore in dir_ignores):
        continue

    process = False

    if file_path.name in include_filenames:
        process = True
        process_tabs = include_filenames[file_path.name]
    elif file_path.suffix in include_extensions:
        process = True
        process_tabs = include_extensions[file_path.suffix]

    if process:
        print("Checking %s" % file_path)
        if process_tabs is not None:
            subprocess.run(['sed', '-i', 's/\\t/%s/g' % tab_spaces, file_name], check=True)
        subprocess.run(['sed', '-i', 's/[ \\t]*$//', file_name], check=True)
