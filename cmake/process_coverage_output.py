
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

import re
import sys

if __name__ == "__main__":
    ok = False

    genhtml_output = open(sys.argv[1]+'/index.html', "r")
    table_cell_counter=0
    match_counter=0;
    for line in genhtml_output:
        #print line
        if table_cell_counter > 0:
            table_cell_counter = table_cell_counter + 1
        if re.match('(.*)headerItem(.*)Lines(.*)',line):
            # This line gives the header of the table in the top right of the genhtml index page.
            # It should only match one line in the entire file.            
            table_cell_counter = 1
            match_counter = match_counter + 1
        if table_cell_counter == 4:
            # If we go down 4 lines we get to the 
            if re.match('(.*)headerCovTableEntryHi(.*)',line):
                ok = True

    if match_counter is not 1:
        print('Did not find a match for coverage summary line, or found more than one match. process_coverage_output.py needs tweaking!')
        sys.exit(2)
        
    if not ok:
        print('Coverage is not 100% - failing coverage test.')
        sys.exit(1)
    else:
        print('Coverage 100% - test passed.')
