
"""Copyright (c) 2005-2023, University of Oxford.
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
This little script makes it slightly easier to maintain Valgrind memory suppression
files.  By removing duplicates and providing a cannonical ordering for suppression files
side-by-side comparison of files is tractible.

Typical usage:
python python/utils/TidySuppressionFile.py chaste.supp chaste-lucid.supp

For each filename in the arguments:
* a backup copy is made (with appended tilde)
* all suppressions are parsed and stored
* suppressions that obvious match known library calls (PetscInitialize, tetgen etc.) are tagged
* tags are added into the user-defined name of the suppression where appropriate
* duplicate suppressions (where the body is character-for-character the same) are removed
* suppressions are reordered: first by tag and then lexicographically by the body of the suppresion
* the original file is overwritten

"""

import re
import shutil
import sys





def TidyFile(filepath):
    """Process, tidy and reorder a suppresion file."""
    print("Tidying file ", filepath)

    # Make a backup of the original file
    shutil.copy2(filepath, filepath+'~')

    # Now open the original and parse
    out = open(filepath, 'w')
    count = 0
    comment = ""
    suppressions = []
    in_suppression = False
    for line in open(filepath+'~', 'r'):
        if (not in_suppression):
            # Look for a pre-suppression comment
            if (line[0] == "#"):
                comment += line
            # Clear pre-suppressions comment (if there was one)
            if (count==0 and line[0] == "\n"):
                # Write out anything at the top of the file
                out.write(comment)
                out.write(line) # Retain blank line
                comment = ""
            if (line.strip()=="{"):
                count += 1
                in_suppression = True
                name = ""
                match = ""
                library_call = ""
        else:
            if (line.strip()=="}"):
                in_suppression = False
                suppressions.append((library_call, body, comment, name, memcheck, match))
                comment = ""
            elif (name==""):
                name = line
            elif (line.strip()[:16] == "match-leak-kinds"):
                assert(name != "")
                assert(memcheck != "")
                assert(match == "")
                match = line
            elif (line.strip()[:8] == "Memcheck"):
                assert (name != "")
                memcheck = line
                body = ""
            else:
                # Add line to body
                for pattern in ["PetscInitialize", "PetscFinalize", "tetgen", "H5F", "vtk", "CVode", "boost", "xerces"]:
                    if (pattern in line):
                        library_call = pattern
                body += line
            #end if
    # Sort the list of tuples by the body (lexicographically)
    suppressions.sort()
    
    # Make note of all suppression bodies seen so far
    seen_set = set()
    
    # Write out - overwrites the original file
    for (library_call, body, comment, name, memcheck, match) in suppressions:
        if (not (match, body) in seen_set):
            seen_set.add((match, body))
            
            # Top of suppression
            out.write(comment)
            out.write("{\n")
            
            # Optional name - add in the type if it's not already there
            if (not library_call in name):
                out.write("   "+library_call+" ")
            out.write(name)
            
            # Main part of suppression
            out.write(memcheck)
            out.write(match)
            out.write(body)

            # Bottom of suppression
            out.write("}\n")
        else:
            print("Found a duplicate:\n"+body)

    out.close()



if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: %s <path_to_release_notes.html>" % sys.argv[0])
    for i in sys.argv[1:]:
        TidyFile(i)
