
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
This little script makes it easier to create the canned copy of the ReleaseNotes
wiki page that is shipped with Chaste releases as docs/ReleaseNotes.html.  It
allows you merely to copy the HTML source from Trac into that file, then you can
run this script to tidy it up for you.

Things tidied are:
 * Remove class="..." from a tags
 * Remove id="..." from heading tags, and add a blank line before each heading which had an id
 * Turn relative URLs into absolute ones by prepending https://chaste.cs.ox.ac.uk
 * Remove <span class="icon">...</span>
"""

import re
import shutil
import sys

# The regular expressions to search for, and their replacements
modifiers = [(re.compile(r'(<a[^>]+)class=".*?"'), r'\1'),                  # Remove class attribute
             (re.compile(r'(<h.) id=".*?"'), r'\n\n\1'),                 # Remove id attribute & add blank line
             (re.compile(r'href="/'), r'href="https://chaste.cs.ox.ac.uk/'), # Make URLs absolute
             (re.compile(r'<span class="icon">.*?</span>'), r'')                # Remove icon spans
             ]

def TidyFile(filepath):
    """Tidy a release notes style HTML file."""
    print "Tidying file", filepath
    # Make a backup of the original file
    shutil.copy2(filepath, filepath+'~')
    # Now open the backup and edit each line, writing to the original location
    out = open(filepath, 'w')
    for line in open(filepath+'~', 'r'):
        for regexp, repl in modifiers:
            line = regexp.sub(repl, line)
        out.write(line)
    out.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print >>sys.stderr, "Usage:", sys.argv[0], "<path_to_release_notes.html>"
    TidyFile(sys.argv[1])
