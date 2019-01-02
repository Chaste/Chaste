
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

import glob
import os
import re
import sys

class ProcessValgrind:

    def StatusColour(self, status):
        """
        Return a colour string indicating whether the given status string
        represents a 'successful' test suite under this build type.
        """
        if status == 'OK':
            return 'green'
        elif status == 'Warn':
            return 'orange'
        else:
            return 'red'

    def DisplayStatus(self, status):
        "Return a (more) human readable version of the given status string."
        if status == 'OK':
            return 'No leaks found'
        elif status == 'Unknown':
            return 'Test output unrecognised (RED)'
        elif status == 'Warn':
            return 'Possible leak found'
        elif status == 'Killed':
            return 'Test exceeded time limit (RED)'
        else:
            return 'Memory leaks found (RED)'


    def EncodeStatus(self, logFile, outputLines=None):
        """
        Encode the output from a test program as a status string.
        The output from valgrind needs to be parsed to check for a leak summary.
        If one is found the status is 'Leaky', otherwise 'OK'.
        Return the encoded status.
        """
        status = 'Unknown'

        # Regexps to check for
        import re
        invalid = re.compile(r'==\d+== Invalid ')
        glibc = re.compile(r'__libc_freeres')
        leaks = re.compile(r'==\d+== LEAK SUMMARY:')
        lost = re.compile(r'==\d+==\s+(definitely|indirectly|possibly) lost: ([0-9,]+) bytes in ([0-9,]+) blocks')
        petsc = re.compile(r'\[0]Total space allocated (\d+) bytes')
        uninit = re.compile(r'==\d+== (Conditional jump or move depends on uninitialised value\(s\)|Use of uninitialised value)')
        open_files = re.compile(r'==(\d+)== Open (?:file descriptor|AF_UNIX socket) (?![012])(\d+): (?!(?:/home/bob/eclipse/lockfile|/dev/urandom))(.*)')
        orte_init = re.compile(r'==(\d+)==    (?:by|at) .*(: orte_init)?.*')
        test_killed = 'Test killed due to exceeding time limit'

        if outputLines is None:
            outputLines = logFile.readlines()
        for lineno in range(len(outputLines)):
            if outputLines[lineno].startswith(test_killed):
                status = 'Killed'
                break
            m = petsc.match(outputLines[lineno])
            if m and int(m.group(1)) > 0:
                # PETSc Vec or Mat allocated and not destroyed
                status = 'Leaky'
                break

            m = uninit.match(outputLines[lineno])
            if m:
                # Uninitialised values problem
                status = 'Uninit'
                break

            m = invalid.match(outputLines[lineno])
            if m:
                # Invalid read/write/free()/etc. found. This is bad, unless it's glibc's fault.
                match = glibc.search(outputLines[lineno+3])
                if not match:
                    status = 'Leaky'
                    break

            m = leaks.match(outputLines[lineno])
            if m:
                # Check we have really lost some memory
                # (i.e. ignore 'still reachable' memory)
                lineno += 1
                match = lost.match(outputLines[lineno])
                while match:
                    blocks = int(match.group(3).replace(',', ''))
                    if blocks > 0:
                        # Indirectly lost memory should only be a warning, unless we also have
                        # directly lost memory, since indirect losses could be due to library
                        # errors that we're suppressing.
                        if match.group(1) == 'indirectly' and status == 'Unknown':
                            status = 'Warn'
                        else:
                            status = 'Leaky'
                            break
                    lineno += 1
                    match = lost.match(outputLines[lineno])
                break

            m = open_files.match(outputLines[lineno])
            if m:
                # There's a file open that shouldn't be.
                # Descriptors 0, 1 and 2 are ok, as are names /dev/urandom
                # and /home/bob/eclipse/lockfile, and the log files.
                # All these OK files are inherited from the parent process.
                if (not outputLines[lineno+1].strip().endswith("<inherited from parent>")
                    and not self._CheckOpenmpiFile(outputLines, lineno+1, orte_init)):
                    status = 'Openfile'
                    break
        if status == 'Unknown':
            status = 'OK'
        return status

    def _CheckOpenmpiFile(self, outputLines, lineno, regexp):
        """Check whether a purported open file is actually something from OpenMPI."""
        result = False
        m = regexp.match(outputLines[lineno])
        while m:
            if m and m.group(1):
                result = True
                break
            if not m:
                break
            lineno += 1
            m = regexp.match(outputLines[lineno])
        return result

if __name__ == "__main__":
    files = glob.glob(sys.argv[1]+'/*_valgrind.out')
    index_file = open(sys.argv[1]+'/index.html','w')
    procVal = ProcessValgrind()
    ok = True
    index_file.write('<!DOCTYPE html>\n')
    index_file.write('<html>\n')
    index_file.write('<body>\n')
    for file in files:
        filename = os.path.basename(file)
        testname = re.match('(.*)_valgrind.out',filename).group(1)
        status = procVal.EncodeStatus(open(file,'r'))
        colour = procVal.StatusColour(status)
        index_file.write('<p> <font color="%s">%s: %s <a href="%s">(test output)</a>\n'%(colour,testname,procVal.DisplayStatus(status),filename))
        if colour == 'red':
            ok = False
    index_file.write('</body>\n')
    index_file.write('</html>\n')
    index_file.close()
    if not ok:
        print('Memory testing not 100% pass rate - failing memory testing.')
        sys.exit(1)
    else:
        print('Memory testing 100% - test passed.')
