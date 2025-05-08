"""Copyright (c) 2005-2025, University of Oxford.
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
import pathlib
import re
import subprocess
import sys


class ProcessValgrind:

    @staticmethod
    def status_colour(_status):
        """
        Return a colour string indicating whether the given status string
        represents a 'successful' test suite under this build type.
        """
        if _status == 'OK':
            return 'green'
        elif _status == 'Warn':
            return 'orange'
        else:
            return 'red'

    @staticmethod
    def display_status(_status):
        """Return a (more) human-readable version of the given status string."""
        if _status == 'OK':
            return 'No leaks found'
        elif _status == 'Unknown':
            return 'Test output unrecognised (RED)'
        elif _status == 'Warn':
            return 'Possible leak found'
        elif _status == 'Killed':
            return 'Test exceeded time limit (RED)'
        else:
            return 'Memory leaks found (RED)'

    def encode_status(self, log_file, output_lines=None):
        """
        Encode the output from a test program as a status string.
        The output from valgrind needs to be parsed to check for a leak summary.
        If one is found the status is 'Leaky', otherwise 'OK'.
        Return the encoded status.
        """
        _status = 'Unknown'

        # Regexps to check for
        import re
        invalid = re.compile(r'==\d+== Invalid ')
        glibc = re.compile(r'__libc_freeres')
        leaks = re.compile(r'==\d+== LEAK SUMMARY:')
        lost = re.compile(r'==\d+==\s+(definitely|indirectly|possibly) lost: ([0-9,]+) bytes in ([0-9,]+) blocks')
        petsc = re.compile(r'\[0]Total space allocated (\d+) bytes')
        uninit = re.compile(
            r'==\d+== (Conditional jump or move depends on uninitialised value\(s\)|Use of uninitialised value)')
        open_files = re.compile(
            r'==(\d+)== Open (?:file descriptor|AF_UNIX socket) (?![012])(\d+): (?!(?:/home/bob/eclipse/lockfile|/dev/urandom))(.*)')
        orte_init = re.compile(r'==(\d+)==    (?:by|at) .*(: orte_init)?.*')
        test_killed = 'Test killed due to exceeding time limit'

        if output_lines is None:
            output_lines = log_file.readlines()
        for lineno in range(len(output_lines)):
            if output_lines[lineno].startswith(test_killed):
                _status = 'Killed'
                break
            m = petsc.match(output_lines[lineno])
            if m and int(m.group(1)) > 0:
                # PETSc Vec or Mat allocated and not destroyed
                _status = 'Leaky'
                break

            m = uninit.match(output_lines[lineno])
            if m:
                # Uninitialised values problem
                _status = 'Uninit'
                break

            m = invalid.match(output_lines[lineno])
            if m:
                # Invalid read/write/free()/etc. found. This is bad, unless it's glibc's fault.
                match = glibc.search(output_lines[lineno + 3])
                if not match:
                    _status = 'Leaky'
                    break

            m = leaks.match(output_lines[lineno])
            if m:
                # Check we have really lost some memory
                # (i.e. ignore 'still reachable' memory)
                lineno += 1
                match = lost.match(output_lines[lineno])
                while match:
                    blocks = int(match.group(3).replace(',', ''))
                    if blocks > 0:
                        # Indirectly lost memory should only be a warning, unless we also have
                        # directly lost memory, since indirect losses could be due to library
                        # errors that we're suppressing.
                        if match.group(1) == 'indirectly' and _status == 'Unknown':
                            _status = 'Warn'
                        else:
                            _status = 'Leaky'
                            break
                    lineno += 1
                    match = lost.match(output_lines[lineno])
                break

            m = open_files.match(output_lines[lineno])
            if m:
                # There's a file open that shouldn't be.
                # Descriptors 0, 1 and 2 are ok, as are names /dev/urandom
                # and /home/bob/eclipse/lockfile, and the log files.
                # All these OK files are inherited from the parent process.
                if (not output_lines[lineno + 1].strip().endswith("<inherited from parent>")
                        and not self._check_openmpi_file(output_lines, lineno + 1, orte_init)):
                    _status = 'Openfile'
                    break
        if _status == 'Unknown':
            _status = 'OK'
        return _status

    @staticmethod
    def _check_openmpi_file(output_lines, lineno, regexp):
        """Check whether a purported open file is actually something from OpenMPI."""
        result = False
        m = regexp.match(output_lines[lineno])
        while m:
            if m and m.group(1):
                result = True
                break
            if not m:
                break
            lineno += 1
            m = regexp.match(output_lines[lineno])
        return result

    @staticmethod
    def get_html_head():
        return """
<head>
  <style>
    .test-green {
      color: green;
    }
    .test-orange {
      color: orange;
    }
    .test-red {
      color: red;
    }
  </style>
  <title>Chaste valgrind memtest output</title>
</head>
"""

    @staticmethod
    def _get_git_info():
        """
        Get branch and commit information either from GitHub Actions environment or local git repository.
        Returns:
            tuple: branch name, commit SHA
        """

        # Most likely, we're on GitHub actions, and we can get the info we need from environment variables that are
        # set to appropriate values in the memory-testing.yml workflow file
        branch = os.environ.get('Chaste_GH_WORKFLOW_BRANCH')
        commit = os.environ.get('Chaste_GH_WORKFLOW_SHA')
        if branch and commit:
            return branch, commit

        # Otherwise, we interrogate the Git repository directly
        chaste_source_dir = pathlib.Path(__file__).parent.parent

        try:
            branch = subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=chaste_source_dir).strip().decode("utf-8")
        except subprocess.CalledProcessError as _:
            branch = "unknown branch"

        try:
            commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=chaste_source_dir).strip().decode("utf-8")
        except subprocess.CalledProcessError as _:
            commit = "unknown commit"

        return branch, commit

    @staticmethod
    def get_index_heading():
        branch, commit = ProcessValgrind._get_git_info()

        if commit == "unknown commit":
            return f'Memtest output for {commit} on {branch}'
        else:
            return f'Memtest output for commit <a href="https://github.com/Chaste/Chaste/commit/{commit}">{commit}</a> on branch {branch}'


if __name__ == "__main__":

    files = glob.glob(sys.argv[1] + '/*_valgrind.txt')
    procVal = ProcessValgrind()

    green_lines = []
    orange_lines = []
    red_lines = []

    for file in files:
        filename = os.path.basename(file)
        testname = re.match('(.*)_valgrind.txt', filename).group(1)
        status = procVal.encode_status(open(file, 'r'))
        colour = procVal.status_colour(status)
        disp_status = procVal.display_status(status)

        if colour == 'green':
            green_lines.append(f'  <p class="test-green">{testname}: {disp_status} <a href="{filename}">(test output)</a></p>')
        elif colour == 'orange':
            orange_lines.append(f'  <p class="test-orange">{testname}: {disp_status} <a href="{filename}">(test output)</a></p>')
        else:  # colour == 'red'
            red_lines.append(f'  <p class="test-red">{testname}: {disp_status} <a href="{filename}">(test output)</a></p>')

    num_tests = len(green_lines) + len(orange_lines) + len(red_lines)

    with open(sys.argv[1] + '/index.html', 'w') as index_file:

        index_file.write('<!DOCTYPE html>\n')
        index_file.write('<html lang="en">\n')
        index_file.write(ProcessValgrind.get_html_head())
        index_file.write('<body>\n')

        # Write heading
        index_file.write(f'  <h2>{ProcessValgrind.get_index_heading()}</h2>\n')
        index_file.write(f'  <br>\n')

        # Write summary
        index_file.write(f'  <p><strong>Summary of {num_tests} tests:</strong></p>\n')
        index_file.write(f'  <p class="test-green"><strong>green: {len(green_lines)}</strong></p>\n')
        index_file.write(f'  <p class="test-orange"><strong>orange: {len(orange_lines)}</strong></p>\n')
        index_file.write(f'  <p class="test-red"><strong>red: {len(red_lines)}</strong></p>\n')
        index_file.write(f'  <br>\n')

        # Write red lines
        if len(red_lines) > 0:
            index_file.write('\n'.join(sorted(red_lines)))
            index_file.write(f'  <br>\n')

        # Write orange lines
        if len(orange_lines) > 0:
            index_file.write('\n'.join(sorted(orange_lines)))
            index_file.write(f'  <br>\n')

        # Write green lines
        if len(green_lines) > 0:
            index_file.write('\n'.join(sorted(green_lines)))
            index_file.write(f'  <br>\n')

        index_file.write('</body>\n')
        index_file.write('</html>\n')

    if len(red_lines) > 0:
        print('Memory testing not 100% pass rate - failing memory testing.')
        sys.exit(1)
    else:
        print('Memory testing 100% - test passed.')
