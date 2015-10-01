import os, sys, inspect
import glob

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

    def TestName(self, logFile):
        """
        Return the name of the test that generated the test file
        """
        import re
        match = re.match('test_summary/(.*)_valgrind.out',logFile)
        if match:
            return match.group(1)
        else:
            return 'unknown test name!'


if __name__ == "__main__":
    files = glob.glob(sys.argv[1]+'/*_valgrind.out')
    procVal = ProcessValgrind()
    print '<!DOCTYPE html>'
    print '<html>'
    print '<body>'
    for file in files:
        status = procVal.EncodeStatus(open(file,'r'))
        print '<p> <font color="',procVal.StatusColour(status),'">',procVal.TestName(file),': ',procVal.DisplayStatus(status),'<a href="',file,'">(test output)</a>'
    print '</body>'
    print '</html>'


