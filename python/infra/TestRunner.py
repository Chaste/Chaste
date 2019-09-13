#!/usr/bin/env python


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


# Script to run a test executable and record the status
# (whether all tests passed or how many failed) as part of
# the output filename.

# The test executable should be run and the output captured,
# and stored to a file.

# This script expects 3 - 5 arguments.
# The first is the executable to run.
# The second is the name of the output .log file, in which the output
#  of the program will be stored.
# The 3rd is the short name for the type of build performed.
# The 4th, if present, is any extra flags to pass to the executable.

# The directory in which to store copies of the log files
# with altered names encoding the status, as determined by the
# BuildTypes module, is determined from the build object.
# It will be created if necessary.
# The .log file basename, without extension, will be prepended to the status.

import copy
import os
import signal
import subprocess
import sys
import time
import threading

try:
    import psutil
    if psutil.version_info[0] > 1:
        # New API :(
        get_prop = lambda proc, prop: getattr(proc, prop)()
    else:
        get_prop = lambda proc, prop: getattr(proc, prop)
except ImportError:
    psutil = None

def usage():
    print "Usage:",sys.argv[0],\
        "<test exe> <.log file> <build type> [run time flags] [--no-stdout]"

def KillTest(pid=None, exe=None):
    """Kill off a running test, specified either by pid or by name, or both.

    First, recursively kill pid and all its children.  Then check all running
    processes to see if they are instances of exe, and recursively kill them
    if so.

    Requires 'easy_install psutil' to function.
    """
    # Recursive kill of pid and its children, killing innermost first
    if pid:
        try:
            for proc in psutil.process_iter():
                if get_prop(proc, 'ppid') == pid:
                    KillTest(proc.pid)
            print "Killing", pid
            os.kill(pid, signal.SIGKILL)
        except (psutil.NoSuchProcess, OSError):
            pass
    # Kill by name - figure out pid and call ourselves with that
    if exe:
        #print "Killing", exe
        for proc in psutil.process_iter():
            try:
                if get_prop(proc, 'exe') == exe:
                    KillTest(proc.pid)
            except (psutil.NoSuchProcess, psutil.AccessDenied, OSError):
                pass


def GetTestNameFromLogFilePath(logFilePath):
    """Figure out what name to display for a test suite based on the results log path.

    Except for special infrastructure tests, the log file will always live at
    <component>/build/<build_dir>/<relative_path>.log.  We return the test name
    <component>/test/<relative_path>, with all occurrences of os.sep replaced by '-'.
    """
    #We used to do: os.path.splitext(os.path.basename(logfile))[0]
    path_no_ext = os.path.splitext(logFilePath)[0]
    parts = path_no_ext.split(os.sep)
    try:
        build_idx = parts.index('build')
        test_name = os.sep.join(parts[:build_idx] + ['test'] + parts[build_idx+2:])
        test_name = test_name.replace(os.sep, '-')
    except ValueError:
        test_name = os.path.basename(path_no_ext)
    return test_name


class TestKillerTimer(object):
    """A specialised timer creator for killing running tests if they take too long."""
    def __init__(self, timeLimit, pid, exe, log):
        """Create a new timer."""
        self.timeLimit = timeLimit
        self.pid = pid
        self.exe = exe
        self.log = log
        self.extraTimer = None
        self.mainTimer = threading.Timer(timeLimit, self.MaybeDoKill)
        self.mainTimer.start()

    def Stop(self):
        """Cancel the timer (if it hasn't already fired) and wait for it to stop fully."""
        self.mainTimer.cancel()
        if self.extraTimer:
            self.extraTimer.cancel()
            # Wait for it to finish
            self.extraTimer.join()
        # Wait for the main timer to finish
        self.mainTimer.join()

    def GetCpuTime(self, process):
        """Recursively compute the total CPU time used by a process and its children."""
        cpu_time = sum(process.get_cpu_times())
        for child in process.get_children():
            cpu_time += self.GetCpuTime(child)
        return cpu_time

    def MaybeDoKill(self):
        """Only kill the test if it has consumed a reasonable amount of CPU, otherwise wait a bit longer."""
        try:
            test_proc = psutil.Process(self.pid)
            cpu_time = self.GetCpuTime(test_proc)
            if cpu_time < 0.5 * self.timeLimit:
                # The machine is probably loaded, so allow some extra grace time.
                self.extraTimer = threading.Timer(self.timeLimit, self.DoKill)
                self.extraTimer.start()
        except:
            self.DoKill()

    def DoKill(self):
        """Really kill off the test and its children."""
        KillTest(self.pid, self.exe)
        msg = '\n\nTest killed due to exceeding time limit of %d seconds\n\n' % self.timeLimit
        self.log.write(msg)
        print msg


def run_test(exefile, logfile, build, run_time_flags='', echo=True, time_limit=0, env=os.environ):
    """Actually run the given test."""
    # Find out how we're supposed to run tests under this build type
    if exefile.startswith("python/infra/Check"):
        command = exefile + ' 2>&1 ' + run_time_flags
    elif exefile.endswith('.py'):
        if build.build_dir == 'line_profile':
            prof = ' --line-profile'
        elif build.is_profile:
            prof = ' --profile'
        else:
            prof = ''
        command = './python/infra/TestPythonCode.py --num-procs ' + str(build.GetNumProcesses()) + ' ' + exefile + prof + ' 2>&1 ' + run_time_flags
    else:
        command = build.GetTestRunnerCommand(exefile, '2>&1 ' + run_time_flags)

    # Run the test program and record output & exit code
    log_fp = file(logfile, 'w')
    if not log_fp:
        raise IOError("Unable to open log file")
    if build.no_store_results:
        stdout = None
        echo = False
    else:
        stdout = subprocess.PIPE
    start_time = time.time()
    test_proc = subprocess.Popen(command, shell=True, stdout=stdout, stderr=subprocess.STDOUT, env=env)
    if time_limit and psutil:
        # Set a Timer to kill the test if it runs too long
        kill_timer = TestKillerTimer(time_limit, test_proc.pid, os.path.abspath(exefile), log_fp)
    for line in test_proc.stdout or []:
        log_fp.write(line)
        if echo:
            print line,
    exit_code = test_proc.wait()
    end_time = time.time()
    if time_limit and psutil:
        kill_timer.Stop()
    log_fp.close()

    #print "Time",end_time,start_time

    if True:
        import glob
        # Get the test status and copy log file
        test_dir = build.output_dir
        if not os.path.exists(test_dir):
            os.makedirs(test_dir)
        if not os.path.isdir(test_dir):
            print test_dir, "is not a directory; unable to copy output."
            sys.exit(1)
        test_name = GetTestNameFromLogFilePath(logfile)
        log_fp = open(logfile, 'r')
        status = build.EncodeStatus(exit_code, log_fp)
        log_fp.close()
        #print test_name, status
        # Remove any old copies of results from this test
        oldfiles = glob.glob(os.path.join(test_dir, test_name+'.*'))
        for oldfile in oldfiles:
            os.remove(oldfile)
        # Copy the new results
        copy_to = build.ResultsFileName(dir=test_dir, testsuite=test_name, status=status,
                                        runtime=end_time-start_time)
        #print copy_to
        os.system("/bin/cp -f " + logfile + " " + copy_to)
        return copy_to


if __name__ == '__main__':
    # Check for valid arguments.
    if '--no-stdout' in sys.argv:
        echo = False
        sys.argv.remove('--no-stdout')
    else:
        echo = True

    if len(sys.argv) < 4:
        print "Syntax error: insufficient arguments."
        usage()
        sys.exit(1)

    exefile, logfile, build_type = sys.argv[1:4]

    # Get any extra command line args for the test
    if len(sys.argv) > 4:
        run_time_flags = sys.argv[4]
    else:
        run_time_flags = ''

    # Get a build object for this build type
    parent_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    sys.path[0:0] = [parent_path]
    import BuildTypes
    build = BuildTypes.GetBuildType(build_type)

    run_test(exefile, logfile, build, run_time_flags, echo)

# Builder function for running via SCons
def get_build_function(build, run_time_flags='', test_time_limit=0):
    """Return a function that can be used as a Builder by SCons."""

    def build_function(target, source, env):
        # Set up the environment from env['ENV']
        # Python tests need PYTHONPATH set from env['ENV']['PYTHONPATH'] and env['PYINCPATH']
        orig_python_path = copy.copy(env['ENV']['PYTHONPATH'])
        new_paths = env.Flatten(env.subst(env.get('PYINCPATH', ''), conv=lambda n: n))
        env.PrependENVPath('PYTHONPATH', new_paths)
        test_env = copy.deepcopy(env['ENV'])
        env['ENV']['PYTHONPATH'] = orig_python_path
        # Functional curation also needs CFLAGS & LDFLAGS for Cython compilation
        test_env['CFLAGS'] = os.environ.get('CHASTE_CYTHON_CFLAGS', '')
        test_env['LDFLAGS'] = os.environ.get('CHASTE_CYTHON_LDFLAGS', '')
        # Run the test
        log = run_test(str(source[0]), str(target[0]), build, run_time_flags, time_limit=test_time_limit, env=test_env)
        # Note the extra dependency of the copied log file
        #env.SideEffect(log, target)
        env.Depends(os.path.join(os.path.dirname(log), 'index.html'), log)
        return None

    return build_function
