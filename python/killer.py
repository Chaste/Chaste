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


# Kill off any processes run by our user in this directory (except us, of course).
# Will use psutil if available, direct reads of /proc if not.
# Note that the latter method only works on Linux.

import os
import signal
import sys
import time

try:
    import psutil
except ImportError:
    psutil = None


sim = '-s' in sys.argv
if '-d' in sys.argv:
    i = sys.argv.index('-d')
    kill_dir = os.path.realpath(sys.argv[i+1])
else:
    kill_dir = os.path.realpath(os.getcwd())

our_uid = os.getuid()
our_pid = os.getpid()

print "Killing processes owned by", our_uid, "in", kill_dir

if psutil:
    Process = psutil.Process
    get_procs = psutil.process_iter
else:
    # Define minimal psutil functionality implemented with /proc reads
    class Process(object):
        """A simple class representing a running process."""
        def __init__(self, pid):
            self._pid = str(pid)
            self.pid = int(pid)

        def cmdline(self):
            try:
                f = open('/proc/' + self._pid + '/cmdline')
                cmdline = f.read().split('\x00')[:-1]
                f.close()
            except OSError:
                cmdline = None
            return cmdline

        def uids(self):
            try:
                s = os.stat('/proc/' + self._pid)
                return (s.st_uid,)
            except OSError:
                return (None,)

        def send_signal(self, sig):
            try:
                os.kill(int(self._pid), sig)
            except OSError:
                pass

        def terminate(self):
            self.send_signal(signal.SIGTERM)

        def kill(self):
            self.send_signal(signal.SIGKILL)

        def getcwd(self):
            return os.path.realpath('/proc/' + self._pid + '/cwd')

    def get_procs():
        """Return an iterator over running processes, yielding a Process instance for each."""
        for pid in filter(lambda pid: pid[0] in "0123456789", os.listdir('/proc')):
            yield Process(pid)


def check_pid(proc):
    """Get key information about the given process.

    Returns a list with the command line elements if it's
    running as our_uid in the kill_dir, or None otherwise.
    """
    result = None
    try:
        if proc.uids()[0] == our_uid:
            if proc.getcwd() == kill_dir and proc.pid != our_pid:
                result = proc.cmdline()
    except:
        # We can't read all our processes; that's ok
        pass
    return result

def try_kill(proc, sig):
    """Try to kill a process, but ignore errors (e.g. because a process is already dead)."""
    try:
        proc.send_signal(sig)
    except:
        pass

# First, try killing off scons
for pid in get_procs():
    cmdline = check_pid(pid)
    if cmdline is not None:
        if len(cmdline) > 1 and 'scons' in cmdline[1]:
            print "SCons is running as PID", pid
            if not sim:
                try_kill(pid, signal.SIGTERM)
                print "  ** Killing (sent SIGTERM)"
                # Now sleep for a bit to let it die
                time.sleep(10) # seconds
                break

# Next, try to make the builder script save any results folders to the right place
for pid in get_procs():
    cmdline = check_pid(pid)
    if cmdline is not None:
        if len(cmdline) > 2 and 'builder' in cmdline[1] and 'no-lock' in cmdline[2]:
            print "Builder is running as PID", pid
            if not sim:
                try_kill(pid, signal.SIGUSR1)
                print "  ** Poking (sent SIGUSR1)"
                # Now sleep for a bit to let it move files
                time.sleep(5) # seconds
                break

# Next, try killing everything still running
for pid in get_procs():
    cmdline = check_pid(pid)
    if cmdline is not None:
        print pid, "is running from our dir as", cmdline[0:3]
        if not sim:
            try_kill(pid, signal.SIGKILL)
            print "  ** Killing (sent SIGKILL)"
