#!/usr/bin/env python


"""Copyright (c) 2005-2014, University of Oxford.
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


# Kill off any processes run by our user in this directory
# (except us, of course)

import os
import signal
import sys
import time

sim = '-s' in sys.argv
if '-d' in sys.argv:
    i = sys.argv.index('-d')
    kill_dir = os.path.realpath(sys.argv[i+1])
else:
    kill_dir = os.path.realpath(os.getcwd())

our_uid = os.getuid()
our_pid = os.getpid()

print "Killing processes owned by", our_uid, "in", kill_dir

def check_pid(pid):
    """Get information about the given process (pid is a string).

    Returns a list with the command line elements if it's
    running as our_uid in the kill_dir, or None otherwise.
    """
    result = None
    try:
        s = os.stat('/proc/' + pid)
        if s.st_uid == our_uid:
            cwd = os.path.realpath('/proc/' + pid + '/cwd')
            if cwd == kill_dir and int(pid) != our_pid:
                f = open('/proc/' + pid + '/cmdline')
                cmdline = f.read().split('\x00')[:-1]
                f.close()
                result = cmdline
    except OSError:
        # We can't read all our processes; that's ok
        pass
    return result

test_pids = filter(lambda pid: pid[0] in "0123456789", os.listdir('/proc'))

# First, try killing off scons
for pid in test_pids:
    cmdline = check_pid(pid)
    if cmdline is not None:
        if len(cmdline) > 1 and 'scons' in cmdline[1]:
            print "SCons is running as PID", pid
            if not sim:
                os.kill(int(pid), signal.SIGTERM)
                print "  ** Killing (sent SIGTERM)"
                # Now sleep for a bit to let it die
                time.sleep(10) # seconds
                # Then re-check running proceses
                test_pids = filter(lambda pid: pid[0] in "0123456789",
                                   os.listdir('/proc'))
                break

# Next, try killing everything still running
for pid in test_pids:
    cmdline = check_pid(pid)
    if cmdline is not None:
        print pid, "is running from our dir as", cmdline[0:3]
        if not sim:
            os.kill(int(pid), signal.SIGKILL)
            print "  ** Killing (sent SIGKILL)"
