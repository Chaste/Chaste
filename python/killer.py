#!/usr/bin/env python


"""Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.
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
            os.kill(int(pid), signal.SIGTERM)
            print "  ** Killing (sent SIGTERM)"
