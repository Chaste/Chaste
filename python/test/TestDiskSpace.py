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
Check that there's enough disk space for regular Chaste building.
"""

import os
import unittest

class TestDiskSpace(unittest.TestCase):
    """Test that there is plenty of disk available for Chaste."""
    def TestDiskSpace(self):
        """Test that there is plenty of disk available for Chaste."""
        try:
            import psutil
            du = psutil.disk_usage
        except:
            # psutil 0.3.0 or newer is needed for this test to work
            return
        gb = 1024*1024*1024
        
        source_free_gb = du(__file__).free / gb
        print "Free space on Chaste source partition: %dGB" % source_free_gb
        self.failIf(source_free_gb < 10,
                    "The disk containing the Chaste source tree has less than 10GB of space left.")

        test_output_dir = os.path.abspath(os.getenv('CHASTE_TEST_OUTPUT', os.curdir))
        if not os.path.exists(test_output_dir):
            try:
                os.makedirs(test_output_dir)
            except:
                self.fail("Unable to create test output folder '%s'" % test_output_dir)
                return
        test_free_gb = du(test_output_dir).free / gb
        print "Free space on Chaste test output partition: %dGB" % test_free_gb 
        self.failIf(test_free_gb < 10,
                    "The disk containing the Chaste test output has less than 10GB of space left.")
