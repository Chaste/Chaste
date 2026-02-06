""""

Copyright (c) 2005-2026, University of Oxford.
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
import unittest
import subprocess

CHASTE_BUILD_DIR = "."

# This Python script accompanies the TestCommandLineArgumentsTutorial.
# This script will execute the TestCommandLineArgumentsTutorial
# multiple times with several command line arguments as test variables.


class TestPythonCommandLineArgumentsTutorial(unittest.TestCase):
    def setUp(self) -> None:
        self.test_exe = f"{CHASTE_BUILD_DIR}/global/test/TestCommandLineArgumentsTutorial"

    def test_command_line_default_tutorial(self) -> None:
        # First we will pass options to TestCommandLineDefaultTutorial without passing
        # a vector to TestCommandLineDoubleTutorial.
        # Here we will declare some values we wish to later pass to a for loop.
        N = 2
        L = 3
        M = 4

        # Here we set up a simple for loop over variables i, j and k based on the values of N, L and M.
        for i in range(0, N+1):
            for j in range(1, L+1):
                for k in range(2, M+1):
                    # Each loop runs an instance of the TestCommandLineArgumentsTutorial with
                    # opt1,opt2 and opt3 taking on the values of i, j and k respectively.
                    subprocess.run([self.test_exe, "-opt1", f"{i}", "-opt2", f"{j}", "-opt3", f"{k}"], check=True)

    def test_command_line_double_tutorial(self) -> None:
        # Now we will pass a vector to TestCommandLineDoubleTutorial without passing any
        # options to TestCommandLineDefaultTutorial.
        # Here we will declare some values we wish to later pass to a for loop.
        N = 2000
        L = 3000
        M = 4000

        # Here we set up a simple for loop over variables i, j and k based on the values of N, L and M.
        for i in range(1, N+1, 1000):
            idouble = i/120.0
            for j in range(1001, L+1, 1000):
                jdouble = j/150.0
                for k in range(2001, M+1, 1000):
                    kdouble = k/180.0
                    # Each loop runs an instance of the TestCommandLineArgumentsTutorial
                    # with a vector containing our variables.
                    subprocess.run([self.test_exe, "--my-vector-of-arguments",
                                   f"{idouble}", f"{jdouble}", f"{kdouble}"], check=True)


if __name__ == "__main__":
    unittest.main()
