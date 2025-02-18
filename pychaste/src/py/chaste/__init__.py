"""PyChaste Module"""

__copyright__ = """Copyright (c) 2005-2025, University of Oxford.
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

import os
import sys

import petsc4py

import chaste.cell_based
import chaste.core
import chaste.mesh
import chaste.ode
import chaste.pde
import chaste.visualization


def init(test_output=None, comm=None):
    """
     Initialize PETSc and set the CHASTE_TEST_OUTPUT environment variable.

    :param test_output: The CHASTE_TEST_OUTPUT directory.
    :param comm: MPI communicator.
    """
    # Set CHASTE_TEST_OUTPUT
    if test_output:
        # Set to user specified value if provided
        os.environ["CHASTE_TEST_OUTPUT"] = test_output
    elif os.environ.get("CHASTE_TEST_OUTPUT") is None:
        # Set to the current working directory if not already set
        os.environ["CHASTE_TEST_OUTPUT"] = os.getcwd()

    # Initialize PETSc
    if comm is None:
        petsc4py.init(sys.argv)
    else:
        petsc4py.init(comm=comm)

    return chaste.core.OutputFileHandler("", False)


init()
