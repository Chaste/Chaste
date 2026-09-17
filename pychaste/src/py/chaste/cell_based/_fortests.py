"""Helper classes for running cell-based tests."""

__copyright__ = """Copyright (c) 2005-2026, University of Oxford.
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

from chaste.cell_based import CellId, SimulationTime
from chaste.core import RandomNumberGenerator, Timer


def SetupNotebookTest():
    """Reset the global simulation state at the start of a tutorial notebook.

    A notebook cannot inherit the setUp of the test-suite classes below, so it
    calls this explicitly to keep notebooks that run in sequence independent.
    """
    simulation_time = SimulationTime.Instance()
    simulation_time.SetStartTime(0.0)
    RandomNumberGenerator.Instance().Reseed(0)
    CellId.ResetMaxCellId()


def TearDownNotebookTest():
    """Destroy the global simulation state at the end of a tutorial notebook.

    The teardown counterpart of :func:`SetupNotebookTest`.
    """
    simulation_time = SimulationTime.Instance()
    simulation_time.Destroy()
    RandomNumberGenerator.Instance().Destroy()


class AbstractCellBasedTestSuite(unittest.TestCase):
    """Base test suite that resets the global simulation state around each test.

    The Python equivalent of the C++ base class of the same name (see
    cell_based/src/fortests/AbstractCellBasedTestSuite.hpp). It is not wrapped
    from C++; it is re-implemented directly in Python as a unittest.TestCase
    subclass so Python test suites get the same per-test setUp/tearDown behaviour.
    """

    def setUp(self):
        """Reset the global simulation state before each test."""
        SetupNotebookTest()

    def tearDown(self):
        """Destroy the global simulation state after each test."""
        TearDownNotebookTest()


class AbstractCellBasedWithTimingsTestSuite(AbstractCellBasedTestSuite):
    """AbstractCellBasedTestSuite that also times each test.

    The Python equivalent of the C++ base class of the same name (see
    cell_based/src/fortests/AbstractCellBasedWithTimingsTestSuite.hpp).
    """

    def setUp(self):
        """Reset the timer, then the global simulation state, before each test."""
        Timer().Reset()
        super().setUp()

    def tearDown(self):
        """Destroy the global simulation state, then report the elapsed time."""
        super().tearDown()
        Timer().Print("Test elapsed")
