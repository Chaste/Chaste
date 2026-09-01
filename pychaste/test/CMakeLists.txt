# Copyright (c) 2005-2026, University of Oxford.
# All rights reserved.
#
# University of Oxford means the Chancellor, Masters and Scholars of the
# University of Oxford, having an administrative office at Wellington
# Square, Oxford OX1 2JD, UK.
#
# This file is part of Chaste.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * Neither the name of the University of Oxford nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if(NOT Chaste_ENABLE_PYCHASTE)
    return()
endif()

project(testpychaste)
chaste_do_test_component(pychaste)

set_tests_properties(TestPyImports PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/pychaste")
set_tests_properties(TestPyWrapperChanges PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/pychaste")

# These tests render a VTK scene off-screen. Under parallel ctest they share a
# single Xvfb display, and concurrent renders collide on the older VTK/mesa stack
# (VTK7 on Ubuntu 22.04), giving an intermittent SEGV. A shared RESOURCE_LOCK
# makes the rendering tests mutually exclusive - no two scene renders overlap -
# while the non-rendering tests keep running in parallel. (RUN_SERIAL would
# instead pause the whole run for each of them.)
set_tests_properties(
    # visualization
    TestPyJupyterNotebookManager
    TestVtkSceneWithCaBased
    TestVtkSceneWithMeshBased
    TestVtkSceneWithPottsBased
    # cell populations
    TestPyCaBasedCellPopulation
    TestPyMeshBasedCellPopulation
    TestPyNodeBasedCellPopulation
    TestPyVertexBasedCellPopulation
    # tutorials
    TestPyCellSortingTutorial
    TestPyImmersedBoundaryTutorial
    TestPyMeshBasedCellSimulationsTutorial
    TestPyNodeBasedCellSimulationsTutorial
    TestPyPottsBasedCellSimulationsTutorial
    TestPyScratchAssayTutorial
    TestPySpheroidTutorial
    TestPyTensileTestTutorial
    TestPyVertexBasedCellSimulationsTutorial
    PROPERTIES RESOURCE_LOCK vtk_offscreen_display
)
