"""Copyright (c) 2005-2026, University of Oxford.
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

# ifndef
# define TRIGGER_WIKI

## ## Introduction
##
## This tutorial demonstrates how to use the `RK4NumericalMethod` and the
## **adaptive timestep** feature with off-lattice cell-based simulations in Chaste.
##
## Chaste's off-lattice simulations advance node positions by solving the
## equations of motion
##
##     dr/dt = F(r) / nu
##
## where `F(r)` is the net force on a cell and `nu` is a damping coefficient.
## Two numerical integration methods are available:
##
## - **`ForwardEulerNumericalMethod`** (default): a first-order explicit scheme.
##   Simple and fast, but can require very small time steps to remain stable in
##   stiff problems.
## - **`RK4NumericalMethod`**: the classic 4th-order Runge-Kutta scheme.
##   Evaluates forces at four intermediate positions per step, giving higher
##   accuracy at the same time-step size or allowing larger stable time steps
##   for smooth problems.
##
## ### Adaptive Timestep
##
## When cells overlap or move too far in a single step, Chaste detects this by
## checking that no node has displaced more than a fraction of its cell radius
## (this check is performed inside `AbstractNumericalMethod::DetectStepSizeExceptions`
## and delegates to `AbstractCellPopulation::CheckForStepSizeException`).
##
## If the displacement is too large, a `StepSizeException` is thrown carrying a
## suggested smaller step size.  The adaptive control loop inside
## `OffLatticeSimulation::UpdateCellLocationsAndTopology` catches this exception,
## **reverts all node positions to their state at the beginning of the attempted
## sub-step**, and retries with the suggested smaller step.
##
## The outer loop advances time in sub-steps that sum to the user-supplied `dt`
## (the macro time step).  On a **successful** sub-step, the counter is reset and
## the sub-step size is gently increased (by 1%) so that the scheme can recover
## toward the full macro step size.  If the sub-step fails more than
## `max_adaptive_steps` times in a row, the simulation terminates with an error
## instead of looping forever.
##
## To enable adaptive time stepping:
##
## 1. Call `numerical_method.SetUseAdaptiveTimestep(True)` on the numerical
##    method object.
## 2. Optionally call `simulator.SetMaxAdaptiveTimeStep(n)` to control how many
##    consecutive failed sub-steps are allowed before the simulation aborts
##    (default: 5).
##
## ## The Test

import unittest

import chaste
import chaste.cell_based
import chaste.mesh
import chaste.visualization

from chaste.cell_based import AbstractCellBasedTestSuite


class TestPyNumericalMethodsTutorial(AbstractCellBasedTestSuite):

    ## ### Test 1 – using `RK4NumericalMethod` for a node-based simulation
    ##
    ## In this test we show how to replace the default Forward Euler integrator
    ## with the 4th-order Runge-Kutta (`RK4`) method.  The simulation is otherwise
    ## identical to the basic node-based tutorial so that the substitution is as
    ## clear as possible.

    def test_rk4_node_based_simulation(self):

        # JUPYTER_SETUP

        ## Set up the output directory and generate a small honeycomb mesh.

        chaste.core.OutputFileHandler("Python/TestNumericalMethodsTutorial/Rk4")
        generator = chaste.mesh.HoneycombMeshGenerator(3, 3)
        generating_mesh = generator.GetMesh()

        ## Build a `NodesOnlyMesh` from the generating mesh, specifying an
        ## interaction (cut-off) radius of 1.5 cell diameters.

        mesh = chaste.mesh.NodesOnlyMesh[2]()
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5)

        ## Create cells — one per node — using a uniform cell-cycle model.

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 2]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), transit_type)

        ## Assemble the cell population.

        cell_population = chaste.cell_based.NodeBasedCellPopulation[2](mesh, cells)

        ## Set up the off-lattice simulator.

        simulator = chaste.cell_based.OffLatticeSimulation[2, 2](cell_population)
        simulator.SetOutputDirectory("Python/TestNumericalMethodsTutorial/Rk4")
        simulator.SetSamplingTimestepMultiple(12)
        simulator.SetEndTime(2.0)

        ## **Switch to RK4.**
        ##
        ## We create an `RK4NumericalMethod` object and pass it to the simulator
        ## with `SetNumericalMethod`.  This replaces the default Forward Euler
        ## method for the entire simulation.
        ##
        ## The RK4 scheme evaluates forces at four intermediate positions per
        ## step (k1 through k4) and combines them as
        ##
        ##     r^(t+1) = r^t + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
        ##
        ## giving 4th-order accuracy in the time step `dt`.

        numerical_method = chaste.cell_based.RK4NumericalMethod[2, 2]()
        simulator.SetNumericalMethod(numerical_method)

        ## Add a spring force between neighbouring cells.

        force = chaste.cell_based.GeneralisedLinearSpringForce[2, 2]()
        simulator.AddForce(force)

        ## Run the simulation.

        simulator.Solve()

        ## Verify the final cell count and simulation time.

        self.assertEqual(cell_population.GetNumRealCells(), 9)
        self.assertAlmostEqual(
            chaste.cell_based.SimulationTime.Instance().GetTime(), 2.0, 6
        )

        # JUPYTER_TEARDOWN

    ## ### Test 2 – adaptive time stepping with `RK4NumericalMethod`
    ##
    ## In this test we demonstrate the **adaptive timestep** feature.  We start
    ## with a configuration that would cause cells to overlap badly at the
    ## default time-step size, then let the adaptive control loop shrink the
    ## sub-step automatically until the step succeeds.
    ##
    ## #### How adaptivity works
    ##
    ## Each macro time step `dt` may be broken into several smaller sub-steps.
    ## The control loop in `OffLatticeSimulation::UpdateCellLocationsAndTopology`
    ## works as follows:
    ##
    ## ```
    ## time_advanced = 0
    ## sub_step = dt            # start optimistically with the full step
    ## adaptive_timer = 0       # counts consecutive failures
    ##
    ## while time_advanced < dt:
    ##     try:
    ##         move all nodes by sub_step             # can throw StepSizeException
    ##         apply boundary conditions
    ##         time_advanced += sub_step              # success
    ##         adaptive_timer = 0                     # reset failure counter
    ##         sub_step = min(1.01 * sub_step,        # gently increase
    ##                        dt - time_advanced)
    ##     except StepSizeException as e:
    ##         if adaptive_timer < max_adaptive_steps:
    ##             revert nodes to start-of-sub-step positions
    ##             sub_step = e.suggested_step        # smaller step from exception
    ##             adaptive_timer += 1
    ##         else:
    ##             raise                              # too many failures → abort
    ## ```
    ##
    ## The `StepSizeException` is raised inside the numerical method when
    ## `AbstractCellPopulation::CheckForStepSizeException` finds that a node
    ## has moved further than its permitted maximum displacement.  The exception
    ## carries a suggested new step size that should keep the displacement within
    ## bounds.
    ##
    ## To enable this behaviour you must call
    ## `numerical_method.SetUseAdaptiveTimestep(True)`.  You can also increase the
    ## number of allowed consecutive failures with
    ## `simulator.SetMaxAdaptiveTimeStep(n)` (default is 5).

    def test_adaptive_timestep_with_rk4(self):

        # JUPYTER_SETUP

        ## Set up the output directory and build a small mesh.

        chaste.core.OutputFileHandler(
            "Python/TestNumericalMethodsTutorial/AdaptiveRk4"
        )
        generator = chaste.mesh.HoneycombMeshGenerator(3, 3)
        generating_mesh = generator.GetMesh()

        mesh = chaste.mesh.NodesOnlyMesh[2]()
        mesh.ConstructNodesWithoutMesh(generating_mesh, 1.5)

        ## Create cells.

        transit_type = chaste.cell_based.TransitCellProliferativeType()
        cell_generator = chaste.cell_based.CellsGenerator["UniformCellCycleModel", 2]()
        cells = cell_generator.GenerateBasicRandom(mesh.GetNumNodes(), transit_type)

        cell_population = chaste.cell_based.NodeBasedCellPopulation[2](mesh, cells)

        ## Set up the simulator.  We use a moderately large time step (`dt = 0.1`)
        ## that would be problematic without adaptivity in certain configurations.

        simulator = chaste.cell_based.OffLatticeSimulation[2, 2](cell_population)
        simulator.SetOutputDirectory(
            "Python/TestNumericalMethodsTutorial/AdaptiveRk4"
        )
        simulator.SetSamplingTimestepMultiple(10)
        simulator.SetEndTime(2.0)
        simulator.SetDt(0.1)

        ## Create an RK4 numerical method and **enable adaptive time stepping**.
        ##
        ## With `SetUseAdaptiveTimestep(True)` the solver will automatically
        ## halve (or further reduce) the sub-step whenever a `StepSizeException`
        ## is detected, then retry the sub-step from the saved node positions.

        numerical_method = chaste.cell_based.RK4NumericalMethod[2, 2]()
        numerical_method.SetUseAdaptiveTimestep(True)
        simulator.SetNumericalMethod(numerical_method)

        ## Optionally increase the maximum number of consecutive adaptive
        ## retries from the default of 5.  Here we allow up to 10 retries so
        ## that the simulation can recover from larger overlaps if needed.

        simulator.SetMaxAdaptiveTimeStep(10)

        ## Add a spring force.

        force = chaste.cell_based.GeneralisedLinearSpringForce[2, 2]()
        simulator.AddForce(force)

        ## Run the simulation.  The adaptive loop will silently reduce the
        ## sub-step whenever cells come too close, and gradually increase it
        ## again after a successful sub-step, so no manual intervention is
        ## required.

        simulator.Solve()

        ## Verify the results.

        self.assertEqual(cell_population.GetNumRealCells(), 9)
        self.assertAlmostEqual(
            chaste.cell_based.SimulationTime.Instance().GetTime(), 2.0, 6
        )

        # JUPYTER_TEARDOWN


if __name__ == "__main__":
    unittest.main(verbosity=2)

# endif END_WIKI
