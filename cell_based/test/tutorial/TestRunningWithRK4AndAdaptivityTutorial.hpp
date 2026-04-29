/*

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

*/
/*
 *
 *  Chaste tutorial - this page gets automatically changed to a wiki page
 *  DO NOT remove the comments below, and if the code has to be changed in
 *  order to run, please check the comments are still accurate
 *
 *
 */

#ifndef TESTRUNNINGWITHRK4ANDADAPTIVITYTUTORIAL_HPP_
#define TESTRUNNINGWITHRK4ANDADAPTIVITYTUTORIAL_HPP_

/*
 * ## Using RK4 and adaptive time-stepping in off-lattice simulations
 *
 * ### Introduction
 *
 * This tutorial demonstrates how to use the `RK4NumericalMethod` and the **adaptive
 * timestep** feature with off-lattice cell-based simulations in Chaste.
 *
 * Off-lattice simulations advance node positions by solving the overdamped equations
 * of motion
 *
 *   dr/dt = F(r) / nu
 *
 * where F(r) is the net force on a cell and nu is a drag coefficient. Two integration
 * schemes are available:
 *
 *  * `ForwardEulerNumericalMethod` (the default): a first-order explicit scheme.
 *  * `RK4NumericalMethod`: the classic 4th-order Runge-Kutta scheme, which evaluates
 *    forces at four intermediate positions per step and combines them as
 *
 *      k1 = F(r^t) / nu
 *      k2 = F(r^t + dt*k1/2) / nu
 *      k3 = F(r^t + dt*k2/2) / nu
 *      k4 = F(r^t + dt*k3) / nu
 *      r^(t+1) = r^t + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
 *
 * #### Adaptive time-stepping
 *
 * When cells overlap or a node moves too far in a single step, the numerical method
 * calls `AbstractOffLatticeCellPopulation::CheckForStepSizeException()`, which throws a
 * `StepSizeException` carrying a suggested smaller step size.
 *
 * The control loop inside `OffLatticeSimulation::UpdateCellLocationsAndTopology()` then
 *
 *  1. reverts all node positions to their state at the start of the sub-step, and
 *  2. retries the step with the suggested smaller sub-step size.
 *
 * Each macro time step `dt` may therefore be completed in several sub-steps that sum to
 * `dt`. If more than `mMaxAdaptiveTimeSteps` consecutive sub-steps fail, the simulation
 * terminates with an error.
 *
 * Adaptive time-stepping is enabled by calling
 *
 *   p_numerical_method->SetUseAdaptiveTimestep(true);
 *
 * The maximum number of consecutive failed sub-steps can be set with
 *
 *   simulator.SetMaxAdaptiveTimeStep(n);
 *
 * ### The test
 *
 * We begin by including the necessary header files.
 */
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

/* The following header is included in all cell-based test suites. */
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

/* Standard cell-based simulation headers. */
#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"

/* The following header gives us the RK4 numerical method used in this tutorial. */
#include "RK4NumericalMethod.hpp"

/*
 * Next, we define the test class.
 */
class TestRunningWithRK4AndAdaptivityTutorial : public AbstractCellBasedTestSuite
{
public:

    /*
     * ### Test 1 - a node-based simulation using `RK4NumericalMethod`
     *
     * In the first test we run a simple node-based simulation and replace the default
     * Forward Euler integrator with the 4th-order Runge-Kutta (`RK4`) method.
     * The simulation is otherwise identical to the basic node-based tutorial so that
     * the substitution is as clear as possible.
     */
    void TestRK4NodeBasedSimulation()
    {
        /** The next line is needed because HoneycombMeshGenerator is not designed to be run in parallel */
        EXIT_IF_PARALLEL;

        /* First we generate a small honeycomb mesh to use as the generating mesh for a
         * `NodesOnlyMesh`. We choose a cut-off radius of 1.5 cell diameters to define
         * which pairs of cells interact via the spring force. */
        HoneycombMeshGenerator generator(3, 3);
        boost::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        /* Create one cell per node using a uniform cell-cycle model. */
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        /* Assemble the cell population. */
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        /* Set up the simulator. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedRK4");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(2.0);

        /* **Switch to RK4.**
         *
         * We create an `RK4NumericalMethod` object using `boost::make_shared` and pass
         * it to the simulator.  This replaces the default Forward Euler method for the
         * entire simulation. */
        boost::shared_ptr<RK4NumericalMethod<2> > p_numerical_method
            = boost::make_shared<RK4NumericalMethod<2> >();
        simulator.SetNumericalMethod(p_numerical_method);

        /* Add a spring force between neighbouring cells. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        /* Run the simulation. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * If different simulation input parameters are being explored the lines should be removed. */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 9u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 2.0, 1e-10);
    }

    /*
     * ### Test 2 - adaptive time-stepping with `RK4NumericalMethod`
     *
     * In this test we demonstrate the **adaptive timestep** feature.  We enable
     * adaptivity so that if a node tries to move further than allowed in a single
     * sub-step, the solver automatically reduces the sub-step, reverts the node
     * positions, and retries.
     *
     * #### How the adaptive loop works
     *
     * Inside `OffLatticeSimulation::UpdateCellLocationsAndTopology()` each macro
     * time step `dt` is completed by a `while` loop:
     *
     * ```
     * time_advanced = 0
     * sub_step = dt              // start with the full macro step
     * adaptive_timer = 0         // counts consecutive failures
     *
     * while (time_advanced < dt):
     *
     *     try:
     *         move nodes by sub_step             // may throw StepSizeException
     *         apply boundary conditions
     *         time_advanced += sub_step
     *         adaptive_timer = 0                 // reset counter after success
     *
     *     catch StepSizeException e:
     *         if adaptive_timer < max_adaptive_steps:
     *             revert nodes to saved positions
     *             sub_step = e.GetSuggestedNewStep()
     *             adaptive_timer += 1
     *         else:
     *             throw                          // too many failures, abort
     * ```
     *
     * The `StepSizeException` is raised inside the numerical method whenever
     * `AbstractOffLatticeCellPopulation::CheckForStepSizeException()` detects that a node
     * would move further than its permitted maximum displacement.  The exception
     * carries a smaller suggested step size.
     */
    void TestAdaptiveTimestepWithRK4()
    {
        /** The next line is needed because HoneycombMeshGenerator is not designed to be run in parallel */
        EXIT_IF_PARALLEL;

        /* Build a small node-based mesh. */
        HoneycombMeshGenerator generator(3, 3);
        boost::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        /* Create one cell per node. */
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        /* Set up the simulator with a moderately large time step. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("NodeBasedRK4Adaptive");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(2.0);
        simulator.SetDt(0.1);

        /* Create an RK4 numerical method and **enable adaptive time-stepping**.
         *
         * With `SetUseAdaptiveTimestep(true)` the solver will automatically
         * reduce the sub-step whenever a `StepSizeException` is detected,
         * revert the node positions to those saved at the start of the sub-step,
         * and retry with the smaller step. */
        boost::shared_ptr<RK4NumericalMethod<2> > p_numerical_method
            = boost::make_shared<RK4NumericalMethod<2> >();
        p_numerical_method->SetUseAdaptiveTimestep(true);
        simulator.SetNumericalMethod(p_numerical_method);

        /* Optionally increase the maximum number of consecutive failed sub-steps
         * from the default of 5.  Here we allow up to 10 retries so that the
         * simulation can recover from larger overlaps if needed. */
        simulator.SetMaxAdaptiveTimeStep(10);

        /* Add a spring force. */
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        /* Run the simulation.  The adaptive loop will silently reduce the sub-step
         * whenever cells come too close and retry until the step succeeds. */
        simulator.Solve();

        /* The next two lines are for test purposes only and are not part of this tutorial.
         * If different simulation input parameters are being explored the lines should be removed. */
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 9u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 2.0, 1e-10);
    }
};

/*
 * To visualize the results, use Paraview. See the [Visualizing With Paraview](../visualizingwithparaview/) tutorial for more information.
 *
 * Load the file `$CHASTE_TEST_OUTPUT/NodeBasedRK4/results_from_time_0/results.pvd`
 * and add glyphs to represent cells.
 */

#endif /* TESTRUNNINGWITHRK4ANDADAPTIVITYTUTORIAL_HPP_ */
