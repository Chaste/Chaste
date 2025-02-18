/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYDEMOTUTORIAL_HPP_
#define TESTIMMERSEDBOUNDARYDEMOTUTORIAL_HPP_

/*
 * ## Example showing how to create and run an immersed boundary simulation in Chaste
 *
 * We create a simple palisade of cells with a basement membrane, and see how to:
 *   * set the initial conditions;
 *   * change the cell-level properties.
 *
 * ### The test
 *
 * We begin by including the necessary header files.
 */

/* Required for the test environment */
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

/* Required for core elements of the Chaste simulation environment */
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"

/* Required for the immersed boundary functionality */
#include "ImmersedBoundaryLinearInteractionForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"

/* Required for setting up the numerical method */
#include "ForwardEulerNumericalMethod.hpp"
#include <boost/make_shared.hpp>

/* This test is never run in parallel */
#include "FakePetscSetup.hpp"

/* In Chaste, every simulation is run as a 'test', and here we define a test class which inherits from
 * `AbstractCellBasedTestSuite`. The class `AbstractCellBasedTestSuite` sets up various parameters for us.  Of
 * particular use is that `RandomNumberGenerator` is re-seeded with zero.  This means any random numbers generated
 * (for instance, random variation in cell size) is reproducible from one simulation to the next.
 */
class TestRunningImmersedBoundarySimulationsTutorial : public AbstractCellBasedTestSuite
{
public:
    /*
     * ## Simulation - a basic immersed boundary simulation
     *
     * Immersed boundary methods simulate two-way coupled cell-fluid interactions. They can be used to model
     * cells within fluids, and fluids within cells. The chaste implementation supports multiple fluid sources,
     * cell wall shape changes, laminas, leaky laminas and biological noise.
     *
     * In this simulation, we create an immersed boundary framework with a palisade of cells. The cells have a slight
     * variation in size, and a basement membrane is included. Each cell in the simulation is assigned a basic
     * cell-cycle model; no proliferation occurs.
     */
    void TestImmersedBoundaryPalisadeSimulation()
    {
        /* The first thing we define is a 2D (specified by the <2,2>) mesh. This holds spatial information of the
         * simulation, including that of the underlying fluid grid as well as the location of cell boundary-points. To
         * create this we use an `ImmersedBoundaryPalisadeMeshGenerator`. The parameters are, in order:
         *   * 5: number of cells in the palisade
         *   * 100: number of boundary points in each cell
         *   * 0.2: 'superellipse' exponent which determines initial cell shape; 1.0 is an ellipse, 0.0 a rectangle
         *   * 2.0: the initial aspect ratio of each cell; they start twice as high as their width
         *   * 0.15: the proportion of random variation in cell heights
         *   * true: whether the mesh should contain a basement membrane
         */

        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();

        /* We now generate a collection of cells. We do this by using a `CellsGenerator` and we specify the
         * proliferative behaviour of the cell by choosing a `CellCycleModel`. Here we choose an
         * `UniformCellCycleModel` which does not allow proliferation. For an immersed boundary
         * simulation we need as may cells as elements in the mesh. */
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        /* We now create a `CellPopulation` object (passing in the mesh and cells) to connect the mesh and the cells
         * together. Here we use an `ImmersedBoundaryCellPopulation` and the dimension is <2>.*/
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        /* We now create an `OffLatticeSimulation` object and pass in the `CellPopulation`. We also set some
         * options for the simulation like output directory, output multiple (so we don't visualize every timestep),
         * and end time.
         * Additionally, we tell the numerical method that we want the cell population to update node locations.*/
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2,2> >());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

        double dt = 0.01;
        simulator.SetOutputDirectory("TestImmersedBoundaryDemoTutorial");
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(100.0 * dt);

        /* All of the machinery for the immersed boundary method is handled in the following `SimulationModifier`.
         * Here, we create a 'shared pointer' to an `ImmersedBoundarySimulationModifier` object and pass it to the
         * `OffLatticeSimulation`.*/
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        /* We now associate an `ImmersedBoundaryLinearMembraneForce` and
         * `ImmersedBoundaryLinearInteractionForce` to the `SimulationModifier` which
         * handles the membrane elasticity forces.  These are created in a similar manner as above.*/
        MAKE_PTR(ImmersedBoundaryLinearMembraneForce<2>, p_boundary_force);
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetElementSpringConst(0.5 * 1e8);
        p_boundary_force->SetLaminaSpringConst(1.0 * 1e8);

        MAKE_PTR(ImmersedBoundaryLinearInteractionForce<2>, p_cell_cell_force);
        p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
        p_cell_cell_force->SetSpringConst(1.0 * 1e6);

        /* Finally we call the `Solve` method on the simulation to run the simulation.*/
        simulator.Solve();
    }
};

#endif /*TESTIMMERSEDBOUNDARYDEMOTUTORIAL_HPP_*/
