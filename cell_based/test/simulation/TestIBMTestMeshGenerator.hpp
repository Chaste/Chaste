/*

Copyright (c) 2005-2015, University of Oxford.
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
 * TestIBMMouseEmbryoSimulation.hpp
 *
 *  Created on: 29 Apr 2016
 *      Author: Bartosz Bartmanski
 */

#ifndef TESTIBMTESTMESHGENERATOR_HPP_
#define TESTIBMTESTMESHGENERATOR_HPP_

/* The following are used in every Chaste test */
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

/* Required for the core elements of the Chaste simulation environment */
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"

/* Required for the Immersed Boundary functionality */
#include "ImmersedBoundaryLinearInteractionForce.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundarySingleCellMigrationForce.hpp"

/* Required for setting up the numerical method */
#include "ForwardEulerNumericalMethod.hpp"
#include <boost/make_shared.hpp>

/* To make sure that the test is not run in parallel */
#include "FakePetscSetup.hpp"

/* Required specifically for this simulation */
#include "ImmersedBoundaryTestMeshGenerator.hpp"
#include "ImmersedBoundaryCellSizeWriter.hpp"

class TestIBMTestMeshGenerator : public AbstractCellBasedTestSuite
{
public:
	/*
	 * == Mouse Embryo simulation - simulation of the Distal half of the Visceral Endoderm cells ==
	 *
	 * We will model the VE cells of a mouse embryo using the Immersed Boundary method,
	 * which will be placed on a membrane to keep these ordered and having the required shape.
	 * As a first approximation, the VE cells will be on a circular basement membrane, with the
	 * lower half representing the VE cells around the epiblast
	 */
	void TestMeshGenerator() throw (Exception)
	{
        /*
         * First use the mesh generator to set up the immersed boundary elements
         */
		ImmersedBoundaryTestMeshGenerator generator(0.5, 0.5, 100, 1.0);
		ImmersedBoundaryMesh<2,2>* p_mesh = generator.GetMesh();

		/*
		 * In the next three lines, we generate cells that will be used in the simulation
		 * We can choose the the proliferative type and Cell Cycle model.
		 * Apparently {{{UniformlyDistributedCellCycleModel}}} doesn't allow proliferation, which doesn't matter for this simulation, as on the timescales we're looking at,
		 * there is very little proliferation
		 */
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        /* We now create a {{{CellPopulation}}} object (passing in the mesh and cells) to connect the mesh and the cells
         * together. Here we use an {{{ImmersedBoundaryCellPopulation}}} and the dimension is <2>.*/
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<ImmersedBoundaryCellSizeWriter>();

        /* We now create an {{{OffLatticeSimulation}}} object and pass in the {{{CellPopulation}}}. We also set some
         * options for the simulation like output directory, output multiple (so we don't visualize every timestep),
         * and end time.
         * Additionally, we tell the numerical method that we want the cell population to update node locations.*/
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2,2> >());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

        double dt = 0.01;
        simulator.SetOutputDirectory("TestIBMTestMeshGenerator");
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(200.0 * dt);

        /* All of the machinery for the immersed boundary method is handled in the following {{{SimulationModifier}}}.
         * Here, we create a 'shared pointer' to an {{{ImmersedBoundarySimulationModifier}}} object and pass it to the
         * {{{OffLatticeSimulation}}}.*/
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        /* We now associate an {{{ImmersedBoundaryLinearMembraneForce}}} and
         * {{{ImmersedBoundaryLinearInteractionForce}}} to the {{{SimulationModifier}}} which
         * handles the membrane elasticity forces.  These are created in a similar manner as above.*/
        MAKE_PTR(ImmersedBoundaryLinearMembraneForce<2>, p_boundary_force);
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetElementSpringConst(0.5 * 1e8);
        p_boundary_force->SetLaminaSpringConst(1.0 * 1e8);

        MAKE_PTR(ImmersedBoundaryLinearInteractionForce<2>, p_cell_cell_force);
        p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
        p_cell_cell_force->SetSpringConst(1.0 * 1e6);

        MAKE_PTR(SingleCellMigrationForce<2>, p_migration_force);
        p_main_modifier->AddImmersedBoundaryForce(p_migration_force);
        p_migration_force->SetElementIndex(0);
        p_migration_force->SetStrength(1.0);

        /* Finally we call the {{{Solve}}} method on the simulation to run the simulation.*/
        simulator.Solve();


	}

};

#endif /* TESTIBMTESTMESHGENERATOR_HPP_ */
