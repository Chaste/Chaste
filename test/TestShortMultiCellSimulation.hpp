/*

Copyright (c) 2005-2014, University of Oxford.
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

// Needed for the test environment
#include <cxxtest/cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// External library - not part of Chaste
#include <fftw3.h>

#include "OffLatticeSimulation.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"
#include "ImmersedBoundaryCellCellInteractionForce.hpp"
#include "Timer.hpp"

#include "Debug.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestShortMultiCellSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestShortMultiCellSim() throw(Exception)
    {
        /*
         * 1: Num cells
         * 2: Num nodes per cell
         * 3: Superellipse exponent
         * 4: Superellipse aspect ratio
         * 5: Random y-variation
         * 6: Include membrane
         */
        ImmersedBoundaryPalisadeMeshGenerator gen(7, 128, 0.1, 2.5, 0.0, true);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

        p_mesh->SetNumGridPtsXAndY(256);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(false);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetConsoleProgressOutput(true);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        MAKE_PTR_ARGS(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force, (cell_population));
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetSpringConstant(1.0 * 1e7);

        MAKE_PTR_ARGS(ImmersedBoundaryCellCellInteractionForce<2>, p_cell_cell_force, (cell_population));
        p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
        p_cell_cell_force->SetSpringConstant(1.0 * 1e6);



        // Get height of basement lamina
        double lamina_height = 0.0;
        for (unsigned node_idx = 0 ; node_idx < p_mesh->GetElement(0)->GetNumNodes() ; node_idx++)
        {
            lamina_height += p_mesh->GetElement(0)->GetNode(node_idx)->rGetModifiableLocation()[1];
        }
        lamina_height /= p_mesh->GetElement(0)->GetNumNodes();


        // Kick the second cell in from the left and set its E-cad level
        double kick = 1.05;
//        unsigned e_cad_location = p_cell_cell_force->rGetProteinNodeAttributeLocations()[0];
//        unsigned p_cad_location = p_cell_cell_force->rGetProteinNodeAttributeLocations()[1];

//        double centroid_x = p_mesh->GetCentroidOfElement(3)[0];

        for (unsigned node_idx = 0 ; node_idx < p_mesh->GetElement(3)->GetNumNodes() ; node_idx++)
        {
            double new_height = lamina_height + kick * (p_mesh->GetElement(3)->GetNode(node_idx)->rGetLocation()[1] - lamina_height);

            p_mesh->GetElement(3)->GetNode(node_idx)->rGetModifiableLocation()[1] = new_height;

//            if (p_mesh->GetElement(3)->GetNode(node_idx)->GetRegion() == 2 && p_mesh->GetElement(3)->GetNode(node_idx)->rGetLocation()[0] > centroid_x)
//            {
//                p_mesh->GetElement(3)->GetNode(node_idx)->rGetNodeAttributes()[e_cad_location] = 0.0;
//            }
//            else
//            {
//                p_mesh->GetElement(3)->GetNode(node_idx)->rGetNodeAttributes()[e_cad_location] = 0.0;
//            }
//
//            p_mesh->GetElement(3)->GetNode(node_idx)->rGetNodeAttributes()[e_cad_location] = 0.0;
//
//            p_mesh->GetElement(3)->GetNode(node_idx)->rGetNodeAttributes()[p_cad_location] = 2.0;

        }
//
//        for (unsigned node_idx = 0 ; node_idx < p_mesh->GetElement(4)->GetNumNodes() ; node_idx++)
//        {
//            if (p_mesh->GetElement(4)->GetNode(node_idx)->rGetLocation()[1] > 0.6)
//            {
//                p_mesh->GetElement(4)->GetNode(node_idx)->rGetNodeAttributes()[p_cad_location] = 2.0;
//            }
//        }

        // Set simulation properties
        double dt = 0.05;
        simulator.SetOutputDirectory("TestShortMultiCellSimulation");
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(5000.0 * dt);
        simulator.Solve();

//        // Add a cell-cell interaction force with the same intrinsic strength as the membrane force
//        MAKE_PTR_ARGS(ImmersedBoundaryCellCellInteractionForce<2>, p_cell_cell_force, (cell_population));
//        p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
//        p_cell_cell_force->SetSpringConstant(2.0 * 1e6);

    }
};
