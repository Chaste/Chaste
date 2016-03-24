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
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#include "OffLatticeSimulation.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"

#include "Debug.hpp"
#include "Timer.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestShortSingleCellSimulation : public AbstractCellBasedTestSuite
{
public:

    void TestShortSingleCellSim() throw(Exception)
    {
        std::vector<c_vector<double, 2> > offsets(9);
        offsets[0][0] = 0.35; offsets[0][1] = 0.3;
        offsets[1][0] = 0.35; offsets[1][1] = 0.6;
        offsets[2][0] = 0.35; offsets[2][1] = 0.9;
        offsets[3][0] = 0.50; offsets[3][1] = 0.15;
        offsets[4][0] = 0.50; offsets[4][1] = 0.45;
        offsets[5][0] = 0.50; offsets[5][1] = 0.75;
        offsets[6][0] = 0.65; offsets[6][1] = 0.3;
        offsets[7][0] = 0.65; offsets[7][1] = 0.6;
        offsets[8][0] = 0.65; offsets[8][1] = 0.9;

        std::vector<c_vector<double, 2> > node_locations(32);
        for (unsigned location = 0 ; location < node_locations.size() ; location++)
        {
            double theta = 2.0 * M_PI * (double)location / (double)node_locations.size();
            node_locations[location][0] = cos(theta);
            node_locations[location][1] = sin(theta);
        }

        std::vector<Node<2>*> nodes;
        std::vector<ImmersedBoundaryElement<2,2>*> elems;

        for (unsigned offset = 0 ; offset < offsets.size() ; offset++)
        {
            std::vector<Node<2>*> nodes_this_elem;

            for (unsigned location = 0 ; location < node_locations.size() ; location++)
            {
                unsigned index = offset * node_locations.size() + location;
                Node<2>* p_node = new Node<2>(index, offsets[offset] + 0.0725 * node_locations[location], true);

                nodes_this_elem.push_back(p_node);
                nodes.push_back(p_node);
            }

            ImmersedBoundaryElement<2,2>* p_elem = new ImmersedBoundaryElement<2,2>(offset, nodes_this_elem);
            p_elem->rGetCornerNodes().push_back(nodes_this_elem[0]);
            p_elem->rGetCornerNodes().push_back(nodes_this_elem[0]);
            p_elem->rGetCornerNodes().push_back(nodes_this_elem[0]);
            p_elem->rGetCornerNodes().push_back(nodes_this_elem[0]);

            elems.push_back(p_elem);
        }
MARK;
        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2,2>(nodes, elems);
MARK;
        p_mesh->SetNumGridPtsXAndY(64);
MARK;
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
MARK;
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        MARK;
        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);
        MARK;
        // Add force laws
        MAKE_PTR_ARGS(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force, (cell_population));
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetSpringConstant(0.5 * 1e7);
        MARK;



        PRINT_2_VARIABLES(p_mesh->GetNumNodes(), p_mesh->GetNumElements());

        // Set simulation properties
        simulator.SetOutputDirectory("TestShortSingleCellSimulation");
        MARK;
        double dt = 0.05;
        MARK;
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(100.0 * dt);
        MARK;
        simulator.Solve();
        MARK;
        delete(p_mesh);
    }
};