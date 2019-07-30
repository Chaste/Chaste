/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef TESTREPRESENTATIVE3DNODEBASEDSIMULATION_HPP_
#define TESTREPRESENTATIVE3DNODEBASEDSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"


#include "OffLatticeSimulation.hpp"
#include "RepulsionForce.hpp"
#include "UniformCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "CellAncestorWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "FakePetscSetup.hpp"

/**
 * This class consists of a single test, in which a 3D model
 * of a cell population is simulated. Cells are restrained to a cuboid and
 * cells are removed on reaching the top of the domain.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class Test3dTissueRepresentativeSimulation : public AbstractCellBasedWithTimingsTestSuite
{
public:

    /*
     * Create and simulate a simple 3D cell population of about 1000 cells within a cuboid box with sloughing on the top edge
     */
    void Test3dNodeBasedInBoxWithSloughing()
    {
        double size_of_box = 8.0;
        unsigned cells_across = 12;
        double scaling = size_of_box/(double(cells_across-1));

        // Create a simple 3D NodeBasedCellPopulation consisting of cells evenly spaced in a regular grid
        std::vector<Node<3>*> nodes;
        unsigned index = 0;
        for (unsigned i=0; i<cells_across; i++)
        {
            for (unsigned j=0; j<cells_across; j++)
            {
                for (unsigned k=0; k<cells_across; k++)
                {
                    nodes.push_back(new Node<3>(index, false,  (double) i * scaling , (double) j * scaling, (double) k * scaling));
                    index++;
                }
            }
        }

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<3> node_based_cell_population(mesh, cells);
        //node_based_cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(node_based_cell_population);
        simulator.SetOutputDirectory("Representative3dNodeBased");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(RepulsionForce<3>, p_linear_force);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Create some plane boundary conditions to restrict to box and pass them to the simulation

        // Restrict to x>0
        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition_1, (&node_based_cell_population, zero_vector<double>(3), -unit_vector<double>(3,0)));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_1);

        // Restrict to x < size_of_box
        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition_2, (&node_based_cell_population, size_of_box*unit_vector<double>(3,0), unit_vector<double>(3,0)));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_2);

        // Restrict to y>0
        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition_3, (&node_based_cell_population, zero_vector<double>(3), -unit_vector<double>(3,1)));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_3);

        // Restrict to y < size_of_box
        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition_4, (&node_based_cell_population, size_of_box*unit_vector<double>(3,1), unit_vector<double>(3,1)));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_4);

        // Restrict to z > 0
        MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, p_boundary_condition_5, (&node_based_cell_population, zero_vector<double>(3), -unit_vector<double>(3,2)));
        simulator.AddCellPopulationBoundaryCondition(p_boundary_condition_5);

        // Create cell killer at z= size_of_box
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer,(&node_based_cell_population,  size_of_box*unit_vector<double>(3,2), unit_vector<double>(3,2)));
        simulator.AddCellKiller(p_cell_killer);


        // Run simulation
        simulator.Solve();

        // Check some results
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 1128u);

        AbstractCellPopulation<3>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();

        for (unsigned i=0; i<100; i++)
        {
            ++cell_iter;
        }
        c_vector<double,3> node_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);
        TS_ASSERT_DELTA(node_location[0],0.9604, 1e-4);
        TS_ASSERT_DELTA(node_location[1],8.0, 1e-4);
        TS_ASSERT_DELTA(node_location[2],3.6539, 1e-4);

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /*TESTREPRESENTATIVE3DNODEBASEDSIMULATION_HPP_*/
