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

#ifndef TEST3DOFFLATTICEPRESENTATIVESIMULATION_HPP_
#define TEST3DOFFLATTICEPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "TrianglesMeshReader.hpp"
#include "OffLatticeSimulation.hpp"
#include "TrianglesMeshWriter.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "CellAncestorWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "FakePetscSetup.hpp"

/**
 * This class consists of a single test, in which a 3D model
 * of a cell population is simulated. There is not yet any division
 * in the model, and all cells have the same proliferative state
 * (differentiated) and mutation state (wildtype).
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class Test3dTissueRepresentativeSimulation : public AbstractCellBasedTestSuite
{
public:

    /*
     * Create and simulate a simple 3D cell population with a cuboid
     * mesh, with ghost nodes around the outside
     */
    void Test3DHoneycombMeshWithGhostNodes()
    {
        /*          _ _ _ _ _
         *        /        /|
         *       /        / |
         *         /_ _ _ _ /  | depth (z-direction)
         *        |         |  |
         *     |         |  |
         *     |         |  /
         *     |         | / height (y-direction)
         *        |_ _ _ _ _|/
         *        width
         *    (x-direction)
         */

        unsigned nodes_across = 12;  // Corresponds to x axis
        unsigned nodes_depth = 12;   // Corresponds to y axis
        unsigned nodes_up = 12;        // Corresponds to z axis
        unsigned ghosts = 2;        // Layers of ghost nodes

        /*
         *    *     *     *     *
         *   / \   / \   / \   / \       (1 layer of ghost nodes = 1 honeycomb layer)
         *  /   \ /   \ /   \ /   \
         * *     *     *     *     *
         *
         */

        std::vector<Node<3>*> nodes;

        double x_coordinate, y_coordinate, z_coordinate;
        std::vector<unsigned> ghost_node_indices, real_node_indices;
        unsigned node_index = 0;

        for (unsigned k=0; k<2*nodes_up-1; k++)        // Each layer going up
        {
            z_coordinate = (double)k/2.0;

            bool is_even_layer = ((int)k % 2 == 0);

            if (is_even_layer)
            {
                // Want the nodes that sit at x=0,2,4,...
                for (unsigned j=0; j<nodes_depth; j++)
                {
                    y_coordinate = (double)j;
                    x_coordinate = 0.0;

                    for (unsigned i=0; i<nodes_across; i++)
                    {
                        nodes.push_back(new Node<3>(node_index,  false,  x_coordinate, y_coordinate, z_coordinate));

                        if (x_coordinate < ghosts || x_coordinate > (double)nodes_across-1.0-ghosts
                                || y_coordinate < ghosts || y_coordinate > (double)nodes_depth-1.0-ghosts
                                || z_coordinate < ghosts|| z_coordinate > (double)nodes_up-1.0-ghosts)
                        {
                            ghost_node_indices.push_back(node_index);
                        }
                        else
                        {
                            real_node_indices.push_back(node_index);
                        }

                        node_index++;
                        x_coordinate += 1.0;
                    }
                }
            }
            else
            {
                // Want the nodes that sit at x=1,3,5,...
                for (unsigned j=0; j<nodes_depth-1; j++)
                {
                    y_coordinate = (double)j + 0.5;
                    x_coordinate = 0.5;

                    for (unsigned i=0; i<nodes_across-1; i++)
                    {
                        nodes.push_back(new Node<3>(node_index,  false,  x_coordinate, y_coordinate, z_coordinate));

                        if (x_coordinate < ghosts || x_coordinate > (double)nodes_across-1.0-ghosts
                                || y_coordinate < ghosts || y_coordinate > (double)nodes_depth-1.0-ghosts
                                || z_coordinate < ghosts|| z_coordinate > (double)nodes_up-1.0-ghosts)
                        {
                            ghost_node_indices.push_back(node_index);
                        }
                        else
                        {
                            real_node_indices.push_back(node_index);
                        }

                        node_index++;
                        x_coordinate += 1.0;
                    }
                }
            }
        }

        MutableMesh<3,3> mesh(nodes);

        unsigned total_nodes_including_ghosts = nodes_up*nodes_across*nodes_depth + (nodes_up-1)*(nodes_across-1)*(nodes_depth-1);
        TS_ASSERT_EQUALS(mesh.GetNumNodes(),total_nodes_including_ghosts);
        TS_ASSERT_LESS_THAN(real_node_indices.size(), mesh.GetNumNodes());

        unsigned across_without_ghosts = nodes_across-2*ghosts;
        unsigned depth_without_ghosts = nodes_depth-2*ghosts;
        unsigned up_without_ghosts = nodes_up-2*ghosts;

        unsigned total_nodes_without_ghosts = up_without_ghosts*across_without_ghosts*depth_without_ghosts + (up_without_ghosts-1)*(across_without_ghosts-1)*(depth_without_ghosts-1);
        TS_ASSERT_EQUALS(total_nodes_without_ghosts, real_node_indices.size());

        // Now only assign cells to real_node_indices

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (std::vector<unsigned>::iterator real_node_iter=real_node_indices.begin();
                                            real_node_iter != real_node_indices.end();
                                            ++real_node_iter)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            cells.push_back(p_cell);
        }

        TS_ASSERT_EQUALS(real_node_indices.size(), cells.size());

        MeshBasedCellPopulationWithGhostNodes<3> cell_population(mesh, cells, real_node_indices);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        TS_ASSERT_EQUALS(ghost_node_indices.size(), cell_population.GetGhostNodeIndices().size());

        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("Test3DHoneycombBoxMeshWithGhostNodes");
        simulator.SetEndTime(1.0);
        simulator.SetSamplingTimestepMultiple(12);

        // Create a force law and pass it to the OffLatticeSimulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        simulator.Solve();
    }
};

#endif /*TEST3DOFFLATTICEPRESENTATIVESIMULATION_HPP_*/
