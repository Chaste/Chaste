/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTVERTEXCRYPTBOUNDARYFORCE_HPP_
#define TESTVERTEXCRYPTBOUNDARYFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexCryptBoundaryForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"


class TestVertexCryptBoundaryForce : public AbstractCellBasedTestSuite
{
public:

    void TestForceOutputParameters()
    {
        std::string output_directory = "TestForcesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with VertexCryptBoundaryForce
        VertexCryptBoundaryForce<2> boundary_force;
        TS_ASSERT_EQUALS(boundary_force.GetIdentifier(), "VertexCryptBoundaryForce-2");

        out_stream boundary_force_parameter_file = output_file_handler.OpenOutputFile("boundary_results.parameters");
        boundary_force.OutputForceParameters(boundary_force_parameter_file);
        boundary_force_parameter_file->close();

        std::string boundary_force_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + boundary_force_results_dir + "boundary_results.parameters crypt/test/data/TestForcesForCrypt/boundary_results.parameters").c_str()), 0);

    }

    void TestVertexCryptBoundaryForceMethods() throw (Exception)
    {
        // Create a simple 2D VertexMesh
        HoneycombVertexMeshGenerator generator(5, 5, false, 0.1, 0.5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Translate mesh so that some points are below y=0
        p_mesh->Translate(0.0, -3.0);

        // Set up cells, one for each VertexElement. Give each cell
        // a birth time of -elem_index, so its age is elem_index
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);

            CellPtr p_cell(new Cell(p_state, p_model));
            double birth_time = 0.0 - elem_index;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force system
        VertexCryptBoundaryForce<2> force(100);

        // Initialise a vector of new node forces
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(cell_population.GetNumNodes());

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        force.AddForceContribution(node_forces, cell_population);

        // Check forces are correct
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(node_forces[i][0], 0.0, 1e-4);

            double y = cell_population.GetNode(i)->rGetLocation()[1];
            if (y >= 0.0)
            {
                // If y > 0, the force contribution should be zero...
                TS_ASSERT_DELTA(node_forces[i][1], 0.0, 1e-4);
            }
            else
            {
                // ...otherwise, the force contribution should be quadratic in y
                double expected_force = force.GetForceStrength()*y*y;
                TS_ASSERT_DELTA(node_forces[i][1], expected_force, 1e-4);
            }
        }
    }

    void TestVertexCryptBoundaryForceForceWithNonVertexCellPopulation() throw (Exception)
    {
        // Create a NodeBasedCellPopulation
        std::vector<Node<2>*> nodes;
        unsigned num_nodes = 10;
        for (unsigned i=0; i<num_nodes; i++)
        {
            double x = (double)(i);
            double y = (double)(i);
            nodes.push_back(new Node<2>(i, true, x, y));
        }
        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, num_nodes);

        NodeBasedCellPopulation<2> non_vertex_cell_population(mesh, cells);

        // Create a vector forces on nodes
        std::vector<c_vector<double, 2> > node_forces;
        node_forces.reserve(num_nodes);
        for (unsigned i=0; i<num_nodes; i++)
        {
            node_forces.push_back(zero_vector<double>(2));
        }

        // Test that VertexCryptBoundaryForce throws the correct exception
        VertexCryptBoundaryForce<2> force(100);
        TS_ASSERT_THROWS_THIS(force.AddForceContribution(node_forces, non_vertex_cell_population),
                "VertexCryptBoundaryForce is to be used with VertexBasedCellPopulations only");

        // When the node-only mesh goes out of scope, then it's a different set of nodes that get destroyed
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /*TESTVERTEXCRYPTBOUNDARYFORCE_HPP_*/
