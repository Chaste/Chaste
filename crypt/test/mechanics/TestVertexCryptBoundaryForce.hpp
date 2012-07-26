/*

Copyright (c) 2005-2012, University of Oxford.
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
#include "FileComparison.hpp"


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
        FileComparison( boundary_force_results_dir + "boundary_results.parameters", "crypt/test/data/TestForcesForCrypt/boundary_results.parameters").CompareFiles();

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
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(DIFFERENTIATED);
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
