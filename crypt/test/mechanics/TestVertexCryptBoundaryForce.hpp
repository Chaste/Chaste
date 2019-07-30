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

#ifndef TESTVERTEXCRYPTBOUNDARYFORCE_HPP_
#define TESTVERTEXCRYPTBOUNDARYFORCE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "FixedG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexCryptBoundaryForce.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "PetscSetupAndFinalize.hpp"

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

    void TestVertexCryptBoundaryForceMethods()
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
        boost::shared_ptr<AbstractCellProliferativeType> p_diff_type(new DifferentiatedCellProliferativeType);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0 - elem_index;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create a force system
        VertexCryptBoundaryForce<2> force(100);

        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            cell_population.GetNode(i)->ClearAppliedForce();
        }

        force.AddForceContribution(cell_population);

        // Check forces are correct
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[0], 0.0, 1e-4);

            double y = cell_population.GetNode(i)->rGetLocation()[1];
            if (y >= 0.0)
            {
                // If y > 0, the force contribution should be zero...
                TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], 0.0, 1e-4);
            }
            else
            {
                // ...otherwise, the force contribution should be quadratic in y
                double expected_force = force.GetForceStrength()*y*y;
                TS_ASSERT_DELTA(cell_population.GetNode(i)->rGetAppliedForce()[1], expected_force, 1e-4);
            }
        }
    }

    void TestVertexCryptBoundaryForceForceWithNonVertexCellPopulation()
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
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<2> non_vertex_cell_population(mesh, cells);

        for (AbstractMesh<2,2>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
                node_iter != mesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        // Test that VertexCryptBoundaryForce throws the correct exception
        VertexCryptBoundaryForce<2> force(100);
        TS_ASSERT_THROWS_THIS(force.AddForceContribution(non_vertex_cell_population),
                "VertexCryptBoundaryForce is to be used with VertexBasedCellPopulations only");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif /*TESTVERTEXCRYPTBOUNDARYFORCE_HPP_*/
