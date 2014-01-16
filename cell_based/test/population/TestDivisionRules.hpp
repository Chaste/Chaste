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

#ifndef TESTDIVISIONRULE_HPP_
#define TESTDIVISIONRULE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "CellsGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellDivisionRule.hpp"
#include "DiagonalDivisionRule.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestDivisionRules : public AbstractCellBasedTestSuite
{
public:

    void TestAddCellwithDiagonalDivisionRule()
    {
        // Make some nodes
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.0, 2.0));
        nodes.push_back(new Node<2>(2, true, -2.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 0.0, -2.0));

        // Make a rectangular element out of nodes 0,1,2,3
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);

        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));

        // Make a vertex mesh
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        TS_ASSERT_EQUALS(vertex_mesh.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(vertex_mesh.GetNumNodes(), 4u);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        unsigned old_num_nodes = vertex_mesh.GetNumNodes();

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_temp_cell(new Cell(p_state, p_model));
        p_temp_cell->SetCellProliferativeType(p_stem_type);
        p_temp_cell->SetBirthTime(-1);

        // Set the division rule for our population to be the diagonal division rule

        boost::shared_ptr<AbstractCellDivisionRule<2> > p_division_rule_to_set(new DiagonalDivisionRule<2>());
        cell_population.SetDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population and add new cell by dividing element 0 along diagonal axis
        boost::shared_ptr<AbstractCellDivisionRule<2> > p_division_rule = cell_population.GetDivisionRule();
        c_vector<double, 2> diagonal_axis = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);

        // Check that the axis is pointing in direction (1,1)
        TS_ASSERT_DELTA(diagonal_axis[0], 1.0, 1e-9);
        TS_ASSERT_DELTA(diagonal_axis[1], 1.0, 1e-9);

        cell_population.AddCell(p_temp_cell, diagonal_axis, p_cell0);

        // Check the location of the new nodes
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes)->rGetLocation()[0], 1.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes)->rGetLocation()[1], 1.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes+1)->rGetLocation()[0], -1.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetNode(old_num_nodes+1)->rGetLocation()[1], -1.0, 1e-12);
    }

    void TestArchiveDiagonalDivisionRule() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "division_rules.arch";

        // Create data structures to store variables to test for equality here
        {
            boost::shared_ptr<AbstractCellDivisionRule<2> > p_division_rule(new DiagonalDivisionRule<2>());

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Record values to test into data structures
            // If necessary you can use static_cast<ConcreteClass*>(p_abstract_class)
            // (if your abstract class doesn't contain the necessary variables and methods)

            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractCellDivisionRule<2> > p_division_rule;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_division_rule;

            // Check things in the data structures with TS_ASSERTS here.
            // If necessary you can use static_cast<ConcreteClass*>(p_abstract_class_2)
            // (if your abstract class doesn't contain the necessary variables and methods)
            // Check that we have got back the right kind of division rule.
            TS_ASSERT(dynamic_cast <DiagonalDivisionRule<2>* > (p_division_rule.get()));
        }
    }

};

#endif /*TESTDIVISIONRULE_HPP_*/
