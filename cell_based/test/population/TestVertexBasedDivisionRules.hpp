/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTVERTEXBASEDDIVISIONRULE_HPP_
#define TESTVERTEXBASEDDIVISIONRULE_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "CellsGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"
#include "DiagonalVertexBasedDivisionRule.hpp"
#include "RandomDirectionVertexBasedDivisionRule.hpp"
#include "HoneycombVertexMeshGenerator.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestVertexBasedDivisionRules : public AbstractCellBasedTestSuite
{
public:

    void TestAddCellwithDiagonalVertexBasedDivisionRule()
    {
        /**
         * In this test we basically test that the AbstractVertexBasedDivisionRule is implemented and joined with the population
         * correctly. We make a new DiagonalVertexBasedDivisionRule, divide a cell with it and check that the new vertices
         * are in the correct position.
         */
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

        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule_to_set(new DiagonalVertexBasedDivisionRule<2>());
        cell_population.SetVertexBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population and add new cell by dividing element 0 along diagonal axis
        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule = cell_population.GetVertexBasedDivisionRule();
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

    void TestRandomDirectionVertexBasedDivisionRule()
    {
        /**
         * This tests the RandomDirectionVertexBasedDivisionRule. We first create a vertex based cell population and check whether we can
         * give the division rule to the population and get it back. Then we create 10000 division vectors and check that they point
         * uniformly in random directions.
         */

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

        // Create a cell population
        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);

        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule_to_set(new RandomDirectionVertexBasedDivisionRule<2>());
        cell_population.SetVertexBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule = cell_population.GetVertexBasedDivisionRule();

        // Get 10000 division vectors, check each length, their mean and their variance.
        c_vector<double, 2> average_axis = zero_vector<double>(2);
        c_vector<double, 2> axis_variance = zero_vector<double>(2);
        double average_angle = 0.0;
        double angle_variance = 0.0;
        for (unsigned iteration = 0; iteration < 10000; iteration++)
        {
            c_vector<double, 2> random_axis = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);
            TS_ASSERT_DELTA(norm_2(random_axis), 1.0,1e-6);
            average_axis(0) += random_axis(0);
            axis_variance(0) += random_axis(0)*random_axis(0);
            average_axis(1) += random_axis(1);
            axis_variance(1) += random_axis(1)*random_axis(1);
            average_angle += asin(random_axis(0));
            angle_variance += asin(random_axis(0))*asin(random_axis(0));
        }
        average_axis(0) /= 10000.0;
        average_axis(1) /= 10000.0;
        axis_variance(0) /= 10000.0;
        axis_variance(1) /= 10000.0;
        average_angle /= 10000.0;
        angle_variance /= 10000.0;

        TS_ASSERT_DELTA(average_axis(0), 0.0, 1e-2);
        TS_ASSERT_DELTA(average_axis(1), 0.0, 1e-2);
        TS_ASSERT_DELTA(axis_variance(0), 0.5, 1e-2);
        TS_ASSERT_DELTA(axis_variance(1), 0.5, 1e-2);
        TS_ASSERT_DELTA(average_angle, 0.0, 1e-2);
        TS_ASSERT_DELTA(angle_variance, M_PI*M_PI/12.0, 1e-2);
    }

    void TestArchiveRandomDirectionVertexBasedDivisionRule() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "division_rules.arch";

        // Create data structures to store variables to test for equality here
        {
            boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule(new RandomDirectionVertexBasedDivisionRule<2>());

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Record values to test into data structures
            // If necessary you can use static_cast<ConcreteClass*>(p_abstract_class)
            // (if your abstract class doesn't contain the necessary variables and methods)

            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_division_rule;

            // Check things in the data structures with TS_ASSERTS here.
            // If necessary you can use static_cast<ConcreteClass*>(p_abstract_class_2)
            // (if your abstract class doesn't contain the necessary variables and methods)
            // Check that we have got back the right kind of division rule.
            TS_ASSERT(dynamic_cast <RandomDirectionVertexBasedDivisionRule<2>* > (p_division_rule.get()));
        }
    }

    void TestArchiveDiagonalVertexBasedDivisionRule() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "division_rules.arch";

        // Create data structures to store variables to test for equality here
        {
            boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule(new DiagonalVertexBasedDivisionRule<2>());

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Record values to test into data structures
            // If necessary you can use static_cast<ConcreteClass*>(p_abstract_class)
            // (if your abstract class doesn't contain the necessary variables and methods)

            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // restore from the archive
            (*p_arch) >> p_division_rule;

            // Check things in the data structures with TS_ASSERTS here.
            // If necessary you can use static_cast<ConcreteClass*>(p_abstract_class_2)
            // (if your abstract class doesn't contain the necessary variables and methods)
            // Check that we have got back the right kind of division rule.
            TS_ASSERT(dynamic_cast <DiagonalVertexBasedDivisionRule<2>* > (p_division_rule.get()));
        }
    }

};

#endif /*TESTVERTEXBASEDDIVISIONRULE_HPP_*/
