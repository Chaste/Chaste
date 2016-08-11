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

#ifndef TESTCRYPTDIVISIONRULES_HPP_
#define TESTCRYPTDIVISIONRULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CryptCellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "CryptCentreBasedDivisionRule.hpp"
#include "CryptVertexBasedDivisionRule.hpp"
#include "SmartPointers.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestCryptDivisionRules : public AbstractCellBasedTestSuite
{
public:

    void TestCryptCentreBasedDivisionRule1d()
    {
        // Create a 2d cell population
        MutableMesh<1,1> mesh;
        mesh.ConstructLinearMesh(21);

        ChastePoint<1> shifted_point;
        shifted_point.rGetLocation()[0] = 10.5;
        mesh.SetNode(10, shifted_point);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0 - node_index;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<1,1> cell_population(mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);

        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_division_rule_to_set(new CryptCentreBasedDivisionRule<1,1>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractCentreBasedDivisionRule<1,1> > p_division_rule = cell_population.GetCentreBasedDivisionRule();

        std::pair<c_vector<double, 1>, c_vector<double, 1> > positions = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 1> parent_position = positions.first;
        c_vector<double, 1> daughter_position = positions.second;

        TS_ASSERT_DELTA(parent_position[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(daughter_position[0], 0.15, 1e-4);
    }

    void TestCryptCentreBasedDivisionRule2d()
    {
        // Create a 2d cell population
        c_vector<double,2> location;
        location[0] = 1.0;
        location[1] = 1.0;
        Node<2>* p_node = new Node<2>(0u,location, false);

        MutableMesh<2,2> conf_mesh;
        conf_mesh.AddNode(p_node);

        // Create cells
        std::vector<CellPtr> conf_cells;
        CryptCellsGenerator<FixedDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(conf_cells, &conf_mesh, std::vector<unsigned>(), true);

        // Create cell population
        MeshBasedCellPopulation<2,2> cell_population(conf_mesh, conf_cells);
        cell_population.SetMeinekeDivisionSeparation(0.1);

        CellPtr p_cell0 = *(cell_population.Begin());

        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule_to_set(new CryptCentreBasedDivisionRule<2,2>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule = cell_population.GetCentreBasedDivisionRule();

        std::pair<c_vector<double, 2>, c_vector<double, 2> > positions = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);

        c_vector<double, 2> parent_position = positions.first;
        c_vector<double, 2> daughter_position = positions.second;

        TS_ASSERT_DELTA(parent_position[0], 1.0417, 1e-4);
        TS_ASSERT_DELTA(parent_position[1], 1.0275, 1e-4);
        TS_ASSERT_DELTA(daughter_position[0], 0.9582, 1e-4);
        TS_ASSERT_DELTA(daughter_position[1], 0.9724, 1e-4);
    }

    void TestCryptCentreBasedDivisionRule3d()
    {
        // Create a 3d cell population
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/3D_Single_tetrahedron_element");
        MutableMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel,3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<3,3> cell_population(mesh, cells);

        CellPtr p_cell0 = *(cell_population.Begin());

        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule_to_set(new CryptCentreBasedDivisionRule<3,3>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3> > p_division_rule = cell_population.GetCentreBasedDivisionRule();

        TS_ASSERT_THROWS_THIS(p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population),
            "CryptCentreBasedDivisionRule is not implemented for SPACE_DIM == 3");
    }

    void TestCryptVertexBasedDivisionRule()
    {
        // Create a vertex cell population
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, true, 2.0, 0.0));
        nodes.push_back(new Node<2>(1, true, 0.0, 2.0));
        nodes.push_back(new Node<2>(2, true, -2.0, 0.0));
        nodes.push_back(new Node<2>(3, true, 0.0, -2.0));
        std::vector<Node<2>*> nodes_elem_1;
        nodes_elem_1.push_back(nodes[0]);
        nodes_elem_1.push_back(nodes[1]);
        nodes_elem_1.push_back(nodes[2]);
        nodes_elem_1.push_back(nodes[3]);
        std::vector<VertexElement<2,2>*> vertex_elements;
        vertex_elements.push_back(new VertexElement<2,2>(0, nodes_elem_1));
        MutableVertexMesh<2,2> vertex_mesh(nodes, vertex_elements);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, vertex_mesh.GetNumElements());

        VertexBasedCellPopulation<2> cell_population(vertex_mesh, cells);

        CellPtr p_cell0 = cell_population.GetCellUsingLocationIndex(0);

        // Set the division rule for our population to be the random direction division rule
        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule_to_set(new CryptVertexBasedDivisionRule<2>());
        cell_population.SetVertexBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population
        boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule = cell_population.GetVertexBasedDivisionRule();

        c_vector<double, 2> division_axis = p_division_rule->CalculateCellDivisionVector(p_cell0, cell_population);
        TS_ASSERT_DELTA(division_axis[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(division_axis[1], 0.0, 1e-4);
    }

    void TestArchiveCryptCentreBasedDivisionRule() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "CryptCentreBasedDivisionRule.arch";

        {
            boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule(new CryptCentreBasedDivisionRule<2,2>());

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> p_division_rule;

            typedef CryptCentreBasedDivisionRule<2,2> CryptRule;
            TS_ASSERT(dynamic_cast<CryptRule*>(p_division_rule.get()));
        }
    }

    void TestArchiveCryptVertexBasedDivisionRule() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "CryptVertexBasedDivisionRule.arch";

        {
            boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule(new CryptVertexBasedDivisionRule<2>());

            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) << p_division_rule;
        }

        {
            boost::shared_ptr<AbstractVertexBasedDivisionRule<2> > p_division_rule;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> p_division_rule;

            typedef CryptVertexBasedDivisionRule<2> CryptRule;
            TS_ASSERT(dynamic_cast<CryptRule*>(p_division_rule.get()));
        }
    }
};

#endif /*TESTCRYPTDIVISIONRULES_HPP_*/
