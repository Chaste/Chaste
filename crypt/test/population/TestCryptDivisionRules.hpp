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

#ifndef TESTCRYPTDIVISIONRULES_HPP_
#define TESTCRYPTDIVISIONRULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "CryptCellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "RandomDirectionCentreBasedDivisionRule.hpp"
#include "CryptCentreBasedDivisionRule.hpp"
#include "CryptVertexBasedDivisionRule.hpp"
#include "CryptShovingCaBasedDivisionRule.hpp"
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
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
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
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
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
        CellsGenerator<FixedG1GenerationalCellCycleModel,3> cells_generator;
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
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
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

    void TestAddCellWithCryptShovingBasedDivisionRule()
    {
        /**
         * In this test we create a new CryptShovingCaBasedDivisionRule, divide a cell with it
         * and check that the new cells are in the correct locations. First, we test where
         * there is space around the cells. This is the default setup.
         */

        // Create a simple Potts mesh
        PottsMeshGenerator<2> generator(3, 0, 0, 4, 0, 0,1,0,0,false, true); // Periodic in x
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create 6 cells in the bottom 2 rows
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<6; index++)
        {
            location_indices.push_back(index);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size()); // Note all cells are stem cells by default.

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Check the cell locations
        unsigned cell_locations[6] = {0,1,2,3,4,5};
        unsigned index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter),cell_locations[index])
            ++index;
        }

        // Make a new cell to add
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetCellProliferativeType(p_transit_type);
        p_new_cell->SetBirthTime(-1);

        // Set the division rule for our population to be the cryot shoving division rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule_to_set(new CryptShovingCaBasedDivisionRule());
        cell_population.SetCaBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population and try to add new cell by dividing cell at site 0;
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule = cell_population.GetCaBasedDivisionRule();

        // Select left cell in bottom row
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);

        // The CryptShovingCaBasedDivisionRule method IsRoomToDivide() always returns true
        TS_ASSERT((p_division_rule->IsRoomToDivide(p_cell_0,cell_population)));

        /*
         * Test adding the new cell in the population; this calls CalculateDaughterNodeIndex().
         * The new cell moves into node 3. This is because stem cells always divide upwards
         */
        cell_population.AddCell(p_new_cell, p_cell_0);

        // Now check the cells are in the correct place
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 7u);

        // Note the cell on node 3 has been shoved to node 6 and the new cell is on node 3
        unsigned new_cell_locations[7] = {0,1,2,6,4,5,3};
        index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter),new_cell_locations[index])
            ++index;
        }

        //Now divide a non stem cell
        // Select transit cell we just added
        CellPtr p_cell_3 = cell_population.GetCellUsingLocationIndex(3);

        FixedG1GenerationalCellCycleModel* p_model_2 = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell_2(new Cell(p_state, p_model_2));
        p_new_cell_2->SetCellProliferativeType(p_transit_type);
        p_new_cell_2->SetBirthTime(-1);

        /*
         * Test adding the new cell in the population; this calls CalculateDaughterNodeIndex().
         * The new cell moves into node 7.
         */
        cell_population.AddCell(p_new_cell_2, p_cell_3);

        // Now check the cells are in the correct place
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 8u);

        // Note the cell on node 4 has been shoved to node 7 and the new cell is on node 4
        unsigned new_cell_locations_2[8] = {0,1,2,6,7,5,3,4};
        index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter),new_cell_locations_2[index])
            ++index;
        }
    }

    void TestAddCellWithCryptShovingBasedDivisionRuleAndShovingRequired()
    {

        /**
         * In this test of CryptShovingCaBasedDivisionRule we check the case where there is
         * no room to divide without the cells being shoved to the edge of the mesh.
         */

        // Create a simple Potts mesh
        PottsMeshGenerator<2> generator(3, 0, 0, 3, 0, 0, 1, 0, 0, false, true); // x periodic
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create 9 cells, one for each node
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<9; index++)
        {
            location_indices.push_back(index);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Make a new cell to add
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetCellProliferativeType(p_stem_type);
        p_new_cell->SetBirthTime(-1);

        // Set the division rule for our population to be the crypt shoving division rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule_to_set(new CryptShovingCaBasedDivisionRule());
        cell_population.SetCaBasedDivisionRule(p_division_rule_to_set);

        // Get the division rule back from the population and try to add new cell by dividing cell at site 0
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule = cell_population.GetCaBasedDivisionRule();

        // Select bottom left cell
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);

        // Can't divide without shoving into top
        TS_ASSERT_THROWS_THIS(p_division_rule->CalculateDaughterNodeIndex(p_new_cell, p_cell_0, cell_population),
            "Cells reaching the top of the crypt need to increase length to at least double the sloughing height.");
    }

    void TestArchiveCryptCentreBasedDivisionRule()
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

    void TestArchiveCryptVertexBasedDivisionRule()
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

    void TestArchivingCryptShovingCaBasedDivisionRule()
    {
        EXIT_IF_PARALLEL; // Beware of processes overwriting the identical archives of other processes
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CryptShovingCaBasedDivisionRule.arch";

        {
            CryptShovingCaBasedDivisionRule division_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer to most abstract class possible
            AbstractCaBasedDivisionRule<2>* const p_division_rule = &division_rule;
            output_arch << p_division_rule;
        }

        {
            AbstractCaBasedDivisionRule<2>* p_division_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_division_rule;

            TS_ASSERT(p_division_rule != NULL);

            // Tidy up
            delete p_division_rule;
        }
    }
};

#endif /*TESTCRYPTDIVISIONRULES_HPP_*/
