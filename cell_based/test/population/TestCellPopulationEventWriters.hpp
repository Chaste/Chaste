/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef TESTCELLPOPULATIONEVENTWRITERS_HPP_
#define TESTCELLPOPULATIONEVENTWRITERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"

// Cell population eventwriters
#include "CellDivisionLocationsWriter.hpp"
#include "CellRemovalLocationsWriter.hpp"

// Files to create populations
#include "HoneycombVertexMeshGenerator.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "PetscSetupAndFinalize.hpp"

// Note that high level tests of all cell writers can be found in the
// TestMeshBasedCellPopulation::TestMeshBasedCellPopulationWriteResultsToFile

class TestCellPopulationWriters : public AbstractCellBasedTestSuite
{
public:

    void TestCellDivisionLocationsWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2, 2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestCellDivisionLocationsWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellDivisionLocationsWriter and test that the correct output is generated
        CellDivisionLocationsWriter<2, 2> division_writer;
        division_writer.OpenOutputFile(output_file_handler);
        division_writer.WriteTimeStamp();
        division_writer.Visit(&cell_population);
        division_writer.WriteNewline();
        division_writer.CloseFile();

        FileComparison(results_dir + "divisions.dat", "cell_based/test/data/TestCellPopulationEventWriters/divisions.dat").CompareFiles();

        // Test that we can append to files
        division_writer.OpenOutputFileForAppend(output_file_handler);
        division_writer.WriteTimeStamp();
        division_writer.Visit(&cell_population);
        division_writer.WriteNewline();
        division_writer.CloseFile();

        FileComparison(results_dir + "divisions.dat", "cell_based/test/data/TestCellPopulationEventWriters/divisions_twice.dat").CompareFiles();

        {
            // Coverage of the Visit() method when called on a MeshBasedCellPopulation
            HoneycombMeshGenerator tet_generator(5, 5, 0);
            MutableMesh<2, 2>* p_tet_mesh = tet_generator.GetMesh();
            std::vector<CellPtr> mesh_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> mesh_based_cells_generator;
            mesh_based_cells_generator.GenerateBasic(mesh_based_cells, p_tet_mesh->GetNumNodes());
            MeshBasedCellPopulation<2> mesh_based_cell_population(*p_tet_mesh, mesh_based_cells);

            TS_ASSERT_THROWS_NOTHING(division_writer.Visit(&mesh_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a CaBasedCellPopulation
            PottsMeshGenerator<2> ca_based_generator(5, 0, 0, 5, 0, 0);
            PottsMesh<2>* p_ca_based_mesh = ca_based_generator.GetMesh();
            std::vector<CellPtr> ca_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> ca_based_cells_generator;
            ca_based_cells_generator.GenerateBasic(ca_based_cells, 5);
            std::vector<unsigned> location_indices;
            location_indices.push_back(7);
            location_indices.push_back(11);
            location_indices.push_back(12);
            location_indices.push_back(13);
            location_indices.push_back(17);
            CaBasedCellPopulation<2> ca_based_cell_population(*p_ca_based_mesh, ca_based_cells, location_indices);

            TS_ASSERT_THROWS_NOTHING(division_writer.Visit(&ca_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a NodeBasedCellPopulation
            std::vector<Node<2>*> node_based_nodes;
            node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
            NodesOnlyMesh<2> node_based_mesh;
            node_based_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
            std::vector<CellPtr> node_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
            node_based_generator.GenerateBasic(node_based_cells, node_based_mesh.GetNumNodes());
            NodeBasedCellPopulation<2> node_based_cell_population(node_based_mesh, node_based_cells);

            TS_ASSERT_THROWS_NOTHING(division_writer.Visit(&node_based_cell_population));

            // Tidy up
            delete node_based_nodes[0];
            delete node_based_nodes[1];
        }

        {
            // Coverage of the Visit() method when called on a PottsBasedCellPopulation
            PottsMeshGenerator<2> potts_based_generator(4, 1, 2, 4, 1, 2);
            PottsMesh<2>* p_potts_based_mesh = potts_based_generator.GetMesh();
            std::vector<CellPtr> potts_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> potts_based_cells_generator;
            potts_based_cells_generator.GenerateBasic(potts_based_cells, p_potts_based_mesh->GetNumElements());
            PottsBasedCellPopulation<2> potts_based_cell_population(*p_potts_based_mesh, potts_based_cells);

            TS_ASSERT_THROWS_NOTHING(division_writer.Visit(&potts_based_cell_population));
        }
    }

    void TestCellDivisionLocationsWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellDivisionLocationsWriter.arch";

        {
            AbstractCellBasedWriter<2, 2>* const p_cell_writer = new CellDivisionLocationsWriter<2, 2>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_cell_writer;
            delete p_cell_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2, 2>* p_cell_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell_writer_2;
            delete p_cell_writer_2;
        }
    }

    void TestCellRemovalLocationsWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2, 2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestCellRemovalLocationsWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellRemovalLocationsWriter and test that the correct output is generated
        CellRemovalLocationsWriter<2, 2> removal_writer;
        removal_writer.OpenOutputFile(output_file_handler);
        removal_writer.WriteTimeStamp();
        removal_writer.Visit(&cell_population);
        removal_writer.WriteNewline();
        removal_writer.CloseFile();

        FileComparison(results_dir + "removals.dat", "cell_based/test/data/TestCellPopulationEventWriters/removals.dat").CompareFiles();

        // Test that we can append to files
        removal_writer.OpenOutputFileForAppend(output_file_handler);
        removal_writer.WriteTimeStamp();
        removal_writer.Visit(&cell_population);
        removal_writer.WriteNewline();
        removal_writer.CloseFile();

        FileComparison(results_dir + "removals.dat", "cell_based/test/data/TestCellPopulationEventWriters/removals_twice.dat").CompareFiles();

        {
            // Coverage of the Visit() method when called on a MeshBasedCellPopulation
            HoneycombMeshGenerator tet_generator(5, 5, 0);
            MutableMesh<2, 2>* p_tet_mesh = tet_generator.GetMesh();
            std::vector<CellPtr> mesh_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> mesh_based_cells_generator;
            mesh_based_cells_generator.GenerateBasic(mesh_based_cells, p_tet_mesh->GetNumNodes());
            MeshBasedCellPopulation<2> mesh_based_cell_population(*p_tet_mesh, mesh_based_cells);

            TS_ASSERT_THROWS_NOTHING(removal_writer.Visit(&mesh_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a CaBasedCellPopulation
            PottsMeshGenerator<2> ca_based_generator(5, 0, 0, 5, 0, 0);
            PottsMesh<2>* p_ca_based_mesh = ca_based_generator.GetMesh();
            std::vector<CellPtr> ca_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> ca_based_cells_generator;
            ca_based_cells_generator.GenerateBasic(ca_based_cells, 5);
            std::vector<unsigned> location_indices;
            location_indices.push_back(7);
            location_indices.push_back(11);
            location_indices.push_back(12);
            location_indices.push_back(13);
            location_indices.push_back(17);
            CaBasedCellPopulation<2> ca_based_cell_population(*p_ca_based_mesh, ca_based_cells, location_indices);

            TS_ASSERT_THROWS_NOTHING(removal_writer.Visit(&ca_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a NodeBasedCellPopulation
            std::vector<Node<2>*> node_based_nodes;
            node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
            NodesOnlyMesh<2> node_based_mesh;
            node_based_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
            std::vector<CellPtr> node_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
            node_based_generator.GenerateBasic(node_based_cells, node_based_mesh.GetNumNodes());
            NodeBasedCellPopulation<2> node_based_cell_population(node_based_mesh, node_based_cells);

            TS_ASSERT_THROWS_NOTHING(removal_writer.Visit(&node_based_cell_population));

            // Tidy up
            delete node_based_nodes[0];
            delete node_based_nodes[1];
        }

        {
            // Coverage of the Visit() method when called on a PottsBasedCellPopulation
            PottsMeshGenerator<2> potts_based_generator(4, 1, 2, 4, 1, 2);
            PottsMesh<2>* p_potts_based_mesh = potts_based_generator.GetMesh();
            std::vector<CellPtr> potts_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> potts_based_cells_generator;
            potts_based_cells_generator.GenerateBasic(potts_based_cells, p_potts_based_mesh->GetNumElements());
            PottsBasedCellPopulation<2> potts_based_cell_population(*p_potts_based_mesh, potts_based_cells);

            TS_ASSERT_THROWS_NOTHING(removal_writer.Visit(&potts_based_cell_population));
        }
    }

    void TestCellRemovalLocationsWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellRemovalLocationsWriter.arch";

        {
            AbstractCellBasedWriter<2, 2>* const p_cell_writer = new CellRemovalLocationsWriter<2, 2>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_cell_writer;
            delete p_cell_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2, 2>* p_cell_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell_writer_2;
            delete p_cell_writer_2;
        }
    }
};

#endif /*TESTCELLPOPULATIONEVENTWRITERS_HPP_*/
