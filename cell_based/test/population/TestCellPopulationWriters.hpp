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

#ifndef TESTCELLPOPULATIONWRITERS_HPP_
#define TESTCELLPOPULATIONWRITERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLabel.hpp"

// Cell population writers
#include "BoundaryNodeWriter.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "HeterotypicBoundaryLengthWriter.hpp"
#include "NodeLocationWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "RadialCellDataDistributionWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"
#include "VertexT2SwapLocationsWriter.hpp"
#include "VertexT3SwapLocationsWriter.hpp"
#include "VoronoiDataWriter.hpp"

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

    void TestBoundaryNodeWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 3D MeshBasedCellPopulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestBoundaryNodeWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a BoundaryNodeWriter and test that the correct output is generated
        BoundaryNodeWriter<3,3> boundary_writer;
        boundary_writer.OpenOutputFile(output_file_handler);
        boundary_writer.WriteTimeStamp();
        boundary_writer.Visit(&cell_population);
        boundary_writer.WriteNewline();
        boundary_writer.CloseFile();

        FileComparison(results_dir + "results.vizboundarynodes", "cell_based/test/data/TestCellPopulationWriters/results.vizboundarynodes").CompareFiles();

        // Test that we can append to files
        boundary_writer.OpenOutputFileForAppend(output_file_handler);
        boundary_writer.WriteTimeStamp();
        boundary_writer.Visit(&cell_population);
        boundary_writer.WriteNewline();
        boundary_writer.CloseFile();

        FileComparison(results_dir + "results.vizboundarynodes", "cell_based/test/data/TestCellPopulationWriters/results.vizboundarynodes_twice").CompareFiles();
    }

    void TestBoundaryNodeWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "BoundaryNodeWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new BoundaryNodeWriter<2,2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_population_writer;
            delete p_population_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_population_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_population_writer_2;
            delete p_population_writer_2;
       }
    }

    void TestCellPopulationAdjacencyMatrixWriter()
    {
        EXIT_IF_PARALLEL;

        // Test with a MeshBasedCellPopulation
        HoneycombMeshGenerator tet_generator(5, 5, 0);
        MutableMesh<2,2>* p_tet_mesh = tet_generator.GetMesh();
        std::vector<CellPtr> mesh_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> mesh_based_cells_generator;
        mesh_based_cells_generator.GenerateBasic(mesh_based_cells, p_tet_mesh->GetNumNodes());
        MeshBasedCellPopulation<2> mesh_based_cell_population(*p_tet_mesh, mesh_based_cells);

        // Label a subset of the cells
        boost::shared_ptr<AbstractCellProperty> p_label(mesh_based_cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
        mesh_based_cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        mesh_based_cell_population.GetCellUsingLocationIndex(1)->AddCellProperty(p_label);

        // Create an output directory for the writer
        std::string output_directory = "TestCellPopulationAdjacencyMatrixWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellPopulationAreaWriter and test that the correct output is generated
        CellPopulationAdjacencyMatrixWriter<2,2> adjacency_writer;
        adjacency_writer.OpenOutputFile(output_file_handler);
        adjacency_writer.WriteTimeStamp();
        adjacency_writer.Visit(&mesh_based_cell_population);
        adjacency_writer.WriteNewline();
        adjacency_writer.CloseFile();

        FileComparison(results_dir + "cellpopulationadjacency.dat", "cell_based/test/data/TestCellPopulationWriters/cellpopulationadjacency.dat").CompareFiles();

        // Test that we can append to files
        adjacency_writer.OpenOutputFileForAppend(output_file_handler);
        adjacency_writer.WriteTimeStamp();
        adjacency_writer.Visit(&mesh_based_cell_population);
        adjacency_writer.WriteNewline();
        adjacency_writer.CloseFile();

        FileComparison(results_dir + "cellpopulationadjacency.dat", "cell_based/test/data/TestCellPopulationWriters/cellpopulationadjacency_twice.dat").CompareFiles();

        // Coverage of the Visit() method when called on a CaBasedCellPopulation
        {
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

            TS_ASSERT_THROWS_NOTHING(adjacency_writer.Visit(&ca_based_cell_population));
        }

        // Coverage of the Visit() method when called on a NodeBasedCellPopulation
        {
            std::vector<Node<2>* > node_based_nodes;
            node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
            NodesOnlyMesh<2> node_based_mesh;
            node_based_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
            std::vector<CellPtr> node_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
            node_based_generator.GenerateBasic(node_based_cells, node_based_mesh.GetNumNodes());
            NodeBasedCellPopulation<2> node_based_cell_population(node_based_mesh, node_based_cells);

            TS_ASSERT_THROWS_NOTHING(adjacency_writer.Visit(&node_based_cell_population));

            // Tidy up
            delete node_based_nodes[0];
            delete node_based_nodes[1];
        }

        // Coverage of the Visit() method when called on a PottsBasedCellPopulation
        {
            PottsMeshGenerator<2> potts_based_generator(4, 1, 2, 4, 1, 2);
            PottsMesh<2>* p_potts_based_mesh = potts_based_generator.GetMesh();
            std::vector<CellPtr> potts_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> potts_based_cells_generator;
            potts_based_cells_generator.GenerateBasic(potts_based_cells, p_potts_based_mesh->GetNumElements());
            PottsBasedCellPopulation<2> potts_based_cell_population(*p_potts_based_mesh, potts_based_cells);

            TS_ASSERT_THROWS_NOTHING(adjacency_writer.Visit(&potts_based_cell_population));
        }

        // Coverage of the Visit() method when called on a VertexBasedCellPopulation
        {
            HoneycombVertexMeshGenerator vertex_based_generator(4, 6);
            MutableVertexMesh<2,2>* p_vertex_based_mesh = vertex_based_generator.GetMesh();
            std::vector<CellPtr> vertex_based_cells;
            boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> vertex_based_cells_generator;
            vertex_based_cells_generator.GenerateBasic(vertex_based_cells, p_vertex_based_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
            VertexBasedCellPopulation<2> vertex_based_cell_population(*p_vertex_based_mesh, vertex_based_cells);

            boost::shared_ptr<AbstractCellProperty> p_new_label(vertex_based_cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
            vertex_based_cell_population.GetCellUsingLocationIndex(3)->AddCellProperty(p_new_label);
            vertex_based_cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_new_label);
            vertex_based_cell_population.GetCellUsingLocationIndex(5)->AddCellProperty(p_new_label);
            vertex_based_cell_population.GetCellUsingLocationIndex(6)->AddCellProperty(p_new_label);

            TS_ASSERT_THROWS_NOTHING(adjacency_writer.Visit(&vertex_based_cell_population));
        }
    }

    void TestCellPopulationAdjacencyMatrixWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellPopulationAdjacencyMatrixWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new CellPopulationAdjacencyMatrixWriter<2,2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_population_writer;
            delete p_population_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_population_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_population_writer_2;
            delete p_population_writer_2;
       }
    }

    void TestCellPopulationAreaWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 3D MeshBasedCellPopulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        cell_population.CreateVoronoiTessellation();

        // Create an output directory for the writer
        std::string output_directory = "TestCellPopulationAreaWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellPopulationAreaWriter and test that the correct output is generated
        CellPopulationAreaWriter<3,3> area_writer;
        area_writer.OpenOutputFile(output_file_handler);
        area_writer.WriteTimeStamp();
        area_writer.Visit(&cell_population);
        area_writer.WriteNewline();
        area_writer.CloseFile();

        FileComparison(results_dir + "cellpopulationareas.dat", "cell_based/test/data/TestCellPopulationWriters/cellpopulationareas.dat").CompareFiles();

        // Test that we can append to files
        area_writer.OpenOutputFileForAppend(output_file_handler);
        area_writer.WriteTimeStamp();
        area_writer.Visit(&cell_population);
        area_writer.WriteNewline();
        area_writer.CloseFile();

        FileComparison(results_dir + "cellpopulationareas.dat", "cell_based/test/data/TestCellPopulationWriters/cellpopulationareas_twice.dat").CompareFiles();

        CellPopulationAreaWriter<2,2> area_writer_2d;

        // Test the correct exception is thrown if using a NodeBasedCellPopulation
        std::vector<Node<2>* > node_based_nodes;
        node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
        NodesOnlyMesh<2> nodes_only_mesh;
        nodes_only_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
        std::vector<CellPtr> node_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
        node_based_generator.GenerateBasic(node_based_cells, nodes_only_mesh.GetNumNodes());
        NodeBasedCellPopulation<2> node_based_cell_population(nodes_only_mesh, node_based_cells);

        TS_ASSERT_THROWS_THIS(area_writer_2d.Visit(&node_based_cell_population),
            "CellPopulationAreaWriter cannot be used with a NodeBasedCellPopulation");

        // Tidy up
        for (unsigned i=0; i<node_based_nodes.size(); i++)
        {
            delete node_based_nodes[i];
        }

        // Test the correct exception is thrown if using a CaBasedCellPopulation
        PottsMeshGenerator<2> ca_based_generator(4, 0, 0, 4, 0, 0);
        PottsMesh<2>* p_ca_based_mesh = ca_based_generator.GetMesh();
        std::vector<CellPtr> ca_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> ca_based_cells_generator;
        ca_based_cells_generator.GenerateBasic(ca_based_cells, 4);
        std::vector<unsigned> location_indices;
        location_indices.push_back(7);
        location_indices.push_back(11);
        location_indices.push_back(12);
        location_indices.push_back(13);
        CaBasedCellPopulation<2> ca_based_cell_population(*p_ca_based_mesh, ca_based_cells, location_indices);

        TS_ASSERT_THROWS_THIS(area_writer_2d.Visit(&ca_based_cell_population),
            "CellPopulationAreaWriter cannot be used with a CaBasedCellPopulation");

        // Test the correct exception is thrown if using a PottsBasedCellPopulation
        PottsMeshGenerator<2> potts_based_generator(4, 1, 2, 4, 1, 2);
        PottsMesh<2>* p_potts_based_mesh = potts_based_generator.GetMesh();
        std::vector<CellPtr> potts_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> potts_based_cells_generator;
        potts_based_cells_generator.GenerateBasic(potts_based_cells, p_potts_based_mesh->GetNumElements());
        PottsBasedCellPopulation<2> potts_based_cell_population(*p_potts_based_mesh, potts_based_cells);

        TS_ASSERT_THROWS_THIS(area_writer_2d.Visit(&potts_based_cell_population),
            "CellPopulationAreaWriter cannot be used with a PottsBasedCellPopulation");

        // Test the correct exception is thrown if using a PottsBasedCellPopulation
        HoneycombVertexMeshGenerator vertex_based_generator(4, 6);
        MutableVertexMesh<2,2>* p_vertex_mesh = vertex_based_generator.GetMesh();
        std::vector<CellPtr> vertex_based_cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> vertex_based_cells_generator;
        vertex_based_cells_generator.GenerateBasic(vertex_based_cells, p_vertex_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> vertex_cell_population(*p_vertex_mesh, vertex_based_cells);

        TS_ASSERT_THROWS_THIS(area_writer_2d.Visit(&vertex_cell_population),
            "CellPopulationAreaWriter cannot be used with a VertexBasedCellPopulation");
    }

    void TestCellPopulationAreaWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellPopulationAreaWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new CellPopulationAreaWriter<2,2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_population_writer;
            delete p_population_writer;
        }
        PetscTools::Barrier(); // Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_population_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_population_writer_2;
            delete p_population_writer_2;
       }
    }

    void TestCellPopulationElementWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestCellPopulationElementWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellPopulationElementWriter and test that the correct output is generated
        CellPopulationElementWriter<2,2> element_writer;
        element_writer.OpenOutputFile(output_file_handler);
        element_writer.WriteTimeStamp();
        element_writer.Visit(&cell_population);
        element_writer.WriteNewline();
        element_writer.CloseFile();

        FileComparison(results_dir + "results.vizelements", "cell_based/test/data/TestCellPopulationWriters/results.vizelements").CompareFiles();

        // Test that we can append to files
        element_writer.OpenOutputFileForAppend(output_file_handler);
        element_writer.WriteTimeStamp();
        element_writer.Visit(&cell_population);
        element_writer.WriteNewline();
        element_writer.CloseFile();

        FileComparison(results_dir + "results.vizelements", "cell_based/test/data/TestCellPopulationWriters/results.vizelements_twice").CompareFiles();

        // Test the correct exception is thrown if using a NodeBasedCellPopulation
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
        NodesOnlyMesh<2> nodes_only_mesh;
        nodes_only_mesh.ConstructNodesWithoutMesh(nodes, 1.5);
        std::vector<CellPtr> node_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
        node_based_generator.GenerateBasic(node_based_cells, nodes_only_mesh.GetNumNodes());
        NodeBasedCellPopulation<2> node_based_cell_population(nodes_only_mesh, node_based_cells);

        TS_ASSERT_THROWS_THIS(element_writer.Visit(&node_based_cell_population),
            "CellPopulationElementWriter cannot be used with a NodeBasedCellPopulation");

        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

        // Test the correct exception is thrown if using a CaBasedCellPopulation
        PottsMeshGenerator<2> ca_based_generator(4, 0, 0, 4, 0, 0);
        PottsMesh<2>* p_ca_based_mesh = ca_based_generator.GetMesh();
        std::vector<CellPtr> ca_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> ca_based_cells_generator;
        ca_based_cells_generator.GenerateBasic(ca_based_cells, 4);
        std::vector<unsigned> location_indices;
        location_indices.push_back(7);
        location_indices.push_back(11);
        location_indices.push_back(12);
        location_indices.push_back(13);
        CaBasedCellPopulation<2> ca_based_cell_population(*p_ca_based_mesh, ca_based_cells, location_indices);

        TS_ASSERT_THROWS_THIS(element_writer.Visit(&ca_based_cell_population),
            "CellPopulationElementWriter cannot be used with a CaBasedCellPopulation");
    }

    void TestCellPopulationElementWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellPopulationElementWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new CellPopulationElementWriter<2,2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_population_writer;
            delete p_population_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_population_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_population_writer_2;
            delete p_population_writer_2;
       }
    }

    void TestHeterotypicBoundaryLengthWriter()
    {
        EXIT_IF_PARALLEL;

        // Test with a MeshBasedCellPopulationWithGhostNodes
        {
            // Create a simple 2D cell population (use ghost nodes to avoid infinite edge lengths in the Voronoi tessellation)
            HoneycombMeshGenerator generator(5, 3, 2);
            MutableMesh<2,2>* p_mesh = generator.GetMesh();
            std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel,2> cells_generator;
            cells_generator.GenerateGivenLocationIndices(cells, location_indices);
            MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
            cell_population.InitialiseCells();

            // Label a subset of the cells
            boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
            cell_population.GetCellUsingLocationIndex(20)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(21)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(22)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(23)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(24)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(32)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(39)->AddCellProperty(p_label);

            /*
             * In this case, the group of labelled cells 0, 1, 2, 3, 4, 8 share 11 edges
             * with unlabelled cells, while the isolated labelled cell 11 shares 4 edges
             * with unlabelled cells.
             *
             * Thus there are 15 edges shared between labelled and unlabelled cells.
             * There are 30 shared edges in total (regardless of label).
             * Each edge has length 1/sqrt(3) = 0.577350.
             *
             * Thus the total length is 30/sqrt(3) and the heterotypic boundary
             * length is 15/sqrt(3) = 8.660254.
             *
             * This can be verified by eyeballing the output file.
             */

            // Create an output directory for the writer
            std::string output_directory = "TestHeterotypicBoundaryLengthWriterMesh";
            OutputFileHandler output_file_handler(output_directory, false);
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

            // Create a BoundaryNodeWriter and test that the correct output is generated
            HeterotypicBoundaryLengthWriter<2,2> labelled_boundary_writer;
            labelled_boundary_writer.OpenOutputFile(output_file_handler);
            labelled_boundary_writer.WriteTimeStamp();
            labelled_boundary_writer.Visit(&cell_population);
            labelled_boundary_writer.WriteNewline();
            labelled_boundary_writer.CloseFile();

            FileComparison(results_dir + "heterotypicboundary.dat", "cell_based/test/data/TestCellPopulationWriters/heterotypicboundary.dat_mesh").CompareFiles();
        }

        // Test with a NodeBasedCellPopulation
        {
            // Create a simple 2D cell population
            std::vector<Node<2>* > nodes;
            for (unsigned j=0; j<4; j++)
            {
                for (unsigned i=0; i<6; i++)
                {
                    unsigned node_index = i + 6*j;
                    nodes.push_back(new Node<2>(node_index, false, (double)i, (double)j));
                }
            }
            NodesOnlyMesh<2> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 1.5);
            for (unsigned index=0; index<mesh.GetNumNodes(); index++)
            {
                mesh.GetNode(index)->SetRadius(0.6);
            }
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());
            NodeBasedCellPopulation<2> cell_population(mesh, cells);
            cell_population.InitialiseCells();

            // Label a subset of the cells
            boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
            cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(6)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(9)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(10)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(11)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(12)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(15)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(16)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(17)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(19)->AddCellProperty(p_label);

            /*
             * In this case, the group of labelled cells 0, 1, 12 share 4 edges with
             * unlabelled cells; the group of labelled cells 9, 10, 11, 15, 16, 17
             * share 8 edges with unlabelled cells; and the isolated labelled cell 19
             * shares 3 edges with unlabelled cells.
             * Thus there are 15 edges shared between labelled and unlabelled cells.
             *
             * In total there are 38 shared edges (regardless of label).
             *
             * Since each cell's radius is set to 0.6 and neighbours are a distance
             * 1.0 apart, the approximate shared edge is given by sqrt(11)/5 = 0.663324.
             *
             * Thus the total length is 38*sqrt(11)/5 = 25.206348 and the heterotypic
             * boundary length is 15*sqrt(11)/5 = 9.949874.
             *
             * This can be verified by eyeballing the output file.
             */

            // Create an output directory for the writer
            std::string output_directory = "TestHeterotypicBoundaryLengthWriterNode";
            OutputFileHandler output_file_handler(output_directory, false);
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

            // Create a BoundaryNodeWriter and test that the correct output is generated
            HeterotypicBoundaryLengthWriter<2,2> labelled_boundary_writer;
            labelled_boundary_writer.OpenOutputFile(output_file_handler);
            labelled_boundary_writer.WriteTimeStamp();
            labelled_boundary_writer.Visit(&cell_population);
            labelled_boundary_writer.WriteNewline();
            labelled_boundary_writer.CloseFile();

            FileComparison(results_dir + "heterotypicboundary.dat", "cell_based/test/data/TestCellPopulationWriters/heterotypicboundary.dat_node").CompareFiles();

            // Avoid memory leak
            for (unsigned i=0; i<nodes.size(); i++)
            {
                delete nodes[i];
            }
        }

        // Test with a PottsBasedCellPopulation
        {
            // Create a simple 2D cell population
            PottsMeshGenerator<2> generator(9, 3, 3, 6, 3, 2);
            PottsMesh<2>* p_mesh = generator.GetMesh();
            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());
            PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.InitialiseCells();

            // Label a subset of the cells
            boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
            cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(1)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(5)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(8)->AddCellProperty(p_label);

            /*
             * In this case, the group of labelled cells 0, 1, 4, 5 share 3 'long' edges
             * and 3 'short' edges with unlabelled cells.
             *
             * Thus there are 6 edges shared between labelled and unlabelled cells.
             *
             * In total there are 6 'long' edges and 6 'short' edges, thus 12 shared
             * edges in total (regardless of label). 'Long' edges have length 3 and
             * 'short' edges have length 2.
             *
             * Thus the total length is 6*3 + 6*2 = 30 and the heterotypic boundary
             * length is 3*3 + 3*2 = 15.
             *
             * Note the number of cell pairs is 12 and the number of heterotypic
             * cell pairs is 6.
             *
             * This can be verified by eyeballing the output file.
             */

            // Create an output directory for the writer
            std::string output_directory = "TestHeterotypicBoundaryLengthWriterPotts";
            OutputFileHandler output_file_handler(output_directory, false);
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

            // Create a BoundaryNodeWriter and test that the correct output is generated
            HeterotypicBoundaryLengthWriter<2,2> labelled_boundary_writer;
            labelled_boundary_writer.OpenOutputFile(output_file_handler);
            labelled_boundary_writer.WriteTimeStamp();
            labelled_boundary_writer.Visit(&cell_population);
            labelled_boundary_writer.WriteNewline();
            labelled_boundary_writer.CloseFile();

            FileComparison(results_dir + "heterotypicboundary.dat", "cell_based/test/data/TestCellPopulationWriters/heterotypicboundary.dat_potts").CompareFiles();
        }

        // Test with a CaBasedCellPopulation
        {
            // Create a simple 2D cell population
            PottsMeshGenerator<2> generator(3, 0, 0, 3, 0, 0);
            PottsMesh<2>* p_mesh = generator.GetMesh();

            std::vector<unsigned> location_indices;
            for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
            {
                location_indices.push_back(i);
            }

            std::vector<CellPtr> cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, location_indices.size());
            CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
            cell_population.InitialiseCells();

            // Label a subset of the cells
            boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
            cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(1)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(5)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(8)->AddCellProperty(p_label);

            /*
             * In this case, the group of labelled cells 0, 1, 4, 5 share 6  edges
             * edges with unlabelled cells.
             *
             * Thus there are 5 edges shared between labelled and unlabelled cells.
             *
             * In total there are 12 shared edges in total (regardless of label).
             * All edges have length 1.
             *
             * Thus the total length is 12 and the heterotypic boundary
             * length is 5. Note the number of associated neighbour connections is the same.
             *
             * This can be verified by eyeballing the output file.
             */

            // Create an output directory for the writer
            std::string output_directory = "TestHeterotypicBoundaryLengthWriterCa";
            OutputFileHandler output_file_handler(output_directory, false);
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

            // Create a BoundaryNodeWriter and test that the correct output is generated
            HeterotypicBoundaryLengthWriter<2,2> labelled_boundary_writer;
            labelled_boundary_writer.OpenOutputFile(output_file_handler);
            labelled_boundary_writer.WriteTimeStamp();
            labelled_boundary_writer.Visit(&cell_population);
            labelled_boundary_writer.WriteNewline();
            labelled_boundary_writer.CloseFile();

            FileComparison(results_dir + "heterotypicboundary.dat", "cell_based/test/data/TestCellPopulationWriters/heterotypicboundary.dat_ca").CompareFiles();
        }

        // Test with a VertexBasedCellPopulation
        {
            // Create a simple 2D cell population
            HoneycombVertexMeshGenerator generator(4, 4);
            MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
            std::vector<CellPtr> cells;
            boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
            cell_population.InitialiseCells();

            // Label a subset of the cells
            boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
            cell_population.GetCellUsingLocationIndex(1)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(7)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(8)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(9)->AddCellProperty(p_label);
            cell_population.GetCellUsingLocationIndex(12)->AddCellProperty(p_label);

            /*
             * In this case, the group of labelled cells 1, 4, 8, 9, 12 share 9 edges
             * with unlabelled cells, while the isolated labelled cell 7 shares 3 edges
             * with unlabelled cells.
             *
             * Thus there are 12 edges shared between labelled and unlabelled cells.
             * There are 33 shared edges in total (regardless of label).
             * Each edge has length 1/sqrt(3) = 0.577350.
             *
             * Thus the total length is 33/sqrt(3) = 19.052558 and the heterotypic boundary
             * length is 12/sqrt(3) = 6.928203.
             *
             * This can be verified by eyeballing the output file.
             */

            // Create an output directory for the writer
            std::string output_directory = "TestHeterotypicBoundaryLengthWriterVertex";
            OutputFileHandler output_file_handler(output_directory, false);
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

            // Create a BoundaryNodeWriter and test that the correct output is generated
            HeterotypicBoundaryLengthWriter<2,2> labelled_boundary_writer;
            labelled_boundary_writer.OpenOutputFile(output_file_handler);
            labelled_boundary_writer.WriteTimeStamp();
            labelled_boundary_writer.Visit(&cell_population);
            labelled_boundary_writer.WriteNewline();
            labelled_boundary_writer.CloseFile();

            FileComparison(results_dir + "heterotypicboundary.dat", "cell_based/test/data/TestCellPopulationWriters/heterotypicboundary.dat_vertex").CompareFiles();

            // Test that we can append to files
            labelled_boundary_writer.OpenOutputFileForAppend(output_file_handler);
            labelled_boundary_writer.WriteTimeStamp();
            labelled_boundary_writer.Visit(&cell_population);
            labelled_boundary_writer.WriteNewline();
            labelled_boundary_writer.CloseFile();

            FileComparison(results_dir + "heterotypicboundary.dat", "cell_based/test/data/TestCellPopulationWriters/heterotypicboundary.dat_vertex_twice").CompareFiles();
        }
    }

    void TestHeterotypicBoundaryLengthWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "HeterotypicBoundaryLengthWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new HeterotypicBoundaryLengthWriter<2,2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_population_writer;
            delete p_population_writer;
        }
        PetscTools::Barrier(); // Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_population_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_population_writer_2;
            delete p_population_writer_2;
       }
    }

    void TestNodeLocationWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 3D NodeBasedCellPopulation
        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false));
        nodes.push_back(new Node<3>(1, false, 1.0, 1.0, 1.0));
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestNodeLocationWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a NodeLocationWriter and test that the correct output is generated
        NodeLocationWriter<3,3> location_writer;
        location_writer.OpenOutputFile(output_file_handler);
        location_writer.WriteTimeStamp();
        location_writer.Visit(&cell_population);
        location_writer.WriteNewline();
        location_writer.CloseFile();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestCellPopulationWriters/results.viznodes").CompareFiles();

        // Test that we can append to files
        location_writer.OpenOutputFileForAppend(output_file_handler);
        location_writer.WriteTimeStamp();
        location_writer.Visit(&cell_population);
        location_writer.WriteNewline();
        location_writer.CloseFile();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestCellPopulationWriters/results.viznodes_twice").CompareFiles();

        // Avoid memory leaks
        delete nodes[0];
        delete nodes[1];
    }

    void TestNodeLocationWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "NodeLocationWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new NodeLocationWriter<2,2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_population_writer;
            delete p_population_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_population_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_population_writer_2;
            delete p_population_writer_2;
       }
    }

    void TestNodeVelocityWriterWithMeshBasedCellPopulation()
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (needed to avoid tripping an assertion when accessing the time step)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple 3D MeshBasedCellPopulation
        std::vector<Node<3>*> mesh_based_nodes;
        mesh_based_nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        mesh_based_nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        mesh_based_nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        mesh_based_nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        mesh_based_nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh_based_mesh(mesh_based_nodes);
        std::vector<CellPtr> mesh_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> mesh_based_cells_generator;
        mesh_based_cells_generator.GenerateBasic(mesh_based_cells, mesh_based_mesh.GetNumNodes());
        MeshBasedCellPopulation<3> mesh_based_cell_population(mesh_based_mesh, mesh_based_cells);

        // Call ClearAppliedForce() on each node (needed to avoid tripping an assertion when accessing node attributes)
        for (AbstractMesh<3,3>::NodeIterator node_iter = mesh_based_cell_population.rGetMesh().GetNodeIteratorBegin();
             node_iter != mesh_based_cell_population.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        // Add a non-zero force to some nodes
        c_vector<double, 3> force_on_node_2;
        force_on_node_2[0] = 11.0;
        force_on_node_2[1] = 1.3;
        force_on_node_2[2] = 3.7;
        mesh_based_cell_population.GetNode(2)->AddAppliedForceContribution(force_on_node_2);

        c_vector<double, 3> force_on_node_3;
        force_on_node_3[0] = 4.5;
        force_on_node_3[1] = 5.9;
        force_on_node_3[2] = 0.6;
        mesh_based_cell_population.GetNode(3)->AddAppliedForceContribution(force_on_node_3);

        // Create an output directory for the writer
        std::string mesh_based_output_directory = "TestNodeVelocityWriterWithMeshBasedCellPopulation";
        OutputFileHandler mesh_based_output_file_handler(mesh_based_output_directory, false);
        std::string mesh_based_results_dir = mesh_based_output_file_handler.GetOutputDirectoryFullPath();

        // Create a NodeNelocityWriter and test that the correct output is generated
        NodeVelocityWriter<3,3> mesh_based_writer;
        mesh_based_writer.OpenOutputFile(mesh_based_output_file_handler);
        mesh_based_writer.WriteTimeStamp();
        mesh_based_writer.Visit(&mesh_based_cell_population);
        mesh_based_writer.WriteNewline();
        mesh_based_writer.CloseFile();

        FileComparison(mesh_based_results_dir + "nodevelocities.dat", "cell_based/test/data/TestCellPopulationWriters/nodevelocities_mesh.dat").CompareFiles();
    }

    void TestNodeVelocityWriterWithNodeBasedCellPopulation()
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (needed to avoid tripping an assertion when accessing the time step)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple 3D NodeBasedCellPopulation
        std::vector<Node<3>* > node_based_nodes;
        node_based_nodes.push_back(new Node<3>(0, false, 0.0, 0.0, 0.0));
        node_based_nodes.push_back(new Node<3>(1, false, 1.0, 1.0, 1.0));
        NodesOnlyMesh<3> node_based_mesh;
        node_based_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
        std::vector<CellPtr> node_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> node_based_generator;
        node_based_generator.GenerateBasic(node_based_cells, node_based_mesh.GetNumNodes());
        NodeBasedCellPopulation<3> node_based_cell_population(node_based_mesh, node_based_cells);

        // Call ClearAppliedForce() on each node (needed to avoid tripping an assertion when accessing node attributes)
        for (AbstractMesh<3,3>::NodeIterator node_iter = node_based_cell_population.rGetMesh().GetNodeIteratorBegin();
             node_iter != node_based_cell_population.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        // Create an output directory for the writer
        std::string node_based_output_directory = "TestNodeVelocityWriterWithNodeBasedCellPopulation";
        OutputFileHandler node_based_output_file_handler(node_based_output_directory, false);
        std::string node_based_results_dir = node_based_output_file_handler.GetOutputDirectoryFullPath();

        // Create a NodeNelocityWriter and test that the correct output is generated
        NodeVelocityWriter<3,3> node_based_writer;
        node_based_writer.OpenOutputFile(node_based_output_file_handler);
        node_based_writer.WriteTimeStamp();
        node_based_writer.Visit(&node_based_cell_population);
        node_based_writer.WriteNewline();
        node_based_writer.CloseFile();

        // At this time step, the node velocity components should all be zero
        FileComparison(node_based_results_dir + "nodevelocities.dat", "cell_based/test/data/TestCellPopulationWriters/nodevelocities_node.dat").CompareFiles();

        // Now increment time and add a non-zero force for each node
        SimulationTime::Instance()->IncrementTimeOneStep();

        c_vector<double, 3> force_on_node_0;
        force_on_node_0[0] = 1.0;
        force_on_node_0[1] = 2.0;
        force_on_node_0[2] = 3.0;
        node_based_cell_population.GetNode(0)->AddAppliedForceContribution(force_on_node_0);

        c_vector<double, 3> force_on_node_1;
        force_on_node_1[0] = 4.0;
        force_on_node_1[1] = 5.0;
        force_on_node_1[2] = 6.0;
        node_based_cell_population.GetNode(1)->AddAppliedForceContribution(force_on_node_1);

        // Test that we can append to files
        node_based_writer.OpenOutputFileForAppend(node_based_output_file_handler);
        node_based_writer.WriteTimeStamp();
        node_based_writer.Visit(&node_based_cell_population);
        node_based_writer.WriteNewline();
        node_based_writer.CloseFile();

        // At the next time step, the node velocity components should be increasing positive integers
        FileComparison(node_based_results_dir + "nodevelocities.dat", "cell_based/test/data/TestCellPopulationWriters/nodevelocities_node_twice.dat").CompareFiles();

        // Avoid memory leaks
        delete node_based_nodes[0];
        delete node_based_nodes[1];
    }

    void TestNodeVelocityWriterWithVertexBasedCellPopulation()
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (needed to avoid tripping an assertion when accessing the time step)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator vertex_based_generator(4, 6);
        MutableVertexMesh<2,2>* p_vertex_based_mesh = vertex_based_generator.GetMesh();
        std::vector<CellPtr> vertex_based_cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> vertex_based_cells_generator;
        vertex_based_cells_generator.GenerateBasic(vertex_based_cells, p_vertex_based_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> vertex_based_cell_population(*p_vertex_based_mesh, vertex_based_cells);

        // Call ClearAppliedForce() on each node (needed to avoid tripping an assertion when accessing node attributes)
        for (AbstractMesh<2,2>::NodeIterator node_iter = vertex_based_cell_population.rGetMesh().GetNodeIteratorBegin();
             node_iter != vertex_based_cell_population.rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            node_iter->ClearAppliedForce();
        }

        // Add a non-zero force to some nodes
        c_vector<double, 2> force_on_node_0;
        force_on_node_0[0] = 0.01;
        force_on_node_0[1] = 14.8;
        vertex_based_cell_population.GetNode(0)->AddAppliedForceContribution(force_on_node_0);

        c_vector<double, 2> force_on_node_2;
        force_on_node_2[0] = 3.96;
        force_on_node_2[1] = 12.6;
        vertex_based_cell_population.GetNode(2)->AddAppliedForceContribution(force_on_node_2);

        // Create an output directory for the writer
        std::string vertex_based_output_directory = "TestNodeVelocityWriterWithVertexBasedCellPopulation";
        OutputFileHandler vertex_based_output_file_handler(vertex_based_output_directory, false);
        std::string vertex_based_results_dir = vertex_based_output_file_handler.GetOutputDirectoryFullPath();

        // Create a NodeNelocityWriter and test that the correct output is generated
        NodeVelocityWriter<2,2> vertex_based_writer;
        vertex_based_writer.OpenOutputFile(vertex_based_output_file_handler);
        vertex_based_writer.WriteTimeStamp();
        vertex_based_writer.Visit(&vertex_based_cell_population);
        vertex_based_writer.WriteNewline();
        vertex_based_writer.CloseFile();

        FileComparison(vertex_based_results_dir + "nodevelocities.dat", "cell_based/test/data/TestCellPopulationWriters/nodevelocities_vertex.dat").CompareFiles();

    }

    void TestNodeVelocityWriterExceptions()
    {
        EXIT_IF_PARALLEL;

        // Test the correct exception is thrown if using a CaBasedCellPopulation
        PottsMeshGenerator<2> ca_based_generator(4, 0, 0, 4, 0, 0);
        PottsMesh<2>* p_ca_based_mesh = ca_based_generator.GetMesh();
        std::vector<CellPtr> ca_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> ca_based_cells_generator;
        ca_based_cells_generator.GenerateBasic(ca_based_cells, 4);
        std::vector<unsigned> location_indices;
        location_indices.push_back(7);
        location_indices.push_back(11);
        location_indices.push_back(12);
        location_indices.push_back(13);
        CaBasedCellPopulation<2> ca_based_cell_population(*p_ca_based_mesh, ca_based_cells, location_indices);

        NodeVelocityWriter<2,2> node_velocity_writer;
        TS_ASSERT_THROWS_THIS(node_velocity_writer.Visit(&ca_based_cell_population),
            "NodeVelocityWriter cannot be used with a CaBasedCellPopulation");

        // Test the correct exception is thrown if using a PottsBasedCellPopulation
        PottsMeshGenerator<2> potts_based_generator(4, 1, 2, 4, 1, 2);
        PottsMesh<2>* p_potts_based_mesh = potts_based_generator.GetMesh();
        std::vector<CellPtr> potts_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> potts_based_cells_generator;
        potts_based_cells_generator.GenerateBasic(potts_based_cells, p_potts_based_mesh->GetNumElements());
        PottsBasedCellPopulation<2> potts_based_cell_population(*p_potts_based_mesh, potts_based_cells);

        TS_ASSERT_THROWS_THIS(node_velocity_writer.Visit(&potts_based_cell_population),
            "NodeVelocityWriter cannot be used with a PottsBasedCellPopulation");
    }

    void TestNodeVelocityWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "NodeVelocityWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new NodeVelocityWriter<2,2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_population_writer;
            delete p_population_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_population_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_population_writer_2;
            delete p_population_writer_2;
       }
    }

    void TestRadialCellDataDistributionWriter()
    {
        EXIT_IF_PARALLEL;

        // Test with a VerexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 4);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        double value = 0.0;
        for (AbstractCellPopulation<2>::Iterator cell_iter=cell_population.Begin();
             cell_iter!=cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("this average", value);
            value += 1.0;
        }
        // Create an output directory for the writer
        std::string output_directory = "TestRadialCellDataDistributionWriterVertex";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a RadialCellDataDistributionWriter and test that the correct output is generated
        RadialCellDataDistributionWriter<2,2> radial_writer;
        radial_writer.SetVariableName("this average");
        radial_writer.SetNumRadialBins(3);
        radial_writer.OpenOutputFile(output_file_handler);
        radial_writer.WriteTimeStamp();
        radial_writer.Visit(&cell_population);
        radial_writer.WriteNewline();
        radial_writer.CloseFile();

        FileComparison(results_dir + "radial_dist.dat", "cell_based/test/data/TestCellPopulationWriters/radial_dist.dat").CompareFiles();

        // Test that we can append to files
        radial_writer.OpenOutputFileForAppend(output_file_handler);
        radial_writer.WriteTimeStamp();
        radial_writer.Visit(&cell_population);
        radial_writer.WriteNewline();
        radial_writer.CloseFile();

        FileComparison(results_dir + "radial_dist.dat", "cell_based/test/data/TestCellPopulationWriters/radial_dist_twice.dat").CompareFiles();

        ///\todo Improve tests below (#2847)

        // Test with a MeshBasedCellPopulation
        {
            HoneycombMeshGenerator tet_generator(5, 5, 0);
            MutableMesh<2,2>* p_tet_mesh = tet_generator.GetMesh();
            std::vector<CellPtr> mesh_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> mesh_based_cells_generator;
            mesh_based_cells_generator.GenerateBasic(mesh_based_cells, p_tet_mesh->GetNumNodes());
            MeshBasedCellPopulation<2> mesh_based_cell_population(*p_tet_mesh, mesh_based_cells);
            for (AbstractCellPopulation<2>::Iterator cell_iter=mesh_based_cell_population.Begin();
                 cell_iter!=mesh_based_cell_population.End();
                 ++cell_iter)
            {
                 cell_iter->GetCellData()->SetItem("this average", 1.0);
            }
            radial_writer.Visit(&mesh_based_cell_population);
        }

        // Test with a CaBasedCellPopulation
        {
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
            for (AbstractCellPopulation<2>::Iterator cell_iter=ca_based_cell_population.Begin();
                 cell_iter!=ca_based_cell_population.End();
                 ++cell_iter)
            {
                 cell_iter->GetCellData()->SetItem("this average", 1.0);
            }
            radial_writer.Visit(&ca_based_cell_population);
        }

        // Test with a NodeBasedCellPopulation
        {
            std::vector<Node<2>* > node_based_nodes;
            node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
            NodesOnlyMesh<2> node_based_mesh;
            node_based_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
            std::vector<CellPtr> node_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
            node_based_generator.GenerateBasic(node_based_cells, node_based_mesh.GetNumNodes());
            NodeBasedCellPopulation<2> node_based_cell_population(node_based_mesh, node_based_cells);
            for (AbstractCellPopulation<2>::Iterator cell_iter=node_based_cell_population.Begin();
                 cell_iter!=node_based_cell_population.End();
                 ++cell_iter)
            {
                 cell_iter->GetCellData()->SetItem("this average", 1.0);
            }
            TS_ASSERT_THROWS_NOTHING(radial_writer.Visit(&node_based_cell_population));

            // Tidy up
            delete node_based_nodes[0];
            delete node_based_nodes[1];
        }

        // Test with a PottsBasedCellPopulation
        {
            PottsMeshGenerator<2> potts_based_generator(4, 1, 2, 4, 1, 2);
            PottsMesh<2>* p_potts_based_mesh = potts_based_generator.GetMesh();
            std::vector<CellPtr> potts_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> potts_based_cells_generator;
            potts_based_cells_generator.GenerateBasic(potts_based_cells, p_potts_based_mesh->GetNumElements());
            PottsBasedCellPopulation<2> potts_based_cell_population(*p_potts_based_mesh, potts_based_cells);
            for (AbstractCellPopulation<2>::Iterator cell_iter=potts_based_cell_population.Begin();
                 cell_iter!=potts_based_cell_population.End();
                 ++cell_iter)
            {
                 cell_iter->GetCellData()->SetItem("this average", 1.0);
            }
            TS_ASSERT_THROWS_NOTHING(radial_writer.Visit(&potts_based_cell_population));
        }
    }

    void TestRadialCellDataDistributionWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "RadialCellDataDistributionWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new RadialCellDataDistributionWriter<2,2>();
            static_cast<RadialCellDataDistributionWriter<2,2>*>(p_cell_writer)->SetVariableName("radial average");
            static_cast<RadialCellDataDistributionWriter<2,2>*>(p_cell_writer)->SetNumRadialBins(5);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_cell_writer;

            delete p_cell_writer;
        }
        PetscTools::Barrier(); // Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell_writer_2;

            typedef RadialCellDataDistributionWriter<2,2> RadialWriter;
            TS_ASSERT_EQUALS(static_cast<RadialWriter*>(p_cell_writer_2)->GetVariableName(), "radial average");
            TS_ASSERT_EQUALS(static_cast<RadialWriter*>(p_cell_writer_2)->GetNumRadialBins(), 5u);
            delete p_cell_writer_2;
       }
    }

    void TestVertexT1SwapLocationsWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestVertexT1SwapLocationsWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a VertexT1SwapLocationsWriter and test that the correct output is generated
        VertexT1SwapLocationsWriter<2,2> t1_swaps_writer;
        t1_swaps_writer.OpenOutputFile(output_file_handler);
        t1_swaps_writer.WriteTimeStamp();
        t1_swaps_writer.Visit(&cell_population);
        t1_swaps_writer.WriteNewline();
        t1_swaps_writer.CloseFile();

        FileComparison(results_dir + "T1SwapLocations.dat", "cell_based/test/data/TestCellPopulationWriters/T1SwapLocations.dat").CompareFiles();

        // Test that we can append to files
        t1_swaps_writer.OpenOutputFileForAppend(output_file_handler);
        t1_swaps_writer.WriteTimeStamp();
        t1_swaps_writer.Visit(&cell_population);
        t1_swaps_writer.WriteNewline();
        t1_swaps_writer.CloseFile();

        FileComparison(results_dir + "T1SwapLocations.dat", "cell_based/test/data/TestCellPopulationWriters/T1SwapLocations_twice.dat").CompareFiles();

        {
            // Coverage of the Visit() method when called on a MeshBasedCellPopulation
            HoneycombMeshGenerator tet_generator(5, 5, 0);
            MutableMesh<2,2>* p_tet_mesh = tet_generator.GetMesh();
            std::vector<CellPtr> mesh_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> mesh_based_cells_generator;
            mesh_based_cells_generator.GenerateBasic(mesh_based_cells, p_tet_mesh->GetNumNodes());
            MeshBasedCellPopulation<2> mesh_based_cell_population(*p_tet_mesh, mesh_based_cells);

            TS_ASSERT_THROWS_NOTHING(t1_swaps_writer.Visit(&mesh_based_cell_population));
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

            TS_ASSERT_THROWS_NOTHING(t1_swaps_writer.Visit(&ca_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a NodeBasedCellPopulation
            std::vector<Node<2>* > node_based_nodes;
            node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
            NodesOnlyMesh<2> node_based_mesh;
            node_based_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
            std::vector<CellPtr> node_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
            node_based_generator.GenerateBasic(node_based_cells, node_based_mesh.GetNumNodes());
            NodeBasedCellPopulation<2> node_based_cell_population(node_based_mesh, node_based_cells);

            TS_ASSERT_THROWS_NOTHING(t1_swaps_writer.Visit(&node_based_cell_population));

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

            TS_ASSERT_THROWS_NOTHING(t1_swaps_writer.Visit(&potts_based_cell_population));
        }
    }

    void TestVertexT1SwapLocationsWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "VertexT1SwapLocationsWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new VertexT1SwapLocationsWriter<2,2>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_cell_writer;
            delete p_cell_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell_writer_2;
            delete p_cell_writer_2;
       }
    }

    void TestVertexT2SwapLocationsWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestVertexT2SwapLocationsWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a VertexT2SwapLocationsWriter and test that the correct output is generated
        VertexT2SwapLocationsWriter<2,2> t2_swaps_writer;
        t2_swaps_writer.OpenOutputFile(output_file_handler);
        t2_swaps_writer.WriteTimeStamp();
        t2_swaps_writer.Visit(&cell_population);
        t2_swaps_writer.WriteNewline();
        t2_swaps_writer.CloseFile();

        FileComparison(results_dir + "T2SwapLocations.dat", "cell_based/test/data/TestCellPopulationWriters/T2SwapLocations.dat").CompareFiles();

        // Test that we can append to files
        t2_swaps_writer.OpenOutputFileForAppend(output_file_handler);
        t2_swaps_writer.WriteTimeStamp();
        t2_swaps_writer.Visit(&cell_population);
        t2_swaps_writer.WriteNewline();
        t2_swaps_writer.CloseFile();

        FileComparison(results_dir + "T2SwapLocations.dat", "cell_based/test/data/TestCellPopulationWriters/T2SwapLocations_twice.dat").CompareFiles();

        {
            // Coverage of the Visit() method when called on a MeshBasedCellPopulation
            HoneycombMeshGenerator tet_generator(5, 5, 0);
            MutableMesh<2,2>* p_tet_mesh = tet_generator.GetMesh();
            std::vector<CellPtr> mesh_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> mesh_based_cells_generator;
            mesh_based_cells_generator.GenerateBasic(mesh_based_cells, p_tet_mesh->GetNumNodes());
            MeshBasedCellPopulation<2> mesh_based_cell_population(*p_tet_mesh, mesh_based_cells);

            TS_ASSERT_THROWS_NOTHING(t2_swaps_writer.Visit(&mesh_based_cell_population));
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

            TS_ASSERT_THROWS_NOTHING(t2_swaps_writer.Visit(&ca_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a NodeBasedCellPopulation
            std::vector<Node<2>* > node_based_nodes;
            node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
            NodesOnlyMesh<2> node_based_mesh;
            node_based_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
            std::vector<CellPtr> node_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
            node_based_generator.GenerateBasic(node_based_cells, node_based_mesh.GetNumNodes());
            NodeBasedCellPopulation<2> node_based_cell_population(node_based_mesh, node_based_cells);

            TS_ASSERT_THROWS_NOTHING(t2_swaps_writer.Visit(&node_based_cell_population));

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

            TS_ASSERT_THROWS_NOTHING(t2_swaps_writer.Visit(&potts_based_cell_population));
        }
    }

    void TestVertexT2SwapLocationsWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "VertexT2SwapLocationsWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new VertexT2SwapLocationsWriter<2,2>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_cell_writer;
            delete p_cell_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell_writer_2;
            delete p_cell_writer_2;
       }
    }

    void TestVertexT3SwapLocationsWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestVertexT2SwapLocationsWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a VertexT1SwapLocationsWriter and test that the correct output is generated
        VertexT3SwapLocationsWriter<2,2> t3_swaps_writer;
        t3_swaps_writer.OpenOutputFile(output_file_handler);
        t3_swaps_writer.WriteTimeStamp();
        t3_swaps_writer.Visit(&cell_population);
        t3_swaps_writer.WriteNewline();
        t3_swaps_writer.CloseFile();

        FileComparison(results_dir + "T3SwapLocations.dat", "cell_based/test/data/TestCellPopulationWriters/T3SwapLocations.dat").CompareFiles();

        // Test that we can append to files
        t3_swaps_writer.OpenOutputFileForAppend(output_file_handler);
        t3_swaps_writer.WriteTimeStamp();
        t3_swaps_writer.Visit(&cell_population);
        t3_swaps_writer.WriteNewline();
        t3_swaps_writer.CloseFile();

        FileComparison(results_dir + "T3SwapLocations.dat", "cell_based/test/data/TestCellPopulationWriters/T3SwapLocations_twice.dat").CompareFiles();

        {
            // Coverage of the Visit() method when called on a MeshBasedCellPopulation
            HoneycombMeshGenerator tet_generator(5, 5, 0);
            MutableMesh<2,2>* p_tet_mesh = tet_generator.GetMesh();
            std::vector<CellPtr> mesh_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> mesh_based_cells_generator;
            mesh_based_cells_generator.GenerateBasic(mesh_based_cells, p_tet_mesh->GetNumNodes());
            MeshBasedCellPopulation<2> mesh_based_cell_population(*p_tet_mesh, mesh_based_cells);

            TS_ASSERT_THROWS_NOTHING(t3_swaps_writer.Visit(&mesh_based_cell_population));
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

            TS_ASSERT_THROWS_NOTHING(t3_swaps_writer.Visit(&ca_based_cell_population));
        }

        {
            // Coverage of the Visit() method when called on a NodeBasedCellPopulation
            std::vector<Node<2>* > node_based_nodes;
            node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
            node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
            NodesOnlyMesh<2> node_based_mesh;
            node_based_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
            std::vector<CellPtr> node_based_cells;
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
            node_based_generator.GenerateBasic(node_based_cells, node_based_mesh.GetNumNodes());
            NodeBasedCellPopulation<2> node_based_cell_population(node_based_mesh, node_based_cells);

            TS_ASSERT_THROWS_NOTHING(t3_swaps_writer.Visit(&node_based_cell_population));

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

            TS_ASSERT_THROWS_NOTHING(t3_swaps_writer.Visit(&potts_based_cell_population));
        }
    }

    void TestVertexT3SwapLocationsWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "VertexT3SwapLocationsWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new VertexT3SwapLocationsWriter<2,2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_cell_writer;
            delete p_cell_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_cell_writer_2;
            delete p_cell_writer_2;
       }
    }

    void TestVoronoiDataWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 3D MeshBasedCellPopulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        cell_population.CreateVoronoiTessellation();

        // Create an output directory for the writer
        std::string output_directory = "TestVoronoiDataWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a VoronoiDataWriter and test that the correct output is generated
        VoronoiDataWriter<3,3> voronoi_writer;
        voronoi_writer.OpenOutputFile(output_file_handler);
        voronoi_writer.WriteTimeStamp();
        voronoi_writer.Visit(&cell_population);
        voronoi_writer.WriteNewline();
        voronoi_writer.CloseFile();

        FileComparison(results_dir + "voronoi.dat", "cell_based/test/data/TestCellPopulationWriters/voronoi.dat").CompareFiles();

        // Test that we can append to files
        voronoi_writer.OpenOutputFileForAppend(output_file_handler);
        voronoi_writer.WriteTimeStamp();
        voronoi_writer.Visit(&cell_population);
        voronoi_writer.WriteNewline();
        voronoi_writer.CloseFile();

        FileComparison(results_dir + "voronoi.dat", "cell_based/test/data/TestCellPopulationWriters/voronoi_twice.dat").CompareFiles();

        VoronoiDataWriter<2,2> voronoi_writer_2d;

        // Test the correct exception is thrown if using a NodeBasedCellPopulation
        std::vector<Node<2>* > node_based_nodes;
        node_based_nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        node_based_nodes.push_back(new Node<2>(1, false, 1.0, 1.0));
        NodesOnlyMesh<2> nodes_only_mesh;
        nodes_only_mesh.ConstructNodesWithoutMesh(node_based_nodes, 1.5);
        std::vector<CellPtr> node_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> node_based_generator;
        node_based_generator.GenerateBasic(node_based_cells, nodes_only_mesh.GetNumNodes());
        NodeBasedCellPopulation<2> node_based_cell_population(nodes_only_mesh, node_based_cells);

        TS_ASSERT_THROWS_THIS(voronoi_writer_2d.Visit(&node_based_cell_population),
            "VoronoiDataWriter cannot be used with a NodeBasedCellPopulation");

        // Tidy up
        for (unsigned i=0; i<node_based_nodes.size(); i++)
        {
            delete node_based_nodes[i];
        }

        // Test the correct exception is thrown if using a CaBasedCellPopulation
        PottsMeshGenerator<2> ca_based_generator(4, 0, 0, 4, 0, 0);
        PottsMesh<2>* p_ca_based_mesh = ca_based_generator.GetMesh();
        std::vector<CellPtr> ca_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> ca_based_cells_generator;
        ca_based_cells_generator.GenerateBasic(ca_based_cells, 4);
        std::vector<unsigned> location_indices;
        location_indices.push_back(7);
        location_indices.push_back(11);
        location_indices.push_back(12);
        location_indices.push_back(13);
        CaBasedCellPopulation<2> ca_based_cell_population(*p_ca_based_mesh, ca_based_cells, location_indices);

        TS_ASSERT_THROWS_THIS(voronoi_writer_2d.Visit(&ca_based_cell_population),
            "VoronoiDataWriter cannot be used with a CaBasedCellPopulation");

        // Test the correct exception is thrown if using a PottsBasedCellPopulation
        PottsMeshGenerator<2> potts_based_generator(4, 1, 2, 4, 1, 2);
        PottsMesh<2>* p_potts_based_mesh = potts_based_generator.GetMesh();
        std::vector<CellPtr> potts_based_cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> potts_based_cells_generator;
        potts_based_cells_generator.GenerateBasic(potts_based_cells, p_potts_based_mesh->GetNumElements());
        PottsBasedCellPopulation<2> potts_based_cell_population(*p_potts_based_mesh, potts_based_cells);

        TS_ASSERT_THROWS_THIS(voronoi_writer_2d.Visit(&potts_based_cell_population),
            "VoronoiDataWriter cannot be used with a PottsBasedCellPopulation");

        // Test the correct exception is thrown if using a PottsBasedCellPopulation
        HoneycombVertexMeshGenerator vertex_based_generator(4, 6);
        MutableVertexMesh<2,2>* p_vertex_mesh = vertex_based_generator.GetMesh();
        std::vector<CellPtr> vertex_based_cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> vertex_based_cells_generator;
        vertex_based_cells_generator.GenerateBasic(vertex_based_cells, p_vertex_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);
        VertexBasedCellPopulation<2> vertex_cell_population(*p_vertex_mesh, vertex_based_cells);

        TS_ASSERT_THROWS_THIS(voronoi_writer_2d.Visit(&vertex_cell_population),
            "VoronoiDataWriter cannot be used with a VertexBasedCellPopulation");
    }

    void TestVoronoiDataWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "VoronoiDataWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new VoronoiDataWriter<2,2>();
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);
            output_arch << p_population_writer;
            delete p_population_writer;
        }
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_population_writer_2;
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);
            input_arch >> p_population_writer_2;
            delete p_population_writer_2;
       }
    }
};

#endif /*TESTCELLPOPULATIONWRITERS_HPP_*/
