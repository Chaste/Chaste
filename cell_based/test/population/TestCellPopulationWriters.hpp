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

#ifndef TESTCELLPOPULATIONWRITERS_HPP_
#define TESTCELLPOPULATIONWRITERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"

// Cell writers
#include "CellIdWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"

// Cell population writers
#include "NodeLocationWriter.hpp"
#include "BoundaryNodeWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellPopulationAreaWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"
#include "VertexT3SwapLocationsWriter.hpp"
#include "VoronoiDataWriter.hpp"

// Files to create populations
#include "HoneycombVertexMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"

#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCellPopulationWriters : public AbstractCellBasedTestSuite
{
public:

    void TestNodeLocationWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 3D NodeBasedCellPopulation
        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false));
        nodes.push_back(new Node<3>(1, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestNodeLocationWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a NodeLocationWriter and test that the correct output is generated
        NodeLocationWriter<3,3> location_writer;
        location_writer.OpenOutputFile(output_directory);
        location_writer.WriteTimeStamp();
        location_writer.Visit(&cell_population);
        location_writer.WriteNewline();
        location_writer.CloseFile();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestNodeLocationWriter/results.viznodes").CompareFiles();

        // Test that we can append to files
        location_writer.OpenOutputFileForAppend(output_directory);
        location_writer.WriteTimeStamp();
        location_writer.Visit(&cell_population);
        location_writer.WriteNewline();
        location_writer.CloseFile();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestNodeLocationWriter/results.viznodes_twice").CompareFiles();

        // Avoid memory leaks
        delete nodes[0];
        delete nodes[1];
    }

    void TestBoundaryNodeWriter() throw (Exception)
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
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestBoundaryNodeWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a BoundaryNodeWriter and test that the correct output is generated
        BoundaryNodeWriter<3,3> boundary_writer;
        boundary_writer.OpenOutputFile(output_directory);
        boundary_writer.WriteTimeStamp();
        boundary_writer.Visit(&cell_population);
        boundary_writer.WriteNewline();
        boundary_writer.CloseFile();

        FileComparison(results_dir + "results.vizboundarynodes", "cell_based/test/data/TestBoundaryNodeWriter/results.vizboundarynodes").CompareFiles();

        // Test that we can append to files
        boundary_writer.OpenOutputFileForAppend(output_directory);
        boundary_writer.WriteTimeStamp();
        boundary_writer.Visit(&cell_population);
        boundary_writer.WriteNewline();
        boundary_writer.CloseFile();

        FileComparison(results_dir + "results.vizboundarynodes", "cell_based/test/data/TestBoundaryNodeWriter/results.vizboundarynodes_twice").CompareFiles();
    }

    void TestCellPopulationElementWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestCellPopulationElementWriter";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellPopulationElementWriter and test that the correct output is generated
        CellPopulationElementWriter<2,2> element_writer;
        element_writer.OpenOutputFile(output_directory);
        element_writer.WriteTimeStamp();
        element_writer.Visit(&cell_population);
        element_writer.WriteNewline();
        element_writer.CloseFile();

        FileComparison(results_dir + "results.vizelements", "cell_based/test/data/TestCellPopulationElementWriter/results.vizelements").CompareFiles();

        // Test that we can append to files
        element_writer.OpenOutputFileForAppend(output_directory);
        element_writer.WriteTimeStamp();
        element_writer.Visit(&cell_population);
        element_writer.WriteNewline();
        element_writer.CloseFile();

        FileComparison(results_dir + "results.vizelements", "cell_based/test/data/TestCellPopulationElementWriter/results.vizelements_twice").CompareFiles();
    }

    void TestCellMutationStatesWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D MultipleCaBasedCellPopulation
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 5u);

        std::vector<unsigned> location_indices;
        location_indices.push_back(7);
        location_indices.push_back(11);
        location_indices.push_back(12);
        location_indices.push_back(13);
        location_indices.push_back(17);

        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddPopulationWriter<CellMutationStatesWriter>();
        cell_population.GenerateCellResults();

        // Create an output directory for the writer
        std::string output_directory = "TestCellMutationStatesWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellMutationStatesWriter and test that the correct output is generated
        CellMutationStatesWriter<2,2> mutation_states_writer;
        mutation_states_writer.OpenOutputFile(output_directory);
        mutation_states_writer.WriteHeader(&cell_population);
        mutation_states_writer.WriteTimeStamp();
        mutation_states_writer.Visit(&cell_population);
        mutation_states_writer.WriteNewline();
        mutation_states_writer.CloseFile();

        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestCellMutationStatesWriter/cellmutationstates.dat").CompareFiles();

        // Test that we can append to files
        mutation_states_writer.OpenOutputFileForAppend(output_directory);
        mutation_states_writer.WriteTimeStamp();
        mutation_states_writer.Visit(&cell_population);
        mutation_states_writer.WriteNewline();
        mutation_states_writer.CloseFile();

        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestCellMutationStatesWriter/cellmutationstates_twice.dat").CompareFiles();
    }

    void TestArchivingOfCellMutationStatesWriter() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellMutationStatesWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new CellMutationStatesWriter<2,2>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_population_writer;

            delete p_population_writer;
        }

        {
            AbstractCellBasedWriter<2,2>* p_population_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_population_writer_2;

            delete p_population_writer_2;
       }
    }

    void TestCellProliferativeTypesAndPhasesCountWriters() throw (Exception)
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
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        cell_population.AddPopulationWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddPopulationWriter<CellProliferativePhasesCountWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        cell_population.GenerateCellResults();

        // Create an output directory for the writer
        std::string output_directory = "TestCellProliferativeTypesAndPhasesCountWriters";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellProliferativeTypesCountWriter and test that the correct output is generated
        CellProliferativeTypesCountWriter<3,3> types_count_writer;
        types_count_writer.OpenOutputFile(output_directory);
        types_count_writer.WriteTimeStamp();
        types_count_writer.Visit(&cell_population);
        types_count_writer.WriteNewline();
        types_count_writer.CloseFile();

        FileComparison(results_dir + "celltypes.dat", "cell_based/test/data/TestCellProliferativeTypesAndPhasesCountWriters/celltypes.dat").CompareFiles();

        // Test that we can append to files
        types_count_writer.OpenOutputFileForAppend(output_directory);
        types_count_writer.WriteTimeStamp();
        types_count_writer.Visit(&cell_population);
        types_count_writer.WriteNewline();
        types_count_writer.CloseFile();

        FileComparison(results_dir + "celltypes.dat", "cell_based/test/data/TestCellProliferativeTypesAndPhasesCountWriters/celltypes_twice.dat").CompareFiles();

        // Create a CellProliferativePhasesCountWriter and test that the correct output is generated
        CellProliferativePhasesCountWriter<3,3> phases_count_writer;
        phases_count_writer.OpenOutputFile(output_directory);
        phases_count_writer.WriteTimeStamp();
        phases_count_writer.Visit(&cell_population);
        phases_count_writer.WriteNewline();
        phases_count_writer.CloseFile();

        FileComparison(results_dir + "cellcyclephases.dat", "cell_based/test/data/TestCellProliferativeTypesAndPhasesCountWriters/cellcyclephases.dat").CompareFiles();

        // Test that we can append to files
        phases_count_writer.OpenOutputFileForAppend(output_directory);
        phases_count_writer.WriteTimeStamp();
        phases_count_writer.Visit(&cell_population);
        phases_count_writer.WriteNewline();
        phases_count_writer.CloseFile();

        FileComparison(results_dir + "cellcyclephases.dat", "cell_based/test/data/TestCellProliferativeTypesAndPhasesCountWriters/cellcyclephases_twice.dat").CompareFiles();
    }

    void TestVoronoiDataWriter() throw (Exception)
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
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.CreateVoronoiTessellation();

        // Create an output directory for the writer
        std::string output_directory = "TestVoronoiDataWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a VoronoiDataWriter and test that the correct output is generated
        VoronoiDataWriter<3,3> voronoi_writer;
        voronoi_writer.OpenOutputFile(output_directory);
        voronoi_writer.WriteTimeStamp();
        voronoi_writer.Visit(&cell_population);
        voronoi_writer.WriteNewline();
        voronoi_writer.CloseFile();

        FileComparison(results_dir + "voronoi.dat", "cell_based/test/data/TestVoronoiDataWriter/voronoi.dat").CompareFiles();

        // Test that we can append to files
        voronoi_writer.OpenOutputFileForAppend(output_directory);
        voronoi_writer.WriteTimeStamp();
        voronoi_writer.Visit(&cell_population);
        voronoi_writer.WriteNewline();
        voronoi_writer.CloseFile();

        FileComparison(results_dir + "voronoi.dat", "cell_based/test/data/TestVoronoiDataWriter/voronoi_twice.dat").CompareFiles();
    }

    void TestVertexT1AndT3SwapLocationsWriters() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create an output directory for the writer
        std::string output_directory = "TestVertexT1AndT3SwapLocationsWriters";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a VertexT1SwapLocationsWriter and test that the correct output is generated
        VertexT1SwapLocationsWriter<2,2> t1_swaps_writer;
        t1_swaps_writer.OpenOutputFile(output_directory);
        t1_swaps_writer.WriteTimeStamp();
        t1_swaps_writer.Visit(&cell_population);
        t1_swaps_writer.WriteNewline();
        t1_swaps_writer.CloseFile();

        FileComparison(results_dir + "T1SwapLocations.dat", "cell_based/test/data/TestVertexT1AndT3SwapLocationsWriters/T1SwapLocations.dat").CompareFiles();

        // Test that we can append to files
        t1_swaps_writer.OpenOutputFileForAppend(output_directory);
        t1_swaps_writer.WriteTimeStamp();
        t1_swaps_writer.Visit(&cell_population);
        t1_swaps_writer.WriteNewline();
        t1_swaps_writer.CloseFile();

        FileComparison(results_dir + "T1SwapLocations.dat", "cell_based/test/data/TestVertexT1AndT3SwapLocationsWriters/T1SwapLocations_twice.dat").CompareFiles();

        // Create a VertexT3SwapLocationsWriter and test that the correct output is generated
        VertexT3SwapLocationsWriter<2,2> t3_swaps_writer;
        t3_swaps_writer.OpenOutputFile(output_directory);
        t3_swaps_writer.WriteTimeStamp();
        t3_swaps_writer.Visit(&cell_population);
        t3_swaps_writer.WriteNewline();
        t3_swaps_writer.CloseFile();

        FileComparison(results_dir + "T3SwapLocations.dat", "cell_based/test/data/TestVertexT1AndT3SwapLocationsWriters/T3SwapLocations.dat").CompareFiles();

        // Test that we can append to files
        t3_swaps_writer.OpenOutputFileForAppend(output_directory);
        t3_swaps_writer.WriteTimeStamp();
        t3_swaps_writer.Visit(&cell_population);
        t3_swaps_writer.WriteNewline();
        t3_swaps_writer.CloseFile();

        FileComparison(results_dir + "T3SwapLocations.dat", "cell_based/test/data/TestVertexT1AndT3SwapLocationsWriters/T3SwapLocations_twice.dat").CompareFiles();
    }

    void TestCellPopulationAreaWriter() throw (Exception)
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
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.CreateVoronoiTessellation();

        // Create an output directory for the writer
        std::string output_directory = "TestCellPopulationAreaWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellMutationStatesWriter and test that the correct output is generated
        CellPopulationAreaWriter<3,3> area_writer;
        area_writer.OpenOutputFile(output_directory);
        area_writer.WriteTimeStamp();
        area_writer.Visit(&cell_population);
        area_writer.WriteNewline();
        area_writer.CloseFile();

        FileComparison(results_dir + "cellpopulationareas.dat", "cell_based/test/data/TestCellPopulationAreaWriter/cellpopulationareas.dat").CompareFiles();

        // Test that we can append to files
        area_writer.OpenOutputFileForAppend(output_directory);
        area_writer.WriteTimeStamp();
        area_writer.Visit(&cell_population);
        area_writer.WriteNewline();
        area_writer.CloseFile();

        FileComparison(results_dir + "cellpopulationareas.dat", "cell_based/test/data/TestCellPopulationAreaWriter/cellpopulationareas_twice.dat").CompareFiles();
    }

    void TestNodeVelocityWriterWithMeshBasedCellPopulation() throw (Exception)
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
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> mesh_based_cells_generator;
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
        mesh_based_writer.OpenOutputFile(mesh_based_output_directory);
        mesh_based_writer.WriteTimeStamp();
        mesh_based_writer.Visit(&mesh_based_cell_population);
        mesh_based_writer.WriteNewline();
        mesh_based_writer.CloseFile();

        FileComparison(mesh_based_results_dir + "nodevelocities.dat", "cell_based/test/data/TestNodeVelocityWriter/nodevelocities_mesh.dat").CompareFiles();
    }

    void TestNodeVelocityWriterWithNodeBasedCellPopulation() throw (Exception)
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
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> node_based_generator;
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
        node_based_writer.OpenOutputFile(node_based_output_directory);
        node_based_writer.WriteTimeStamp();
        node_based_writer.Visit(&node_based_cell_population);
        node_based_writer.WriteNewline();
        node_based_writer.CloseFile();

        // At this time step, the node velocity components should all be zero
        FileComparison(node_based_results_dir + "nodevelocities.dat", "cell_based/test/data/TestNodeVelocityWriter/nodevelocities_node.dat").CompareFiles();

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
        node_based_writer.OpenOutputFileForAppend(node_based_output_directory);
        node_based_writer.WriteTimeStamp();
        node_based_writer.Visit(&node_based_cell_population);
        node_based_writer.WriteNewline();
        node_based_writer.CloseFile();

        // At the next time step, the node velocity components should be increasing positive integers
        FileComparison(node_based_results_dir + "nodevelocities.dat", "cell_based/test/data/TestNodeVelocityWriter/nodevelocities_node_twice.dat").CompareFiles();

        // Avoid memory leaks
        delete node_based_nodes[0];
        delete node_based_nodes[1];
    }

    void TestNodeVelocityWriterWithVertexBasedCellPopulation() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (needed to avoid tripping an assertion when accessing the time step)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple 2D VertexBasedCellPopulation
        HoneycombVertexMeshGenerator vertex_based_generator(4, 6);
        MutableVertexMesh<2,2>* p_vertex_based_mesh = vertex_based_generator.GetMesh();
        std::vector<CellPtr> vertex_based_cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> vertex_based_cells_generator;
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
        vertex_based_writer.OpenOutputFile(vertex_based_output_directory);
        vertex_based_writer.WriteTimeStamp();
        vertex_based_writer.Visit(&vertex_based_cell_population);
        vertex_based_writer.WriteNewline();
        vertex_based_writer.CloseFile();

        FileComparison(vertex_based_results_dir + "nodevelocities.dat", "cell_based/test/data/TestNodeVelocityWriter/nodevelocities_vertex.dat").CompareFiles();
    }

    void TestAddWritersToAPopulation() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a 3D NodeBasedCellPopulation
        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false));
        nodes.push_back(new Node<3>(1, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Create a writer and test that it is correctly added to the cell population
        cell_population.AddPopulationWriter<CellPopulationElementWriter>();

        ///\todo test something here (#2404, #2441)

        // Create another writer and test that it is correctly added to the cell population
        cell_population.AddCellWriter<CellIdWriter>();

        ///\todo test something here (#2404, #2441)

        // Avoid memory leaks (note that the writers are deleted by the population destructor)
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestArchivingOfCellPopulationAreaWriter() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellPopulationAreaWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellPopulationAreaWriter<2,2>();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_cell_writer;

            delete p_cell_writer;
        }

        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;

            delete p_cell_writer_2;
       }
    }

    void TestArchivingOfVertexT1SwapLocationsWriter() throw (Exception)
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

        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;

            delete p_cell_writer_2;
       }
    }

    void TestArchivingOfVertexT3SwapLocationsWriter() throw (Exception)
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

        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;

            delete p_cell_writer_2;
       }
    }
};

#endif /*TESTCELLPOPULATIONWRITERS_HPP_*/
