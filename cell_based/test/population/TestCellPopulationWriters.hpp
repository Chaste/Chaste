/*

Copyright (c) 2005-2013, University of Oxford.
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

#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"

// Writers
#include "NodeLocationWriter.hpp"
#include "BoundaryNodeWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "VertexSwapWriters.hpp"
#include "CellWriters.hpp"

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

        // NodeBasedCellPopulation
        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false));
        nodes.push_back(new Node<3>(1, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        std::string output_directory = "TestWriteNodeLocations";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a node location writer and write the files.
        NodeLocationWriter<3,3> location_writer(output_directory);

        location_writer.OpenOutputFile();

        location_writer.Visit(&cell_population);

        location_writer.CloseFile();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestNodeLocationWriter/results.viznodes").CompareFiles();

        // Make sure we can append to files.
        location_writer.OpenOutputFileForAppend();

        location_writer.Visit(&cell_population);

        location_writer.CloseFile();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestNodeLocationWriter/results.viznodes_twice").CompareFiles();

        delete nodes[0];
        delete nodes[1];
    }

    void TestBoundaryNodeWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);

        std::string output_directory = "TestWriteBoundaryNodes";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a node location writer and write the files.
        BoundaryNodeWriter<3,3> boundary_writer(output_directory);

        boundary_writer.OpenOutputFile();

        boundary_writer.Visit(&cell_population);

        boundary_writer.CloseFile();

        FileComparison(results_dir + "results.vizboundarynodes", "cell_based/test/data/TestBoundaryNodeWriter/results.vizboundarynodes").CompareFiles();

        // Make sure we can append to files.
        boundary_writer.OpenOutputFileForAppend();

        boundary_writer.Visit(&cell_population);

        boundary_writer.CloseFile();

        FileComparison(results_dir + "results.vizboundarynodes", "cell_based/test/data/TestBoundaryNodeWriter/results.vizboundarynodes_twice").CompareFiles();
    }

    void TestCellPopulationElementWriter()    throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        std::string output_directory = "TestWriteCellPopulationElements";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        CellPopulationElementWriter<2,2> element_writer(output_directory);

        element_writer.OpenOutputFile();

        element_writer.Visit(&cell_population);

        element_writer.CloseFile();

        FileComparison(results_dir + "results.vizelements", "cell_based/test/data/TestCellPopulationElementWriter/results.vizelements").CompareFiles();

        // Make sure we can append to files.
        element_writer.OpenOutputFileForAppend();

        element_writer.Visit(&cell_population);

        element_writer.CloseFile();

        FileComparison(results_dir + "results.vizelements", "cell_based/test/data/TestCellPopulationElementWriter/results.vizelements_twice").CompareFiles();
    }

    void TestWriteCellMutationStates()    throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 5u);

        std::vector<unsigned> location_indices;
        location_indices.push_back(7);
        location_indices.push_back(11);
        location_indices.push_back(12);
        location_indices.push_back(13);
        location_indices.push_back(17);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.SetOutputCellMutationStates(true);

        std::string output_directory = "TestWriteCellMutationStates";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        CellMutationStatesWriter<2,2> mutation_states_writer(output_directory);

        mutation_states_writer.OpenOutputFile();

        mutation_states_writer.Visit(&cell_population);

        mutation_states_writer.CloseFile();

        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestWriteCellMutationStates/cellmutationstates.dat").CompareFiles();

        // Make sure we can append to files.
        mutation_states_writer.OpenOutputFileForAppend();

        mutation_states_writer.Visit(&cell_population);

        mutation_states_writer.CloseFile();

        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestWriteCellMutationStates/cellmutationstates_twice.dat").CompareFiles();
    }

    void TestWritePopulationCellProliferativeTypesCount()         throw (Exception)
	{
    	EXIT_IF_PARALLEL;

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);

        std::string output_directory = "TestWriteCellProlifertiveTypesCount";
		OutputFileHandler output_file_handler(output_directory, false);

		std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

		/*
		 * Write cell proliferative types count
		 */

		// Create a cell proliferative types writer and write the files.
		CellProliferativeTypesCountWriter<3,3> types_count_writer(output_directory);

		types_count_writer.OpenOutputFile();

		types_count_writer.Visit(&cell_population);

		types_count_writer.CloseFile();

		FileComparison(results_dir + "celltypes.dat", "cell_based/test/data/TestWriteCellProlifertiveTypesCount/celltypes.dat").CompareFiles();

		// Make sure we can append to files.
		types_count_writer.OpenOutputFileForAppend();

		types_count_writer.Visit(&cell_population);

		types_count_writer.CloseFile();

		FileComparison(results_dir + "celltypes.dat", "cell_based/test/data/TestWriteCellProlifertiveTypesCount/celltypes_twice.dat").CompareFiles();

		cell_population.CreateOutputFiles(output_directory, true);

		for (AbstractCellPopulation<3,3>::Iterator cell_iter = cell_population.Begin();
						cell_iter != cell_population.End();
						++cell_iter)
		{
				cell_population.GenerateCellResults(*cell_iter);
		}

		/*
		 * Write cell cycle phases count file.
		 */
		// Create a node location writer and write the files.
		CellProliferativePhasesCountWriter<3,3> phases_count_writer(output_directory);

		phases_count_writer.OpenOutputFile();

		phases_count_writer.Visit(&cell_population);

		phases_count_writer.CloseFile();

		FileComparison(results_dir + "cellcyclephases.dat", "cell_based/test/data/TestWriteCellProlifertiveTypesCount/cellcyclephases.dat").CompareFiles();

		// Make sure we can append to files.
		phases_count_writer.OpenOutputFileForAppend();

		phases_count_writer.Visit(&cell_population);

		phases_count_writer.CloseFile();

		FileComparison(results_dir + "cellcyclephases.dat", "cell_based/test/data/TestWriteCellProlifertiveTypesCount/cellcyclephases_twice.dat").CompareFiles();
	}

    void TestWriteVoronoiData() throw (Exception)
	{
        EXIT_IF_PARALLEL;

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true,  0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true,  1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(2, true,  1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(3, true,  0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(4, false, 0.5, 0.5, 0.5));
        MutableMesh<3,3> mesh(nodes);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create cell population
        MeshBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.CreateVoronoiTessellation();

        std::string output_directory = "TestWriteVoronoiFile";
		OutputFileHandler output_file_handler(output_directory, false);

		std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

		VoronoiDataWriter<3,3> voronoi_writer(output_directory);

		voronoi_writer.OpenOutputFile();

		voronoi_writer.Visit(&cell_population);

		voronoi_writer.CloseFile();

		FileComparison(results_dir + "voronoi.dat", "cell_based/test/data/TestWriteVoronoiFile/voronoi.dat").CompareFiles();

		// Make sure we can append to files.
		voronoi_writer.OpenOutputFileForAppend();

		voronoi_writer.Visit(&cell_population);

		voronoi_writer.CloseFile();

		FileComparison(results_dir + "voronoi.dat", "cell_based/test/data/TestWriteVoronoiFile/voronoi_twice.dat").CompareFiles();
	}

    void TestVertexSwapsWriter()    throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple vertex-based mesh
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        std::string output_directory = "TestWriteVertexSwaps";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        /**
         * T3 Swaps
         */
        VertexT1SwapLocationsWriter<2,2> t1_swaps_writer(output_directory);

        t1_swaps_writer.OpenOutputFile();

        t1_swaps_writer.Visit(&cell_population);

        t1_swaps_writer.CloseFile();

        FileComparison(results_dir + "T1SwapLocations.dat", "cell_based/test/data/TestWriteVertexSwaps/T1SwapLocations.dat").CompareFiles();

        // Make sure we can append to files.
        t1_swaps_writer.OpenOutputFileForAppend();

        t1_swaps_writer.Visit(&cell_population);

        t1_swaps_writer.CloseFile();

        FileComparison(results_dir + "T1SwapLocations.dat", "cell_based/test/data/TestWriteVertexSwaps/T1SwapLocations_twice.dat").CompareFiles();

        /**
         * T3 Swaps
         */
        VertexT3SwapLocationsWriter<2,2> t3_swaps_writer(output_directory);

        t3_swaps_writer.OpenOutputFile();

        t3_swaps_writer.Visit(&cell_population);

        t3_swaps_writer.CloseFile();

        FileComparison(results_dir + "T3SwapLocations.dat", "cell_based/test/data/TestWriteVertexSwaps/T3SwapLocations.dat").CompareFiles();

        // Make sure we can append to files.
        t3_swaps_writer.OpenOutputFileForAppend();

        t3_swaps_writer.Visit(&cell_population);

        t3_swaps_writer.CloseFile();

        FileComparison(results_dir + "T3SwapLocations.dat", "cell_based/test/data/TestWriteVertexSwaps/T3SwapLocations_twice.dat").CompareFiles();
    }

    void TestAddWritersToAPopulation() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // NodeBasedCellPopulation
        std::vector<Node<3>* > nodes;
        nodes.push_back(new Node<3>(0, false));
        nodes.push_back(new Node<3>(1, false, 1.0, 1.0, 1.0));

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> generator;
        generator.GenerateBasic(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        CellPopulationElementWriter<3,3> writer("output");
        cell_population.AddPopulationWriter(&writer);

        CellIdWriter<3,3> id_writer("output_directory");
        cell_population.AddCellWriter(&id_writer);
    }
};

#endif /*TESTCELLPOPULATIONWRITERS_HPP_*/
