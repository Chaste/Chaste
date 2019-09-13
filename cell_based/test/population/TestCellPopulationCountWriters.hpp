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

#ifndef TESTCELLPOPULATIONCOUNTWRITERS_HPP_
#define TESTCELLPOPULATIONCOUNTWRITERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"

// Cell population count writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

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

class TestCellPopulationCountWriters : public AbstractCellBasedTestSuite
{
public:

    void TestCellMutationStatesCountWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D CaBasedCellPopulation
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 5u);
        std::vector<unsigned> location_indices;
        location_indices.push_back(7);
        location_indices.push_back(11);
        location_indices.push_back(12);
        location_indices.push_back(13);
        location_indices.push_back(17);
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        cell_population.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        // Create an output directory for the writer
        std::string output_directory = "TestCellMutationStatesCountWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellMutationStatesCountWriter and test that the correct output is generated
        CellMutationStatesCountWriter<2,2> mutation_states_writer;
        mutation_states_writer.OpenOutputFile(output_file_handler);
        mutation_states_writer.WriteHeader(&cell_population);
        mutation_states_writer.WriteTimeStamp();
        mutation_states_writer.Visit(&cell_population);
        mutation_states_writer.WriteNewline();
        mutation_states_writer.CloseFile();

        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestCellPopulationCountWriters/cellmutationstates.dat").CompareFiles();

        // Test that we can append to files
        mutation_states_writer.OpenOutputFileForAppend(output_file_handler);
        mutation_states_writer.WriteTimeStamp();
        mutation_states_writer.Visit(&cell_population);
        mutation_states_writer.WriteNewline();
        mutation_states_writer.CloseFile();

        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestCellPopulationCountWriters/cellmutationstates_twice.dat").CompareFiles();
    }

    void TestCellMutationStatesCountWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellMutationStatesCountWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new CellMutationStatesCountWriter<2,2>();
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

    void TestCellProliferativePhasesCountWriters()
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

        cell_population.AddCellPopulationCountWriter<CellProliferativePhasesCountWriter>();

        // Create an output directory for the writer
        std::string output_directory = "TestCellProliferativePhasesCountWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellProliferativePhasesCountWriter and test that the correct output is generated
        CellProliferativePhasesCountWriter<3,3> phases_count_writer;
        phases_count_writer.OpenOutputFile(output_file_handler);
        phases_count_writer.WriteTimeStamp();
        phases_count_writer.Visit(&cell_population);
        phases_count_writer.WriteNewline();
        phases_count_writer.CloseFile();

        FileComparison(results_dir + "cellcyclephases.dat", "cell_based/test/data/TestCellPopulationCountWriters/cellcyclephases.dat").CompareFiles();

        // Test that we can append to files
        phases_count_writer.OpenOutputFileForAppend(output_file_handler);
        phases_count_writer.WriteTimeStamp();
        phases_count_writer.Visit(&cell_population);
        phases_count_writer.WriteNewline();
        phases_count_writer.CloseFile();

        FileComparison(results_dir + "cellcyclephases.dat", "cell_based/test/data/TestCellPopulationCountWriters/cellcyclephases_twice.dat").CompareFiles();
    }

    void TestCellProliferativePhasesCountWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellProliferativePhasesCountWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new CellProliferativePhasesCountWriter<2,2>();
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

    void TestCellProliferativeTypesCountWriters()
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

        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // An ordering must be specified for cell mutation states and cell proliferative types
        cell_population.SetDefaultCellMutationStateAndProliferativeTypeOrdering();

        // Create an output directory for the writer
        std::string output_directory = "TestCellProliferativeTypesCountWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create a CellProliferativeTypesCountWriter and test that the correct output is generated
        CellProliferativeTypesCountWriter<3,3> types_count_writer;
        types_count_writer.OpenOutputFile(output_file_handler);
        types_count_writer.WriteTimeStamp();
        types_count_writer.Visit(&cell_population);
        types_count_writer.WriteNewline();
        types_count_writer.CloseFile();

        FileComparison(results_dir + "celltypes.dat", "cell_based/test/data/TestCellPopulationCountWriters/celltypes.dat").CompareFiles();

        // Test that we can append to files
        types_count_writer.OpenOutputFileForAppend(output_file_handler);
        types_count_writer.WriteTimeStamp();
        types_count_writer.Visit(&cell_population);
        types_count_writer.WriteNewline();
        types_count_writer.CloseFile();

        FileComparison(results_dir + "celltypes.dat", "cell_based/test/data/TestCellPopulationCountWriters/celltypes_twice.dat").CompareFiles();
    }

    void TestCellProliferativeTypesCountWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellProliferativeTypesCountWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_population_writer = new CellProliferativeTypesCountWriter<2,2>();
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

#endif /*TESTCELLPOPULATIONCOUNTWRITERS_HPP_*/
