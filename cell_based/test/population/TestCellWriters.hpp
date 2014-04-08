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

#ifndef TESTCELLWRITERS_HPP_
#define TESTCELLWRITERS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "ArchiveOpener.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "Cell.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "CellsGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellAncestor.hpp"
#include "SimulationTime.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLocationWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVariablesWriter.hpp"
#include "CellVolumesWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCellWriters : public AbstractCellBasedTestSuite
{
public:

    void TestCellAgesWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (this is usually done by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        // Create a simple node-based cell population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false,  1.4));
        nodes.push_back(new Node<2>(1, false,  2.3));
        nodes.push_back(new Node<2>(2, false, -6.1));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<3; i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create output directory
        std::string output_directory = "TestCellAgesWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellAgesWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_directory);
        cell_writer.WriteTimeStamp();
		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			 cell_iter != cell_population.End();
			 ++cell_iter)
		{
		    cell_writer.VisitCell(*cell_iter, &cell_population);
		}
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "cellages.dat", "cell_based/test/data/TestCellWriters/cellages.dat").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 0.7, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Ages");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellAgesWriterArchiving() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellAgesWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellAgesWriter<2,2>();

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

    void TestCellAncestorWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (this is usually done by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        // Create a simple node-based cell population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false,  1.4));
        nodes.push_back(new Node<2>(1, false,  2.3));
        nodes.push_back(new Node<2>(2, false, -6.1));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        std::vector<CellPtr> cells;
        for (unsigned i=0; i<3; i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Initialise each cell's ancestor
        cell_population.SetCellAncestorsToLocationIndices();

        // Create output directory
        std::string output_directory = "TestCellAncestorWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellAncestorWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_directory);
        cell_writer.WriteTimeStamp();
		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			 cell_iter != cell_population.End();
			 ++cell_iter)
		{
		    cell_writer.VisitCell(*cell_iter, &cell_population);
		}
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "results.vizancestors", "cell_based/test/data/TestCellWriters/results.vizancestors").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Ancestors");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellAncestorWriterArchiving() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellAncestorWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellAncestorWriter<2,2>();

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

    void TestCellIdWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (this is usually done by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        // Create a simple node-based cell population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false,  1.4));
        nodes.push_back(new Node<2>(1, false,  2.3));
        nodes.push_back(new Node<2>(2, false, -6.1));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        std::vector<CellPtr> cells;
        for (unsigned i=0; i<3; i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Initialise each cell's ancestor
        cell_population.SetCellAncestorsToLocationIndices();

        // Create output directory
        std::string output_directory = "TestCellIdWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellIdWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_directory);
        cell_writer.WriteTimeStamp();
		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			 cell_iter != cell_population.End();
			 ++cell_iter)
		{
		    cell_writer.VisitCell(*cell_iter, &cell_population);
		}
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "loggedcell.dat", "cell_based/test/data/TestCellWriters/loggedcell.dat").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        // (there have been six cells created since the start of this test suite)
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 6.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell IDs");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellIdWriterArchiving() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellIdWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellIdWriter<2,2>();

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

    void TestCellLocationWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple CA-based cell population
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

        // Create output directory
        std::string output_directory = "TestCellLocationWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellLocationWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_directory);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2,2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }

        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "results.vizlocations", "cell_based/test/data/TestCellWriters/results.vizlocations").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 0.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell locations");
    }

    void TestCellLocationWriterArchiving() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellLocationWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellLocationWriter<2,2>();

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

    void TestCellProliferativePhasesWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (this is usually done by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        // Create a simple node-based cell population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false,  1.4));
        nodes.push_back(new Node<2>(1, false,  2.3));
        nodes.push_back(new Node<2>(2, false, -6.1));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<3; i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create output directory
        std::string output_directory = "TestCellProliferativePhasesWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellProliferativePhasesWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_directory);
        cell_writer.WriteTimeStamp();
		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			 cell_iter != cell_population.End();
			 ++cell_iter)
		{
		    cell_writer.VisitCell(*cell_iter, &cell_population);
		}
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "results.vizcellphases", "cell_based/test/data/TestCellWriters/results.vizcellphases").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 4.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cycle phases");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellProliferativePhasesWriterArchiving() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellProliferativePhasesWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellProliferativePhasesWriter<2,2>();

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

    void TestCellProliferativeTypesWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (this is usually done by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        // Create a simple node-based cell population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false,  1.4));
        nodes.push_back(new Node<2>(1, false,  2.3));
        nodes.push_back(new Node<2>(2, false, -6.1));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_mutant_state(CellPropertyRegistry::Instance()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());

        std::vector<CellPtr> cells;
        for (unsigned i=0; i<2; i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }
        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_cell(new Cell(p_mutant_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->SetBirthTime(-0.1);
        cells.push_back(p_cell);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // For coverage of GetCellDataForVtkOutput() label a cell and set a cell to be apoptotic
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        cell_iter->AddCellProperty(p_label);
        ++cell_iter;
        cell_iter->AddCellProperty(p_apoptotic_state);

        // Create output directory
        std::string output_directory = "TestCellProliferativeTypesWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellProliferativeTypesWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_directory);
        cell_writer.WriteTimeStamp();
		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			 cell_iter != cell_population.End();
			 ++cell_iter)
		{
		    cell_writer.VisitCell(*cell_iter, &cell_population);
		}
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "results.vizcelltypes", "cell_based/test/data/TestCellWriters/results.vizcelltypes").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 5.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell types");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellProliferativeTypesWriterArchiving() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellProliferativeTypesWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellProliferativeTypesWriter<2,2>();

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

    void TestCellVariablesWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (this is usually done by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        // Create a simple node-based cell population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false,  1.4));
        nodes.push_back(new Node<2>(1, false,  2.3));
        nodes.push_back(new Node<2>(2, false, -6.1));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver> >
            p_solver(CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver>::Instance());
        p_solver->SetSizeOfOdeSystem(6);
        p_solver->Initialise();

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<3; i++)
        {
            TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel(p_solver);
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create output directory
        std::string output_directory = "TestCellVariablesWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellVariablesWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_directory);
        cell_writer.WriteTimeStamp();
		for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
			 cell_iter != cell_population.End();
			 ++cell_iter)
		{
		    cell_writer.VisitCell(*cell_iter, &cell_population);
		}
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "cellvariables.dat", "cell_based/test/data/TestCellWriters/cellvariables.dat").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 0.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell variables");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellVariablesWriterArchiving() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellVariablesWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellVariablesWriter<2,2>();

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

    void TestCellVolumesWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple vertex-based cell population
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create output directory
        std::string output_directory = "TestCellVolumesWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellVolumesWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_directory);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2,2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestCellWriters/cellareas.dat").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 0.8660254, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell volumes");
    }

    void TestCellVolumesWriterArchiving() throw (Exception)
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellVolumesWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellVolumesWriter<2,2>();

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

#endif /*TESTCELLWRITERS_HPP_*/
