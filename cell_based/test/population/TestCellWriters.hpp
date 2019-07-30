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

#ifndef TESTCELLWRITERS_HPP_
#define TESTCELLWRITERS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "ArchiveOpener.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "ApoptoticCellProperty.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CaBasedCellPopulation.hpp"
#include "Cell.hpp"
#include "CellAncestor.hpp"
#include "CellLabel.hpp"
#include "CellsGenerator.hpp"
#include "CellsGenerator.hpp"
#include "DeltaNotchSrnModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "SimulationTime.hpp"
#include "SmartPointers.hpp"
#include "StemCellProliferativeType.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "UblasCustomFunctions.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "WildTypeCellMutationState.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellAppliedForceWriter.hpp"
#include "CellCycleModelProteinConcentrationsWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "CellDeltaNotchWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellRadiusWriter.hpp"
#include "CellRosetteRankWriter.hpp"
#include "CellVolumesWriter.hpp"

// Boost
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "PetscSetupAndFinalize.hpp"

// Note that high level tests of all cell writers can be found in the
// TestMeshBasedCellPopulation::TestMeshBasedCellPopulationWriteResultsToFile

class TestCellWriters : public AbstractCellBasedTestSuite
{
public:

    void TestCellAgesWriter()
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
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
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

        // Create cell writer
        CellAgesWriter<2,2> cell_writer;

        // Test get and set methods
        TS_ASSERT_EQUALS(cell_writer.GetFileName(), "cellages.dat");
        cell_writer.SetFileName("new_name.txt");
        TS_ASSERT_EQUALS(cell_writer.GetFileName(), "new_name.txt");
        cell_writer.SetFileName("cellages.dat");

        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Ages");
        cell_writer.SetVtkCellDataName("Names");
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Names");
        cell_writer.SetVtkCellDataName("Ages");

        // Output data for each cell to file
        cell_writer.OpenOutputFile(output_file_handler);
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

    void TestCellAgesWriterArchiving()
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

    void TestCellAncestorWriter()
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
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
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
        cell_writer.OpenOutputFile(output_file_handler);
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

    void TestCellAncestorWriterArchiving()
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
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;

            delete p_cell_writer_2;
        }
    }

    void TestCellDeltaNotchWriter()
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (this is usually done by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        // Create a regular vertex mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create some cells, each with a cell-cycle model that incorporates a delta-notch ODE system
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();

            UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell-based population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.SetDataOnAllCells("delta", 1.56);
        cell_population.SetDataOnAllCells("notch", 9.54);
        cell_population.SetDataOnAllCells("mean delta", 87.3);

        // Create output directory
        std::string output_directory = "TestCellDeltaNotchWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellDeltaNotchWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "celldeltanotch.dat", "cell_based/test/data/TestCellWriters/celldeltanotch.dat").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 1.56, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell delta");
    }

    void TestCellDeltaNotchWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellDeltaNotchWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellDeltaNotchWriter<2,2>();

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

    void TestCellDataItemWriter()
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (this is usually done by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        // Create a regular vertex mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create some cells
        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        std::vector<CellPtr> cells;
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }

        // Create cell-based population object
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDataOnAllCells("test_variable", 1.56);

        // Create output directory
        std::string output_directory = "TestCellDataItemWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellDataItemWriter<2,2> cell_writer("test_variable");
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "celldata_test_variable.dat", "cell_based/test/data/TestCellWriters/celldataitem.dat").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 1.56, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "CellData test_variable");
    }

    void TestCellDataItemWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellDataItemWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellDataItemWriter<2,2>("test_variable");

            CellDataItemWriter<2,2>* p_static_cast_cell_writer = static_cast <CellDataItemWriter<2,2>* >(p_cell_writer);
            TS_ASSERT_EQUALS(p_static_cast_cell_writer->GetVtkCellDataName(), "CellData test_variable");

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

            // Check the member variables have been updated properly
            CellDataItemWriter<2,2>* p_static_cast_cell_writer = static_cast <CellDataItemWriter<2,2>* >(p_cell_writer_2);
            TS_ASSERT_EQUALS(p_static_cast_cell_writer->GetVtkCellDataName(), "CellData test_variable");

            delete p_cell_writer_2;
        }
    }

    void TestCellIdWriter()
    {

        EXIT_IF_PARALLEL;

        // Resetting the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

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
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
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
        cell_writer.OpenOutputFile(output_file_handler);
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
        TS_ASSERT_DELTA(vtk_data, 0.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell IDs");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellIdWriterArchiving()
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
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;

            delete p_cell_writer_2;
        }
    }

    void TestCellLabelWriter()
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
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Begin()->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellLabel>());

        // Create output directory
        std::string output_directory = "TestCellLabelWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellLabelWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "results.vizlabels", "cell_based/test/data/TestCellWriters/results.vizlabels").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 5.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell labels");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellLabelWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellLabelWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellLabelWriter<2,2>();

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

    void TestCellLocationIndexWriter()
    {
        EXIT_IF_PARALLEL;

        // Create a simple CA-based cell population
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

        // Create output directory
        std::string output_directory = "TestCellLocationIndexWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellLocationIndexWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2,2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }

        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "results.vizlocationindices", "cell_based/test/data/TestCellWriters/results.vizlocationindices").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 0.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Location indices");
    }

    void TestCellLocationIndexWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellLocationIndexWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellLocationIndexWriter<2,2>();

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

    void TestCellMutationStatesWriter()
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
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<2; i++)
        {
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }
        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        CellPtr p_cell(new Cell(p_mutant_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->SetBirthTime(-0.1);
        cells.push_back(p_cell);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create output directory
        std::string output_directory = "TestCellMutationStatesWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellMutationStatesWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "results.vizmutationstates", "cell_based/test/data/TestCellWriters/results.vizmutationstates").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 0.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Mutation states");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellMutationStatesWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellMutationStatesWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellMutationStatesWriter<2,2>();

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

    void TestCellProliferativePhasesWriter()
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
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
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
        cell_writer.OpenOutputFile(output_file_handler);
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

    void TestCellProliferativePhasesWriterArchiving()
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
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;

            delete p_cell_writer_2;
        }
    }

    void TestCellProliferativeTypesWriter()
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
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }
        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
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
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2>::Iterator other_cell_iter = cell_population.Begin();
             other_cell_iter != cell_population.End();
             ++other_cell_iter)
        {
            cell_writer.VisitCell(*other_cell_iter, &cell_population);
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

    void TestCellProliferativeTypesWriterArchiving()
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
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;

            delete p_cell_writer_2;
        }
    }

    void TestCellCycleModelProteinConcentrationsWriter()
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
        std::string output_directory = "TestCellCycleModelProteinConcentrationsWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellCycleModelProteinConcentrationsWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "proteinconcentrations.dat", "cell_based/test/data/TestCellWriters/proteinconcentrations.dat").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 0.0, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Protein concentrations");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellCycleModelProteinConcentrationsWriterException()
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
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->InitialiseCellCycleModel();
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create output directory
        std::string output_directory = "TestCellCycleModelProteinConcentrationsWriterException";
        OutputFileHandler output_file_handler(output_directory, false);

        // Create cell writer and output data for each cell to file
        CellCycleModelProteinConcentrationsWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();

        TS_ASSERT_THROWS_THIS(cell_writer.VisitCell(*cell_iter, &cell_population),
            "CellCycleModelProteinConcentrationsWriter cannot be used with a cell-cycle model that does not inherit from CellCycleModelOdeHandler");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellCycleModelProteinConcentrationsWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellCycleModelProteinConcentrationsWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellCycleModelProteinConcentrationsWriter<2,2>();

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

    void TestCellVolumesWriter()
    {
        EXIT_IF_PARALLEL;

        // Resetting the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple vertex-based cell population
        HoneycombVertexMeshGenerator generator(4, 6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_diff_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create output directory
        std::string output_directory = "TestCellVolumesWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellVolumesWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_file_handler);
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

    void TestCellVolumesWriterArchiving()
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
        PetscTools::Barrier(); //Processes read after last process has (over-)written archive
        {
            AbstractCellBasedWriter<2,2>* p_cell_writer_2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_cell_writer_2;

            delete p_cell_writer_2;
        }
    }

    void TestCellRosetteRankWriter()
    {
        EXIT_IF_PARALLEL;

        /*
         * Test regular functionality
         */
        {
            // Resetting the maximum cell ID to zero (to account for previous tests)
            CellId::ResetMaxCellId();

            // Create a simple vertex-based cell population
            HoneycombVertexMeshGenerator generator(4, 6);
            MutableVertexMesh<2, 2> *p_mesh = generator.GetMesh();

            std::vector<CellPtr> cells;
            boost::shared_ptr<AbstractCellProperty> p_diff_type(
                    CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
            CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_diff_type);

            VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

            // Create output directory
            std::string output_directory = "TestCellRosetteRankWriter";
            OutputFileHandler output_file_handler(output_directory, false);
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

            // Create cell writer and output data for each cell to file
            CellRosetteRankWriter<2, 2> cell_writer;
            cell_writer.OpenOutputFile(output_file_handler);
            cell_writer.WriteTimeStamp();
            for (AbstractCellPopulation<2, 2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                cell_writer.VisitCell(*cell_iter, &cell_population);
            }
            cell_writer.CloseFile();

            // Test that the data are output correctly
            FileComparison(results_dir + "cellrosetterank.dat",
                           "cell_based/test/data/TestCellWriters/cellrosetterank.dat").CompareFiles();

            // Test the correct data are returned for VTK output for the first cell
            double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
            TS_ASSERT_DELTA(vtk_data, 3.0, 1e-6);

            // Test GetVtkCellDataName() method
            TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell rosette rank");
        }

        /*
         * Test exception for non-vertex-based cell populations
         */
        {
            // Create a simple CA-based cell population
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

            // Create output directory
            std::string output_directory = "TestCellRosetteRankWriter";
            OutputFileHandler output_file_handler(output_directory, false);
            std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

            // Create cell writer and output data for each cell to file
            CellRosetteRankWriter<2, 2> cell_writer;
            cell_writer.OpenOutputFile(output_file_handler);
            cell_writer.WriteTimeStamp();
            for (AbstractCellPopulation<2, 2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                TS_ASSERT_THROWS_THIS(cell_writer.VisitCell(*cell_iter, &cell_population),
                                      "Rosettte rank is only associated with vertex-based cell populations");
            }
            cell_writer.CloseFile();
        }
    }

    void TestCellRosetteRankWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellRosetteRankWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellRosetteRankWriter<2,2>();

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

    void TestCellRadiusWriter()
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
            FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(-0.7 - i*0.5);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create output directory
        std::string output_directory = "TestCellRadiusWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellRadiusWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        for (AbstractCellPopulation<2,2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "cellradii.dat", "cell_based/test/data/TestCellWriters/cellradii.dat").CompareFiles();

        // Test the correct data are returned for VTK output for the first cell
        double vtk_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_DELTA(vtk_data, 0.5, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "Cell radii");

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestCellRadiusWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellRadiusWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellRadiusWriter<2,2>();

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

            delete p_cell_writer_2;
       }
    }

    void TestCellAppliedForceWriter()
    {
        EXIT_IF_PARALLEL;

        // Set up SimulationTime (this is usually done by a simulation object)
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        // Create a simple node-based cell population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0u));
        nodes.push_back(new Node<2>(1u));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Add an applied force to the nodes
        c_vector<double, 2> force_0 = Create_c_vector(1.23, 2.34);
        c_vector<double, 2> force_1 = Create_c_vector(3.45, 4.56);
        mesh.GetNode(0u)->AddAppliedForceContribution(force_0);
        mesh.GetNode(1u)->AddAppliedForceContribution(force_1);

        // Create output directory
        std::string output_directory = "TestCellAppliedForceWriter";
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Create cell writer and output data for each cell to file
        CellAppliedForceWriter<2,2> cell_writer;
        cell_writer.OpenOutputFile(output_file_handler);
        cell_writer.WriteTimeStamp();
        for (auto cell_iter = cell_population.Begin(); cell_iter != cell_population.End(); ++cell_iter)
        {
            cell_writer.VisitCell(*cell_iter, &cell_population);
        }
        cell_writer.CloseFile();

        // Test that the data are output correctly
        FileComparison(results_dir + "cellappliedforce.dat", "cell_based/test/data/TestCellWriters/cellappliedforce.dat").CompareFiles();

        // Test the correct data are returned for VTK vector output for the first cell
        c_vector<double, 2> vtk_data_0 = cell_writer.GetVectorCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        c_vector<double, 2> vtk_data_1 = cell_writer.GetVectorCellDataForVtkOutput(*(++cell_population.Begin()), &cell_population);

        TS_ASSERT_DELTA(vtk_data_0[0], 1.23, 1e-6);
        TS_ASSERT_DELTA(vtk_data_0[1], 2.34, 1e-6);
        TS_ASSERT_DELTA(vtk_data_1[0], 3.45, 1e-6);
        TS_ASSERT_DELTA(vtk_data_1[1], 4.56, 1e-6);

        // Test GetVtkCellDataName() method
        TS_ASSERT_EQUALS(cell_writer.GetVtkVectorCellDataName(), "Cell applied force");

        cell_writer.SetVtkVectorCellDataName("New name");
        TS_ASSERT_EQUALS(cell_writer.GetVtkVectorCellDataName(), "New name");

        // Avoid memory leak
        for (auto& p_node : nodes)
        {
            delete p_node;
        }
    }

    void TestCellAppliedForceWriterArchiving()
    {
        // The purpose of this test is to check that archiving can be done for this class
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "CellRadiusWriter.arch";

        {
            AbstractCellBasedWriter<2,2>* const p_cell_writer = new CellAppliedForceWriter<2,2>();

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

            delete p_cell_writer_2;
        }
    }

    void TestDefaultVecBehaviourWhenWritingScalars()
    {
        // We test here that a writer designed to only output scalar data has the correct default behaviour for
        // outputting vectors
        EXIT_IF_PARALLEL;

        CellRadiusWriter<2,2> cell_writer;

        TS_ASSERT_EQUALS(cell_writer.GetVtkVectorCellDataName(), "DefaultVtkVectorCellDataName");
        TS_ASSERT(cell_writer.GetOutputScalarData());
        TS_ASSERT(!cell_writer.GetOutputVectorData());

        // Create a simple node-based cell population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0u));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        c_vector<double, 2> vec_data = cell_writer.GetVectorCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        for(auto& component : vec_data)
        {
            TS_ASSERT_EQUALS(component, DOUBLE_UNSET);
        }

        // Avoid memory leak
        for (auto& p_node : nodes)
        {
            delete p_node;
        }
    }

    void TestDefaultScalarBehaviourWhenWritingVectors()
    {
        // We test here that a writer designed to only output vector data has the correct default behaviour for
        // outputting scalars
        EXIT_IF_PARALLEL;

        CellAppliedForceWriter<2,2> cell_writer;

        TS_ASSERT_EQUALS(cell_writer.GetVtkCellDataName(), "DefaultVtkCellDataName");
        TS_ASSERT(!cell_writer.GetOutputScalarData());
        TS_ASSERT(cell_writer.GetOutputVectorData());

        // Create a simple node-based cell population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0u));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        auto p_diff_type = boost::make_shared<DifferentiatedCellProliferativeType>();
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_diff_type);

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        double scalar_data = cell_writer.GetCellDataForVtkOutput(*(cell_population.Begin()), &cell_population);
        TS_ASSERT_EQUALS(scalar_data, DOUBLE_UNSET);

        // Avoid memory leak
        for (auto& p_node : nodes)
        {
            delete p_node;
        }
    }
};

#endif /*TESTCELLWRITERS_HPP_*/
