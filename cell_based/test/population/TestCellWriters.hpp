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
#include "AbstractCellBasedTestSuite.hpp"
#include "FileComparison.hpp"
#include "CellWriters.hpp"

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


#include "PetscSetupAndFinalize.hpp"

class TestCellWriters : public AbstractCellBasedTestSuite
{
public:

    void TestCellAncestorWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();
        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (UNSIGNED_UNSET));
        p_cell->SetAncestor(p_cell_ancestor);

        // Create a simple population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        cells.push_back(p_cell);

        NodeBasedCellPopulation<2>* p_cell_population = new NodeBasedCellPopulation<2>(mesh, cells);

        std::string output_directory = "TestCellAncestorWriter";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        CellAncestorWriter<2,2> ancestor_writer(output_directory);

        ancestor_writer.OpenOutputFile();

        ancestor_writer.WriteTimeStamp();

        ancestor_writer.VisitCell(p_cell, p_cell_population);

        ancestor_writer.CloseFile();

        FileComparison(results_dir + "results.vizancestors", "cell_based/test/data/TestCellAncestorWriter/results.vizancestors").CompareFiles();

        ancestor_writer.OpenOutputFileForAppend();

        ancestor_writer.CloseFile();

        delete p_cell_population;
        delete nodes[0];

    }

    void TestWriteCellProliferativeTypesAndPhases() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);

        // Create a simple population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        cells.push_back(p_cell);

        NodeBasedCellPopulation<2>* p_cell_population = new NodeBasedCellPopulation<2>(mesh, cells);

        std::string output_directory = "TestCellProliferativeTypesAndPhasesWriter";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        CellProliferativeTypesWriter<2,2> types_writer(output_directory);

        types_writer.OpenOutputFile();

        types_writer.WriteTimeStamp();

        types_writer.VisitCell(p_cell, p_cell_population);

        types_writer.CloseFile();

        FileComparison(results_dir + "results.vizcelltypes", "cell_based/test/data/TestCellProliferativeTypesWriter/results.vizcelltypes").CompareFiles();

        CellProliferativePhasesWriter<2,2> phases_writer(output_directory);

        phases_writer.OpenOutputFile();

        phases_writer.WriteTimeStamp();

        phases_writer.VisitCell(p_cell, p_cell_population);

        phases_writer.CloseFile();

        FileComparison(results_dir + "results.vizcellphases", "cell_based/test/data/TestCellProliferativeTypesWriter/results.vizcellphases").CompareFiles();

        delete p_cell_population;
        delete nodes[0];
    }

    void TestWriteCellVariables() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        boost::shared_ptr<CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver> >
            p_solver(CellCycleModelOdeSolver<TysonNovakCellCycleModel, BackwardEulerIvpOdeSolver>::Instance());
        p_solver->SetSizeOfOdeSystem(6);
        p_solver->Initialise();

        TysonNovakCellCycleModel* p_cell_model = new TysonNovakCellCycleModel(p_solver);
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->InitialiseCellCycleModel();

        // Create a simple population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        cells.push_back(p_cell);

        NodeBasedCellPopulation<2>* p_cell_population = new NodeBasedCellPopulation<2>(mesh, cells);

        std::string output_directory = "TestCellVariablesWriter";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        CellVariablesWriter<2,2> variables_writer(output_directory);

        variables_writer.OpenOutputFile();

        variables_writer.WriteTimeStamp();

        variables_writer.VisitCell(p_cell, p_cell_population);

        variables_writer.CloseFile();

        FileComparison(results_dir + "cellvariables.dat", "cell_based/test/data/TestCellVariablesWriter/cellvariables.dat").CompareFiles();

        delete p_cell_population;
        delete nodes[0];
    }

    void TestWriteCellAges() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);

        // Create a simple population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        cells.push_back(p_cell);

        NodeBasedCellPopulation<2>* p_cell_population = new NodeBasedCellPopulation<2>(mesh, cells);

        std::string output_directory = "TestCellAgesWriter";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        CellAgesWriter<2,2> ages_writer(output_directory);

        ages_writer.OpenOutputFile();

        ages_writer.WriteTimeStamp();

        ages_writer.VisitCell(p_cell, p_cell_population);

        ages_writer.CloseFile();

        FileComparison(results_dir + "cellages.dat", "cell_based/test/data/TestCellAgesWriter/cellages.dat").CompareFiles();

        delete p_cell_population;
        delete nodes[0];
    }

    void TestCellIdWriter() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(25, 2);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

        FixedDurationGenerationBasedCellCycleModel* p_cell_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_cell(new Cell(p_healthy_state, p_cell_model));
        p_cell->SetCellProliferativeType(p_type);

        // Create a simple population
        std::vector<Node<2>* > nodes;
        nodes.push_back(new Node<2>(0, false, 0.0));

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        cells.push_back(p_cell);

        NodeBasedCellPopulation<2>* p_cell_population = new NodeBasedCellPopulation<2>(mesh, cells);

        std::string output_directory = "TestCellIdWriter";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        CellIdWriter<2,2> id_writer(output_directory);

        id_writer.OpenOutputFile();

        id_writer.WriteTimeStamp();

        id_writer.VisitCell(p_cell, p_cell_population);

        id_writer.CloseFile();

        FileComparison(results_dir + "loggedcell.dat", "cell_based/test/data/TestCellIdWriter/loggedcell.dat").CompareFiles();

        delete nodes[0];
        delete p_cell_population;
    }

    void TestWriteCellVolumesToFile() throw (Exception)
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

        std::string output_directory = "TestCellVolumesWriter";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        CellVolumesWriter<2,2> volumes_writer(output_directory);

        volumes_writer.OpenOutputFile();

        volumes_writer.WriteTimeStamp();

        for (AbstractCellPopulation<2,2>::Iterator cell_iter = cell_population.Begin();
                        cell_iter != cell_population.End();
                        ++cell_iter)
        {
                volumes_writer.VisitCell(*cell_iter, &cell_population);
        }

        volumes_writer.CloseFile();

        FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestCellVolumesWriter/cellareas.dat").CompareFiles();
    }

    void TestWriteLocations() throw (Exception)
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

        std::string output_directory = "TestCellLocationsWriter";
        OutputFileHandler output_file_handler(output_directory, false);

        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        CellLocationWriter<2,2> location_writer(output_directory);

        location_writer.OpenOutputFile();

        location_writer.WriteTimeStamp();

        for (AbstractCellPopulation<2,2>::Iterator cell_iter = cell_population.Begin();
                        cell_iter != cell_population.End();
                        ++cell_iter)
        {
                location_writer.VisitCell(*cell_iter, &cell_population);
        }

        location_writer.CloseFile();

        FileComparison(results_dir + "results.vizlocations", "cell_based/test/data/TestCellLocationsWriter/results.vizlocations").CompareFiles();
    }
};

#endif /*TESTCELLWRITERS_HPP_*/
