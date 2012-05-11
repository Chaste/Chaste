/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef TESTCELLWISEDATA_HPP_
#define TESTCELLWISEDATA_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "CellsGenerator.hpp"
#include "CellwiseData.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "ArchiveOpener.hpp"
#include "ArchiveLocationInfo.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"

/**
 * This class contains tests for methods on the class CellwiseData.
 */
class TestCellwiseData : public AbstractCellBasedTestSuite
{
public:

    void TestCellwiseDataSimple() throw(Exception)
    {
        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index,
        // so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);

        // One variable tests

        CellwiseData<2>* p_data = CellwiseData<2>::Instance();

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);

        p_data->SetPopulationAndNumVars(&cell_population, 1);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->GetNumVariables(), 1u);

        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        cell_iter->GetCellData()->SetItem(0, 1.23);
        TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem(0), 1.23, 1e-12);

        ++cell_iter;
        cell_iter->GetCellData()->SetItem(0, 2.23);
        TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem(0), 2.23, 1e-12);

        // Test that CellData objects are copied correctly on cell division.
        MAKE_PTR(WildTypeCellMutationState, p_state);

        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);

        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetBirthTime(-1);

        c_vector<double,2> new_cell_location;
        new_cell_location[0] = 0.2;
        new_cell_location[1] = 0.3;
        cell_population.AddCell(p_new_cell, new_cell_location, cells[0] /*random choice of parent*/);

        p_data->Destroy();

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), false);

        // Two variable tests

        p_data = CellwiseData<2>::Instance();

        p_data->SetPopulationAndNumVars(&cell_population, 2);

        TS_ASSERT_THROWS_THIS(p_data->SetPopulationAndNumVars(&cell_population, 1),"Can't call SetPopulationAndNumVars() once CellwiseData is setup.");

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);

        TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->GetNumVariables(), 2u);

        AbstractCellPopulation<2>::Iterator cell_iter2 = cell_population.Begin();
        cell_iter2->GetCellData()->SetItem(1, 3.23);       
        TS_ASSERT_DELTA(cell_iter2->GetCellData()->GetItem(1), 3.23, 1e-12);

        TS_ASSERT_THROWS_THIS(cell_iter2->GetCellData()->GetItem(2),
                    "Request for variable above the number of variables stored.");



        ++cell_iter2;
        cell_iter2->GetCellData()->SetItem(1, 4.23);       

        TS_ASSERT_DELTA(cell_iter2->GetCellData()->GetItem(1), 4.23, 1e-12);

        // Other values should have been initialised to zero
        ++cell_iter2;
        TS_ASSERT_THROWS_THIS(cell_iter2->GetCellData()->GetItem(0),"SetItem must be called before using GetItem");

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestCellwiseDataExceptions() throw(Exception)
    {
        // Create a simple mesh with 2 layers of ghost nodes
        HoneycombMeshGenerator generator(2, 2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Setup some cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), TRANSIT);

        // Create a cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        TS_ASSERT_THROWS_THIS(p_data->SetPopulationAndNumVars(&cell_population, 1),
                "CellwiseData does not work with ghost nodes.");

        // Tidy up
        CellwiseData<2>::Destroy();
    }


    void TestArchiveCellwiseData()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index, so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Work out where to put the archive
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "cellwise_data.arch";
        ArchiveLocationInfo::SetMeshFilename("cellwise_data_mesh");

        {
            // Set up the data store
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();
            p_data->SetPopulationAndNumVars(&cell_population, 1);

            // Put some data in
            unsigned i = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem(0, (double) i);
                i++;
            }

            TS_ASSERT_EQUALS(p_data->IsSetUp(), true);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write to the archive
            SerializableSingleton<CellwiseData<2> >* const p_wrapper = p_data->GetSerializationWrapper();
            (*p_arch) << p_wrapper;

            CellwiseData<2>::Destroy();
        }

        {
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            SerializableSingleton<CellwiseData<2> >* p_wrapper;
            (*p_arch) >> p_wrapper;

            // Check the data
            TS_ASSERT_EQUALS(CellwiseData<2>::Instance(), p_data);
            TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);
            TS_ASSERT_EQUALS(p_data->IsSetUp(), true);
            TS_ASSERT_EQUALS(p_data->mUseConstantDataForTesting, false);
            TS_ASSERT_EQUALS(p_data->GetNumVariables(), 1u);

            // We will have constructed a new cell population on load, so use the new cell population
            AbstractCellPopulation<2>& cell_population = p_data->rGetCellPopulation();

            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem(0), (double) cell_population.GetLocationIndexUsingCell(*cell_iter), 1e-12);
            }

            // Tidy up
            CellwiseData<2>::Destroy();
            delete (&cell_population);
        }
    }

    void TestArchiveCellwiseDataWithVertexBasedCellPopulation()
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "vertex_cellwise.arch";
        // The following line is required because the loading of a cell population
        // is usually called by the method CellBasedSimulation::Load()
        ArchiveLocationInfo::SetMeshFilename("vertex_cellwise");

        // Create mesh
        VertexMeshReader<2,2> mesh_reader("mesh/test/data/TestVertexMeshWriter/vertex_mesh_2d");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Need to set up time
        unsigned num_steps = 10;
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<2>* const p_cell_population = new VertexBasedCellPopulation<2>(mesh, cells);

        // Cells have been given birth times of 0 and -1.
        // Loop over them to run to time 0.0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
             cell_iter != p_cell_population->End();
             ++cell_iter)
        {
            cell_iter->ReadyToDivide();
        }

        {
            // Set up the data store
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();
            p_data->SetPopulationAndNumVars(p_cell_population, 1);

            // Put some data in
            unsigned i = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem(0, (double) i);
                i++;
            }

            TS_ASSERT_EQUALS(p_data->IsSetUp(), true);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write to the archive
            SerializableSingleton<CellwiseData<2> >* const p_wrapper = p_data->GetSerializationWrapper();
            (*p_arch) << p_wrapper;

            CellwiseData<2>::Destroy();
            delete p_cell_population;
        }

        {
            CellwiseData<2>* p_data = CellwiseData<2>::Instance();

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            SerializableSingleton<CellwiseData<2> >* p_wrapper;
            (*p_arch) >> p_wrapper;

            // Check the data
            TS_ASSERT_EQUALS(CellwiseData<2>::Instance(), p_data);
            TS_ASSERT_EQUALS(CellwiseData<2>::Instance()->IsSetUp(), true);
            TS_ASSERT_EQUALS(p_data->IsSetUp(), true);

            // We will have constructed a new cell population on load, so use the new cell population
            AbstractCellPopulation<2>& cell_population = p_data->rGetCellPopulation();

            for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem(0), (double) cell_population.GetLocationIndexUsingCell(*cell_iter), 1e-12);
            }

            // Tidy up
            CellwiseData<2>::Destroy();
            delete (&cell_population);
        }
    }
    
    ///\todo #1515 This test doesn't use CellwiseData and should move into another test suite when CellwiseData is redundant
    void TestAddCellData()
    {
        // Set up simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Set up cells, one for each node. Give each a birth time of -node_index, so the age = node_index
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(mesh, cells);
        
        MAKE_PTR_ARGS(CellData, p_cell_data, (1)); 
        p_cell_data->SetItem(0, 100.0);
        cell_population.AddClonedDataToAllCells(p_cell_data);
        
        //Check that the data made it there and that copies of the data are independent
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                 cell_iter != cell_population.End();
                 ++cell_iter)
        {   
            TS_ASSERT_EQUALS(cell_iter->GetCellData()->GetItem(0), 100.0);
            cell_iter->GetCellData()->SetItem(0, 1.0);
        }
        
        //Try it again
        TS_ASSERT_THROWS_THIS(cell_population.AddClonedDataToAllCells(p_cell_data),
            "AddClonedDataToAllCells() assumes that cells have no data");
    }
};

#endif /*TESTCELLWISEDATA_HPP_*/
