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

#ifndef TESTCELLDATA_HPP_
#define TESTCELLDATA_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "ArchiveOpener.hpp"


/**
 * This class contains tests for methods on the cell property CellData.
 */
class TestCellData : public AbstractCellBasedTestSuite
{
public:

    void TestCellDataSimple() throw(Exception)
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


        // One variable tests
        MAKE_PTR_ARGS(CellData, p_cell_data, (2));
        p_cell_data->SetItem(0, DOUBLE_UNSET);
        p_cell_data->SetItem(1, DOUBLE_UNSET);
        cell_population.AddClonedDataToAllCells(p_cell_data);

        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        cell_iter->GetCellData()->SetItem(0, 1.23);
        TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem(0), 1.23, 1e-12);

        ++cell_iter;
        cell_iter->GetCellData()->SetItem(0, 2.23);
        TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem(0), 2.23, 1e-12);


        cell_iter->GetCellData()->SetItem(1, 3.23);
        TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem(1), 3.23, 1e-12);

        TS_ASSERT_THROWS_THIS(cell_iter->GetCellData()->GetItem(2),
                    "Request for variable above the number of variables stored.");

        // Other values should have been initialised to zero
        ++cell_iter;
        TS_ASSERT_THROWS_THIS(cell_iter->GetCellData()->GetItem(0),"SetItem must be called before using GetItem");
    }


    void TestArchiveCellData()
    {
        // Work out where to put the archive
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "cellwise_data.arch";
        ArchiveLocationInfo::SetMeshFilename("cellwise_data_mesh");
        {
            // Create a simple mesh
            TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_4_elements");
            MutableMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);

            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumNodes());

            // Create a cell population
            MeshBasedCellPopulation<2>* const p_cell_population = new MeshBasedCellPopulation<2>(mesh, cells);

            MAKE_PTR_ARGS(CellData, p_cell_data, (1));
            p_cell_data->SetItem(0, DOUBLE_UNSET);
            p_cell_population->AddClonedDataToAllCells(p_cell_data);

            // Put some non-constant data in
            unsigned i = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->GetCellData()->SetItem(0, (double) i);
                i++;
            }

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << p_cell_population;


            delete p_cell_population;
        }

        {

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            MeshBasedCellPopulation<2>* p_cell_population;
            (*p_arch) >> p_cell_population;
            // Check the data


            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem(0), (double) p_cell_population->GetLocationIndexUsingCell(*cell_iter), 1e-12);
            }

            // Tidy up
            delete (p_cell_population);
        }
    }

    
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
