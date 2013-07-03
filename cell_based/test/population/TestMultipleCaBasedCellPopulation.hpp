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

#ifndef TESTMULTIPLECABASEDCELLPOPULATION_HPP_
#define TESTMULTIPLECABASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "PottsMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "DiffusionMultipleCaUpdateRule.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"
#include "FileComparison.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestMultipleCaBasedCellPopulation : public AbstractCellBasedTestSuite
{
public:

    void TestConstructor() throw(Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 1);

        std::vector<unsigned> location_indices;
        location_indices.push_back(12);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Test that the mesh and cells are correctly assigned
        TS_ASSERT_EQUALS(&(cell_population.rGetMesh()), p_mesh);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), p_mesh->GetNumNodes());
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), location_indices.size());

        ///\todo this doesnt do anything as there aren't any cells!
        TS_ASSERT_EQUALS(cells.size(), 0u);
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        for (unsigned i=0; i<cells.size(); i++)
        {
            TS_ASSERT_EQUALS(*cell_iter, cells[i]);
            TS_ASSERT_EQUALS(*cell_iter, cell_population.GetCellUsingLocationIndex(location_indices[i]));
            TS_ASSERT_EQUALS(location_indices[i], cell_population.GetLocationIndexUsingCell(*cell_iter)+1);
            ++cell_iter;
        }

        // Test that we do not have any update rules present
        TS_ASSERT_EQUALS(cell_population.rGetUpdateRuleCollection().empty(), true);

        // Test that the other member variables of this object are initialised correctly
        TS_ASSERT_EQUALS(cell_population.GetUpdateNodesInRandomOrder(), true);
        TS_ASSERT_EQUALS(cell_population.GetIterateRandomlyOverUpdateRuleCollection(), false);

        // For coverage of GetVolumeOfCell()
        for (AbstractCellPopulation<2>::Iterator cell_iter2 = cell_population.Begin();
             cell_iter2 != cell_population.End();
             ++cell_iter2)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter2), 1.0, 1e-6);
        }

        // For Ccverage of GetWidth() method (note that this returns the size of the underlying mesh, not the population of cells)
        TS_ASSERT_DELTA(cell_population.GetWidth(0), 4.0, 1e-12);
        TS_ASSERT_DELTA(cell_population.GetWidth(1), 4.0, 1e-12);

        // For coverage of GetNeighbouringNodeIndices() method
        TS_ASSERT_THROWS_THIS(cell_population.GetNeighbouringNodeIndices(0),
            "Cannot call GetNeighbouringNodeIndices() on a MultipleCaBasedCellPopulation, need to go through the PottsMesh instead");
    }

    void TestConstructorWithMultipleCellsPerSite() throw(Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(2, 0, 0, 1, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_diff_type);

        std::vector<unsigned> location_indices;

        // Test that the correct exception is thrown when we try to create a cell population without location indices
        TS_ASSERT_THROWS_THIS(MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 2),
            "No location indices being passed. Specify where cells lie before creating the cell population.");

        cells_generator.GenerateBasicRandom(cells, 3, p_diff_type);

        // Specify where cells lie
        location_indices.push_back(0);
        location_indices.push_back(0);
        location_indices.push_back(0);

        // Test that the correct exception is thrown when we try to add more cells than the carrying capacity
        TS_ASSERT_THROWS_THIS(MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 2),
            "One of the lattice sites has more cells than the carrying capacity. Check the initial cell locations.");

        cells_generator.GenerateBasicRandom(cells, 3, p_diff_type);

           // Change the initial cell location to avoid the above exception
        location_indices[2] = 1u;

        // Test that the correct exception is thrown when we try to create a cell population with validate = true (there is no validation in the case of a MultipleCaBasedCellPopulation)
        TS_ASSERT_THROWS_THIS(MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 2, false, true),
            "There is no validation for MultipleCaBasedCellPopulation.");

        cells_generator.GenerateBasicRandom(cells, 3, p_diff_type);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 2);

        // Test that cells are in the correct locations
        TS_ASSERT_EQUALS(cell_population.GetCellsUsingLocationIndex(0).size(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetCellsUsingLocationIndex(1).size(), 1u);

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(),3u);
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 0u);
        ++cell_iter;
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 0u);
        ++cell_iter;
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 1u);

        // Test that the mesh and cells are correctly assigned
        TS_ASSERT_EQUALS(&(cell_population.rGetMesh()), p_mesh);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), p_mesh->GetNumNodes());
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), location_indices.size());
    }

    /*
     * This test checks that cell populations with multiple cells per lattice site are dealt with correctly.
     */
    void TestMultipleCellExceptions() throw (Exception)
    {
        EXIT_IF_PARALLEL;   // Cell IDs are set differently in parallel

        // Reset the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple 2D PottsMesh with 4 nodes
        PottsMeshGenerator<2> generator(2, 0, 0, 2, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 3);

        std::vector<unsigned> location_indices;
        location_indices.push_back(0);
        location_indices.push_back(3);
        location_indices.push_back(1);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 2);

        // Check cells are in the correct location
        TS_ASSERT(cell_population.IsCellAttachedToLocationIndex(0));
        TS_ASSERT(cell_population.IsCellAttachedToLocationIndex(1));
        TS_ASSERT(!cell_population.IsCellAttachedToLocationIndex(2));
        TS_ASSERT(cell_population.IsCellAttachedToLocationIndex(3));

        cell_population.GetCellUsingLocationIndex(0);
        cell_population.GetCellUsingLocationIndex(1);
        TS_ASSERT_THROWS_THIS(cell_population.GetCellUsingLocationIndex(2),
            "Location index input argument does not correspond to a Cell");
        TS_ASSERT_THROWS_NOTHING(cell_population.GetCellUsingLocationIndex(3));

        // Now remove first cell from lattice 0 and move it to lattice 3
        cell_population.RemoveCellUsingLocationIndex(0,cells[0]);
        cell_population.AddCellUsingLocationIndex(3,cells[0]);

        // Make sure the maps have been set properly
        TS_ASSERT_EQUALS(cell_population.GetCellsUsingLocationIndex(0).size(), 0u);
        TS_ASSERT_EQUALS(cell_population.GetCellsUsingLocationIndex(3).size(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(cells[0]),  3u);

        // Try moving one of the cells
        cell_population.MoveCellInLocationMap(cells[1], 3, 0);
        TS_ASSERT_EQUALS(cell_population.GetCellsUsingLocationIndex(0).size(), 1u);
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(cells[1]), 0u);

        // Now move it back
        cell_population.MoveCellInLocationMap(cells[1], 0, 3);

        // Coverage as cell is no longer there
        TS_ASSERT_THROWS_THIS(cell_population.RemoveCellUsingLocationIndex(0,cells[0]),
            "Tried to remove a cell which is not attached to the given location index");

        // Check cells are in the correct locations
        TS_ASSERT(!cell_population.IsCellAttachedToLocationIndex(0));
        TS_ASSERT(cell_population.IsCellAttachedToLocationIndex(1));
        TS_ASSERT(!cell_population.IsCellAttachedToLocationIndex(2));
        TS_ASSERT(cell_population.IsCellAttachedToLocationIndex(3));

        TS_ASSERT_THROWS_THIS(cell_population.GetCellUsingLocationIndex(0),
            "Location index input argument does not correspond to a Cell");
        cell_population.GetCellUsingLocationIndex(1);
        TS_ASSERT_THROWS_THIS(cell_population.GetCellUsingLocationIndex(2),
            "Location index input argument does not correspond to a Cell");
        TS_ASSERT_THROWS_THIS(cell_population.GetCellUsingLocationIndex(3),
            "Multiple cells are attached to a single location index.");

        // Now remove first cell from lattice 0 and move it to lattice 3
        cell_population.RemoveCellUsingLocationIndex(1, cells[2]);
        TS_ASSERT_THROWS_THIS(cell_population.AddCellUsingLocationIndex(3, cells[2]),
            "No available spaces at location index 3.");

        // Check GetCellsUsingLocationIndex
        std::set<CellPtr> cells_on_lattice = cell_population.GetCellsUsingLocationIndex(3);
        TS_ASSERT_EQUALS(cells_on_lattice.size(), 2u);

        bool found0 = false;
        bool found1 = false;
        bool foundother = false;

        // This set iterator is iterator in pointer order (the order of addresses in memory); this order can be arbitrary
        for (std::set<CellPtr>::iterator iter = cells_on_lattice.begin();
             iter != cells_on_lattice.end();
             iter++)
        {
            if ((*iter)->GetCellId() == 0)
            {
                found0 = true;
            }
            else if ((*iter)->GetCellId() == 1)
            {
                found1 = true;
            }
            else
            {
                foundother = true;
            }
        }
        TS_ASSERT(found0 && found1);
        TS_ASSERT_EQUALS(foundother, false);

        TS_ASSERT_EQUALS(cells[0]->GetCellId(), 0u);
        TS_ASSERT_EQUALS(cells[1]->GetCellId(), 1u);
    }

    void TestWriteResultsToFileAndOutputCellPopulationParameters()
    {
        EXIT_IF_PARALLEL;    // We cannot currently write to file in parallel.

        // Reset the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        std::string output_directory = "TestMultipleCaBasedCellPopulationWriters";
        OutputFileHandler output_file_handler(output_directory, false);

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

        // For coverage, label one cell
        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());
        cell_population.GetCellUsingLocationIndex(location_indices[0])->AddCellProperty(p_label);

        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "MultipleCaBasedCellPopulation-2");

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellVariables(true);
        cell_population.SetOutputCellVolumes(true);

        // VTK writing needs a simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
        TS_ASSERT_THROWS_NOTHING(cell_population.CreateOutputFiles(output_directory, false));
        cell_population.OpenWritersFiles();
        cell_population.WriteResultsToFiles();
        TS_ASSERT_THROWS_NOTHING(cell_population.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        FileComparison(results_dir + "results.viznodes", "cell_based/test/data/TestMultipleCaBasedCellPopulationWriters/results.viznodes").CompareFiles();
        FileComparison(results_dir + "results.vizlocations", "cell_based/test/data/TestMultipleCaBasedCellPopulationWriters/results.vizlocations").CompareFiles();
        FileComparison(results_dir + "results.vizcelltypes", "cell_based/test/data/TestMultipleCaBasedCellPopulationWriters/results.vizcelltypes").CompareFiles();
        FileComparison(results_dir + "results.vizancestors", "cell_based/test/data/TestMultipleCaBasedCellPopulationWriters/results.vizancestors").CompareFiles();
        FileComparison(results_dir + "cellmutationstates.dat", "cell_based/test/data/TestMultipleCaBasedCellPopulationWriters/cellmutationstates.dat").CompareFiles();
        FileComparison(results_dir + "cellages.dat", "cell_based/test/data/TestMultipleCaBasedCellPopulationWriters/cellages.dat").CompareFiles();
        FileComparison(results_dir + "cellareas.dat", "cell_based/test/data/TestMultipleCaBasedCellPopulationWriters/cellareas.dat").CompareFiles();

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Write cell population parameters to file
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        // Compare output with saved files of what they should look like
        FileComparison( results_dir + "results.parameters", "cell_based/test/data/TestMultipleCaBasedCellPopulationWriters/results.parameters").CompareFiles();
#ifdef CHASTE_VTK
        //Test that VTK writer has produced a file
        FileFinder vtk_file(results_dir + "results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());
 #endif //CHASTE_VTK
    }

    void TestRemoveDeadCellsAndUpdate() throw(Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 2);

        std::vector<unsigned> location_indices;
        location_indices.push_back(12);
        location_indices.push_back(13);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 2u);

        // Test RemoveDeadCells() method
        cell_population.Begin()->Kill();

        TS_ASSERT_EQUALS(cell_population.RemoveDeadCells(), 1u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1u);
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*(cell_population.Begin())),13u);

        // Test that Update() throws no errors (note that for this cell population, Update() doesn't actually do anything)
        TS_ASSERT_THROWS_NOTHING(cell_population.Update());
    }

    void TestAddCell() throw(Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 1);

        std::vector<unsigned> location_indices;
        unsigned initial_cell_index = 12;
        location_indices.push_back(initial_cell_index);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 1u);

        // Create a new cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetCellProliferativeType(p_stem_type);

        // Add new cell to the cell population by dividing the cell
        AbstractCellPopulation<2>::Iterator cell_iter_1 = cell_population.Begin();
        cell_population.AddCell(p_new_cell, zero_vector<double>(2), *cell_iter_1);

        // Test that the initial cell is still in the same place
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter_1), initial_cell_index);

        // Test the location of the new cell - should go to 13 with this random seed.
        unsigned lattice_site = 13u;
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_new_cell), lattice_site);

        // Test that the population size is now 2
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 2u);

        // Check locations of parent and daughter cell
        AbstractCellPopulation<2>::Iterator cell_iter_2 = cell_population.Begin();
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter_2), initial_cell_index);
        ++cell_iter_2;
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter_2), lattice_site);

        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetNumElements(),0u);

        ///\todo Check elements are as expected
//        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(0)->GetNodeGlobalIndex(0), 0u);
//        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(0)->GetNodeGlobalIndex(1), 1u);
//        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(1)->GetNodeGlobalIndex(0), 2u);
//        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(1)->GetNodeGlobalIndex(1), 3u);
    }

///\todo reinstate or remove test (#2031)
    void DONOTTestAddCellToManyCells() throw(Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 25);

        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<25; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 25u);

        // Create a new cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        CellPtr p_new_cell(new Cell(p_state, p_model));
        p_new_cell->SetCellProliferativeType(p_stem_type);

        // Test thatAdd new cell to the cell population by dividing the cell
        AbstractCellPopulation<2>::Iterator cell_iter_1 = cell_population.Begin();
        TS_ASSERT_THROWS_THIS(cell_population.AddCell(p_new_cell, zero_vector<double>(2), *cell_iter_1),
            "No free space to divide.");
    }

    void TestUpdateCellLocationsExceptions()
    {
        // Create a simple 2D PottsMesh with one cells
        PottsMeshGenerator<2> generator(3, 0, 0, 3, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 1u);

        std::vector<unsigned> location_indices;
        unsigned initial_cell_index = 0;
        location_indices.push_back(initial_cell_index);

        // Create cell population
        MultipleCaBasedCellPopulation<2u> cell_population(*p_mesh, cells, location_indices);

        // Set node selection to non-random lattice sweeping: this will loop over the nodes in index order
        TS_ASSERT_EQUALS(true, cell_population.GetUpdateNodesInRandomOrder());
        cell_population.SetUpdateNodesInRandomOrder(false);

        // Create a multiple CA update rule and pass to the population
        MAKE_PTR(DiffusionMultipleCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(1.0);
        cell_population.AddUpdateRule(p_diffusion_update_rule);

        // Commence lattice sweeping, updating where necessary
        TS_ASSERT_THROWS_THIS(cell_population.UpdateCellLocations(-1.0),
            "The probability of cellular movement is smaller than zero. In order to prevent it from happening you should change your time step and parameters");
        TS_ASSERT_THROWS_THIS(cell_population.UpdateCellLocations(5.0),
            "The probability of the cellular movement is bigger than one. In order to prevent it from happening you should change your time step and parameters");
        TS_ASSERT_THROWS_THIS(cell_population.UpdateCellLocations(1.0),
            "The probability of the cell not moving is smaller than zero. In order to prevent it from happening you should change your time step and parameters");
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 1u);
    }

    // For coverage try move when all neighbours are occupied
    void TestUpdateCellLocationsWhenFull()
    {
        // Create a simple 2D PottsMesh entirely populated with cells
        PottsMeshGenerator<2> generator(3, 0, 0, 3, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 9);

        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<9; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        MultipleCaBasedCellPopulation<2u> cell_population(*p_mesh, cells, location_indices);

        unsigned location_index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), location_index);
            location_index++;
        }

        // Commence lattice sweeping, no cells will move
        TS_ASSERT_THROWS_NOTHING(cell_population.UpdateCellLocations(0.1));

        location_index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), location_index);
            location_index++;
        }
    }

    void TestUpdateCellLocationsRandomlyExceptions()
    {
        // Create a simple 2D PottsMesh with two cells
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 1);

        std::vector<unsigned> location_indices;
        unsigned initial_cell_index = 12;
        location_indices.push_back(initial_cell_index);

        // Create cell population
        MultipleCaBasedCellPopulation<2u> cell_population(*p_mesh, cells, location_indices);

        // Set node selection to random lattice sweeping and (for coverage) random iteration over update rules
        cell_population.SetUpdateNodesInRandomOrder(true);
        cell_population.SetIterateRandomlyOverUpdateRuleCollection(true);

        // Create a multiple CA update rule and pass to the population
        MAKE_PTR(DiffusionMultipleCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(1.0);
        cell_population.AddUpdateRule(p_diffusion_update_rule);

        // Commence lattice sweeping, updating where necessary
        TS_ASSERT_THROWS_THIS(cell_population.UpdateCellLocations(-1.0),
            "The probability of cellular movement is smaller than zero. In order to prevent it from happening you should change your time step and parameters");
        TS_ASSERT_THROWS_THIS(cell_population.UpdateCellLocations(5.0),
            "The probability of the cellular movement is bigger than one. In order to prevent it from happening you should change your time step and parameters");
        TS_ASSERT_THROWS_THIS(cell_population.UpdateCellLocations(1.0),
            "The probability of the cell not moving is smaller than zero. In order to prevent it from happening you should change your time step and parameters");
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 1u);
    }

    void TestArchiving() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "MultipleCaBasedCellPopulation-2.arch";

        // The following line is required because the loading of a cell population
        // is usually called by the method CellBasedSimulation::Load()
        ArchiveLocationInfo::SetMeshFilename("potts_mesh_2d");

        // Archive cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a multiple CA cell population object
            PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
            PottsMesh<2>* p_mesh = generator.GetMesh();

            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, 40, p_diff_type);

            std::vector<unsigned> location_indices;
            for (unsigned index=0; index<10; index++)
            {
                location_indices.push_back(index);
                location_indices.push_back(index);
                location_indices.push_back(index);
                location_indices.push_back(index);
            }

            AbstractCellPopulation<2>* const p_cell_population =
                new MultipleCaBasedCellPopulation<2>(*p_mesh, cells, location_indices, 4);

            // Run each cell to time 0
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Create an update rule and pass to the population
            MAKE_PTR(DiffusionMultipleCaUpdateRule<2u>, p_diffusion_update_rule);
            p_diffusion_update_rule->SetDiffusionParameter(1.0);
            static_cast<MultipleCaBasedCellPopulation<2>*>(p_cell_population)->AddUpdateRule(p_diffusion_update_rule);

            // Set member variables in order to test that they are archived correctly
            static_cast<MultipleCaBasedCellPopulation<2>*>(p_cell_population)->SetUpdateNodesInRandomOrder(false);
            static_cast<MultipleCaBasedCellPopulation<2>*>(p_cell_population)->SetIterateRandomlyOverUpdateRuleCollection(true);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Archive the cell population
            (*p_arch) << static_cast<const SimulationTime&>(*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Tidy up
            SimulationTime::Destroy();
            delete p_cell_population;
        }

        // Restore cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            AbstractCellPopulation<2>* p_cell_population;

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore the cell population
            (*p_arch) >> *p_simulation_time;
            (*p_arch) >> p_cell_population;

            // Test that the member variables have been archived correctly
            MultipleCaBasedCellPopulation<2>* p_static_population = static_cast<MultipleCaBasedCellPopulation<2>*>(p_cell_population);

            TS_ASSERT_EQUALS(p_static_population->GetNumNodes(), 100u);

            std::vector<unsigned> available_space = p_static_population->rGetAvailableSpaces();
            for (unsigned i=0; i<10; i++)
            {
                TS_ASSERT_EQUALS(available_space[i], 0u);
            }
            for (unsigned i=10; i<100; i++)
            {
                TS_ASSERT_EQUALS(available_space[i], 4u);
            }

            c_vector<double, 2> cell_location = p_static_population->GetLocationOfCellCentre(*(p_cell_population->Begin()));
            TS_ASSERT_DELTA(cell_location[0], 0.0, 1e-4);
            TS_ASSERT_DELTA(cell_location[1], 0.0, 1e-4);

            TS_ASSERT_EQUALS(p_static_population->GetUpdateNodesInRandomOrder(), false);
            TS_ASSERT_EQUALS(p_static_population->GetIterateRandomlyOverUpdateRuleCollection(), true);

            // Test that the update rule has been archived correctly
            std::vector<boost::shared_ptr<AbstractMultipleCaUpdateRule<2> > > update_rule_collection = p_static_population->rGetUpdateRuleCollection();
            TS_ASSERT_EQUALS(update_rule_collection.size(), 1u);
            TS_ASSERT_EQUALS((*update_rule_collection[0]).GetIdentifier(), "DiffusionMultipleCaUpdateRule-2");

            // Tidy up
            delete p_cell_population;
        }
    }
};

#endif /*TESTMULTIPLECABASEDCELLPOPULATION_HPP_*/
