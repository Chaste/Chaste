/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTPOTTSBASEDCELLPOPULATION_HPP_
#define TESTPOTTSBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellwiseData.hpp"
#include "CellsGenerator.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "PottsMeshGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ArchiveOpener.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"

class TestPottsBasedCellPopulation : public AbstractCellBasedTestSuite
{
public:

    // Test construction, accessors and iterator
    void TestCreateSmallPottsBasedCellPopulationAndGetWidth() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        unsigned num_cells = cells.size();
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 4u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 4u);

        // Coverage
        TS_ASSERT_DELTA(cell_population.GetTemperature(), 0.1, 1e-6);
        cell_population.SetTemperature(10.0);
        TS_ASSERT_DELTA(cell_population.GetTemperature(), 10.0, 1e-6);
        cell_population.SetTemperature(0.1);
        TS_ASSERT_EQUALS(cell_population.GetUpdateNodesInRandomOrder(), true);
        cell_population.SetUpdateNodesInRandomOrder(false);
        TS_ASSERT_EQUALS(cell_population.GetUpdateNodesInRandomOrder(), false);
        cell_population.SetUpdateNodesInRandomOrder(true);

        ChastePoint<2> location;
        TS_ASSERT_THROWS_THIS(cell_population.AddNode(NULL), "Cannot call AddNode on a PottsBasedCellPopulation");
        TS_ASSERT_THROWS_THIS(cell_population.SetNode(0, location), "Cannot call SetNode on a PottsBasedCellPopulation");
        TS_ASSERT_THROWS_THIS(cell_population.GetDampingConstant(0), "Cannot call GetDampingConstant on a PottsBasedCellPopulation");

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), p_mesh->GetNumElements());
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), num_cells);

        // Test GetNeighbouringNodeIndices() method
        TS_ASSERT_THROWS_THIS(cell_population.GetNeighbouringNodeIndices(10), "Cannot call GetNeighbouringNodeIndices() on a PottsBasedCellPopulation, need to go through the PottsMesh instead");

        // Test PottsBasedCellPopulation<2>::Iterator
        unsigned counter = 0;
        for (PottsBasedCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            // Test operator* and that cells are in sync
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), counter);

            // Test operator-> and that cells are in sync
            TS_ASSERT_DELTA(cell_iter->GetAge(), (double)counter, 1e-12);

            TS_ASSERT_EQUALS(counter, p_mesh->GetElement(counter)->GetIndex());

            counter++;
        }

        // Test we have gone through all cells in the for loop
        TS_ASSERT_EQUALS(counter, cell_population.GetNumRealCells());

        // Test GetNumNodes() method
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), p_mesh->GetNumNodes());

        // Test GetWidth() method
        double width_x = cell_population.GetWidth(0);
        TS_ASSERT_DELTA(width_x, 3.0, 1e-4);

        double width_y = cell_population.GetWidth(1);
        TS_ASSERT_DELTA(width_y, 3.0, 1e-4);

        // Test RemoveDeadCells() method
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 4u);
        cell_population.Begin()->Kill();
        TS_ASSERT_EQUALS(cell_population.RemoveDeadCells(), 1u);
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 3u);
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 3u);
        TS_ASSERT_EQUALS(p_mesh->GetNumAllElements(), 4u);
    }

    void TestValidate() throw (Exception)
    {
        // Create a simple potts-based mesh
    	PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements()-1);

        std::vector<unsigned> cell_location_indices;
        for (unsigned i=0; i<cells.size(); i++)
        {
            cell_location_indices.push_back(i);
        }

        // This should throw an exception as the number of cells does not equal the number of elements
        std::vector<CellPtr> cells_copy(cells);
        TS_ASSERT_THROWS_THIS(PottsBasedCellPopulation<2> cell_population(*p_mesh, cells_copy),
                "Element 3 does not appear to have a cell associated with it");

        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_cell(new Cell(p_state, p_model));

        double birth_time = 0.0 - p_mesh->GetNumElements()-1;
        p_cell->SetBirthTime(birth_time);

        cells.push_back(p_cell);
        cell_location_indices.push_back(p_mesh->GetNumElements()-1);

        // This should pass as the number of cells equals the number of elements
        std::vector<CellPtr> cells_copy2(cells);
        TS_ASSERT_THROWS_NOTHING(PottsBasedCellPopulation<2> cell_population(*p_mesh, cells_copy2));

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Check correspondence between elements and cells
        for (PottsMesh<2>::PottsElementIterator iter = p_mesh->GetElementIteratorBegin();
             iter != p_mesh->GetElementIteratorEnd();
             ++iter)
        {
            std::set<unsigned> expected_node_indices;
            unsigned expected_index = iter->GetIndex();

            for (unsigned i=0; i<iter->GetNumNodes(); i++)
            {
                expected_node_indices.insert(iter->GetNodeGlobalIndex(i));
            }

            std::set<unsigned> actual_node_indices;
            unsigned elem_index = iter->GetIndex();
            CellPtr p_cell = cell_population.GetCellUsingLocationIndex(elem_index);
            PottsElement<2>* p_actual_element = cell_population.GetElementCorrespondingToCell(p_cell);
            unsigned actual_index = p_actual_element->GetIndex();

            for (unsigned i=0; i<p_actual_element->GetNumNodes(); i++)
            {
                actual_node_indices.insert(p_actual_element->GetNodeGlobalIndex(i));
            }

            TS_ASSERT_EQUALS(actual_index, expected_index);
            TS_ASSERT_EQUALS(actual_node_indices, expected_node_indices);
        }

        // Create another simple potts-based mesh
        PottsMeshGenerator<2> generator2(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh2 = generator2.GetMesh();

        // Create cells
        std::vector<CellPtr> cells2;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator2;
        cells_generator2.GenerateBasic(cells2, p_mesh2->GetNumElements()+1);

        std::vector<unsigned> cell_location_indices2;
        for (unsigned i=0; i<cells2.size(); i++)
        {
            cell_location_indices2.push_back(i%p_mesh2->GetNumElements()); // Element 0 will have 2 cells
        }

        // This should throw an exception as the number of cells
        // does not equal the number of elements
        TS_ASSERT_THROWS_THIS(PottsBasedCellPopulation<2> cell_population2(*p_mesh2, cells2, false, true, cell_location_indices2),
                "Element 0 appears to have 2 cells associated with it");
    }

    void TestCellDivision() throw(Exception)
    {
        // Create a simple 2D PottsMesh with one cell
    	PottsMeshGenerator<2> generator(2, 1, 2, 2, 1, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 1u);

        // Create a new cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_new_cell(new Cell(p_state, p_model));

        // Add new cell to the cell population by dividing the cell
        AbstractCellPopulation<2>::Iterator cell_iter_1 = cell_population.Begin();
        cell_population.AddCell(p_new_cell, zero_vector<double>(2), *cell_iter_1);

        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_new_cell), 1u);

        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 2u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 2u);

        // Check locations of parent and daughter cell
        AbstractCellPopulation<2>::Iterator cell_iter_2 = cell_population.Begin();
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter_2), 0u);
        ++cell_iter_2;
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter_2), 1u);

        // Check Elements are as expected
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(0)->GetNodeGlobalIndex(0), 0u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(0)->GetNodeGlobalIndex(1), 1u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(1)->GetNodeGlobalIndex(0), 2u);
        TS_ASSERT_EQUALS(cell_population.rGetMesh().GetElement(1)->GetNodeGlobalIndex(1), 3u);
    }

///\todo implement this test (#1666)
//    void TestVoronoiMethods()
//    {
//        // Create a simple 2D PottsMesh
//        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
//        PottsMesh<2>* p_mesh = generator.GetMesh();
//
//        // Create cells
//        std::vector<CellPtr> cells;
//        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());
//
//        // Create cell population
//        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
//
//        // Create element tessellation
//        cell_population.CreateElementTessellation();
//
//        VertexMesh<2,2>* p_tessellation = cell_population.GetElementTessellation();
//
//        TS_ASSERT(p_tessellation != NULL);
//
//        TS_ASSERT_EQUALS(p_tessellation->GetNumNodes(), 25u);
//        TS_ASSERT_EQUALS(p_tessellation->GetNumElements(), 16u);
//    }

    void TestPottsBasedCellPopulationWriters()
    {
        std::string output_directory = "TestPottsBasedCellPopulationWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(4, 2, 2, 4, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "PottsBasedCellPopulation-2");

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellVariables(true);
        cell_population.SetOutputCellVolumes(true);

        //VTK writing needs a simulation time
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        TS_ASSERT_THROWS_NOTHING(cell_population.CreateOutputFiles(output_directory, false));

        cell_population.WriteResultsToFiles();

        TS_ASSERT_THROWS_NOTHING(cell_population.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     notforrelease_cell_based/test/data/TestPottsBasedCellPopulationWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     notforrelease_cell_based/test/data/TestPottsBasedCellPopulationWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     notforrelease_cell_based/test/data/TestPottsBasedCellPopulationWriters/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     notforrelease_cell_based/test/data/TestPottsBasedCellPopulationWriters/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellages.dat     notforrelease_cell_based/test/data/TestPottsBasedCellPopulationWriters/cellages.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat     notforrelease_cell_based/test/data/TestPottsBasedCellPopulationWriters/cellareas.dat").c_str()), 0);

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Write cell population parameters to file
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        // Compare output with saved files of what they should look like
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.parameters    notforrelease_cell_based/test/data/TestPottsBasedCellPopulationWriters/results.parameters").c_str()), 0);
    }

    void TestAddingUpdateRules() throw(Exception)
    {
        // Create a simple 2D PottsMesh with one cell
        PottsMeshGenerator<2> generator(2, 1, 2, 2, 1, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Test we have the correct number of cells and elements
        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1u);
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 1u);

        // Create an update rule and pass to the population
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        cell_population.AddUpdateRule(p_volume_constraint_update_rule);

        // Check the update rules are correct
        std::vector<boost::shared_ptr<AbstractPottsUpdateRule<2> > > update_rule_collection = cell_population.rGetUpdateRuleCollection();
        TS_ASSERT_EQUALS(update_rule_collection.size(),1u);
        TS_ASSERT_EQUALS((*update_rule_collection[0]).GetIdentifier(), "VolumeConstraintPottsUpdateRule-2");
    }

    void TestArchiving2dPottsBasedCellPopulation() throw(Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "potts_cell_population_2d.arch";

        // The following line is required because the loading of a cell population
        // is usually called by the method CellBasedSimulation::Load()
        ArchiveLocationInfo::SetMeshFilename("potts_mesh_2d");

        // Create mesh
        PottsMeshReader<2> mesh_reader("notforrelease_cell_based/test/data/TestPottsMeshWriter/potts_mesh_2d");
        PottsMesh<2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        // Archive cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create cells
            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, mesh.GetNumElements());

            // Create cell population
            PottsBasedCellPopulation<2>* const p_cell_population = new PottsBasedCellPopulation<2>(mesh, cells);

            // Cells have been given birth times of 0 and -1.
            // Loop over them to run to time 0.0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Create an update rule and pass to the population
            MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
            p_cell_population->AddUpdateRule(p_volume_constraint_update_rule);

            // Create output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Write the cell population to the archive
            (*p_arch) << static_cast<const SimulationTime&> (*p_simulation_time);
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

            // Create an input archive
            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            // Restore the cell population
            (*p_arch) >> *p_simulation_time;
            PottsBasedCellPopulation<2>* p_cell_population;
            (*p_arch) >> p_cell_population;

            // Check the cell population has been restored correctly
            TS_ASSERT_EQUALS(p_cell_population->rGetCells().size(), 2u);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // this checks that individual cells and their models are archived.
            unsigned counter = 0;
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                TS_ASSERT_DELTA(cell_iter->GetAge(), (double)(counter), 1e-7);
                counter++;
            }

            // Check the simulation time has been restored (through the cell)
            TS_ASSERT_EQUALS(p_simulation_time->GetTime(), 0.0);

            // Check the mesh has been restored correctly
            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumNodes(), 6u);
            TS_ASSERT_EQUALS(p_cell_population->rGetMesh().GetNumElements(), 2u);

            // Compare the loaded mesh against the original
            PottsMesh<2>& loaded_mesh = p_cell_population->rGetMesh();

            TS_ASSERT_EQUALS(mesh.GetNumNodes(), loaded_mesh.GetNumNodes());

            for (unsigned node_index=0; node_index<mesh.GetNumNodes(); node_index++)
            {
                Node<2>* p_node = mesh.GetNode(node_index);
                Node<2>* p_node2 = loaded_mesh.GetNode(node_index);

                TS_ASSERT_EQUALS(p_node->IsDeleted(), p_node2->IsDeleted());
                TS_ASSERT_EQUALS(p_node->GetIndex(), p_node2->GetIndex());

                for (unsigned dimension=0; dimension<2; dimension++)
                {
                    TS_ASSERT_DELTA(p_node->rGetLocation()[dimension], p_node2->rGetLocation()[dimension], 1e-4);
                }
            }

            TS_ASSERT_EQUALS(mesh.GetNumElements(), loaded_mesh.GetNumElements());

            for (unsigned elem_index=0; elem_index < mesh.GetNumElements(); elem_index++)
            {
                TS_ASSERT_EQUALS(mesh.GetElement(elem_index)->GetNumNodes(),
                                 loaded_mesh.GetElement(elem_index)->GetNumNodes());

                for (unsigned local_index=0; local_index<mesh.GetElement(elem_index)->GetNumNodes(); local_index++)
                {
                    unsigned this_index = mesh.GetElement(elem_index)->GetNodeGlobalIndex(local_index);
                    unsigned loaded_index = loaded_mesh.GetElement(elem_index)->GetNodeGlobalIndex(local_index);

                    TS_ASSERT_EQUALS(this_index, loaded_index);
                }
            }

            // Check Update rules are restored correctly
            std::vector<boost::shared_ptr<AbstractPottsUpdateRule<2> > > update_rule_collection = p_cell_population->rGetUpdateRuleCollection();
            TS_ASSERT_EQUALS(update_rule_collection.size(),1u);
            TS_ASSERT_EQUALS((*update_rule_collection[0]).GetIdentifier(), "VolumeConstraintPottsUpdateRule-2");

            // Tidy up
            delete p_cell_population;
        }
    }
};

#endif /*TESTPOTTSBASEDCELLPOPULATION_HPP_*/
