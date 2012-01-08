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

#ifndef TESTLATTICEBASEDCELLPOPULATION_HPP_
#define TESTLATTICEBASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "ArchiveOpener.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "AdvectionCaUpdateRule.hpp"

class TestCaBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    /**
     * Helper method. Create a single, wild type, differentiated cell and
     * return as a vector for passing into a cell population constructor.
     *
     * @param cellProliferativeType  the proliferative type of the cell (defaults to DIFFERENTIATED)
     */
    std::vector<CellPtr> CreateSingleCellPtr(CellProliferativeType cellProliferativeType=DIFFERENTIATED)
    {
        std::vector<CellPtr> cells;

        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(cellProliferativeType);
        CellPtr p_cell(new Cell(p_state, p_model));

        cells.push_back(p_cell);
        return cells;
    }

    /**
     * Helper method. Return a given location index corresponding to a single
     * cell as a vector for passing into a cell population constructor.
     *
     * @param locationIndex the location index
     */
    std::vector<unsigned> CreateSingleLocationIndex(unsigned locationIndex)
    {
        std::vector<unsigned> location_indices;
        location_indices.push_back(locationIndex);
        return location_indices;
    }

public:

    /*
     * Note that this also tests the method SetEmptySites(), which is
     * called within the constructor.
     */
    void TestConstructorWithLocationIndices() throw(Exception)
    {
        // Create a mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create a vector of two real node indices
        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(4);
        real_node_indices.push_back(1);

        // Create two cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices);

        // Create a CA-based cell population object
        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Test that the mesh and cells are correctly assigned
        TS_ASSERT_EQUALS(&(cell_population.rGetMesh()), &mesh);

        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();

        TS_ASSERT_EQUALS(cell_population.IsCellAssociatedWithADeletedLocation(*cell_iter), false);

        TS_ASSERT_EQUALS(*cell_iter, cells[0]);
        ++cell_iter;
        TS_ASSERT_EQUALS(*cell_iter, cells[1]);

        // Test that we do not have any update rules present
        TS_ASSERT_EQUALS(cell_population.rGetUpdateRuleCollection().empty(), true);

        // Create a vector of expected values of whether each site is empty
        unsigned num_nodes = mesh.GetNumNodes();
        std::vector<bool> expected_empty_sites = std::vector<bool>(num_nodes, true);
        expected_empty_sites[1] = false;
        expected_empty_sites[4] = false;

        // Test that the vector mIsEmptySite has been initialised correctly
        TS_ASSERT_EQUALS(cell_population.rGetEmptySites(), expected_empty_sites);

        // Test this another way, for coverage
        for (unsigned i=0; i<num_nodes; i++)
        {
            bool expected_empty_site = (i!=1 && i!=4);
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(i), expected_empty_site);
        }

        // Test this yet another way, for coverage
        std::set<unsigned> empty_site_indices = cell_population.GetEmptySiteIndices();
        std::set<unsigned>::iterator not_found_iter = empty_site_indices.end();
        for (unsigned i=0; i<num_nodes; i++)
        {
            if (i==1 || i==4)
            {
                TS_ASSERT_EQUALS(empty_site_indices.find(i), not_found_iter);
            }
            else
            {
                TS_ASSERT_DIFFERS(empty_site_indices.find(i), not_found_iter);
            }
        }

        // Test that the other member variables of this object are initialised correctly
        TS_ASSERT_EQUALS(cell_population.GetOnlyUseNearestNeighboursForDivision(), false);
        TS_ASSERT_EQUALS(cell_population.GetUseVonNeumannNeighbourhoods(), false);
        TS_ASSERT_EQUALS(cell_population.GetUpdateNodesInRandomOrder(), true);
        TS_ASSERT_EQUALS(cell_population.GetIterateRandomlyOverUpdateRuleCollection(), false);
    }

    void TestConstructorWithoutLocationIndices() throw(Exception)
    {
        // Create a fully populated CA-based cell population object
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true);

        unsigned num_nodes = mesh.GetNumNodes();

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, num_nodes);

        CaBasedCellPopulation<2> cell_population(mesh, cells);

        // Test that the vector mIsEmptySite has been initialised correctly
        for (unsigned i=0; i<num_nodes; i++)
        {
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(i), false);
        }
    }

    void TestValidateCaBasedCellPopulation()
    {
        // Create a CA-based cell population object with five cells
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true);

        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<5; i++)
        {
            real_node_indices.push_back(i);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices);

        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Test Validate() when node 2 has been set as an empty site
        cell_population.mIsEmptySite[2] = true;

        TS_ASSERT_THROWS_THIS(cell_population.Validate(),
                              "Node 2 is labelled as an empty site and has a cell attached");

        // Test Validate() when node 7 has been set as an occupied site
        cell_population.mIsEmptySite[2] = false;
        cell_population.mIsEmptySite[7] = false;

        TS_ASSERT_THROWS_THIS(cell_population.Validate(),
                              "Node 7 does not appear to be an empty site or have a cell associated with it");
    }

    void TestAddUpdateRule() throw(Exception)
    {
        // Create a CA-based cell population object with five cells
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true);

        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<5; i++)
        {
            real_node_indices.push_back(i);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices);

        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        TS_ASSERT_EQUALS(cell_population.rGetUpdateRuleCollection().empty(), true);

        // Add an update rule
        MAKE_PTR_ARGS(AdvectionCaUpdateRule<2>, p_update_rule, (1, 1.0));
        cell_population.AddUpdateRule(p_update_rule);

        TS_ASSERT_EQUALS(cell_population.rGetUpdateRuleCollection().size(), 1u);
    }

    void TestNodeAndMeshMethods() throw(Exception)
    {
        // Create a CA-based cell population object
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true);

        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(4);
        real_node_indices.push_back(1);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Test GetNumNodes()
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 9u);

        // Test GetNode()
        for (unsigned index=0; index<cell_population.GetNumNodes(); index++)
        {
            Node<2>* p_node = cell_population.GetNode(index);
            TS_ASSERT_EQUALS(p_node->GetIndex(), index);

            c_vector<double, 2> node_location = p_node->rGetLocation();
            double expected_x = (double)(index%3);
            double expected_y = (double)(index>2) + (double)(index>5);
            TS_ASSERT_DELTA(node_location[0], expected_x, 1e-3);
            TS_ASSERT_DELTA(node_location[1], expected_y, 1e-3);
        }

        // Test SetNode()
        ChastePoint<2> unused_point;
        TS_ASSERT_THROWS_THIS(cell_population.SetNode(0, unused_point),
                              "SetNode() cannot be called on a subclass of AbstractOnLatticeCellPopulation.");

        // Test GetWidth() method
        double width_x = cell_population.GetWidth(0);
        TS_ASSERT_DELTA(width_x, 2.0, 1e-4);

        double width_y = cell_population.GetWidth(1);
        TS_ASSERT_DELTA(width_y, 2.0, 1e-4);
    }

    void TestGetLocationOfCellCentre() throw (Exception)
    {
        // Create a CA-based cell population object
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true);

        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(4);
        real_node_indices.push_back(1);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Test GetLocationOfCellCentre()
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        c_vector<double, 2> cell_1_location = cell_population.GetLocationOfCellCentre(*cell_iter);
        TS_ASSERT_DELTA(cell_1_location[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(cell_1_location[1], 1.0, 1e-6);

        ++cell_iter;

        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 1u);
        c_vector<double, 2> cell_2_location = cell_population.GetLocationOfCellCentre(*cell_iter);
        TS_ASSERT_DELTA(cell_2_location[0], 1.0, 1e-6);
        TS_ASSERT_DELTA(cell_2_location[1], 0.0, 1e-6);
    }

    /*
     * Note that this also tests the methods WriteCellVolumeResultsToFile(),
     * GenerateCellResultsAndWriteToFiles() and GenerateCellResults(), which are
     * called within the overridden method  WriteResultsToFiles().
     */
    void TestWriteResultsToFileAndOutputCellPopulationParameters()
    {
        // Resetting the maximum cell ID to zero (to account for previous tests)
        Cell::ResetMaxCellId();

        std::string output_directory = "TestCaBasedCellPopulationWriters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Create a cell population
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true);

        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<5; i++)
        {
            real_node_indices.push_back(i);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        MAKE_PTR(WildTypeCellMutationState, p_wt);
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(0);
            cells[i]->SetMutationState(p_wt);
        }

        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Test this object has the correct identifier
        TS_ASSERT_EQUALS(cell_population.GetIdentifier(), "CaBasedCellPopulation-2");

        cell_population.SetCellAncestorsToLocationIndices();
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);
        cell_population.SetOutputCellVariables(true);
        cell_population.SetOutputCellVolumes(true);

        // For coverage of GetVolumeOfCell()
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_population.GetVolumeOfCell(*cell_iter), 1.0, 1e-6);
        }

        // For coverage of WriteResultsToFiles()
        boost::shared_ptr<AbstractCellProperty> p_state(cell_population.GetCellPropertyRegistry()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc1(cell_population.GetCellPropertyRegistry()->Get<ApcOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apc2(cell_population.GetCellPropertyRegistry()->Get<ApcTwoHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_bcat1(cell_population.GetCellPropertyRegistry()->Get<BetaCateninOneHitCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_state(cell_population.GetCellPropertyRegistry()->Get<ApoptoticCellProperty>());
        boost::shared_ptr<AbstractCellProperty> p_label(cell_population.GetCellPropertyRegistry()->Get<CellLabel>());

        cell_population.GetCellUsingLocationIndex(0)->GetCellCycleModel()->SetCellProliferativeType(TRANSIT);
        cell_population.GetCellUsingLocationIndex(0)->AddCellProperty(p_label);
        cell_population.GetCellUsingLocationIndex(1)->GetCellCycleModel()->SetCellProliferativeType(DIFFERENTIATED);
        cell_population.GetCellUsingLocationIndex(1)->SetMutationState(p_apc1);
        cell_population.GetCellUsingLocationIndex(2)->SetMutationState(p_apc2);
        cell_population.GetCellUsingLocationIndex(3)->SetMutationState(p_bcat1);
        cell_population.GetCellUsingLocationIndex(4)->AddCellProperty(p_apoptotic_state);

        TS_ASSERT_THROWS_NOTHING(cell_population.CreateOutputFiles(output_directory, false));

        cell_population.WriteResultsToFiles();

        TS_ASSERT_THROWS_NOTHING(cell_population.CloseOutputFiles());

        // Compare output with saved files of what they should look like
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.viznodes     cell_based/test/data/TestCaBasedCellPopulationWriters/results.viznodes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizcelltypes     cell_based/test/data/TestCaBasedCellPopulationWriters/results.vizcelltypes").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.vizancestors     cell_based/test/data/TestCaBasedCellPopulationWriters/results.vizancestors").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellmutationstates.dat     cell_based/test/data/TestCaBasedCellPopulationWriters/cellmutationstates.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellages.dat     cell_based/test/data/TestCaBasedCellPopulationWriters/cellages.dat").c_str()), 0);
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "cellareas.dat     cell_based/test/data/TestCaBasedCellPopulationWriters/cellareas.dat").c_str()), 0);

        // Test the GetCellMutationStateCount function
        std::vector<unsigned> cell_mutation_states = cell_population.GetCellMutationStateCount();
        TS_ASSERT_EQUALS(cell_mutation_states.size(), 4u);
        TS_ASSERT_EQUALS(cell_mutation_states[0], 0u);
        for (unsigned i=1; i<4; i++)
        {
            TS_ASSERT_EQUALS(cell_mutation_states[i], 1u);
        }

        // Test the GetCellProliferativeTypeCount function
        std::vector<unsigned> cell_types = cell_population.rGetCellProliferativeTypeCount();
        TS_ASSERT_EQUALS(cell_types.size(), 3u);
        TS_ASSERT_EQUALS(cell_types[0], 0u);
        TS_ASSERT_EQUALS(cell_types[1], 1u);
        TS_ASSERT_EQUALS(cell_types[2], 4u);

        // For coverage
        TS_ASSERT_THROWS_NOTHING(cell_population.WriteResultsToFiles());

        // Test that the cell population parameters are output correctly
        out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

        // Write cell population parameters to file
        cell_population.OutputCellPopulationParameters(parameter_file);
        parameter_file->close();

        // Compare output with saved files of what they should look like
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "results.parameters    cell_based/test/data/TestCaBasedCellPopulationWriters/results.parameters").c_str()), 0);
    }

    /*
     * Note that this also tests the methods SetOnlyUseNearestNeighboursForDivision()
     * and SetVonNeumannNeighbourhoods(), as well as the AbstractOnLatticeCellPopulation
     * methods SetUpdateNodesInRandomOrder(), GetUpdateNodesInRandomOrder(),
     * GetIterateRandomlyOverUpdateRuleCollection() and SetIterateRandomlyOverUpdateRuleCollection().
     */
    void TestArchiving() throw (Exception)
    {
        FileFinder archive_dir("archive", RelativeTo::ChasteTestOutput);
        std::string archive_file = "CaBasedCellPopulation.arch";
        ArchiveLocationInfo::SetMeshFilename("CaBasedCellPopulation");

        // Archive a cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);

            // Create a CA-based cell population object
            TetrahedralMesh<2,2> mesh;
            mesh.ConstructRectangularMesh(6, 6, true);

            std::vector<unsigned> real_node_indices;
            for (unsigned i=0; i<10; i++)
            {
                real_node_indices.push_back(2*i);
            }

            std::vector<CellPtr> cells;
            CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices);

            AbstractCellPopulation<2>* const p_cell_population = new CaBasedCellPopulation<2>(mesh, cells, real_node_indices);

            // Cells have been given birth times of 0, -1, -2, -3, -4.
            // Loop over them to run to time 0.
            for (AbstractCellPopulation<2>::Iterator cell_iter = p_cell_population->Begin();
                 cell_iter != p_cell_population->End();
                 ++cell_iter)
            {
                cell_iter->ReadyToDivide();
            }

            // Set member variables in order to test that they are archived correctly
            static_cast<CaBasedCellPopulation<2>*>(p_cell_population)->SetOnlyUseNearestNeighboursForDivision(true);
            static_cast<CaBasedCellPopulation<2>*>(p_cell_population)->SetUseVonNeumannNeighbourhoods(true);
            static_cast<CaBasedCellPopulation<2>*>(p_cell_population)->SetUpdateNodesInRandomOrder(false);
            static_cast<CaBasedCellPopulation<2>*>(p_cell_population)->SetIterateRandomlyOverUpdateRuleCollection(true);

            // Add an update rule
            MAKE_PTR_ARGS(AdvectionCaUpdateRule<2>, p_update_rule, (1, 1.0));
            static_cast<CaBasedCellPopulation<2>*>(p_cell_population)->AddUpdateRule(p_update_rule);

            // Create an output archive
            ArchiveOpener<boost::archive::text_oarchive, std::ofstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_oarchive* p_arch = arch_opener.GetCommonArchive();

            // Archive the cell population
            (*p_arch) << static_cast<const SimulationTime&>(*p_simulation_time);
            (*p_arch) << p_cell_population;

            // Tidy up
            SimulationTime::Destroy();
            delete p_cell_population;
        }

        // Restore the cell population
        {
            // Need to set up time
            unsigned num_steps = 10;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, num_steps+1);
            p_simulation_time->IncrementTimeOneStep();

            // Restore the cell population
            AbstractCellPopulation<2>* p_cell_population;

            ArchiveOpener<boost::archive::text_iarchive, std::ifstream> arch_opener(archive_dir, archive_file);
            boost::archive::text_iarchive* p_arch = arch_opener.GetCommonArchive();

            (*p_arch) >> *p_simulation_time;
            (*p_arch) >> p_cell_population;

            // Test that the member variables have been archived correctly
            CaBasedCellPopulation<2>* p_static_population = static_cast<CaBasedCellPopulation<2>*>(p_cell_population);
            TS_ASSERT_EQUALS(p_static_population->GetOnlyUseNearestNeighboursForDivision(), true);
            TS_ASSERT_EQUALS(p_static_population->GetUseVonNeumannNeighbourhoods(), true);
            TS_ASSERT_EQUALS(p_static_population->GetUpdateNodesInRandomOrder(), false);
            TS_ASSERT_EQUALS(p_static_population->GetIterateRandomlyOverUpdateRuleCollection(), true);

            // Test that the update rule has been archived correctly
            std::vector<boost::shared_ptr<AbstractCaUpdateRule<2> > > update_rule_collection = p_static_population->rGetUpdateRuleCollection();
            TS_ASSERT_EQUALS(update_rule_collection.size(), 1u);
            TS_ASSERT_EQUALS((*update_rule_collection[0]).GetIdentifier(), "AdvectionCaUpdateRule-2");

            // Tidy up
            delete p_cell_population;
        }
    }

    void TestMoveCell() throw(Exception)
    {
        // Create a CA-based cell population object
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true);

        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(4);
        real_node_indices.push_back(1);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        CellPtr p_cell = *(cell_population.Begin());
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_cell), 4u);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(4), p_cell);

        // Move the cell at node 4 to node 2
        cell_population.MoveCell(p_cell, 2);

        // Test that the cell was moved correctly
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_cell), 2u);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(2), p_cell);

        // Test that the member variable mIsEmptySite was updates correctly
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            bool expected_empty_site = (i!=1 && i!=2);
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(i), expected_empty_site);
        }
    }

    void TestAddCell() throw(Exception)
    {
        c_vector<double,2> division_vector = zero_vector<double>(2);

        // Create a fully populated CA-based cell population object
        TetrahedralMesh<2,2> full_mesh;
        full_mesh.ConstructRectangularMesh(2, 2, true);

        std::vector<CellPtr> full_cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> full_cells_generator;
        full_cells_generator.GenerateBasic(full_cells, full_mesh.GetNumNodes());

        CaBasedCellPopulation<2> full_cell_population(full_mesh, full_cells);

        // Create another cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
        p_model->SetCellProliferativeType(STEM);
        CellPtr p_cell(new Cell(p_state, p_model));

        // Test that we cannot add the cell to the population
        CellPtr p_full_parent_cell = *(full_cell_population.Begin());
        TS_ASSERT_THROWS_THIS(full_cell_population.AddCell(p_cell, division_vector, p_full_parent_cell),
                              "Cell can not divide as there are no free neighbours at maximum degree in any direction.");

        // Now create a sparsely populated CA-based cell population object
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true);

        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(4);
        real_node_indices.push_back(1);

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // When adding a cell, only try to find empty sites among the parent cell's nearest neighbours
        cell_population.SetOnlyUseNearestNeighboursForDivision(true);

        // Test that the correct exception is thrown when attempting to add the cell without providing a parent cell
        TS_ASSERT_THROWS_THIS(cell_population.AddCell(p_cell, division_vector),
                              "A parent cell must be provided when calling AddCell() on a CaBasedCellPopulation.");

        // Add this cell to the population
        CellPtr p_parent_cell = *(cell_population.Begin());
        cell_population.AddCell(p_cell, division_vector, p_parent_cell);

        // Test that the cell was added correctly
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
        ++cell_iter;
        ++cell_iter;
        TS_ASSERT_EQUALS(*cell_iter, p_cell);

        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_cell), 2u);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(2), p_cell);

        // Test that the member variable mIsEmptySite was updates correctly
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            bool expected_empty_site = (i!=1 && i!=2 && i!=4);
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(i), expected_empty_site);
        }

        // Create another cell
        FixedDurationGenerationBasedCellCycleModel* p_model2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model2->SetCellProliferativeType(STEM);
        CellPtr p_cell2(new Cell(p_state, p_model2));

        // When adding a cell, now try to find empty sites anywhere
        cell_population.SetOnlyUseNearestNeighboursForDivision(false);

        // Add this cell to the population
        cell_population.AddCell(p_cell2, division_vector, p_parent_cell);

        // Test that the cell was added correctly
        ++cell_iter;
        TS_ASSERT_EQUALS(*cell_iter, p_cell2);

        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(p_cell2), 7u);
        TS_ASSERT_EQUALS(cell_population.GetCellUsingLocationIndex(7), p_cell2);

        // Test that the member variable mIsEmptySite was updates correctly
        for (unsigned i=0; i<cell_population.GetNumNodes(); i++)
        {
            bool expected_empty_site = (i!=1 && i!=2 && i!=4 && i!=7);
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(i), expected_empty_site);
        }
    }

    void TestAddCellWithOneEmptySite() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(6, 6, true); // 7*7 nodes

        // Initially we create a cell for every node except node 45
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (i != 45)
            {
                real_node_indices.push_back(i);
            }
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        // Create a cell population
        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        TS_ASSERT(cell_population.IsEmptySite(45));
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 49u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 48u);

        // Test GetNeighbouringNodeIndices() method
        std::set<unsigned> node_20_neighbours = cell_population.GetNeighbouringNodeIndices(20);

        std::set<unsigned> expected_node_20_neighbours;
        expected_node_20_neighbours.insert(12);
        expected_node_20_neighbours.insert(13);
        expected_node_20_neighbours.insert(19);
        expected_node_20_neighbours.insert(26);
        expected_node_20_neighbours.insert(27);

        TS_ASSERT_EQUALS(node_20_neighbours, expected_node_20_neighbours);

        // Create a new cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model_2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_2->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_new_cell(new Cell(p_state, p_model_2));

        // Add new cell to the cell population by dividing the cell at node 24
        CellPtr p_parent_cell = cell_population.GetCellUsingLocationIndex(24);
        cell_population.AddCell(p_new_cell, zero_vector<double>(2), p_parent_cell);

        // Check the number of cells is correct
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 49u);

        // Check each cell corresponds to the correct location index
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();

        // The indices of cells 0 to 30 should be unchanged
        for (unsigned i=0; i<31; i++)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
            ++cell_iter;
        }

        // Cell 31 should have been pushed up to index 38
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 38u);
        ++cell_iter;

        // The indices of cells 32 to 37 should be unchanged
        for (unsigned i=32; i<38; i++)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
            ++cell_iter;
        }

        // Cell 38 should have been pushed up to index 45
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 45u);
        ++cell_iter;

        // The indices of cells 39 to 44 and 46 to 48 should be unchanged
        for (unsigned i=39; i<49; i++)
        {
            if (i != 45)
            {
                TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
                ++cell_iter;
            }
        }

        // Lastly the new cell should be at index 31
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 31u);
    }

    void TestGetFreeNeighbouringNodeIndices1d()
    {
        // Create mesh
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructLinearMesh(4); // 5 nodes

        /*
         * Numbering the nodes as follows:
         *
         *  0------1------2------3------4
         */

        // Create one differentiated cell, initially corresponding to the far left node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        CaBasedCellPopulation<1> lattice_based_cell_population(mesh, cells, location_indices);

        // Now find the neighbouring available sites
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 1u);

        std::set<unsigned> expected_free_neighbouring_sites;
        expected_free_neighbouring_sites.insert(1);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-boundary node 2
        lattice_based_cell_population.MoveCell(*lattice_based_cell_population.Begin(), 2);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(2);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(3);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test eastern end
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 4);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(4);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 1u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(3);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    // Checks that this method correctly returns the first degree, i.e. nearest, neighbours
    void TestGetFirstDegreeNeighbouringNodeIndices() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true); // 3*3 nodes

        // Create a single cell and corresponding location index
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> real_node_indices = CreateSingleLocationIndex(4);

        // Create a cell population
        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        std::set<unsigned> nearest_neighbours = cell_population.GetNthDegreeNeighbouringNodeIndices(4, 1);

        unsigned indices[8] = {7, 6, 3, 0, 1, 2, 5, 8};
        std::set<unsigned> expected_neighbours(indices, indices + 8);

        TS_ASSERT_EQUALS(nearest_neighbours, expected_neighbours);

        std::set<unsigned>::iterator neighbour_iter = nearest_neighbours.begin();
        std::set<unsigned>::iterator expected_iter = expected_neighbours.begin();

        for (unsigned i=0; i<nearest_neighbours.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }

        // Now going to check for a corner node
        std::set<unsigned> nearest_neighbours_2 = cell_population.GetNthDegreeNeighbouringNodeIndices(2, 1);

        unsigned indices_2[3] = {5, 4, 1};
        std::set<unsigned> expected_neighbours_2(indices_2, indices_2 + 3);

        TS_ASSERT_EQUALS(nearest_neighbours_2, expected_neighbours_2);

        neighbour_iter = nearest_neighbours_2.begin();
        expected_iter = expected_neighbours_2.begin();

        for (unsigned i=0; i<nearest_neighbours_2.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }

        // Now checking a boundary node
        std::set<unsigned> nearest_neighbours_3 = cell_population.GetNthDegreeNeighbouringNodeIndices(3, 1);

        unsigned indices_3[5] = {6, 0, 1, 4, 7};
        std::set<unsigned> expected_neighbours_3(indices_3, indices_3 + 5);

        TS_ASSERT_EQUALS(nearest_neighbours_3, expected_neighbours_3);

        neighbour_iter = nearest_neighbours_3.begin();
        expected_iter = expected_neighbours_3.begin();

        for (unsigned i=0; i<nearest_neighbours_3.size(); i++)
        {
            TS_ASSERT_EQUALS(*neighbour_iter, *expected_iter);
            expected_iter++;
            neighbour_iter++;
        }
    }

    void TestGetNthDegreeNeighbouringNodeIndices() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(6, 6, true); // 7*7 nodes

        // Create a single cell and corresponding location index
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> real_node_indices = CreateSingleLocationIndex(29);

        // Create a cell population
        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Check the maximum degree possible in each direction
        std::vector<unsigned> max_degrees = cell_population.GetMaximumDegreeInEachDirection(29);
        TS_ASSERT_EQUALS(max_degrees.size(), 8u);

        unsigned expected_max_degrees[8] = {2, 1, 1, 1, 4, 4, 5, 2};
        for (unsigned i=0; i<max_degrees.size(); i++)
        {
            TS_ASSERT_EQUALS(max_degrees[i], expected_max_degrees[i]);
        }

        // Check the second degree neighbours
        std::set<unsigned> nearest_neighbours_2 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 2);

        unsigned indices_2[5] = {43, 15, 17, 31, 45};
        std::set<unsigned> expected_neighbours_2(indices_2, indices_2 + 5);

        TS_ASSERT_EQUALS(nearest_neighbours_2, expected_neighbours_2);

        // Check the third degree neighbours
        std::set<unsigned> nearest_neighbours_3 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 3);

        unsigned indices_3[3] = {8, 11, 32};
        std::set<unsigned> expected_neighbours_3(indices_3, indices_3 + 3);

        TS_ASSERT_EQUALS(nearest_neighbours_3, expected_neighbours_3);

        // Check the fourth degree neighbours
        std::set<unsigned> nearest_neighbours_4 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 4);

        unsigned indices_4[3] = {1, 5, 33};
        std::set<unsigned> expected_neighbours_4(indices_4, indices_4 + 3);

        TS_ASSERT_EQUALS(nearest_neighbours_4, expected_neighbours_4);

        // Check the fifth degree neighbours
        std::set<unsigned> nearest_neighbours_5 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 5);

        std::set<unsigned> expected_neighbours_5;
        expected_neighbours_5.insert(34);

        TS_ASSERT_EQUALS(nearest_neighbours_5, expected_neighbours_5);

        // Check the sixth degree neighbours
        std::set<unsigned> nearest_neighbours_6 = cell_population.GetNthDegreeNeighbouringNodeIndices(29, 6);
        TS_ASSERT_EQUALS(nearest_neighbours_6.empty(), true);
    }

    void TestCellDivisionWithOneEmptySiteOnlySearchingNearestNeighbours() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(3, 3, true); // 4*4 nodes

        // Initially we create a cell for every node except node 13
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (i != 13)
            {
                real_node_indices.push_back(i);
            }
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        // Create a cell population
        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices, true);
        cell_population.SetOnlyUseNearestNeighboursForDivision(true);

        TS_ASSERT(cell_population.IsEmptySite(13));
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 15u);

        // Create a new cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model_2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_2->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_new_cell(new Cell(p_state, p_model_2));

        // Add new cell to the cell population by dividing the cell at node 5
        CellPtr p_parent_cell = cell_population.GetCellUsingLocationIndex(5);

        // Try to divide the parent - should not be possible as there are no free nearest neighbours

        TS_ASSERT_THROWS_THIS(cell_population.AddCell(p_new_cell, zero_vector<double>(2), p_parent_cell),
                              "Cell can not divide as there are no free neighbours at maximum degree in any direction.");
    }

    void TestCellDivisionWithOneNearestNeighbourEmptySiteOnlySearchingNearestNeighbours() throw(Exception)
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(3, 3, true); // 4*4 nodes

        // Initially we create a cell for every node except node 9
        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            if (i!=9)
            {
                real_node_indices.push_back(i);
            }
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);

        // Create a cell population
        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices, true); // Only searching the nearest neighbours (degree=1)

        TS_ASSERT(cell_population.IsEmptySite(9));
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 15u);

        // Create a new cell
        MAKE_PTR(WildTypeCellMutationState, p_state);
        FixedDurationGenerationBasedCellCycleModel* p_model_2 = new FixedDurationGenerationBasedCellCycleModel();
        p_model_2->SetCellProliferativeType(DIFFERENTIATED);
        CellPtr p_new_cell(new Cell(p_state, p_model_2));

        // Add new cell to the cell population by dividing the cell at node 5
        CellPtr p_parent_cell = cell_population.GetCellUsingLocationIndex(5);
        cell_population.AddCell(p_new_cell, zero_vector<double>(2), p_parent_cell);

        // Check the number of cells is correct
        TS_ASSERT_EQUALS(cell_population.rGetCells().size(), 16u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 16u);

        // Check each cell corresponds to the correct location index
        AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();

        // The indices of cells 0 to 8 and 10 to 15 should be unchanged
        for (unsigned i=0; i<16; i++)
        {
            if (i != 9)
            {
                TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), i);
                ++cell_iter;
            }
        }

        // The new cell should be at node index 9
        TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 9u);
    }

    void TestGetFreeNeighbouringNodeIndices2d()
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 3, true); // 3*4 nodes

        /*
         * Numbering the nodes as follows:
         *
         *     9---10---11
         *     |    |    |
         *     6----7----8
         *     |    |    |
         *     3----4----5
         *     |    |    |
         *     0----1----2
         */

        // Create one cell, initially corresponding to the corner node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        CaBasedCellPopulation<2> lattice_based_cell_population(mesh, cells, location_indices);

        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        unsigned indices[3] = {1, 3, 4};
        std::set<unsigned> expected_free_neighbouring_sites(indices, indices + 3);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test east side node (node 5)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 5);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(5);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();
        expected_free_neighbouring_sites.insert(1);
        expected_free_neighbouring_sites.insert(2);
        expected_free_neighbouring_sites.insert(4);
        expected_free_neighbouring_sites.insert(7);
        expected_free_neighbouring_sites.insert(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-boundary node (node 7)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 7);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(7);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 8u);

        expected_free_neighbouring_sites.clear();
        for (unsigned i=0; i<12; i++)
        {
            if (i > 2 && i != 7)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test north boundary node (node 10)

        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 10);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(10);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        expected_free_neighbouring_sites.clear();

        for (unsigned i=0; i<12; i++)
        {
            if (i > 5 && i != 10)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestGetFreeNeighbouringNodeIndices3d()
    {
        // Create mesh
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2, 3, 4); // 3*4*5 nodes

        /*
         *  Numbering the nodes as follows:
         *
         *      Bottom layer                layer 2                      layer 3
         *     9---10---11               21---22---23                  33---34---35
         *     |    |    |               |    |    |                   |    |    |
         *     6----7----8               18---19---20                  30---31---32
         *     |    |    |               |    |    |                   |    |    |
         *     3----4----5               15---16---17                  27---28---29
         *     |    |    |               |    |    |                   |    |    |
         *     0----1----2               12---13---14                  24---25---26
         *
         *        layer 4                 Top layer
         *     45---46---47             57---58---59
         *     |     |    |             |     |    |
         *     42---43---44             54---55---56
         *     |     |    |             |     |    |
         *     39---40---41             51---52---53
         *     |     |    |             |     |    |
         *     36---37---38             48---49---50
         */

        // Create one cell, initially corresponding to the corner node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        CaBasedCellPopulation<3> lattice_based_cell_population(mesh, cells, location_indices);

        // Test bottom left corner node (node 0)
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 7u);

        std::set<unsigned> expected_free_neighbouring_sites;
        for (unsigned i=0; i<17; i++)
        {
            if (i==1 || i==3 || i==4 || i==12 || i==13 || i==15 || i==16)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test west side node (node 30)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 30);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(30);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 17u);

        expected_free_neighbouring_sites.clear();

        for (unsigned i=0; i<47; i++)
        {
            if (i==15 || i==16 || i==18 || i==19 || i==21 || i==22 || i==27 || i==28 || i==31 || i==33 || i==34 || i==39 || i==40 || i==42 || i==43 || i==45 || i==46)
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test east side node (node 44)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 44);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(44);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 17u);

        expected_free_neighbouring_sites.clear();
        for (unsigned i=0; i<60; i++)
        {
            if (   i==28 || i==29 || i==31 || i==32 || i==34 || i==35 || i==40 || i==41 || i==43 || i==46 || i==47
                || i==52 || i==53 || i==55 || i==56 || i==58 || i==59 )
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test top layer node (node 52)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 52);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(52);

        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 17u);

        expected_free_neighbouring_sites.clear();
        for (unsigned i=0; i<57; i++)
        {
            if (   i==36 || i==37 || i==38 || i==39 || i==40 || i==41 || i==42 || i==43 || i==44 || i==48 || i==49
                || i==50 || i==51 || i==53 || i==54 || i==55 || i==56 )
            {
                expected_free_neighbouring_sites.insert(i);
            }
        }
        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);
    }

    void TestGetFreeVonNeumannNeighbouringNodeIndices2d()
    {
        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 3, true); // 3*4 nodes

        // Create one cell, initially corresponding to the bottom left node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> cell_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        CaBasedCellPopulation<2> lattice_based_cell_population(mesh, cells, cell_indices);

        // Set the neighbourhoods to be of von Neumann type
        lattice_based_cell_population.SetUseVonNeumannNeighbourhoods(true);

        /*
         * Numbering the nodes as follows:
         *
         *     9---10---11
         *     |    |    |
         *     6----7----8
         *     |    |    |
         *     3----4----5
         *     |    |    |
         *     0----1----2
         */

        // Test bottom left node
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        unsigned indices[2] = {1, 3};
        std::set<unsigned> expected_free_neighbouring_sites(indices, indices + 2);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test non-corner bottom node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 1);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(1);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        unsigned indices_2[3] = {0, 2, 4};
        std::set<unsigned> expected_free_neighbouring_sites_2(indices_2, indices_2 + 3);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_2);

        // Test bottom right node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 2);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(2);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        unsigned indices_3[2] = {1, 5};
        std::set<unsigned> expected_free_neighbouring_sites_3(indices_3, indices_3 + 2);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_3);

        // Test non-corner left nodes
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 3);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(3);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        unsigned indices_4[3] = {0, 4, 6};
        std::set<unsigned> expected_free_neighbouring_sites_4(indices_4, indices_4 + 3);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_4);

        // Test interior node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 4);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(4);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 4u);

        unsigned indices_5[4] = {1, 3, 5, 7};
        std::set<unsigned> expected_free_neighbouring_sites_5(indices_5, indices_5 + 4);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_5);

        // Test non-corner right nodes
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 5);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(5);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        unsigned indices_6[3] = {2, 4, 8};
        std::set<unsigned> expected_free_neighbouring_sites_6(indices_6, indices_6 + 3);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_6);

        // Test top left node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 9);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(9);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        unsigned indices_7[2] = {6, 10};
        std::set<unsigned> expected_free_neighbouring_sites_7(indices_7, indices_7 + 2);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_7);

        // Test non-corner top node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 10);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(10);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        unsigned indices_8[3] = {7, 9, 11};
        std::set<unsigned> expected_free_neighbouring_sites_8(indices_8, indices_8 + 3);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_8);

        // Test top right node
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 8);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(8);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        unsigned indices_9[2] = {5, 7};
        std::set<unsigned> expected_free_neighbouring_sites_9(indices_9, indices_9 + 2);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_9);

        // Now check the case where there is more than one cell in the cell population
        TetrahedralMesh<2,2> mesh2;
        mesh2.ConstructRectangularMesh(2, 2, true); // 3*3 nodes
        std::vector<CellPtr> cells2;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        for (unsigned i=0; i<3; i++)
        {
            FixedDurationGenerationBasedCellCycleModel* p_model = new FixedDurationGenerationBasedCellCycleModel();
            p_model->SetCellProliferativeType(DIFFERENTIATED);
            CellPtr p_cell(new Cell(p_state, p_model));
            cells2.push_back(p_cell);
        }

        std::vector<unsigned> location_indices;
        location_indices.push_back(1);
        location_indices.push_back(4);
        location_indices.push_back(5);

        // Create a cell population
        CaBasedCellPopulation<2> lattice_based_cell_population2(mesh2, cells2, location_indices);

        // Set the neighbourhoods to be of von Neumann type
        lattice_based_cell_population2.SetUseVonNeumannNeighbourhoods(true);

        free_neighbouring_sites = lattice_based_cell_population2.GetFreeNeighbouringNodeIndices(1);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 2u);

        unsigned indices_10[2] = {0, 2};
        std::set<unsigned> expected_free_neighbouring_sites_10(indices_10, indices_10 + 2);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_10);
    }

    void TestGetFreeVonNeumannNeighbouringNodeIndices3d()
    {
        // Create mesh
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructCuboid(2, 3, 4); // 3*4*5 nodes

        /*
         *  Numbering the nodes as follows:
         *
         *     Bottom layer                layer 2                       layer 3
         *     9---10---11               21---22---23                 33---34---35
         *     |    |    |               |     |    |                  |    |    |
         *     6----7----8               18---19---20                 30---31---32
         *     |    |    |               |     |    |                  |    |    |
         *     3----4----5               15---16---17                 27---28---29
         *     |    |    |               |     |    |                  |    |    |
         *     0----1----2               12---13---14                 24---25---26
         *
         *       layer 4                  Top layer
         *     45---46---47              57---58---59
         *     |     |    |              |     |    |
         *     42---43---44              54---55---56
         *     |     |    |              |     |    |
         *     39---40---41              51---52---53
         *     |     |    |              |     |    |
         *     36---37---38              48---49---50
         */

        // Create one cell, initially corresponding to the corner node
        std::vector<CellPtr> cells = CreateSingleCellPtr();
        std::vector<unsigned> location_indices = CreateSingleLocationIndex(0);

        // Create a cell population
        CaBasedCellPopulation<3> lattice_based_cell_population(mesh, cells, location_indices);

        // Set von Neumann neighbourhoods
        lattice_based_cell_population.SetUseVonNeumannNeighbourhoods(true);

        // Test bottom left corner node (node 0)
        std::set<unsigned> free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(0);

        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 3u);

        unsigned indices[3] = {1, 3, 12};
        std::set<unsigned> expected_free_neighbouring_sites(indices, indices + 3);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites);

        // Test west side node (node 30)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 30);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(30);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        unsigned indices_2[5] = {27, 31, 33, 18, 42};
        std::set<unsigned> expected_free_neighbouring_sites_2(indices_2, indices_2 + 5);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_2);

        // Test east side node (node 44)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 44);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(44);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        unsigned indices_3[5] = {41, 43, 47, 32, 56};
        std::set<unsigned> expected_free_neighbouring_sites_3(indices_3, indices_3 + 5);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_3);

        // Test top layer node (node 52)
        lattice_based_cell_population.MoveCell(*(lattice_based_cell_population.Begin()), 52);
        free_neighbouring_sites = lattice_based_cell_population.GetFreeNeighbouringNodeIndices(52);
        TS_ASSERT_EQUALS(free_neighbouring_sites.size(), 5u);

        unsigned indices_4[5] = {40, 49, 51, 53, 55};
        std::set<unsigned> expected_free_neighbouring_sites_4(indices_4, indices_4 + 5);

        TS_ASSERT_EQUALS(free_neighbouring_sites, expected_free_neighbouring_sites_4);
    }

    void TestRemoveDeadCellsAndUpdate()
    {
        // Set up SimulationTime singleton
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(10.0, 1);

        // Create mesh
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(3, 3, true); // 4*4 nodes

        // Create location indices
        std::vector<unsigned> real_node_indices;
        real_node_indices.push_back(0);
        real_node_indices.push_back(4);
        real_node_indices.push_back(7);
        real_node_indices.push_back(8);
        real_node_indices.push_back(12);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices, DIFFERENTIATED);
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->SetBirthTime(-1.0);
        }

        // Set the last cell to start apoptosis
        cells[2]->StartApoptosis();

        // Create a cell population
        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        // Test we have the right numbers of nodes and cells
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 5u);

        p_simulation_time->IncrementTimeOneStep();

        unsigned num_removed = cell_population.RemoveDeadCells();

        // Test that one cell has been removed
        TS_ASSERT_EQUALS(num_removed, 1u);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 4u);

        // Test that no nodes have been removed
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 16u);

        // Test that each cell's node index has been correctly updated
        unsigned index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_population.GetLocationIndexUsingCell(*cell_iter), 4*index);
            index++;
        }

        // Test that Update() throws no errors
        TS_ASSERT_THROWS_NOTHING(cell_population.Update());
    }

    void TestUpdateCellLocations() throw(Exception)
    {
        // Create a CA-based cell population object with five
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructRectangularMesh(2, 2, true);
        unsigned num_nodes = mesh.GetNumNodes();

        std::vector<unsigned> real_node_indices;
        for (unsigned i=0; i<5; i++)
        {
            real_node_indices.push_back(i);
        }

        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, real_node_indices.size(), real_node_indices);

        CaBasedCellPopulation<2> cell_population(mesh, cells, real_node_indices);

        /*
         * Numbering the nodes as follows:
         *
         *     6----7----8
         *     |    |    |
         *     3----4----5
         *     |    |    |
         *     0----1----2
         *
         * At this point cells should be located at sites 0, 1, 2, 3, 4
         */
        for (unsigned i=0; i<num_nodes; i++)
        {
            bool expected_empty_site = (i > 4);
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(i), expected_empty_site);
        }

        // Add an advection update rule with a northerly direction and unit speed
        MAKE_PTR_ARGS(AdvectionCaUpdateRule<2>, p_update_rule, (0, 1.0));
        cell_population.AddUpdateRule(p_update_rule);

        // For coverage, iterate randomly over the update rule collection
        cell_population.SetIterateRandomlyOverUpdateRuleCollection(true);

        // Update the cell positions once, using unit timestep
        cell_population.UpdateCellLocations(1.0);

        // At this point cells should be located at sites 3, 4, 5, 6, 7
        for (unsigned i=0; i<num_nodes; i++)
        {
            bool expected_empty_site = (i < 3 || i > 7);
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(i), expected_empty_site);
        }

        // For coverage, update nodes in a random order
        cell_population.SetUpdateNodesInRandomOrder(false);

        // Update the cell positions again, using unit timestep
        cell_population.UpdateCellLocations(1.0);

        // At this point cells should be located at sites 3, 4, 6, 7, 8
        for (unsigned i=0; i<num_nodes; i++)
        {
            bool expected_empty_site = (i < 3 || i == 5 || i > 8);
            TS_ASSERT_EQUALS(cell_population.IsEmptySite(i), expected_empty_site);
        }
    }
};

#endif /*TESTLATTICEBASEDCELLPOPULATION_HPP_*/
