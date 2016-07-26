/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef TESTONLATTICESIMULATIONWITHCABASEDCELLPOPULATION_HPP_
#define TESTONLATTICESIMULATIONWITHCABASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AveragedSourcePde.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "OnLatticeSimulation.hpp"
#include "Warnings.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "SmartPointers.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "RandomCaSwitchingUpdateRule.hpp"
#include "CellIdWriter.hpp"
#include "ApoptoticCellKiller.hpp"
#include "ApoptoticCellProperty.hpp"

#include "CellProliferativeTypesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"


class TestOnLatticeSimulationWithCaBasedCellPopulation : public AbstractCellBasedWithTimingsTestSuite
{
    void RandomlyLabelCells(std::vector<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
    {
        for (unsigned i=0; i<rCells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
            {
                rCells[i]->AddCellProperty(pLabel);
            }
        }
    }

public:

    void TestOnLatticeSimulationExceptions()
    {
        EXIT_IF_PARALLEL;

        // Create a simple tetrahedral mesh
        HoneycombMeshGenerator generator(3, 3, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        TS_ASSERT_THROWS_THIS(OnLatticeSimulation<2> simulator(node_based_cell_population),
            "OnLatticeSimulations require a subclass of AbstractOnLatticeCellPopulation.");
    }

    void TestCaSingleCellRandomMovement() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Specify the timestep and size of domain to let us calculate the probabilities of movement
        double delta_t = 1;
        double diffusion_parameter = 0.1;
        unsigned num_runs = 2000;
        unsigned location_of_cell[9] = {0,0,0,0,0,0,0,0,0};

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(3, 0, 0, 3, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 1, p_diff_type);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        location_indices.push_back(4);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        std::string output_directory = "TestCaSingleCellRandomMovement";
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(delta_t);
        simulator.SetEndTime(delta_t);

        // Add update rule
        MAKE_PTR(DiffusionCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(diffusion_parameter);
        simulator.AddCaUpdateRule(p_diffusion_update_rule);

        for (unsigned i=1; i<=num_runs; i++)
        {
            simulator.SetEndTime(delta_t*i);

            // Run simulation
            simulator.Solve();

            TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 1u);
            AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();

            unsigned cell_location = simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);
            TS_ASSERT_LESS_THAN(cell_location, 9u);

            location_of_cell[cell_location]++;

            // Reset the position of the cell
            simulator.rGetCellPopulation().MoveCellInLocationMap(*cell_iter, cell_location, 4u);

            TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter), 4u);
        }

        // Check that we still have only one cell
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 1u);
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

        ///\todo Check that the cell is moving correctly
        double probability_of_occupation[9];
        for (unsigned i=0; i<9; i++)
        {
            probability_of_occupation[i] = (double) location_of_cell[i]/(double) num_runs;
        }

        // Note that these simulations are stochastic and so the tolerances are relatively loose
        TS_ASSERT_DELTA(probability_of_occupation[0], diffusion_parameter*delta_t/4.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[1], diffusion_parameter*delta_t/2.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[2], diffusion_parameter*delta_t/4.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[3], diffusion_parameter*delta_t/2.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[4], 1.0 - 3.0*diffusion_parameter*delta_t, 2e-2);
        TS_ASSERT_DELTA(probability_of_occupation[5], diffusion_parameter*delta_t/2.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[6], diffusion_parameter*delta_t/4.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[7], diffusion_parameter*delta_t/2.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[8], diffusion_parameter*delta_t/4.0, 1e-2);

        // For coverage
        simulator.RemoveAllCaUpdateRules();
    }

    void TestCaMonolayerWithBirth() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 2, p_stem_type);

        // Specify where the cells lie
        std::vector<unsigned> location_indices;
        location_indices.push_back(50);
        location_indices.push_back(51);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        std::string output_directory = "TestCaMonolayerWithBirth";
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(0.1);
        simulator.SetEndTime(40);

        // Adding update rule(s)
        MAKE_PTR(DiffusionCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddCaUpdateRule(p_diffusion_update_rule);
        simulator.SetOutputDivisionLocations(true);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 12u); ///\todo #2066 Check this!

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 10u); ///\todo #2066 Check this!
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

        // Now remove the update rules and check that only birth happens when the simulator runs again
        simulator.RemoveAllPottsUpdateRules();
        simulator.SetEndTime(50);
        simulator.Solve();

        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 20u);
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 22u);

#ifdef CHASTE_VTK
        // Test that the VTK writer has produced a file
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Initial condition file
        FileFinder vtk_file(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_from_time_0/results_400.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());
 #endif //CHASTE_VTK
    }

    void TestCaMonolayerWithDeath() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Reset the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_diff_type);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<p_mesh->GetNumNodes(); index++)
        {
            location_indices.push_back(index);
        }
        TS_ASSERT_EQUALS(location_indices.size(),p_mesh->GetNumNodes());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaMonolayerWithDeath");
        simulator.SetDt(0.1);
        simulator.SetEndTime(0.1); // only one step as we only care about cells being killed

        // No movement rule as only care about cell death

        // Add a cell killer that will kill all cells in the top half of the domain
        c_vector<double,2> point = zero_vector<double>(2);
        point[1] = 4.5;
        c_vector<double,2> normal = unit_vector<double>(2,1);
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point, normal)); // v>4.5
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 50u);

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 50u);

        // Check that cells above y=5.5 (i.e. above index 50) have been killed and removed
        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            if (i < 50)
            {
                TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetCellUsingLocationIndex(i)->GetCellId(),i);
            }
            else
            {
                TS_ASSERT_THROWS_THIS(simulator.rGetCellPopulation().GetCellUsingLocationIndex(i),
                    "Location index input argument does not correspond to a Cell");
            }
        }
    }

    void  TestCaMonolayerWithApoptoticCellKiller() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Reset the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_diff_type);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<p_mesh->GetNumNodes(); index++)
        {
            location_indices.push_back(index);
        }
        TS_ASSERT_EQUALS(location_indices.size(),p_mesh->GetNumNodes());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.SetOutputResultsForChasteVisualizer(false);
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellProliferativePhasesWriter>();

        // update all cells in population to prescribe an initial oxygen concentration and apoptosis time
        std::list<CellPtr> cells2 = cell_population.rGetCells();
        std::list<CellPtr>::iterator it;
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        for (it = cells2.begin(); it != cells2.end(); ++it)
        {
            (*it)->SetApoptosisTime(3);
            double random_number = p_gen->ranf();
            if (random_number < 0.5)
            {
                (*it)->AddCellProperty(CellPropertyRegistry::Instance()->Get<ApoptoticCellProperty>());
            }
        }

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaMonolayerWithApoptoticCellKiller");
        simulator.SetDt(1);
        simulator.SetEndTime(5); // only one step as we only care about cells being killed

        // No movement rule as only care about cell death

        // Add a cell killer that will kill all apoptotic cells
        MAKE_PTR_ARGS(ApoptoticCellKiller<2>, p_killer, (&cell_population)); // v>4.5
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT(simulator.GetNumDeaths() > 0);

    }

    void TestCaMonolayerWithRandomSwitching() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Reset the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(4, 0, 0, 4, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 8, p_diff_type);

        // Specify some cells on the bottom row
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<8u; index++)
        {
            location_indices.push_back(index);
        }

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellIdWriter>();

        // Coverage of the case where AbstractOnLatticeCellPopulation::mIterateRandomlyOverUpdateRuleCollection is true
        cell_population.SetIterateRandomlyOverUpdateRuleCollection(true);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaMonolayerWithRandomSwitching");
        simulator.SetDt(0.1);
        simulator.SetEndTime(0.4); // 4 steps so only 7 switches happens (nodes 7-10, 2-5,1-6,5-9,9-14,11-14,10-6)

        // Add switching Update Rule
        MAKE_PTR(RandomCaSwitchingUpdateRule<2u>, p_switching_update_rule);
        p_switching_update_rule->SetSwitchingParameter(1.0);
        simulator.AddCaSwitchingUpdateRule(p_switching_update_rule);

        // Coverage of remove method
        simulator.RemoveAllCaSwitchingUpdateRules();
        simulator.AddCaSwitchingUpdateRule(p_switching_update_rule);

        // Run simulation
        simulator.Solve();

        // Check the number of cells is still the same
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 8u);

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

        // Loop over the cells and check their new positions
        unsigned cell_locations[8]= {0u, 10u, 11u, 3u, 4u, 2u, 1u, 6u};
        unsigned index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
             cell_iter != simulator.rGetCellPopulation().End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter), cell_locations[index]);

            index++;
        }
    }

    ///\todo #2713 check this test makes sense
    void TestCaMonolayerWithRandomSwitchingAndUpdateNodesInNonRandomOrder() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Set up a cell population identically to the previous test
        CellId::ResetMaxCellId();

        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_diff_type);

        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<p_mesh->GetNumNodes(); index++)
        {
            location_indices.push_back(index);
        }
        TS_ASSERT_EQUALS(location_indices.size(),p_mesh->GetNumNodes());

        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellWriter<CellIdWriter>();

        // Coverage of the case where AbstractOnLatticeCellPopulation::mUpdateNodesInRandomOrder is false
        cell_population.SetUpdateNodesInRandomOrder(false);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaMonolayerWithRandomSwitching");
        simulator.SetDt(0.1);
        simulator.SetEndTime(0.2);

        // Add switching update rule
        MAKE_PTR(RandomCaSwitchingUpdateRule<2u>, p_switching_update_rule);
        p_switching_update_rule->SetSwitchingParameter(0.5);
        simulator.AddCaSwitchingUpdateRule(p_switching_update_rule);

        // Run simulation
        simulator.Solve();

        // Check the number of cells is still the same
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 25u);

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

        // Loop over the cells and check their new positions
        unsigned node_index = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();
             cell_iter != simulator.rGetCellPopulation().End();
             ++cell_iter)
        {
            if (node_index == 3)
            {
                TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter), 7u);
            }
            else if (node_index == 7)
            {
                TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter), 3u);
            }
            else if (node_index == 17)
            {
                TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter), 21u);
            }
            else if (node_index == 21)
            {
                TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter), 17u);
            }
            else
            {
                TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter), node_index);
            }
            node_index++;
        }
    }

    /*
     * RandomMovement has been tested in TestCaSingleCellRandomMovement for one cell per lattice site.
     * This test is just to ensure that the above test works when there are multiple cells per lattice site.
     */
    void TestCaMultipleCellsRandomMovement() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 40, p_diff_type);

        // Specify where cells lie: four cells in each of the first ten sites
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<10; index++)
        {
            location_indices.push_back(index);
            location_indices.push_back(index);
            location_indices.push_back(index);
            location_indices.push_back(index);
        }

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 4);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        std::string output_directory = "TestCaMultipleCellRandomMovement";
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1);
        simulator.SetEndTime(10);

        // Add update rule
        MAKE_PTR(DiffusionCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.1);
        simulator.AddCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 40u);
    }

    void TestCaMultipleCellsRandomMovementIn3d() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 3D PottsMesh
        PottsMeshGenerator<3> generator(10, 0, 0, 10, 0, 0, 10, 0, 0);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 40, p_diff_type);

        // Specify where cells lie: four cells in each of the first ten sites
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<10u; index++)
        {
            location_indices.push_back(index);
            location_indices.push_back(index);
            location_indices.push_back(index);
            location_indices.push_back(index);
        }

        // Create cell population
        CaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices, 4u);

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        std::string output_directory = "TestCaMultipleCellRandomMovementIn3d";
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1);
        simulator.SetEndTime(100);

        // Add update rule
        MAKE_PTR(DiffusionCaUpdateRule<3>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.1);
        simulator.AddCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 40u);
    }

    /*
     * Cellular birth has been tested in TestCaSingleCellWithBirth for one cell per lattice site.
     * This test adds to the above by further testing cellular birth considering multiple cells per lattice site.
     * A  two-lattice mesh was created and only one lattice had free space to add one daughter cell.
     */
    void TestMultipleCellsPerLatticeSiteWithBirth() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(2, 0, 0, 1, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 3, p_stem_type);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        location_indices.push_back(0);
        location_indices.push_back(0);
        location_indices.push_back(1);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 2);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        std::string output_directory = "TestMultipleCellsPerLatticeSiteWithBirth";
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(0.1);
        simulator.SetEndTime(40);

        // Add update rule
        MAKE_PTR(DiffusionCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 1u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

#ifdef CHASTE_VTK
        // Test that VTK writer has produced a file
        OutputFileHandler output_file_handler(output_directory, false);
        std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

        // Initial condition file
        FileFinder vtk_file(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file.Exists());

        // Final file
        FileFinder vtk_file2(results_dir + "results_from_time_0/results_400.vtu", RelativeTo::Absolute);
        TS_ASSERT(vtk_file2.Exists());
#endif //CHASTE_VTK
    }

    /*
     * Cellular death has been tested in TestCaMonolayerWithDeath for one cell per lattice site.
     * This test is just to ensure that the above test works when there are multiple cells per lattice site.
     */
    void TestMultipleCellsPerLatticeSiteWithDeath() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Reset the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 2*p_mesh->GetNumNodes(), p_diff_type);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<p_mesh->GetNumNodes(); index++)
        {
          //adding two cells per lattice site
          location_indices.push_back(index);
          location_indices.push_back(index);
        }
        TS_ASSERT_EQUALS(location_indices.size(),2*p_mesh->GetNumNodes());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 2);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCellsPerLatticeSiteWithDeath");
        simulator.SetDt(0.1);
        simulator.SetEndTime(0.1); //only one step as only care about cells being killed

        // No movement rule, as we only care about cell death

        // Add a cell killer that will kill all cells in the top half of the domain
        c_vector<double,2> point = zero_vector<double>(2);
        point[1] = 4.5;
        c_vector<double,2> normal = unit_vector<double>(2,1);
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point, normal)); // v>4.5
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 100u);

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 100u);

        // Check cells above y=5.5 (i.e. above index 50) have been killed and removed.
        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            if (i < 50)
            {
                TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetCellsUsingLocationIndex(i).size(), 2u);
            }
            else
            {
                TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetCellsUsingLocationIndex(i).size(), 0u);
            }
        }
    }

    void TestStandardResultForArchivingTestsBelow() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 10, p_diff_type);

        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<10; index++)
        {
            location_indices.push_back(index);
        }

        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 4);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOnLatticeSimulationWithCaBasedCellPopulationStandardResult");
        simulator.SetDt(1);
        simulator.SetEndTime(20);

        // Add update rule
        MAKE_PTR(DiffusionCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.1);
        simulator.AddCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check some results
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 10u);

        CellPtr p_cell = *(simulator.rGetCellPopulation().Begin());
        c_vector<double, 2> cell_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(p_cell);
        TS_ASSERT_DELTA(cell_location[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(cell_location[1], 3.0, 1e-4);
    }

    void TestSave() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 10, p_diff_type);

        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<10; index++)
        {
            location_indices.push_back(index);
        }

        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 4);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOnLatticeSimulationWithCaBasedCellPopulationSaveAndLoad");
        simulator.SetDt(1);
        simulator.SetEndTime(10);

        // Add update rule
        MAKE_PTR(DiffusionCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.1);
        simulator.AddCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();

        // Save the results
        CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Save(&simulator);
    }

    ///\todo make the following test pass - see #2066 and comment in OnLatticeSimulation::UpdateCellPopulation()
    void DONOTTestLoad() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Load the simulation from the TestSave() method above and run it from time 10 to 15
        OnLatticeSimulation<2>* p_simulator1;
        p_simulator1 = CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Load("TestOnLatticeSimulationWithCaBasedCellPopulationSaveAndLoad", 10.0);

        p_simulator1->SetEndTime(15);
        p_simulator1->Solve();

        // Save, then reload and run from time 15 to 20
        CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Save(p_simulator1);
        OnLatticeSimulation<2>* p_simulator2
            = CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Load("TestOnLatticeSimulationWithCaBasedCellPopulationSaveAndLoad", 15.0);

        p_simulator2->SetEndTime(20);
        p_simulator2->Solve();

        // These results are from time 20 in TestStandardResultForArchivingTestsBelow()
        TS_ASSERT_EQUALS(p_simulator2->rGetCellPopulation().GetNumRealCells(), 10u);

        CellPtr p_cell = *(p_simulator2->rGetCellPopulation().Begin());
        c_vector<double, 2> cell_location = static_cast <CaBasedCellPopulation<2>*>(&p_simulator2->rGetCellPopulation())->GetLocationOfCellCentre(p_cell);
        TS_ASSERT_DELTA(cell_location[0], 1.0, 1e-4);
        TS_ASSERT_DELTA(cell_location[1], 0.0, 1e-4);

        // Tidy up
        delete p_simulator1;
        delete p_simulator2;
    }
};

#endif /*TESTONLATTICESIMULATIONWITHCABASEDCELLPOPULATION_HPP_*/
