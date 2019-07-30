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

#ifndef TESTCRYPTSTATISTICS_HPP_
#define TESTCRYPTSTATISTICS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"
#include "SimpleDataWriter.hpp"
#include "CryptStatistics.hpp"
#include "CryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "SloughingCellKiller.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "NumericFileComparison.hpp"
#include "CellMutationStatesCountWriter.hpp"
// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * Note that all these tests call setUp() and tearDown() before running,
 * so if you copy them into a new test suite be sure to copy these methods
 * too.
 */
class TestCryptStatistics : public CxxTest::TestSuite
{
private:

    void setUp()
    {
        // Initialise singleton classes
        SimulationTime::Instance()->SetStartTime(0.0);
    }
    void tearDown()
    {
        // Clear up singleton classes
        SimulationTime::Destroy();
    }

public:

    void TestGetSection()
    {
        double crypt_length = 22.0;

        // Create mesh
        unsigned cells_across = 3;
        unsigned cells_up = 3;
        unsigned thickness_of_ghost_layer = 0;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);// true = mature cells

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);
        crypt.InitialiseCells(); // must be called explicitly as there is no simulation

        CryptStatistics crypt_statistics(crypt);

        std::vector<CellPtr> test_section = crypt_statistics.GetCryptSection(sqrt(3.0),0.5,1.5);

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section.size(), 6u);

        unsigned expected_indices[6] = {0,1,3,4,7,8};

        for (unsigned i=0; i<test_section.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section[i]), expected_indices[i]);
        }

        // Test that we get a valid section when the x-values are the same
        std::vector<CellPtr> test_section_vertical = crypt_statistics.GetCryptSection(sqrt(3.0),0.5,0.5);

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section_vertical.size(), 5u);

        unsigned expected_indices_vertical[6] = {0,1,3,6,7};

        for (unsigned i=0; i<test_section_vertical.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section_vertical[i]), expected_indices_vertical[i]);
        }

        std::vector<CellPtr> test_section_periodic = crypt_statistics.GetCryptSectionPeriodic(sqrt(3.0),0.5,2.5);

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic.size(), 6u);

        unsigned expected_indices_periodic[6] = {0,1,3,5,6,8};

        for (unsigned i=0; i<test_section_periodic.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section_periodic[i]), expected_indices_periodic[i]);
        }

        std::vector<CellPtr> test_section_periodic_2 = crypt_statistics.GetCryptSectionPeriodic(sqrt(3.0),2.5,0.5);

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic_2.size(), 6u);

        unsigned expected_indices_periodic_2[6] = {0,2,3,5,6,7};

        for (unsigned i=0; i<test_section_periodic_2.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section_periodic_2[i]), expected_indices_periodic_2[i]);
        }

        // Test an overwritten method
        std::vector<CellPtr> test_section_periodic_3 = crypt_statistics.GetCryptSectionPeriodic(crypt_length + 2.0);

        // Test the cells are correct
        TS_ASSERT_EQUALS(test_section_periodic_3.size(), 3u);
        unsigned expected_indices_periodic_3[6] = {2,4,8};

        for (unsigned i=0; i<test_section_periodic_3.size(); i++)
        {
            TS_ASSERT_EQUALS(crypt.GetLocationIndexUsingCell(test_section_periodic_3[i]), expected_indices_periodic_3[i]);
        }
    }

    void TestMakeMeinekeGraphs()
    {
        // Specify output directory
        std::string output_directory = "MakeMeinekeGraphs";

        // Create mesh
        double crypt_length = 22.0;
        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;
        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> temp_cells;
        CryptCellsGenerator<UniformG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(temp_cells, p_mesh, std::vector<unsigned>(), true, 0.3, 2.0, 3.0, 4.0, true);

        // This awkward way of setting up the cells is a result of #430
        std::vector<CellPtr> cells;
        for (unsigned i=0; i<location_indices.size(); i++)
        {
            cells.push_back(temp_cells[location_indices[i]]);
        }

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices, false, 30.0); // Last parameter adjusts Ghost spring stiffness in line with the linear_force later on

        // Set cell population to output cell types
        crypt.AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

        // Create and configure simulation
        CryptSimulation2d simulator(crypt, false, false);
        simulator.SetOutputDirectory(output_directory);

        double time_of_each_run = simulator.GetDt(); // for each run
        simulator.SetEndTime(time_of_each_run);

        // Create a force law and cell killer and pass then to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(30.0); // normally 15.0;
        simulator.AddForce(p_linear_force);
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&simulator.rGetCellPopulation(), crypt_length));
        simulator.AddCellKiller(p_killer);

        // Specify unusual set up
        simulator.UseJiggledBottomCells();

        // Test CryptStatistics::GetCryptSectionPeriodic() by labelling a column of cells...
        CryptStatistics crypt_statistics(crypt);
        std::vector<CellPtr> test_section = crypt_statistics.GetCryptSectionPeriodic(crypt_length + 2.0, 8.0, 8.0);
        for (unsigned i=0; i<test_section.size(); i++)
        {
            test_section[i]->AddCellProperty(crypt.GetCellPropertyRegistry()->Get<CellLabel>());
        }

        simulator.Solve();
        CellBasedSimulationArchiver<2, CryptSimulation2d>::Save(&simulator);

        // ... and checking visualization of labelled cells against previous run
        OutputFileHandler handler(output_directory, false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.viznodes";
        TS_ASSERT(NumericFileComparison( results_file, "crypt/test/data/MakeMeinekeGraphs/results.viznodes").CompareFiles());

        std::string results_file2 = handler.GetOutputDirectoryFullPath() + "results_from_time_0/results.vizcelltypes";
        FileComparison( results_file2, "crypt/test/data/MakeMeinekeGraphs/results.vizcelltypes").CompareFiles();

        // Next, test LabelSPhaseCells()
        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            cell_iter->RemoveCellProperty<CellLabel>();
        }
        crypt_statistics.LabelSPhaseCells();

        // Iterate over cells checking for correct labels
        unsigned counter = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            bool is_labelled = cell_iter->HasCellProperty<CellLabel>();
            bool in_s_phase = (static_cast <AbstractPhaseBasedCellCycleModel*>(cell_iter->GetCellCycleModel())->GetCurrentCellCyclePhase() == S_PHASE);

            TS_ASSERT_EQUALS(is_labelled, in_s_phase);

            if (in_s_phase)
            {
                counter++;
            }
        }

        TS_ASSERT_EQUALS(counter, 28u);

        // Test that LabelAllCellsAsHealthy sets all cells back to be UNLABELLED wild type cells
        crypt_statistics.LabelAllCellsAsHealthy();

        // Iterate over cells checking for correct labels
        counter = 0;
        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS(cell_iter->GetMutationState()->IsType<WildTypeCellMutationState>(), true);
            TS_ASSERT_EQUALS(cell_iter->HasCellProperty<CellLabel>(), false);
            counter++;
        }

        TS_ASSERT_EQUALS(counter, simulator.rGetCellPopulation().GetNumRealCells());

        crypt_statistics.LabelSPhaseCells();

        simulator.SetEndTime(2*time_of_each_run);
        simulator.Solve();

        // TEST CryptStatistics::AreCryptSectionCellsLabelled

        // Set cells which are not in the crypt section to be in state APC +/-, so that we can
        // see the section
        test_section = crypt_statistics.GetCryptSectionPeriodic(crypt_length + 2.0, 8.0, 8.0);

        MAKE_PTR(ApcOneHitCellMutationState, p_apc1);

        for (AbstractCellPopulation<2>::Iterator cell_iter = crypt.Begin();
             cell_iter != crypt.End();
             ++cell_iter)
        {
            bool in_section = false;
            for (unsigned vector_index=0; vector_index<test_section.size(); vector_index++)
            {
                if (test_section[vector_index] == *cell_iter)
                {
                    in_section = true;
                }
            }
            if (!in_section)
            {
                cell_iter->SetMutationState(p_apc1);
            }
        }
        simulator.SetEndTime(3*time_of_each_run);
        simulator.Solve();

        std::vector<CellPtr> crypt_section = crypt_statistics.GetCryptSection(crypt_length + 2.0, 8.0, 8.0);
        std::vector<bool> labelled = crypt_statistics.AreCryptSectionCellsLabelled(crypt_section);

        // Test that the vector of booleans corresponds with a visualisation of the data -
        // only the first few cells at the base of the crypt have been labelled
        for (unsigned vector_index=0; vector_index<labelled.size(); vector_index++)
        {
            if (vector_index == 2u || vector_index == 3u || vector_index == 4u)
            {
                TS_ASSERT(labelled[vector_index]);
            }
            else
            {
                TS_ASSERT(!labelled[vector_index]);
            }
        }
        RandomNumberGenerator::Destroy();
    }

    /**
     * This test runs multiple crypt simulations and records whether
     * or not labelled cells are in a randomly chosen crypt section.
     */
    void TestMultipleCryptSimulations()
    {
        std::string output_directory = "MakeMoreMeinekeGraphs";

        double crypt_length = 22.0;

        // Create mesh
        unsigned cells_across = 13;
        unsigned cells_up = 25;
        double crypt_width = 12.1;
        unsigned thickness_of_ghost_layer = 3;

        unsigned num_simulations = 2;

        // Guess of maximum number of cells a crypt section may contain
        unsigned max_length_of_crypt_section = 5 * (unsigned)sqrt(pow(cells_across/2.0+1,2.0) + pow((double)cells_up,2.0));

        std::vector<unsigned> labelled_cells_counter(max_length_of_crypt_section);

        for (unsigned i=0; i<max_length_of_crypt_section; i++)
        {
            labelled_cells_counter[i] = 0;
        }

        double time_of_each_run;
        std::vector<bool> labelled;

        CryptStatistics* p_crypt_statistics;

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2>* p_crypt;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        std::vector<unsigned> location_indices;

        Cylindrical2dMesh* p_mesh;
        SimulationTime* p_simulation_time;

        // Loop over the number of simulations
        for (unsigned simulation_index=0; simulation_index<num_simulations; simulation_index++)
        {
            // Create new structures for each simulation
            p_mesh = generator.GetCylindricalMesh();
            location_indices = generator.GetCellLocationIndices();

            // Reset start time
            SimulationTime::Destroy();
            p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);

            // Set up cells
            std::vector<CellPtr> temp_cells;
            CryptCellsGenerator<UniformG1GenerationalCellCycleModel> cells_generator;
            cells_generator.Generate(temp_cells, p_mesh, std::vector<unsigned>(), true, 0.3, 2.0, 3.0, 4.0, true);

            // This awkward way of setting up the cells is a result of #430
            std::vector<CellPtr> cells;
            for (unsigned i=0; i<location_indices.size(); i++)
            {
                cells.push_back(temp_cells[location_indices[i]]);
            }

            // Set up crypt
            p_crypt = new MeshBasedCellPopulationWithGhostNodes<2>(*p_mesh, cells, location_indices);

            // Set cell population to output cell types
            p_crypt->AddCellPopulationCountWriter<CellMutationStatesCountWriter>();

            // Set up crypt simulation
            CryptSimulation2d simulator(*p_crypt, false, false);
            simulator.SetOutputDirectory(output_directory);

            // Set length of simulation here
            time_of_each_run = 10.0*simulator.GetDt(); // for each run
            simulator.SetEndTime(time_of_each_run);

            // Create a force laws and pass it to the simulation
            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
            p_linear_force->SetMeinekeSpringStiffness(30.0); // normally 15.0;
            simulator.AddForce(p_linear_force);

            // Set up cell killer
            MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&simulator.rGetCellPopulation(), crypt_length));
            simulator.AddCellKiller(p_killer);

            simulator.UseJiggledBottomCells();

            // Set up crypt statistics
            p_crypt_statistics = new CryptStatistics(*p_crypt);

            // Run for a bit
            simulator.Solve();
            p_crypt_statistics->LabelSPhaseCells();

            simulator.SetEndTime(2.0*time_of_each_run);
            simulator.Solve();

            std::vector<CellPtr> crypt_section = p_crypt_statistics->GetCryptSection(crypt_length + 2.0, 8.0, 8.0);
            labelled = p_crypt_statistics->AreCryptSectionCellsLabelled(crypt_section);

            // Store information from this simulation in a global vector
            for (unsigned cell_index=0; cell_index < labelled.size(); cell_index++)
            {
                TS_ASSERT(cell_index < labelled_cells_counter.size());

                if (labelled[cell_index])
                {
                    labelled_cells_counter[cell_index]++;
                }
            }

            // Tidy up
            cells.clear();
            labelled.clear();
            WntConcentration<2>::Destroy();

            delete p_crypt_statistics;
            delete p_crypt;
        }

        // Calculate percentage of labelled cells at each position in 'labelled_cells_counter'
        std::vector<double> percentage_of_labelled_cells(max_length_of_crypt_section);
        for (unsigned index=0; index < max_length_of_crypt_section; index ++)
        {
            percentage_of_labelled_cells[index] = (double) labelled_cells_counter[index]/num_simulations;
        }

        // Write data to file
        SimpleDataWriter writer1(output_directory, "percentage_of_labelled_cells.dat", percentage_of_labelled_cells, false);

        // Test against previous run
        // ... and checking visualization of labelled cells against previous run
        OutputFileHandler handler(output_directory, false);
        std::string results_file = handler.GetOutputDirectoryFullPath() + "percentage_of_labelled_cells.dat";
        FileComparison( results_file, "crypt/test/data/MakeMoreMeinekeGraphs/percentage_of_labelled_cells.dat").CompareFiles();

        RandomNumberGenerator::Destroy();
    }
};

#endif /*TESTCRYPTSTATISTICS_HPP_*/
