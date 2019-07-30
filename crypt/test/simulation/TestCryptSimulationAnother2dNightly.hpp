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

#ifndef TESTCRYPTSIMULATIONANOTHER2DNIGHTLY_HPP_
#define TESTCRYPTSIMULATIONANOTHER2DNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"
#include "CryptSimulation2d.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "CryptCellsGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "WntConcentration.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "CellBasedEventHandler.hpp"
#include "NumericFileComparison.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "WntCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestCryptSimulationAnother2dNightly : public AbstractCellBasedWithTimingsTestSuite
{
private:
    void setUp()
    {
        CellBasedEventHandler::Disable(); // these tests fail with event-handling on
        AbstractCellBasedWithTimingsTestSuite::setUp();
    }

public:

    /**
     * Sloughing with a sloughing cell killer and not turning
     * into ghost nodes on a non-periodic mesh.
     */
    void TestSloughingCellKillerOnNonPeriodicCrypt()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        unsigned cells_across = 6;
        double domain_width = (double)cells_across-1.0;
        unsigned cells_up = 12;
        double domain_height = ((double)cells_up-1.0)*sqrt(3.0)/2.0;
        unsigned thickness_of_ghost_layer = 4;

        HoneycombMeshGenerator generator(cells_across, cells_up,thickness_of_ghost_layer);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DSloughingDeathNonPeriodic");
        simulator.SetEndTime(4.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, domain_height, true, domain_width));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();
    }

    void TestSloughingDeathWithPeriodicMesh()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        unsigned cells_across = 7;
        unsigned cells_up = 12;
        double crypt_width = 6.0;
        unsigned thickness_of_ghost_layer = 4;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        double crypt_length = cells_up*(sqrt(3.0)/2)*crypt_width/cells_across;

        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("Crypt2DSloughingDeathPeriodic");
        simulator.SetEndTime(4.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Create cell killer and pass in to crypt simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        std::vector<bool> ghost_node_indices_after = (static_cast<MeshBasedCellPopulationWithGhostNodes<2>* >(&(simulator.rGetCellPopulation())))->rGetGhostNodes();
        unsigned num_ghosts = 0;
        for (unsigned i=0; i<ghost_node_indices_after.size(); i++)
        {
            if (ghost_node_indices_after[i])
            {
                num_ghosts++;
            }
        }

        // Check no new ghost nodes have been created
        TS_ASSERT_EQUALS(num_ghosts, p_mesh->GetNumNodes() - crypt.GetNumRealCells());

        // There should be this number of cells left after this amount of time
        // (we have lost two rows of 7 but had a bit of birth too)
        TS_ASSERT_EQUALS(crypt.GetNumRealCells(), 84u);
    }

    void TestMonolayerWithCutoffPointAndNoGhosts()
    {
        EXIT_IF_PARALLEL;    // HoneycombMeshGenerator doesn't work in parallel.

        unsigned num_cells_depth = 11;
        unsigned num_cells_width = 6;
        //double crypt_length = num_cells_depth - 1.0;

        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<FixedG1GenerationalCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true, -1.0);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        CryptSimulation2d simulator(crypt);
        simulator.SetOutputDirectory("MonolayerCutoffPointNoGhosts");
        simulator.SetEndTime(12.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Run simulation
        simulator.Solve();

        ///\todo there are no tests here!
    }

    /*
    * This tests that the results files are correct (added because of #1130).
    */
    void TestResultsFileForLongerCryptSimulation()
    {
        EXIT_IF_PARALLEL; // HoneycombMeshGenerator doesn't work in parallel

        // Set output directory
        std::string output_directory = "TestResultsFileForLongerCryptSimulation";

        // Create cylindrical mesh
        unsigned cells_across = 16;
        unsigned cells_up = 19;
        unsigned thickness_of_ghost_layer = 0;
        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();
        /*
         * Shift all left-hand cells.
         * The first column of cells are on x=0 ( == x=16) which is the periodic boundary.
         * These might flip to x=16 on the first timestep (depending on the re-mesh implementation).
         * The second column of cells are at x=1/2.
         */
        for (Cylindrical2dMesh::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                        node_iter != p_mesh->GetNodeIteratorEnd();
                        ++node_iter)
        {
            if (node_iter->rGetLocation()[0] <= DBL_EPSILON) // In first column
            {
                node_iter->rGetModifiableLocation()[0] += 1e-8; // Shift right
            }
        }

        double crypt_length = cells_up*(sqrt(3.0)/2);

        // Get location indices corresponding to real cells in mesh
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up each cell with a simple Wnt-based cell-cycle model
        std::vector<CellPtr> cells;
        CryptCellsGenerator<SimpleWntCellCycleModel> cell_generator;
        cell_generator.Generate(cells, p_mesh, location_indices, true);

        // Set some model parameters for the cell-cycle model
        for (unsigned index=0; index < cells.size(); index++)
        {
           static_cast<SimpleWntCellCycleModel*>(cells[index]->GetCellCycleModel())->SetSDuration(7.4);
           static_cast<SimpleWntCellCycleModel*>(cells[index]->GetCellCycleModel())->SetG2Duration(1.4);
           static_cast<SimpleWntCellCycleModel*>(cells[index]->GetCellCycleModel())->SetMDuration(0.72);
           static_cast<SimpleWntCellCycleModel*>(cells[index]->GetCellCycleModel())->SetTransitCellG1Duration(9.4);
           static_cast<SimpleWntCellCycleModel*>(cells[index]->GetCellCycleModel())->SetStemCellG1Duration(9.4);
        }

        // Create cell population
        MeshBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Set cell population to output cell types
        crypt.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // Set up instance of WntConcentration singleton and associate it with crypt
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation
        CryptSimulation2d simulator(crypt);

        // Set where to output simulation results
        simulator.SetOutputDirectory(output_directory);

        // Set length of simulation
        simulator.SetEndTime(20.0);

        // Only save results every tenth time step
        simulator.SetSamplingTimestepMultiple(10);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        // Unusual set-up here (corresponds to the Meineke crypt model parameters)
        p_linear_force->SetMeinekeSpringStiffness(30.0);
        // Sets the MeinekeSpringGrowthDuration to be the default MPhase duration
        p_linear_force->SetMeinekeSpringGrowthDuration(static_cast<SimpleWntCellCycleModel*>(crypt.rGetCells().front()->GetCellCycleModel())->GetMDuration());
        simulator.AddForce(p_linear_force);

        // Set up sloughing cell killer and pass in to simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&(simulator.rGetCellPopulation()), crypt_length));
        simulator.AddCellKiller(p_killer);

        // Unusual set-up here (corresponds to the Meineke crypt model parameters)
        simulator.UseJiggledBottomCells();

        // Run simulation
        simulator.Solve();

        // Test that results files are correct
        OutputFileHandler handler(output_directory, false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_ele(results_dir + "/results.vizelements", "crypt/test/data/TestResultsFileForLongerCryptSimulation/results.vizelements");
        TS_ASSERT(comp_ele.CompareFiles(1e-15));
        FileComparison( results_dir + "/results.vizelements", "crypt/test/data/TestResultsFileForLongerCryptSimulation/results.vizelements").CompareFiles();

        NumericFileComparison comp_nodes(results_dir + "/results.viznodes", "crypt/test/data/TestResultsFileForLongerCryptSimulation/results.viznodes");
        TS_ASSERT(comp_nodes.CompareFiles(1e-15));

        NumericFileComparison comp_celltypes(results_dir + "/results.vizcelltypes", "crypt/test/data/TestResultsFileForLongerCryptSimulation/results.vizcelltypes");
        TS_ASSERT(comp_celltypes.CompareFiles(1e-15));

        FileComparison( results_dir + "/results.vizsetup", "crypt/test/data/TestResultsFileForLongerCryptSimulation/results.vizsetup").CompareFiles();
        FileComparison( results_dir + "/results.parameters", "crypt/test/data/TestResultsFileForLongerCryptSimulation/results.parameters").CompareFiles();

        // Tidy up
        WntConcentration<2>::Destroy();
    }
};

#endif /*TESTCRYPTSIMULATIONANOTHER2DNIGHTLY_HPP_*/
