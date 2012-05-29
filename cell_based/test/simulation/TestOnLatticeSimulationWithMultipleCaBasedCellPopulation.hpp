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

#ifndef TESTONLATTICESIMULATIONWITHMULTIPLECABASEDCELLPOPULATION_HPP_
#define TESTONLATTICESIMULATIONWITHMULTIPLECABASEDCELLPOPULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
//#include "ReplicatableVector.hpp"
#include "NumericFileComparison.hpp"
//#include "FunctionalBoundaryCondition.hpp"
#include "AveragedSourcePde.hpp"


#include "PottsBasedCellPopulation.hpp"

#include "AbstractCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
//#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "OnLatticeSimulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "Warnings.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "DiffusionMultipleCaUpdateRule.hpp"

//  class SimplePdeForTesting : public AbstractLinearEllipticPde<2,2>
//    {
//    public:
//        double ComputeConstantInUSourceTerm(const ChastePoint<2>&, Element<2,2>* pElement)
//        {
//            return 0;
//        }
//
//        double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>*)
//        {
//            return 1;
//        }
//
//        c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
//        {
//            return identity_matrix<double>(2u);
//        }
//    };

    /**
     * For use in TestOnLatticeSimulationWithPdes::TestWithBoundaryConditionVaryingInTime.
     */
    double bc_func(const ChastePoint<2>& p)
    {
        double value = SimulationTime::Instance()->GetTime();
        return value;
    }

class TestOnLatticeSimulationWithMultipleCaBasedCellPopulation : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

    void RandomlyLabelCells(std::vector<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
    {
        for (unsigned i = 0; i<rCells.size(); i++)
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
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes());

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);
        node_based_cell_population.SetMechanicsCutOffLength(1.5);

        TS_ASSERT_THROWS_THIS(OnLatticeSimulation<2> simulator(node_based_cell_population),
            "OnLatticeSimulations require a subclass of AbstractOnLatticeCellPopulation.");
    }

    void TestMoreOnLatticeSimulationExceptions()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Create cell population
        PottsBasedCellPopulation<2> potts_based_cell_population(*p_mesh, cells);

        // Try to set up off lattice simulation
        TS_ASSERT_THROWS_THIS(OffLatticeSimulation<2> simulator(potts_based_cell_population),
            "OffLatticeSimulations require a subclass of AbstractOffLatticeCellPopulation.");
    }

    void TestMultipleCaSingleCellRandomMovement() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // timestep and size of domain to let us calculate the probabilities of movement.
        double delta_t = 1;
        double diffusion_parameter = 0.1;
        unsigned num_runs = 2000;
        unsigned location_of_cell[9] = {0,0,0,0,0,0,0,0,0};

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(3, 0, 0, 3, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 1u, DIFFERENTIATED);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        location_indices.push_back(4u);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        std::string output_directory = "TestMultipleCaSingleCellRandomMovement";
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(delta_t);
        simulator.SetEndTime(delta_t);

        /*
         * Adding update rule(s).
         */
        MAKE_PTR(DiffusionMultipleCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(diffusion_parameter);
        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        for (unsigned i=1; i <= num_runs; i++)
        {
        	simulator.SetEndTime(delta_t*i);
        	// Run simulation
        	simulator.Solve();

        	TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 1u);
        	AbstractCellPopulation<2>::Iterator cell_iter = simulator.rGetCellPopulation().Begin();

			unsigned cell_location = simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter);
        	assert(cell_location<9);
        	location_of_cell[cell_location]++;

        	// Reset the position of the cell
        	simulator.rGetCellPopulation().MoveCellInLocationMap(*cell_iter,cell_location,4u);

			assert(simulator.rGetCellPopulation().GetLocationIndexUsingCell(*cell_iter)==4u);

        }

        // Check that still have only one cell
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 1u);
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

        // TODO Check that its moving correctly
        double probability_of_occupation[9];

        for (unsigned i=0; i<9; i++)
        {
        	probability_of_occupation[i] = (double) location_of_cell[i]/(double) num_runs;
        }

        // Note that these simulations are stochastic and so the tolerances are relatively loose
        TS_ASSERT_DELTA(probability_of_occupation[0],diffusion_parameter*delta_t/4.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[1],diffusion_parameter*delta_t/2.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[2],diffusion_parameter*delta_t/4.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[3],diffusion_parameter*delta_t/2.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[4],1.0 - 3.0 * diffusion_parameter*delta_t, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[5],diffusion_parameter*delta_t/2.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[6],diffusion_parameter*delta_t/4.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[7],diffusion_parameter*delta_t/2.0, 1e-2);
        TS_ASSERT_DELTA(probability_of_occupation[8],diffusion_parameter*delta_t/4.0, 1e-2);
        
        
        //For coverage
        simulator.RemoveAllMultipleCaUpdateRules();

    }

    void TestCaMonolayerWithBirth() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 2u, STEM);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        location_indices.push_back(50u);
        location_indices.push_back(51u);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        std::string output_directory = "TestMultipleCaMonolayerWithBirth";
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(0.1);
        simulator.SetEndTime(40);

        // Adding update rule(s).
        MAKE_PTR(DiffusionMultipleCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);

        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 17u);///\todo #2066 Check this!

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 15u);///\todo #2066 Check this!
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);


        // Now remove the update rules and check that only birth happens when the simulator runs again
        simulator.RemoveAllPottsUpdateRules();
        simulator.SetEndTime(50);
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 17u);



#ifdef CHASTE_VTK
        //Test that VTK writer has produced a file
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

        // Resetting the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), DIFFERENTIATED);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<p_mesh->GetNumNodes(); index++)
        {
        	location_indices.push_back(index);
        }
        TS_ASSERT_EQUALS(location_indices.size(),p_mesh->GetNumNodes());

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCaMonolayerWithDeath");
        simulator.SetDt(0.1);
        simulator.SetEndTime(0.1); //only one step as only care about cells being killed

        // No movement rule as only care about cell death

        // Add a cell Killer that will kill all cells in the top half of the domain
        c_vector<double,2> point = zero_vector<double>(2);
        point[1] = 4.5;
        c_vector<double,2> normal = unit_vector<double>(2,1); 
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point, normal)); //v>4.5
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 50u);

        // Test no deaths and some births
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 50u);

        // Check cells above y=5.5 (i.e. above index 50) have been killed and removed.
        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
        	if(i<50)
        	{
        		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetCellUsingLocationIndex(i)->GetCellId(),i);
        	}
        	else
        	{
        		TS_ASSERT_THROWS_THIS(simulator.rGetCellPopulation().GetCellUsingLocationIndex(i),"Location index input argument does not correspond to a Cell");
        	}
        }
    }

    void TestMultipleCaMultipleCellsRandomMovement() throw (Exception)
    {
        /*
         * RandomMovement has been tested in TestMultipleCaSingleCellRandomMovement for one cell
         * per lattice site.
         * This test is just to ensure that the above test works when there are multiple cells per lattice site.
         */

        EXIT_IF_PARALLEL;

         // Create a simple 2D PottsMesh
         PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
         PottsMesh<2>* p_mesh = generator.GetMesh();

         // Create cells
         std::vector<CellPtr> cells;
         CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
         cells_generator.GenerateBasicRandom(cells, 40u, DIFFERENTIATED);

         // Specify where cells lie 4 cells in the first ten sites
         std::vector<unsigned> location_indices;
         for (unsigned index=0; index<10u; index++)
         {
            location_indices.push_back(index);
            location_indices.push_back(index);
            location_indices.push_back(index);
            location_indices.push_back(index);
         }

         // Create cell population
         MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 4u);

         // Set up cell-based simulation
         OnLatticeSimulation<2> simulator(cell_population);
         std::string output_directory = "TestMultipleCaMultipleCellRandomMovement";
         simulator.SetOutputDirectory(output_directory);
         simulator.SetDt(1);
         simulator.SetEndTime(10);

         /*
          * Adding update rule(s).
          */
         MAKE_PTR(DiffusionMultipleCaUpdateRule<2>, p_diffusion_update_rule);
         p_diffusion_update_rule->SetDiffusionParameter(0.1);
         simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

         // Run simulation
         simulator.Solve();

         TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 40u);

    }

    void TestMultipleCaMultipleCellsRandomMovementIn3d() throw (Exception)
    {
       EXIT_IF_PARALLEL;

        // Create a simple 3D PottsMesh
        PottsMeshGenerator<3> generator(10, 0, 0, 10, 0, 0, 10, 0, 0);
        PottsMesh<3>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 40u, DIFFERENTIATED);

        // Specify where cells lie 4 cells in the first ten sites
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<10u; index++)
        {
           location_indices.push_back(index);
           location_indices.push_back(index);
           location_indices.push_back(index);
           location_indices.push_back(index);
        }

        // Create cell population
        MultipleCaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices, 4u);

        // Set up cell-based simulation
        OnLatticeSimulation<3> simulator(cell_population);
        std::string output_directory = "TestMultipleCaMultipleCellRandomMovementIn3d";
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1);
        simulator.SetEndTime(100);

        /*
         * Adding update rule(s).
         */
        MAKE_PTR(DiffusionMultipleCaUpdateRule<3>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.1);
        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();

        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 40u);

    }


    void TestMultipleCellsPerLatticeSiteWithBirth() throw (Exception)
    {

        /*
         * Cellular birth has been tested in TestMultipleCaSingleCellWithBirth for one cell per lattice site.
         * This test  adds to the above by further testing cellular birth considering multiple cells per lattice site.
         * A  two-lattice mesh was created and only one lattice had free space to add one daughter cell.
         */

          EXIT_IF_PARALLEL;

          // Create a simple 2D PottsMesh
          PottsMeshGenerator<2> generator(2, 0, 0, 1, 0, 0);
          PottsMesh<2>* p_mesh = generator.GetMesh();

          // Create cells
          std::vector<CellPtr> cells;
          CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
          cells_generator.GenerateBasicRandom(cells, 3u, STEM);

          // Specify where cells lie
          std::vector<unsigned> location_indices;
          location_indices.push_back(0u);
          location_indices.push_back(0u);
          location_indices.push_back(1u);


          // Create cell population
          MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 2);

          // Set up cell-based simulation
          OnLatticeSimulation<2> simulator(cell_population);
          std::string output_directory = "TestMultipleCellsPerLatticeSiteWithBirth";
          simulator.SetOutputDirectory(output_directory);
          simulator.SetDt(0.1);
          simulator.SetEndTime(40);

          // Adding update rule(s).
          MAKE_PTR(DiffusionMultipleCaUpdateRule<2u>, p_diffusion_update_rule);
          p_diffusion_update_rule->SetDiffusionParameter(0.5);

          simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

          // Run simulation
          simulator.Solve();

          // Check the number of cells
          TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 4u);

          // Test no deaths and some births
          TS_ASSERT_EQUALS(simulator.GetNumBirths(), 1u);
          TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);


    #ifdef CHASTE_VTK
          //Test that VTK writer has produced a file
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

    void TestMultipleCellsPerLatticeSiteWithDeath() throw (Exception)
    {
        /*
         * Cellular death has been tested in TestCaMonolayerWithDeath for one cell per lattice site.
         * This test is just to ensure that the above test works when there are multiple cells per lattice site.
         */

        EXIT_IF_PARALLEL;

        // Resetting the maximum cell ID to zero (to account for previous tests)
        CellId::ResetMaxCellId();

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 2*p_mesh->GetNumNodes(), DIFFERENTIATED);

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
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 2);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCellsPerLatticeSiteWithDeath");
        simulator.SetDt(0.1);
        simulator.SetEndTime(0.1); //only one step as only care about cells being killed

        // No movement rule as only care about cell death

        // Add a cell Killer that will kill all cells in the top half of the domain
        c_vector<double,2> point = zero_vector<double>(2);
        point[1] = 4.5;
        c_vector<double,2> normal = unit_vector<double>(2,1); 
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, point, normal)); //v>4.5
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
         if(i<50)
         {
             TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetCellsUsingLocationIndex(i).size(),2u);
         }
         else
         {
             TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetCellsUsingLocationIndex(i).size(),0u);
         }
        }
}

//
//    void NoTestStandardResultForArchivingTestsBelow() throw (Exception)
//    {
//        // Create a simple 2D PottsMesh
//        PottsMeshGenerator<2> generator(10, 1, 4, 10, 1, 4);
//        PottsMesh<2>* p_mesh = generator.GetMesh();
//
//        // Create cells
//        std::vector<CellPtr> cells;
//        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), STEM);
//
//        // Create cell population
//        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
//
//        // Set up cell-based simulation
//        OnLatticeSimulation<2> simulator(cell_population);
//        simulator.SetOutputDirectory("TestOnLatticeSimulationWithPottsBasedCellPopulationStandardResult");
//        simulator.SetDt(0.1);
//        simulator.SetEndTime(20);
//        simulator.SetSamplingTimestepMultiple(10);
//
//        // Create update rules and pass to the simulation
//        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
//        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
//        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
//        simulator.AddPottsUpdateRule(p_adhesion_update_rule);
//
//        // Run simulation
//        simulator.Solve();
//
//        // Check some results
//        PottsElement<2>* element_0 = static_cast <PottsBasedCellPopulation<2>*>(&simulator.rGetCellPopulation())->GetElement(0u);
//        TS_ASSERT_EQUALS(element_0->GetNumNodes(), 16u);
//        TS_ASSERT_EQUALS(element_0->GetNode(0)->GetIndex(), 34u);
//        TS_ASSERT_EQUALS(element_0->GetNode(8)->GetIndex(), 24u);
//        TS_ASSERT_EQUALS(element_0->GetNode(15)->GetIndex(), 32u);
//
//        PottsElement<2>* element_1 = static_cast <PottsBasedCellPopulation<2>*>(&simulator.rGetCellPopulation())->GetElement(1u);
//        TS_ASSERT_EQUALS(element_1->GetNumNodes(), 16u);
//        TS_ASSERT_EQUALS(element_1->GetNode(0)->GetIndex(), 46u);
//        TS_ASSERT_EQUALS(element_1->GetNode(8)->GetIndex(), 69u);
//        TS_ASSERT_EQUALS(element_1->GetNode(15)->GetIndex(), 25u);
//    }
//
//    void NoTestSave() throw (Exception)
//    {
//        // Create a simple 2D PottsMesh
//        PottsMeshGenerator<2> generator(10, 1, 4, 10, 1, 4);
//        PottsMesh<2>* p_mesh = generator.GetMesh();
//
//        // Create cells
//        std::vector<CellPtr> cells;
//        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), STEM);
//
//        // Create cell population
//        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
//
//        // Set up cell-based simulation
//        OnLatticeSimulation<2> simulator(cell_population);
//        simulator.SetOutputDirectory("TestOnLatticeSimulationWithPottsBasedCellPopulationSaveAndLoad");
//        simulator.SetDt(0.1);
//        simulator.SetEndTime(10);
//        simulator.SetSamplingTimestepMultiple(10);
//
//        // Create update rules and pass to the simulation
//        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
//        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
//        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
//        simulator.AddPottsUpdateRule(p_adhesion_update_rule);
//
//        // Run simulation
//        simulator.Solve();
//
//        // Save the results
//        CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Save(&simulator);
//    }
//
//    void NoTestLoad() throw (Exception)
//    {
//        // Load the simulation from the TestSave method above and
//        // run it from 10.0 to 15.0
//        OnLatticeSimulation<2>* p_simulator1;
//        p_simulator1 = CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Load("TestOnLatticeSimulationWithPottsBasedCellPopulationSaveAndLoad", 10.0);
//
//        p_simulator1->SetEndTime(15.0);
//        p_simulator1->Solve();
//
//        // Save, then reload and run from 15.0 to 20.0
//        CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Save(p_simulator1);
//        OnLatticeSimulation<2>* p_simulator2
//            = CellBasedSimulationArchiver<2, OnLatticeSimulation<2> >::Load("TestOnLatticeSimulationWithPottsBasedCellPopulationSaveAndLoad", 15.0);
//
//        p_simulator2->SetEndTime(20.0);
//        p_simulator2->Solve();
//
//        // These results are from time 20.0 in TestStandardResultForArchivingTestsBelow()
//        PottsElement<2>* element_0 = static_cast <PottsBasedCellPopulation<2>*>(&p_simulator2->rGetCellPopulation())->GetElement(0u);
//        TS_ASSERT_EQUALS(element_0->GetNumNodes(), 16u);
//        TS_ASSERT_EQUALS(element_0->GetNode(0)->GetIndex(), 34u);
//        TS_ASSERT_EQUALS(element_0->GetNode(8)->GetIndex(), 24u);
//        TS_ASSERT_EQUALS(element_0->GetNode(15)->GetIndex(), 32u);
//
//        PottsElement<2>* element_1 = static_cast <PottsBasedCellPopulation<2>*>(&p_simulator2->rGetCellPopulation())->GetElement(1u);
//        TS_ASSERT_EQUALS(element_1->GetNumNodes(), 16u);
//        TS_ASSERT_EQUALS(element_1->GetNode(0)->GetIndex(), 46u);
//        TS_ASSERT_EQUALS(element_1->GetNode(8)->GetIndex(), 69u);
//        TS_ASSERT_EQUALS(element_1->GetNode(15)->GetIndex(), 25u);
//
//        // Tidy up
//        delete p_simulator1;
//        delete p_simulator2;
//    }
};

#endif /*TESTONLATTICESIMULATIONWITHMULTIPLECABASEDCELLPOPULATION_HPP_*/
