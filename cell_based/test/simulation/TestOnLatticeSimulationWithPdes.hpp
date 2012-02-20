/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef TESTOFFLATTICESIMULATIONWITHPDES_HPP_
#define TESTOFFLATTICESIMULATIONWITHPDES_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ReplicatableVector.hpp"
#include "NumericFileComparison.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FunctionalBoundaryCondition.hpp"
#include "AveragedSourcePde.hpp"
#include "SmartPointers.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "ChemotaxisPottsUpdateRule.hpp"
#include "OnLatticeSimulation.hpp"
#include "Warnings.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellPopulation.hpp"

class SimplePdeForTesting : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>&, Element<2,2>* pElement)
    {
        return -1.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>*)
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
};

/**
 * For use in TestOnLatticeSimulationWithPdes::TestWithBoundaryConditionVaryingInTime.
 */
double bc_func(const ChastePoint<2>& p)
{
    double value = SimulationTime::Instance()->GetTime();
    return value;
}

class TestOnLatticeSimulationWithPdes : public AbstractCellBasedTestSuite
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

public:

    void TestPottsBasedWithoutCoarseMeshThrowsException() throw(Exception)
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
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        /*
         * Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
         * we need to pass it some initial conditions to avoid memory errors.
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(cell_population.GetLocationIndexUsingCell(*cell_iter), cell_population.GetLocationIndexUsingCell(*cell_iter), 0);
        }

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsSimulationWithPdesThrowsException");
        simulator.SetEndTime(0.1);


        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        AveragedSourcePde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.SetImposeBcsOnCoarseBoundary(true);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Solve the system
        TS_ASSERT_THROWS_THIS(simulator.Solve(), "Trying to solve a PDE on a cell population that doesn't have a mesh. Try calling UseCoarsePdeMesh().");

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestPottsBasedWithCoarseMesh() throw(Exception)
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
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 1);
        p_data->SetCellPopulation(&cell_population);

        /*
         * Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
         * we need to pass it some initial conditions to avoid memory errors.
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(cell_population.GetLocationIndexUsingCell(*cell_iter), cell_population.GetLocationIndexUsingCell(*cell_iter), 0);
        }

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsBasedCellPopulationWithPdes");
        simulator.SetEndTime(0.1);


        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        AveragedSourcePde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.UseCoarsePdeMesh(10.0, 50.0);
        pde_handler.SetImposeBcsOnCoarseBoundary(true);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Solve the system
        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {

            double analytic_solution = 1.0;
            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(p_data->GetValue(*cell_iter), analytic_solution, 1e-2);
        }

        // Find centre of cell population
        c_vector<double,2> centre_of_cell_population = zero_vector<double>(2);

        for (unsigned i=0; i<simulator.rGetCellPopulation().GetNumNodes(); i++)
        {
            centre_of_cell_population += simulator.rGetCellPopulation().GetNode(i)->rGetLocation();
        }
        centre_of_cell_population /= simulator.rGetCellPopulation().GetNumNodes();

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);
        TetrahedralMesh<2,2>* p_coarse_mesh = simulator.GetCellBasedPdeHandler()->GetCoarsePdeMesh();
        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();

        // Test that the two centres match
        c_vector<double,2> centre_diff = centre_of_cell_population - centre_of_coarse_pde_mesh;
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);

        // Test FindCoarseElementContainingCell() and initialisation of mCellPdeElementMap
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            unsigned containing_element_index = simulator.GetCellBasedPdeHandler()->mCellPdeElementMap[*cell_iter];
            TS_ASSERT_LESS_THAN(containing_element_index, p_coarse_mesh->GetNumElements());
            TS_ASSERT_EQUALS(containing_element_index, simulator.GetCellBasedPdeHandler()->FindCoarseElementContainingCell(*cell_iter));
        }

        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestPottsBasedWithCoarseMeshTwoEquations() throw(Exception)
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
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up CellwiseData and associate it with the cell population
        CellwiseData<2>* p_data = CellwiseData<2>::Instance();
        p_data->SetNumCellsAndVars(cell_population.GetNumRealCells(), 2);
        p_data->SetCellPopulation(&cell_population);

        /*
         * Since values are first passed in to CellwiseData before it is updated in UpdateAtEndOfTimeStep(),
         * we need to pass it some initial conditions to avoid memory errors.
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_data->SetValue(1.0 , cell_population.GetLocationIndexUsingCell(*cell_iter), 0);
            p_data->SetValue(1.0 , cell_population.GetLocationIndexUsingCell(*cell_iter), 1);
        }

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsBasedCellPopulationWithTwoPdes");
        simulator.SetEndTime(0.1);


        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        AveragedSourcePde<2> pde_1(cell_population, 0.0);
        ConstBoundaryCondition<2> bc_1(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc_1(&pde_1, &bc_1, false);

        AveragedSourcePde<2> pde_2(cell_population, 0.0);
        ConstBoundaryCondition<2> bc_2(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc_2(&pde_2, &bc_2, false);

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc_1);
        pde_handler.AddPdeAndBc(&pde_and_bc_2);
        pde_handler.UseCoarsePdeMesh(10.0, 50.0);
        pde_handler.SetImposeBcsOnCoarseBoundary(true);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Solve the system
        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;

            // Test that PDE solver is working correctly on both pdes
            TS_ASSERT_DELTA(p_data->GetValue(*cell_iter, 0), analytic_solution, 1e-2);
            TS_ASSERT_DELTA(p_data->GetValue(*cell_iter, 1), analytic_solution, 1e-2);
        }
        // Tidy up
        CellwiseData<2>::Destroy();
    }

    void TestPdeOutput() throw(Exception)
    {
    	EXIT_IF_PARALLEL;

        OutputFileHandler handler("TestPottsBasedCellPopulationWithPdes", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_nut(results_dir + "/results.vizcoarsepdesolution", "cell_based/test/data/TestPottsBasedCellPopulationWithPdes/results.vizcoarsepdesolution");
        TS_ASSERT(comp_nut.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizcoarsepdesolution cell_based/test/data/TestPottsBasedCellPopulationWithPdes/results.vizcoarsepdesolution").c_str()), 0);
    }


};

#endif /*TESTOFFLATTICESIMULATIONWITHPDES_HPP_*/
