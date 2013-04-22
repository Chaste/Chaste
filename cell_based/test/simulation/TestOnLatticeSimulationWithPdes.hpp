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
#include "CellwiseSourcePde.hpp"
#include "SmartPointers.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "DiffusionMultipleCaUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"
#include "ChemotaxisPottsUpdateRule.hpp"
#include "OnLatticeSimulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Warnings.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellBasedPdeHandlerOnCuboid.hpp"
#include "CellBasedPdeHandler.hpp"
#include "FileComparison.hpp"


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
        mLastStartTime = (double) std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = (double) std::clock();
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
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsSimulationWithPdesThrowsException");
        simulator.SetEndTime(0.1);


        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        AveragedSourcePde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("nutrient");

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
    }

    void TestPottsBasedWithCoarseMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsBasedCellPopulationWithPdes");
        simulator.SetEndTime(0.1);


        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        AveragedSourcePde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("nutrient");

        CellBasedPdeHandler<2> pde_handler(&cell_population);


        pde_handler.AddPdeAndBc(&pde_and_bc);
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(50.0, 50.0);
        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(10.0, cuboid, true);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

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
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
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
    }

    void TestMultipleCaBasedWithoutCoarseMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 100);

        std::vector<unsigned> location_indices;
        for(unsigned i=0; i<100; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCaBasedCellPopulationWithPdesOnNaturalMesh");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);


        // Set up 2 PDEs (one with dirichlet and one with Neumann BCS) and pass to simulation via handler (zero uptake to check analytic solution)
        // Note even with uptake the PDE has U=0 as solution with zero neumann conditions.
        AveragedSourcePde<2> pde_1(cell_population, 0.0);
        ConstBoundaryCondition<2> bc_1(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc_1(&pde_1, &bc_1, false);
        pde_and_bc_1.SetDependentVariableName("nutrient_dirichlet");

        AveragedSourcePde<2> pde_2(cell_population, 0.0);
        ConstBoundaryCondition<2> bc_2(0.0);
        PdeAndBoundaryConditions<2> pde_and_bc_2(&pde_2, &bc_2, true);
        pde_and_bc_2.SetDependentVariableName("nutrient_neumann");

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc_1);
        pde_handler.AddPdeAndBc(&pde_and_bc_2);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);


        // Create update rules and pass to the simulation
        MAKE_PTR(DiffusionMultipleCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Solve the system
        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;
            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient_dirichlet"), analytic_solution, 1e-2);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient_neumann"), 0.0, 1e-2);
        }

        //Test coarse mesh has the same nodes as the PottsMesh
        TetrahedralMesh<2,2>* p_coarse_mesh = simulator.GetCellBasedPdeHandler()->GetCoarsePdeMesh();

        TS_ASSERT_EQUALS(p_coarse_mesh->GetNumNodes(),p_mesh->GetNumNodes());
        TS_ASSERT_DELTA(p_coarse_mesh->GetWidth(0),p_mesh->GetWidth(0),1e-8);
        TS_ASSERT_DELTA(p_coarse_mesh->GetWidth(1),p_mesh->GetWidth(1),1e-8);

        for(unsigned i=0; i< p_coarse_mesh->GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(p_coarse_mesh->GetNode(i)->rGetLocation()[0],p_mesh->GetNode(i)->rGetLocation()[0],1e-8);
            TS_ASSERT_DELTA(p_coarse_mesh->GetNode(i)->rGetLocation()[1],p_mesh->GetNode(i)->rGetLocation()[1],1e-8);
        }
    }

    void TestMultipleCaBasedWithCellwiseSourcePde() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 100);

        std::vector<unsigned> location_indices;
        for(unsigned i=0; i<100; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCaBasedCellPopulationWithPdesOnNaturalMesh");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);

        CellwiseSourcePde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("nutrient");

        double pde_coefficient = 0.0;
        TS_ASSERT_DELTA(pde.ComputeLinearInUCoeffInSourceTermAtNode(*p_mesh->GetNode(0)), pde_coefficient, 1e-3);

        CellwiseSourcePde<2> non_trivial_pde(cell_population, 1.0);

        double non_trivial_pde_coefficient = 1.0;
        for (PottsMesh<2>::NodeIterator node_iter = p_mesh->GetNodeIteratorBegin();
                node_iter != p_mesh->GetNodeIteratorEnd();
                ++node_iter)
        {
            TS_ASSERT_DELTA(non_trivial_pde.ComputeLinearInUCoeffInSourceTermAtNode(*node_iter), non_trivial_pde_coefficient, 1e-3);
        }

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;
            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
        }
    }


    void TestMultipleCaBasedWithoutCoarseMeshUsingPdeHandlerOnCuboid() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 100);

        std::vector<unsigned> location_indices;
        for(unsigned i=0; i<100; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCaBasedCellPopulationWithPdesOnCuboid");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);

        // Set up 2 PDEs (one with dirichlet and one with Neumann BCS) and pass to simulation via handler (zero uptake to check analytic solution)
        // Note even with uptake the PDE has U=0 as solution with zero neumann conditions.
        AveragedSourcePde<2> pde_1(cell_population, 0.0);
        ConstBoundaryCondition<2> bc_1(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc_1(&pde_1, &bc_1, false);
        pde_and_bc_1.SetDependentVariableName("nutrient_dirichlet");

        AveragedSourcePde<2> pde_2(cell_population, 0.0);
        ConstBoundaryCondition<2> bc_2(0.5);
        PdeAndBoundaryConditions<2> pde_and_bc_2(&pde_2, &bc_2, true);
        pde_and_bc_2.SetDependentVariableName("nutrient_neumann");

        CellBasedPdeHandlerOnCuboid<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc_1);
        pde_handler.AddPdeAndBc(&pde_and_bc_2);
        pde_handler.SetImposeBcsOnCoarseBoundary(false);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Solve the system
        simulator.Solve();

        // Test that PDE solver is working correctly
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
        	c_vector<double, 2> cell_location = simulator.rGetCellPopulation().GetLocationOfCellCentre(*cell_iter);

        	if (cell_location[1] < 1e-6 || cell_location[1] > 9 - 1e-6)
        	{
        		TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient_dirichlet"),1.0, 1e-2);
        		TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient_neumann"), 0.5, 1e-2);
        	}
        	else
        	{
        		TS_ASSERT_LESS_THAN(1.0,cell_iter->GetCellData()->GetItem("nutrient_dirichlet"));
        		TS_ASSERT_LESS_THAN(0.5,cell_iter->GetCellData()->GetItem("nutrient_neumann"));
        	}
        }
    }

    /*
     * This tests that a sensible error is thrown if the coarse mesh is too small/
     */
    void TestMultipleCaBasedCellsOutsideMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 100);

        std::vector<unsigned> location_indices;
        for(unsigned i=0; i<100; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCaBasedCellPopulationWithPdesWithCellsOutsideMesh");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);


        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        AveragedSourcePde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("nutrient");

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);
        pde_handler.SetImposeBcsOnCoarseBoundary(true);

        // Create coarse mesh smaller than the PottsMesh
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(5.0, 5.0);

        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(1.0, cuboid);


        simulator.SetCellBasedPdeHandler(&pde_handler);


        // Solve the system
        TS_ASSERT_THROWS_THIS(simulator.Solve(), "Point [6,0] is not in mesh - all elements tested");
    }


    void TestMultipleCaBasedWithCoarseMesh() throw(Exception)
    {
        EXIT_IF_PARALLEL;

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

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCaBasedCellPopulationWithPdes");
        simulator.SetEndTime(0.1);


        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        AveragedSourcePde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("nutrient");

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Create coarse mesh and centre it on the centre of the Potts Mesh
        c_vector<double,2> centre_of_potts_mesh = zero_vector<double>(2);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            centre_of_potts_mesh += p_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_potts_mesh /= p_mesh->GetNumNodes();
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(5.0, 5.0);

        c_vector<double, 2> translation = 0.5*(lower.rGetLocation() + upper.rGetLocation()) - centre_of_potts_mesh;
        lower.rGetLocation() -= translation;
        upper.rGetLocation() -= translation;
        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(1.0, cuboid);
        pde_handler.SetImposeBcsOnCoarseBoundary(true);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create update rules and pass to the simulation
        MAKE_PTR(DiffusionMultipleCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Solve the system
        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;
            // Test that PDE solver is working correctly
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
        }

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);
        TetrahedralMesh<2,2>* p_coarse_mesh = simulator.GetCellBasedPdeHandler()->GetCoarsePdeMesh();
        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();
        c_vector<double,2> centre_diff = centre_of_coarse_pde_mesh - centre_of_potts_mesh;

        // Test that the two centres match
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

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;

            // Test that PDE solver is working correctly on both pdes
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
        }
    }

    void TestPottsBasedWithCoarseMeshTwoEquations() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsBasedCellPopulationWithTwoPdes");
        simulator.SetEndTime(0.1);

        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        AveragedSourcePde<2> pde_1(cell_population, 0.0);
        ConstBoundaryCondition<2> bc_1(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc_1(&pde_1, &bc_1, false);
        pde_and_bc_1.SetDependentVariableName("quantity 1");

        AveragedSourcePde<2> pde_2(cell_population, 0.0);
        ConstBoundaryCondition<2> bc_2(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc_2(&pde_2, &bc_2, false);
        pde_and_bc_2.SetDependentVariableName("quantity 2");

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc_1);
        pde_handler.AddPdeAndBc(&pde_and_bc_2);
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(50.0, 50.0);
        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(10.0, cuboid, true);
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
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("quantity 1"), analytic_solution, 1e-2);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("quantity 2"), analytic_solution, 1e-2);
        }
#ifdef CHASTE_VTK
        //First file exists
        FileFinder vtk_file("TestPottsBasedCellPopulationWithTwoPdes/results_from_time_0/pde_results_1.vtu", RelativeTo::ChasteTestOutput);
        TS_ASSERT(vtk_file.Exists());
        // Check that the second VTK file for the solution has the dependent quantities
        OutputFileHandler handler("TestPottsBasedCellPopulationWithTwoPdes", false);
        VtkMeshReader<3,3> vtk_reader(handler.GetOutputDirectoryFullPath()+"results_from_time_0/pde_results_2.vtu");
        std::vector<double> data1;
        //There is no Oxygen
        TS_ASSERT_THROWS_CONTAINS(vtk_reader.GetPointData("Oxygen", data1), "No point data");
        TS_ASSERT(data1.empty());
        vtk_reader.GetPointData("quantity 1", data1);
        TS_ASSERT_EQUALS(data1.size(), 6u*6u);
        std::vector<double> data2;
        vtk_reader.GetPointData("quantity 2", data2);
        TS_ASSERT_EQUALS(data1.size(), data2.size());
#endif //CHASTE_VTK
    }

    // Under construction: Test growth of a population of cells that consumes nutrient
    void TestOnLatticeSpheroidWithNutrient() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 4);

        std::vector<unsigned> location_indices;

        location_indices.push_back(55);
        location_indices.push_back(56);
        location_indices.push_back(65);
        location_indices.push_back(66);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOnLatticeSpheroidWithNutrient");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(10);

        // Set up PDE and pass to simulation via handler
        double nutrient_uptake_rate=-0.1;
        AveragedSourcePde<2> pde(cell_population, nutrient_uptake_rate);
        ConstBoundaryCondition<2> bc(1.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("nutrient");

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        pde_handler.SetImposeBcsOnCoarseBoundary(true);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        // Create update rules and pass to the simulation
        MAKE_PTR(DiffusionMultipleCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Solve the system
        simulator.Solve();

        //Test coarse mesh has the same nodes as the PottsMesh
        TetrahedralMesh<2,2>* p_coarse_mesh = simulator.GetCellBasedPdeHandler()->GetCoarsePdeMesh();

        TS_ASSERT_EQUALS(p_coarse_mesh->GetNumNodes(),p_mesh->GetNumNodes());
        TS_ASSERT_DELTA(p_coarse_mesh->GetWidth(0),p_mesh->GetWidth(0),1e-8);
        TS_ASSERT_DELTA(p_coarse_mesh->GetWidth(1),p_mesh->GetWidth(1),1e-8);

        for(unsigned i=0; i< p_coarse_mesh->GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(p_coarse_mesh->GetNode(i)->rGetLocation()[0],p_mesh->GetNode(i)->rGetLocation()[0],1e-8);
            TS_ASSERT_DELTA(p_coarse_mesh->GetNode(i)->rGetLocation()[1],p_mesh->GetNode(i)->rGetLocation()[1],1e-8);
        }
    }

    void TestPdeOutput() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("TestOnLatticeSpheroidWithNutrient", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_nut(results_dir + "/results.vizcoarsepdesolution", "cell_based/test/data/TestOnLatticeSpheroidWithNutrient/results.vizcoarsepdesolution");
        TS_ASSERT(comp_nut.CompareFiles());
        FileComparison( results_dir + "/results.vizcoarsepdesolution", "cell_based/test/data/TestOnLatticeSpheroidWithNutrient/results.vizcoarsepdesolution").CompareFiles();
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHPDES_HPP_*/
