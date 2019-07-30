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

#ifndef TESTOFFLATTICESIMULATIONWITHPDES_HPP_
#define TESTOFFLATTICESIMULATIONWITHPDES_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "ReplicatableVector.hpp"
#include "NumericFileComparison.hpp"
#include "WildTypeCellMutationState.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "PottsMeshGenerator.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "Warnings.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"

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

///\todo move into cell_based/test/cell_based_pde
///\todo merge content into TestSimulationsWith*DomainPdeModifier and remove this test suite
class TestOnLatticeSimulationWithPdes : public AbstractCellBasedWithTimingsTestSuite
{
  public:

    void TestPottsBasedWithCoarseMesh()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsBasedCellPopulationWithPdes");
        simulator.SetEndTime(0.1);

        // Create PDE and boundary condition objects (zero uptake to check analytic solution)
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        c_vector<double,2> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<2> lower(centroid(0)-25.0, centroid(1)-25.0);
        ChastePoint<2> upper(centroid(0)+25.0, centroid(1)+25.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

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
        TetrahedralMesh<2,2>* p_coarse_mesh = p_pde_modifier->GetFeMesh();
        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();

        // Test that the two centres match
        c_vector<double,2> centre_diff = centre_of_cell_population - centre_of_coarse_pde_mesh;
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);
    }

    void TestCaBasedWithoutCoarseMesh()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 100);

        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<100; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaBasedCellPopulationWithPdesOnNaturalMesh");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);

        // Create PDE and boundary condition objects (zero uptake to check analytic solution)
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower = cell_population.rGetMesh().CalculateBoundingBox().rGetLowerCorner();
        ChastePoint<2> upper = cell_population.rGetMesh().CalculateBoundingBox().rGetUpperCorner();
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create update rules and pass to the simulation
        MAKE_PTR(DiffusionCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddUpdateRule(p_diffusion_update_rule);

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

        // Test coarse mesh has the same nodes as the PottsMesh
        TetrahedralMesh<2,2>* p_coarse_mesh = p_pde_modifier->GetFeMesh();

        TS_ASSERT_EQUALS(p_coarse_mesh->GetNumNodes(),p_mesh->GetNumNodes());
        TS_ASSERT_DELTA(p_coarse_mesh->GetWidth(0), p_mesh->GetWidth(0), 1e-8);
        TS_ASSERT_DELTA(p_coarse_mesh->GetWidth(1), p_mesh->GetWidth(1), 1e-8);

        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(p_coarse_mesh->GetNode(i)->rGetLocation()[0], p_mesh->GetNode(i)->rGetLocation()[0], 1e-8);
            TS_ASSERT_DELTA(p_coarse_mesh->GetNode(i)->rGetLocation()[1], p_mesh->GetNode(i)->rGetLocation()[1], 1e-8);
        }
    }

    // In this test there are only 9 cells but 100 lattice sites
    void TestCaBasedWithCellwiseSourceEllipticPde()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<unsigned> location_indices;
        for (unsigned i=3; i<6; i++)
        {
            location_indices.push_back(i+30);
            location_indices.push_back(i+40);
            location_indices.push_back(i+50);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaBasedCellPopulationWithPdeOnNaturalMesh");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
        }
    }

    // In this test there are only 50 cells but 100 lattice sites
    void TestCaBasedWithCellwiseSourceEllipticPdeWithTightBoundaries()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<unsigned> location_indices;
        for (unsigned i=2; i<8; i++)
        {
            location_indices.push_back(i+20);
            location_indices.push_back(i+30);
            location_indices.push_back(i+40);
            location_indices.push_back(i+50);
            location_indices.push_back(i+60);
            location_indices.push_back(i+70);
        }

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, location_indices.size());

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaBasedCellPopulationWithPdeOnNaturalMeshTightBoundaries");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

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

    void TestCaBasedWithoutCoarseMeshUsingPdeHandlerOnCuboid()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 100);

        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<100; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaBasedCellPopulationWithPdesOnCuboid");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);

        // Create PDE and boundary condition objects (use zero uptake to check analytic solution)
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower = cell_population.rGetMesh().CalculateBoundingBox().rGetLowerCorner();
        ChastePoint<2> upper = cell_population.rGetMesh().CalculateBoundingBox().rGetUpperCorner();
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Solve the system
        simulator.Solve();

        // Test that PDE solver is working correctly
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"),1.0, 1e-2);
        }
    }

    /*
     * This tests that a sensible error is thrown if the coarse mesh is too small.
     */
    void TestCaBasedCellsOutsideMesh()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 100);

        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<100; i++)
        {
            location_indices.push_back(i);
        }

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaBasedCellPopulationWithPdesWithCellsOutsideMesh");
        simulator.SetDt(0.1);
        simulator.SetEndTime(1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(5.0, 5.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Solve the system
        TS_ASSERT_THROWS_THIS(simulator.Solve(), "Point [6,0] is not in mesh - all elements tested");
    }

    void TestCaBasedWithCoarseMesh()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 1);

        std::vector<unsigned> location_indices;
        location_indices.push_back(12);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCaBasedCellPopulationWithPdes");
        simulator.SetEndTime(0.1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
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
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("nutrient");

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create update rules and pass to the simulation
        MAKE_PTR(DiffusionCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddUpdateRule(p_diffusion_update_rule);

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
        TetrahedralMesh<2,2>* p_coarse_mesh = p_pde_modifier->GetFeMesh();
        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();
        c_vector<double,2> centre_diff = centre_of_coarse_pde_mesh - centre_of_potts_mesh;

        // Test that the two centres match
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);

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

    void TestPottsBasedWithCoarseMeshTwoEquations()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6, 2, 2, 6, 2, 2);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestPottsBasedCellPopulationWithTwoPdes");
        simulator.SetEndTime(0.1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde_1, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc_1, (1.0));

        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde_2, (cell_population, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc_2, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(50.0, 50.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier_1, (p_pde_1, p_bc_1, false, p_cuboid));
        p_pde_modifier_1->SetDependentVariableName("quantity_1");

        simulator.AddSimulationModifier(p_pde_modifier_1);

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier_2, (p_pde_2, p_bc_2, false, p_cuboid));
        p_pde_modifier_2->SetDependentVariableName("quantity_2");

        simulator.AddSimulationModifier(p_pde_modifier_2);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddUpdateRule(p_adhesion_update_rule);

        // Solve the system
        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;

            // Test that PDE solver is working correctly on both pdes
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("quantity_1"), analytic_solution, 1e-2);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("quantity_2"), analytic_solution, 1e-2);
        }
#ifdef CHASTE_VTK
        // First file exists
        FileFinder vtk_file("TestPottsBasedCellPopulationWithTwoPdes/results_from_time_0/pde_results_quantity_1_0.vtu", RelativeTo::ChasteTestOutput);
        TS_ASSERT(vtk_file.Exists());
        // Check that the second VTK file for the solution has the dependent quantities
        OutputFileHandler handler("TestPottsBasedCellPopulationWithTwoPdes", false);
        VtkMeshReader<3,3> vtk_reader(handler.GetOutputDirectoryFullPath()+"results_from_time_0/pde_results_quantity_2_0.vtu");
        std::vector<double> data1;
        // There is no Oxygen
        TS_ASSERT_THROWS_CONTAINS(vtk_reader.GetPointData("Oxygen", data1), "No point data");
        TS_ASSERT(data1.empty());
#endif //CHASTE_VTK
    }

    // Under construction: Test growth of a population of cells that consumes nutrient
    void TestOnLatticeSpheroidWithNutrient()
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 4);

        std::vector<unsigned> location_indices;
        location_indices.push_back(55);
        location_indices.push_back(56);
        location_indices.push_back(65);
        location_indices.push_back(66);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestOnLatticeSpheroidWithNutrient");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(10);

        // Create PDE and boundary condition objects
        double nutrient_uptake_rate = -0.1;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, nutrient_uptake_rate));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower = cell_population.rGetMesh().CalculateBoundingBox().rGetLowerCorner();
        ChastePoint<2> upper = cell_population.rGetMesh().CalculateBoundingBox().rGetUpperCorner();
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0));
        p_pde_modifier->SetDependentVariableName("nutrient");
        p_pde_modifier->SetOutputSolutionAtPdeNodes(true);

        simulator.AddSimulationModifier(p_pde_modifier);

        // Create update rules and pass to the simulation
        MAKE_PTR(DiffusionCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddUpdateRule(p_diffusion_update_rule);

        // Solve the system
        simulator.Solve();

        // Test coarse mesh has the same nodes as the PottsMesh
        TetrahedralMesh<2,2>* p_coarse_mesh = p_pde_modifier->GetFeMesh();

        TS_ASSERT_EQUALS(p_coarse_mesh->GetNumNodes(),p_mesh->GetNumNodes());
        TS_ASSERT_DELTA(p_coarse_mesh->GetWidth(0), p_mesh->GetWidth(0), 1e-8);
        TS_ASSERT_DELTA(p_coarse_mesh->GetWidth(1), p_mesh->GetWidth(1), 1e-8);

        for (unsigned i=0; i< p_coarse_mesh->GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(p_coarse_mesh->GetNode(i)->rGetLocation()[0], p_mesh->GetNode(i)->rGetLocation()[0], 1e-8);
            TS_ASSERT_DELTA(p_coarse_mesh->GetNode(i)->rGetLocation()[1], p_mesh->GetNode(i)->rGetLocation()[1], 1e-8);
        }
    }

    void TestPdeOutput()
    {
        EXIT_IF_PARALLEL;

        OutputFileHandler handler("TestOnLatticeSpheroidWithNutrient", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_nut(results_dir + "/results.vizpdesolution", "cell_based/test/data/TestOnLatticeSpheroidWithNutrient/results.vizpdesolution");
        TS_ASSERT(comp_nut.CompareFiles());
        FileComparison(results_dir + "/results.vizpdesolution", "cell_based/test/data/TestOnLatticeSpheroidWithNutrient/results.vizpdesolution").CompareFiles();
    }
};

#endif /*TESTOFFLATTICESIMULATIONWITHPDES_HPP_*/
