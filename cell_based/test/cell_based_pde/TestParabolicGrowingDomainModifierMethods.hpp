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

#ifndef TESTPARABOLICGROWINGDOMAINMODIFIERMETHODS_HPP_
#define TESTPARABOLICGROWINGDOMAINMODIFIERMETHODS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/math/special_functions/bessel.hpp>
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "UniformCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "CaBasedCellPopulation.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * In this test suite we check the solution of the CellwiseParabolicPdes
 * against exact solutions.
 *
 * In each case we are solving du/dt = Laplacian U + f where f is constant in different regions.
 *
 * We solve unit disc where the steady state solutions are Bessels functions and logs.
 */
class TestParabolicGrowingDomainModifierMethods : public AbstractCellBasedWithTimingsTestSuite
{
public:

    /*
     * Here the exact steady state solution is
     *
     * u = J0(r)/J0(1)
     * where J0 is the zeroth order bessel fn
     */

    void TestMeshBasedMonolayerWithParabolicPde() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Set initial condition for pde
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.0, 10);

        // Make the PDE and BCs
        CellwiseSourceParabolicPde<2> pde(cell_population, 1, 1, 1);
        ConstBoundaryCondition<2> bc(1.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithMeshOnDisk");

        // Run for 100 timesteps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
         }

        /*
         * Test the solution against the exact solution for the steady state.
         */

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {

            c_vector<double,2> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double r = sqrt(cell_location(0)*cell_location(0) + cell_location(1)*cell_location(1));
            double u_exact = boost::math::cyl_bessel_j(0,r) / boost::math::cyl_bessel_j(0,1);

            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), u_exact, 1e-2);
        }
    }

    /*
     * Here the outer cells (r>1/2) are apoptotic so the  steady state solution is
     *
     * u = C*J0(r)  for r in [0,0.5]
     *     A*ln(r) + 1 for r in [0.5,1]
     *
     *  where J0 is the zeroth order bessel fn and C and A are constants
     *
     */
    void TestMeshBasedHeterogeneousMonolayerWithParabolicPde() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Make cells with r<1/2 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location = mesh.GetNode(i)->rGetLocation();
            double r = sqrt(cell_location(0)*cell_location(0) + cell_location(1)*cell_location(1));
            if (r>0.5)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }

            // Set initial condition for pde
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.0, 10);

        // Make the PDE and BCs
        CellwiseSourceParabolicPde<2> pde(cell_population, 1, 1, 1);
        ConstBoundaryCondition<2> bc(1.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithMeshOnHeterogeneousDisk");

        // Run for 10 timesteps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
         }

        /*
         * Test the solution against the exact steady state solution for the steady state
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {

            c_vector<double,2> cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double r = sqrt(cell_location(0)*cell_location(0) + cell_location(1)*cell_location(1));

            double J005 = boost::math::cyl_bessel_j(0,0.5);
            double J105 = boost::math::cyl_bessel_j(1,0.5);

            double A = -1.0/(2.0*J005/J105 +log(0.5));
            double C = -2*A/J105;

            double u_exact = C*boost::math::cyl_bessel_j(0,r);
            if ( r> 0.5 )
            {
                u_exact = A*log(r)+1.0;
            }

            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), u_exact, 1e-2);
        }
    }


    void TestMeshBasedMonolayerWithParabolicPdeAndNeumannBcs() throw (Exception)
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_984_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Set initial condition for pde
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10.0, 10);

        // Make the PDE and BCs
        CellwiseSourceParabolicPde<2> pde(cell_population, 1, 1, 1);
        ConstBoundaryCondition<2> bc(1.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, true);
        pde_and_bc.SetDependentVariableName("variable");

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithMeshNeumanBcs");

        // Run for 10 timesteps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
         }

        /*
         * Test the solution against the an approximate constant solution.
         */

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double u_approx = 3071;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), u_approx, 1.0);
        }
    }
    // Now test on a square with half apoptotic cells to compare all the population types

    void TestMeshBasedSquareMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(20,20,0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                       cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
           c_vector<double,2> cell_location = p_mesh->GetNode(i)->rGetLocation();
        if (cell_location(0)<10.0)
        {
            cells[i]->AddCellProperty(p_apoptotic_property);
        }
        // Set initial condition for pde
        cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),200u);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Make the PDE and BCs
        CellwiseSourceParabolicPde<2> pde(cell_population, 0.1, 1, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithMeshOnSquare");

        // Run for 10 timesteps
        for (unsigned i=0; i<10; i++)
           {
            SimulationTime::Instance()->IncrementTimeOneStep();
               p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
               p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA( p_cell_210->GetCellData()->GetItem("variable"), 0.6309, 1e-4);
    }

    void TestNodeBasedSquareMonolayer() throw (Exception)
    {

        HoneycombMeshGenerator generator(20,20,0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0)<10.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for pde
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),200u);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Make the PDE and BCs
        CellwiseSourceParabolicPde<2> pde(cell_population, 0.1, 1, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithNodeOnSquare");

        // Run for 10 timesteps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA( p_cell_210->GetCellData()->GetItem("variable"), 0.6309, 1e-2);
        //Checking it doesn't change for this cell population
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.6296, 1e-4);

        // Clear memory
        delete p_mesh;
    }

    void TestVertexBasedSquareMonolayer() throw (Exception)
    {
        HoneycombVertexMeshGenerator generator(20,20);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        p_mesh->Translate(-0.5,-sqrt(3.0)/3); // Shift so cells are on top of those in the above centre based tests.

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
              c_vector<double,2> cell_location = p_mesh->GetCentroidOfElement(i);
              if (cell_location(0)<10.0)
              {
                cells[i]->AddCellProperty(p_apoptotic_property);
              }
            // Set initial condition for pde
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),200u);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Make the PDE and BCs
        CellwiseSourceParabolicPde<2> pde(cell_population, 0.1, 1, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithVertexOnSquare");

        // Run for 10 timesteps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA( p_cell_210->GetCellData()->GetItem("variable"), 0.6309, 1e-1); //low error as mesh is slightlty larger than for centre based models.
          //Checking it doesn't change for this cell population
          TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.6618, 1e-4);
    }

    void TestPottsBasedSquareMonolayer() throw (Exception)
    {
        PottsMeshGenerator<2> generator(100, 20, 4, 100, 20, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Translate and scale so cells are on top of those in the above centre based tests.
        p_mesh->Translate(-11.5,-11.5);
        p_mesh->Scale(0.25,0.25 *sqrt(3.0)*0.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0)<10.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for pde
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),200u);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Make the PDE and BCs
        CellwiseSourceParabolicPde<2> pde(cell_population, 0.1, 1, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithPottsOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA( p_cell_210->GetCellData()->GetItem("variable"), 0.6309, 2e-1); // Low error as mesh is slightly larger than for centre-based models

        // Checking it doesn't change for this cell population
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.6086, 1e-4);
    }

    void TestCaBasedSquareMonolayer() throw (Exception)
    {
        PottsMeshGenerator<2> generator(20, 0, 0, 20, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Scale so cells are on top of those in the above centre-based tests
        p_mesh->Scale(1.0,sqrt(3.0)*0.5);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<400; i++)
        {
            location_indices.push_back(i);
        }

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

        // Make cells with x<10.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 10.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for pde
            cells[i]->GetCellData()->SetItem("variable", 1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 200u);

        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Make the PDE and BCs
        CellwiseSourceParabolicPde<2> pde(cell_population, 0.1, 1, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, false);
        pde_and_bc.SetDependentVariableName("variable");

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithCaOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_210 = cell_population.GetCellUsingLocationIndex(210);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[0], 10, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_210)[1], 5.0*sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA( p_cell_210->GetCellData()->GetItem("variable"), 0.6309, 2e-1); // low error as mesh is slightly larger than for centre-based models

        // Checking it doesn't change for this cell population
        TS_ASSERT_DELTA(p_cell_210->GetCellData()->GetItem("variable"), 0.6086, 2e-3);// Low threshold as slightly different on Intel
    }
};

#endif /*TESTPARABOLICGROWINGDOMAINMODIFIERMETHODS_HPP_*/
