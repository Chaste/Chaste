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

#ifndef TESTELLIPTICBOXDOMAINMODIFIERMETHODS_HPP_
#define TESTELLIPTICBOXDOMAINMODIFIERMETHODS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/math/special_functions/bessel.hpp>

#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "UniformCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellsGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
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
 * In this test suite we check the solution of the AveragedPdes for each population type.
 * In each case we are solving Laplacian U = f where f is constant in different regions.
 * We test on a square with half apoptotic cells and the PDE mesh is twice the size.
 */
class TestEllipticBoxDomainModifierMethods : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestMeshBasedSquareMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(10,10,0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0)<5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Make the PDE and BCs
        AveragedSourceEllipticPde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        MAKE_PTR_ARGS(PdeAndBoundaryConditions<2>, p_pde_and_bc, (&pde, &bc, false));
        p_pde_and_bc->SetDependentVariableName("variable");

        // Create domain
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        ChasteCuboid<2> cuboid(lower, upper);

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde_and_bc, cuboid));

        // For coverage output the solution gradient
        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestAveragedBoxEllipticPdeWithMeshOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[0], 0, .01e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[1], 0.0, 1e-4);
        TS_ASSERT_DELTA( p_cell_0->GetCellData()->GetItem("variable"), 0.8605, 1e-4);

        TS_ASSERT_DELTA( p_cell_0->GetCellData()->GetItem("variable_grad_x"), -0.0323, 1e-4);
        TS_ASSERT_DELTA( p_cell_0->GetCellData()->GetItem("variable_grad_y"), -0.0179, 1e-4);
    }

    void TestNodeBasedSquareMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(10,10,0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0)<5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Make the PDE and BCs
        AveragedSourceEllipticPde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        MAKE_PTR_ARGS(PdeAndBoundaryConditions<2>, p_pde_and_bc, (&pde, &bc, false));
        p_pde_and_bc->SetDependentVariableName("variable");

        // Make domain
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        ChasteCuboid<2> cuboid(lower, upper);

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde_and_bc,cuboid));
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedBoxEllipticPdeWithNodeOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8605, 1e-4);

        // Clear memory
        delete p_mesh;
    }

    void TestVertexBasedSquareMonolayer() throw (Exception)
    {
        HoneycombVertexMeshGenerator generator(10,10);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        p_mesh->Translate(-0.5,-sqrt(3.0)/3); // Shift so cells are on top of those in the above centre based tests

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0)<5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Make the PDE and BCs
        AveragedSourceEllipticPde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        MAKE_PTR_ARGS(PdeAndBoundaryConditions<2>, p_pde_and_bc, (&pde, &bc, false));
        p_pde_and_bc->SetDependentVariableName("variable");

        // Make domain
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        ChasteCuboid<2> cuboid(lower, upper);

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde_and_bc, cuboid));
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedBoxEllipticPdeWithVertexOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[0], 0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[1], 0, 1e-4);
        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8605, 1e-4);
    }

    void TestPottsBasedSquareMonolayer() throw (Exception)
    {
        PottsMeshGenerator<2> generator(50, 10, 4, 50, 10, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Translate and scale so cells are on top of those in the above centre based tests
        p_mesh->Translate(-6.5,-6.5);
        p_mesh->Scale(0.25,0.25 *sqrt(3.0)*0.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0)<5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Make the PDE and BCs
        AveragedSourceEllipticPde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        MAKE_PTR_ARGS(PdeAndBoundaryConditions<2>, p_pde_and_bc, (&pde, &bc, false));
        p_pde_and_bc->SetDependentVariableName("variable");

        // Make domain
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        ChasteCuboid<2> cuboid(lower, upper);

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde_and_bc, cuboid));
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedBoxEllipticPdeWithPottsOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8605, 1e-1); // Testing against off-lattice models

        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8513, 1e-4); // Testing against on-lattice models
    }

    void TestCaBasedSquareMonolayer() throw (Exception)
    {
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Scale so cells are on top of those in the above centre based tests
        p_mesh->Scale(1.0,sqrt(3.0)*0.5);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<100; i++)
        {
            location_indices.push_back(i);
        }

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i =0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0)<5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Make the PDE and BCs
        AveragedSourceEllipticPde<2> pde(cell_population, -0.1);
        ConstBoundaryCondition<2> bc(1.0);
        MAKE_PTR_ARGS(PdeAndBoundaryConditions<2>, p_pde_and_bc, (&pde, &bc, false));
        p_pde_and_bc->SetDependentVariableName("variable");

        // Make domain
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        ChasteCuboid<2> cuboid(lower, upper);

        // Create a PDE modifier object using this PDE and BCs object
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde_and_bc, cuboid));
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedBoxEllipticPdeWithCaOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8605, 1e-1); // Testing against off-lattice models

        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8513, 1e-4); // Testing against on-lattice models
    }
};

#endif /*TESTELLIPTICBOXDOMAINMODIFIERMETHODS_HPP_*/
