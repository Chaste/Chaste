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

#ifndef TESTELLIPTICBOXDOMAINPDEMODIFIER_HPP_
#define TESTELLIPTICBOXDOMAINPDEMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "ApoptoticCellProperty.hpp"
#include "ArchiveOpener.hpp"
#include "AveragedSourceEllipticPde.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "ReplicatableVector.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "VertexBasedCellPopulation.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * In this test suite we check the solution of the AveragedPdes for each population type.
 * In each case we are solving Laplacian U = f where f is constant in different regions.
 * We test on a square with half apoptotic cells and the PDE mesh is twice the size.
 */
class TestEllipticBoxDomainPdeModifier : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestEllipticConstructor()
    {
        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-10.0, -10.0);
        ChastePoint<2> upper(10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 2.0));
        p_pde_modifier->SetDependentVariableName("averaged quantity");

        // Test that member variables are initialised correctly
        TS_ASSERT_EQUALS(p_pde_modifier->rGetDependentVariableName(), "averaged quantity");
        TS_ASSERT_DELTA(p_pde_modifier->GetStepSize(), 2.0, 1e-5);
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(), true);

        // Coverage of some set and methods
        p_pde_modifier->SetBcsOnBoxBoundary(false);
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(), false);

        // Check that the finite element mesh is correct
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumNodes(), 121u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryNodes(), 40u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumElements(), 200u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryElements(), 40u);

        ChasteCuboid<2> bounding_box = p_pde_modifier->mpFeMesh->CalculateBoundingBox();
        TS_ASSERT_DELTA(bounding_box.rGetUpperCorner()[0],  10.0, 1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetUpperCorner()[1],  10.0, 1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetLowerCorner()[0], -10.0, 1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetLowerCorner()[1], -10.0, 1e-5);

        // Coverage
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),false); // Defaults to false
        p_pde_modifier->SetOutputGradient(true);
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),true);
    }

    void TestArchiveEllipticBoxDomainPdeModifier()
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "EllipticBoxDomainPdeModifier.arch";

        // Separate scope to write the archive
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(UniformSourceEllipticPde<2>, p_pde, (-0.1));
            MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

            // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
            ChastePoint<2> lower(-10.0, -10.0);
            ChastePoint<2> upper(10.0, 10.0);
            MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

            // Create a PDE modifier and set the name of the dependent variable in the PDE
            std::vector<double> data(10);
            for (unsigned i=0; i<10; i++)
            {
                data[i] = i + 0.45;
            }
            Vec vector = PetscTools::CreateVec(data);
            EllipticBoxDomainPdeModifier<2> modifier(p_pde, p_bc, false, p_cuboid, 2.0, vector);
            modifier.SetDependentVariableName("averaged quantity");

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize via pointer
            AbstractCellBasedSimulationModifier<2,2>* const p_modifier = &modifier;
            output_arch << p_modifier;
        }

        // Separate scope to read the archive
        {
            AbstractCellBasedSimulationModifier<2,2>* p_modifier2;

            // Restore the modifier
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_modifier2;

            // Test that member variables are correct
            TS_ASSERT_EQUALS((static_cast<EllipticBoxDomainPdeModifier<2>*>(p_modifier2))->rGetDependentVariableName(), "averaged quantity");
            TS_ASSERT_DELTA((static_cast<EllipticBoxDomainPdeModifier<2>*>(p_modifier2))->GetStepSize(), 2.0, 1e-5);
            TS_ASSERT_EQUALS((static_cast<EllipticBoxDomainPdeModifier<2>*>(p_modifier2))->AreBcsSetOnBoxBoundary(), true);

            Vec solution = (static_cast<EllipticBoxDomainPdeModifier<2>*>(p_modifier2))->GetSolution();
            ReplicatableVector solution_repl(solution);

            TS_ASSERT_EQUALS(solution_repl.GetSize(), 10u);
            for (unsigned i=0; i<10; i++)
            {
                TS_ASSERT_DELTA(solution_repl[i], i + 0.45, 1e-6);
            }

            delete p_modifier2;
        }
    }

    void TestMeshBasedSquareMonolayer()
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
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

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

    void TestNodeBasedSquareMonolayer()
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
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestAveragedBoxEllipticPdeWithNodeOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8605, 1e-4);

        // Clear memory
        delete p_mesh;
    }

    void TestVertexBasedSquareMonolayer()
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
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestAveragedBoxEllipticPdeWithVertexOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[0], 0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[1], 0, 1e-4);
        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8605, 1e-4);
    }

    void TestPottsBasedSquareMonolayer()
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
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestAveragedBoxEllipticPdeWithPottsOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8605, 1e-1); // Testing against off-lattice models

        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8513, 1e-4); // Testing against on-lattice models
    }

    void TestCaBasedSquareMonolayer()
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
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(),50u);

        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<2>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestAveragedBoxEllipticPdeWithCaOnSquare");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_0 = cell_population.GetCellUsingLocationIndex(0);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[0], 0.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_0)[1], 0.0, 1e-4);
        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8605, 1e-1); // Testing against off-lattice models

        TS_ASSERT_DELTA(p_cell_0->GetCellData()->GetItem("variable"), 0.8513, 1e-4); // Testing against on-lattice models
    }

    void TestEllipticBoxDomainPdeModifierIn1d()
    {
        // Create mesh
        std::vector<Node<1>*> nodes;
        nodes.push_back(new Node<1>(0, true,  0.0));
        nodes.push_back(new Node<1>(1, false, 1.0));
        nodes.push_back(new Node<1>(2, false, 2.0));
        nodes.push_back(new Node<1>(3, false, 3.0));
        nodes.push_back(new Node<1>(4, true,  4.0));

        NodesOnlyMesh<1> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 1> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Set up cell population
        NodeBasedCellPopulation<1> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<1>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<1>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<1> lower(-5.0);
        ChastePoint<1> upper(15.0);
        MAKE_PTR_ARGS(ChasteCuboid<1>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<1>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestEllipticBoxDomainPdeModifierIn1d");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_2 = cell_population.GetCellUsingLocationIndex(2);
        TS_ASSERT_DELTA(p_cell_2->GetCellData()->GetItem("variable"), 0.3172, 1e-2);
        TS_ASSERT_DELTA(p_cell_2->GetCellData()->GetItem("variable_grad_x"), -0.0250, 1e-2);

        // Clear memory
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestEllipticBoxDomainPdeModifierIn3d()
    {
        // Create a simple mesh
        TetrahedralMesh<3,3> temp_mesh;
        temp_mesh.ConstructCuboid(3,3,3);
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Set up cell population
        NodeBasedCellPopulation<3> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<3>, p_pde, (cell_population, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<3>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<3> lower(-5.0, -5.0, -5.0);
        ChastePoint<3> upper(15.0, 15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<3>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier<3>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestEllipticBoxDomainPdeModifierIn3d");

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_62 = cell_population.GetCellUsingLocationIndex(13);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable"), 0.8841, 1e-2);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable_grad_x"), -0.0117, 1e-2);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable_grad_y"), 0.0510, 1e-2);
        TS_ASSERT_DELTA(p_cell_62->GetCellData()->GetItem("variable_grad_z"), -0.0435, 1e-2);
    }
};

#endif /*TESTELLIPTICBOXDOMAINPDEMODIFIER_HPP_*/
