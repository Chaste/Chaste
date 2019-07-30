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

#ifndef TESTPARABOLICGROWINGDOMAINPDEMODIFIER_HPP_
#define TESTPARABOLICGROWINGDOMAINPDEMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

// This macro prevents errors with GCC 4.8 of form "unable to find numeric literal operator 'operator"" Q'"
// when compiling with -std=gnu++11 (see #2929). \todo: remove when GCC 4.8 is no longer supported.
#define BOOST_MATH_DISABLE_FLOAT128
#include <boost/math/special_functions/bessel.hpp>

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "ApoptoticCellProperty.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "ReplicatableVector.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "VertexBasedCellPopulation.hpp"

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
class TestParabolicGrowingDomainPdeModifier : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestParabolicConstructor()
    {
        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(UniformSourceParabolicPde<2>, p_pde, (-0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("averaged quantity");

        // Test that member variables are initialised correctly
        TS_ASSERT_EQUALS(p_pde_modifier->rGetDependentVariableName(), "averaged quantity");
    }

    void TestGrowingDomainPdeModifierExceptions()
    {
        EXIT_IF_PARALLEL;

        // Create a simple mesh
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        TetrahedralMesh<2,2> temp_mesh;
        temp_mesh.ConstructFromMeshReader(mesh_reader);
        temp_mesh.Scale(5.0,1.0);

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(temp_mesh, 1.5);

        // Set up cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
            p_model->SetDimension(2);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = -RandomNumberGenerator::Instance()->ranf()*18.0;
            p_cell->SetBirthTime(birth_time);

            cells.push_back(p_cell);
        }

        // Set up cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, -1.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("nutrient");

        TS_ASSERT_THROWS_THIS(p_pde_modifier->SetupSolve(cell_population, "output_directory"),
            "ParabolicGrowingDomainPdeModifier cannot be used with an AveragedSourceParabolicPde. Use a ParabolicBoxDomainPdeModifier instead.");
    }

    void TestArchiveParabolicGrowingDomainPdeModifier()
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ParabolicGrowingDomainPdeModifier.arch";

        // Separate scope to write the archive
        {
            // Create PDE and boundary condition objects
            MAKE_PTR_ARGS(UniformSourceParabolicPde<2>, p_pde, (-0.1));
            MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

            // Create a PDE modifier and set the name of the dependent variable in the PDE
            std::vector<double> data(10);
            for (unsigned i=0; i<10; i++)
            {
                data[i] = i + 0.45;
            }
            Vec vector = PetscTools::CreateVec(data);
            ParabolicGrowingDomainPdeModifier<2> modifier(p_pde, p_bc, false, vector);
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

            // See whether we read out the correct variable name
            std::string variable_name = (static_cast<ParabolicGrowingDomainPdeModifier<2>*>(p_modifier2))->rGetDependentVariableName();
            TS_ASSERT_EQUALS(variable_name, "averaged quantity");

            Vec solution = (static_cast<ParabolicGrowingDomainPdeModifier<2>*>(p_modifier2))->GetSolution();
            ReplicatableVector solution_repl(solution);

            TS_ASSERT_EQUALS(solution_repl.GetSize(), 10u);
            for (unsigned i=0; i<10; i++)
            {
                TS_ASSERT_DELTA(solution_repl[i], i + 0.45, 1e-6);
            }

            delete p_modifier2;
        }
    }

    void TestMeshBasedMonolayerWithParabolicPde()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Set initial condition for PDE
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.0, 4);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, 1, 1, 1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithMeshOnDisk");

        // Run for 4 time steps
        for (unsigned i=0; i<4; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
         }

        /*
         * Test the solution against the exact solution,
         *
         * u(r,t) = J0(r)/J0(1),
         *
         * where J0 is the zeroth order Bessel function.
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double,2> cell_location;
            cell_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            double r = sqrt(cell_location(0)*cell_location(0) + cell_location(1)*cell_location(1));
            double u_exact = boost::math::cyl_bessel_j(0,r) / boost::math::cyl_bessel_j(0,1);

            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), u_exact, 1e-2);
        }
    }

    void TestMeshBasedHeterogeneousMonolayerWithParabolicPde()
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
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = mesh.GetNode(i)->rGetLocation();
            double r = sqrt(cell_location(0)*cell_location(0) + cell_location(1)*cell_location(1));
            if (r > 0.5)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }

            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable", 1.0);
        }

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(100.0, 2);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, 1, 1, 1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithMeshOnHeterogeneousDisk");

        // Run for 5 time steps
        for (unsigned i=0; i<2; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
         }

        /*
         * Test the solution against the exact steady-state solution. Here the
         * outer cells (r > 1/2) are apoptotic, so the steady-state solution is
         *
         * u(r) = C*J0(r)  for r in [0,0.5],
         *        A*ln(r) + 1 for r in [0.5,1],
         *
         *  where J0 is the zeroth order Bessel function and C and A are constants.
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
            if (r > 0.5)
            {
                u_exact = A*log(r) + 1.0;
            }

            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), u_exact, 1e-2);
        }
    }

    void TestMeshBasedMonolayerWithParabolicPdeAndNeumannBcs()
    {
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/disk_522_elements");
        MutableMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Set initial condition for PDE
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->GetCellData()->SetItem("variable", 1.0);
        }

        MeshBasedCellPopulation<2> cell_population(mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, 1, 1, 1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithMeshNeumanBcs");

        // Run for 5 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
         }

        // Test the solution against the approximate spatially uniform solution
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double u_approx = 5.8;
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("variable"), u_approx, 0.5);
        }
    }

    // Now test on a square with half apoptotic cells to compare all the population types
    void TestMeshBasedSquareMonolayer()
    {
        HoneycombMeshGenerator generator(6,6,0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x < 3.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                       cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 3.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable", 1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 18u);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10.0, 10);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, 0.1, 1, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithMeshOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_14 = cell_population.GetCellUsingLocationIndex(14);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[0], 2.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[1], sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(p_cell_14->GetCellData()->GetItem("variable"), 0.9604, 1e-4);
    }

    void TestNodeBasedSquareMonolayer()
    {
        HoneycombMeshGenerator generator(6,6,0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x < 3.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 3.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable", 1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 18u);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10.0, 10);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, 0.1, 1, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithNodeOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_14 = cell_population.GetCellUsingLocationIndex(14);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[0], 2.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[1], sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(p_cell_14->GetCellData()->GetItem("variable"), 0.9550, 1e-4);

        // Compare with the mesh-based cell population result
        TS_ASSERT_DELTA(p_cell_14->GetCellData()->GetItem("variable"), 0.9604, 1e-1);

        // Clear memory
        delete p_mesh;
    }

    void TestVertexBasedSquareMonolayer()
    {
        HoneycombVertexMeshGenerator generator(6,6);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        p_mesh->Translate(-0.5,-sqrt(3.0)/3); // Shift so cells are on top of those in the above centre-based tests

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x < 3.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0) < 3.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }

            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 18u);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10.0, 10);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, 0.1, 1, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithVertexOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_14 = cell_population.GetCellUsingLocationIndex(14);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[0], 2.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[1], sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(p_cell_14->GetCellData()->GetItem("variable"), 0.9567, 1e-4);

        // Compare with the mesh-based cell population result
        TS_ASSERT_DELTA(p_cell_14->GetCellData()->GetItem("variable"), 0.9604, 1e-1);
    }

    void TestPottsBasedSquareMonolayer()
    {
        PottsMeshGenerator<2> generator(24, 6, 4, 24, 6, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Scale so cells are on top of those in the above centre-based tests
        p_mesh->Translate(-1.5, -1.5);
        p_mesh->Scale(0.25, sqrt(3.0)*0.5*0.25);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x < 3.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0) < 3.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 18u);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10.0, 10);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, 0.1, 1, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithPottsOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_14 = cell_population.GetCellUsingLocationIndex(14);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[0], 2.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[1], sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(p_cell_14->GetCellData()->GetItem("variable"), 0.9503, 1.2e-3); // Intel optimised build has a slightly different result here

        // Compare with the mesh-based cell population result
        TS_ASSERT_DELTA(p_cell_14->GetCellData()->GetItem("variable"), 0.9604, 1e-1);
    }

    void TestCaBasedSquareMonolayer()
    {
        PottsMeshGenerator<2> generator(6, 0, 0, 6, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Scale so cells are on top of those in the above centre-based tests
        p_mesh->Scale(1.0, 1.0*sqrt(3.0)*0.5);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<36; i++)
        {
            location_indices.push_back(i);
        }

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

        // Make cells with x < 3.0 apoptotic (so no source term)
        boost::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 3.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable", 1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 18u);

        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(10.0, 10);

        // Create PDE and boundary condition objects
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, 0.1, 1, -0.1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false));
        p_pde_modifier->SetDependentVariableName("variable");

        p_pde_modifier->SetupSolve(cell_population,"TestCellwiseParabolicPdeWithCaOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_14 = cell_population.GetCellUsingLocationIndex(14);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[0], 2.0, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_14)[1], sqrt(3.0), 1e-4);
        TS_ASSERT_DELTA(p_cell_14->GetCellData()->GetItem("variable"), 0.9503, 4e-3); // Low threshold as slightly different on Intel

        // Compare with the mesh-based cell population result (low error as mesh is slightly larger than for centre-based models)
        TS_ASSERT_DELTA(p_cell_14->GetCellData()->GetItem("variable"), 0.9604, 1e-1);
    }
};

#endif /*TESTPARABOLICGROWINGDOMAINPDEMODIFIER_HPP_*/
