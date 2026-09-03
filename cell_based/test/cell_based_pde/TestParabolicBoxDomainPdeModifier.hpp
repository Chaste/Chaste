/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef TESTPARABOLICBOXDOMAINPDEMODIFIER_HPP_
#define TESTPARABOLICBOXDOMAINPDEMODIFIER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "ApoptoticCellProperty.hpp"
#include "ArchiveOpener.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "CaBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "ReplicatableVector.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VtkMeshWriter.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/*
 * In this test suite we check the solution of the AveragedParabolicPdes for each population type.
 * In each case we are solving Laplacian U = f where f is constant in different regions.
 * We test on a square with half apoptotic cells and the PDE mesh is twice the size.
 * Note all off-lattice results are the same, and all on-lattice ones are the same as each other.
 */
class TestParabolicBoxDomainPdeModifier : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestParabolicConstructor()
    {
        // Create PDE and boundary condition objects
        double constant_coefficient = -0.1;
        double linear_coefficient = -0.2;
        double diffusion_coefficient = 0.1;
        double rate_coefficient = 0.1;
        MAKE_PTR_ARGS(UniformSourceParabolicPde<2>, p_pde, (constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-10.0, -10.0);
        ChastePoint<2> upper(10.0, 10.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 2.0));
        p_pde_modifier->SetDependentVariableName("averaged quantity");

        // Test that member variables are initialised correctly
        TS_ASSERT_EQUALS(p_pde_modifier->rGetDependentVariableName(), "averaged quantity");
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(), true); //default
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoundingSphere(), false); //default
        TS_ASSERT_EQUALS(p_pde_modifier->GetUseVoronoiCellsForInterpolation(), false); //default
        TS_ASSERT_EQUALS(p_pde_modifier->GetMoveSolutionWithCells(), false); //default

        // Check mesh
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumNodes(),121u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryNodes(),40u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumElements(),200u);
        TS_ASSERT_EQUALS(p_pde_modifier->mpFeMesh->GetNumBoundaryElements(),40u);

        ChasteCuboid<2> bounding_box = p_pde_modifier->mpFeMesh->CalculateBoundingBox();
        TS_ASSERT_DELTA(bounding_box.rGetUpperCorner()[0],10,1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetUpperCorner()[1],10,1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetLowerCorner()[0],-10,1e-5);
        TS_ASSERT_DELTA(bounding_box.rGetLowerCorner()[1],-10,1e-5);

        // Coverage of member varible methods
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),false); // Defaults to false
        p_pde_modifier->SetOutputGradient(true);
        TS_ASSERT_EQUALS(p_pde_modifier->GetOutputGradient(),true);

        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(), true); // Defaults to true
        p_pde_modifier->SetBcsOnBoxBoundary(false);
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(),false);

        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoundingSphere(),false); // Defaults to false
        p_pde_modifier->SetBcsOnBoundingSphere(true);
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoundingSphere(),true);

        TS_ASSERT_EQUALS(p_pde_modifier->GetUseVoronoiCellsForInterpolation(),false); // Defaults to false
        p_pde_modifier->SetUseVoronoiCellsForInterpolation(true);
        TS_ASSERT_EQUALS(p_pde_modifier->GetUseVoronoiCellsForInterpolation(),true);

        TS_ASSERT_EQUALS(p_pde_modifier->GetMoveSolutionWithCells(),false); // Defaults to false
        p_pde_modifier->SetMoveSolutionWithCells(true);
        TS_ASSERT_EQUALS(p_pde_modifier->GetMoveSolutionWithCells(),true);
    }

    void TestArchiveParabolicBoxDomainPdeModifier()
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ParabolicBoxDomainPdeModifier.arch";

        // Separate scope to write the archive
        {

            // Create PDE and boundary condition objects
            double constant_coefficient = -0.1;
            double linear_coefficient = -0.2;
            double diffusion_coefficient = 0.1;
            double rate_coefficient = 0.1;
            MAKE_PTR_ARGS(UniformSourceParabolicPde<2>, p_pde, (constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient));
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
            ParabolicBoxDomainPdeModifier<2> modifier(p_pde, p_bc, false, p_cuboid, 2.0, vector);
            modifier.SetDependentVariableName("averaged quantity");
            modifier.SetMoveSolutionWithCells(true);
            modifier.SetBcsOnBoundingSphere(true);
            modifier.SetUseVoronoiCellsForInterpolation(true);

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
            TS_ASSERT_EQUALS((static_cast<ParabolicBoxDomainPdeModifier<2>*>(p_modifier2))->rGetDependentVariableName(), "averaged quantity");
            TS_ASSERT_DELTA((static_cast<ParabolicBoxDomainPdeModifier<2>*>(p_modifier2))->GetStepSize(), 2.0, 1e-5);
            TS_ASSERT_EQUALS((static_cast<ParabolicBoxDomainPdeModifier<2>*>(p_modifier2))->AreBcsSetOnBoxBoundary(), true);
            TS_ASSERT_EQUALS((static_cast<ParabolicBoxDomainPdeModifier<2>*>(p_modifier2))->AreBcsSetOnBoundingSphere(), true);
            TS_ASSERT_EQUALS((static_cast<ParabolicBoxDomainPdeModifier<2>*>(p_modifier2))->GetUseVoronoiCellsForInterpolation(), true);
            TS_ASSERT_EQUALS((static_cast<ParabolicBoxDomainPdeModifier<2>*>(p_modifier2))->GetMoveSolutionWithCells(),true);

            Vec solution = (static_cast<ParabolicBoxDomainPdeModifier<2>*>(p_modifier2))->GetSolution();
            ReplicatableVector solution_repl(solution);

            TS_ASSERT_EQUALS(solution_repl.GetSize(), 10u);
            for (unsigned i=0; i<10; i++)
            {
                TS_ASSERT_DELTA(solution_repl[i], i + 0.45, 1e-6);
            }

            delete p_modifier2;
        }
    }

    void TestSetupSolveStoresOldCellLocationsWhenMovingSolutionWithCells()
    {
        // Coverage: SetupSolve() only stores mOldCellLocations when mMoveSolutionWithCells is set;
        // elsewhere in this suite the flag is only used via the getter/setter or an archiving
        // round-trip, never through a real call to SetupSolve().
        HoneycombMeshGenerator generator(3, 3, 0);
        boost::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes());
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->GetCellData()->SetItem("variable", 1.0);
        }

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        double constant_coefficient = -0.1;
        double linear_coefficient = -0.2;
        double diffusion_coefficient = 0.1;
        double rate_coefficient = 0.1;
        MAKE_PTR_ARGS(UniformSourceParabolicPde<2>, p_pde, (constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");
        p_pde_modifier->SetMoveSolutionWithCells(true);

        TS_ASSERT_EQUALS(p_pde_modifier->mOldCellLocations.size(), 0u);
        p_pde_modifier->SetupSolve(cell_population, "TestParabolicBoxDomainPdeModifierStoresOldCellLocations");
        TS_ASSERT_EQUALS(p_pde_modifier->mOldCellLocations.size(), cell_population.GetNumRealCells());

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> stored_location = p_pde_modifier->mOldCellLocations[*cell_iter];
            c_vector<double, 2> actual_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            TS_ASSERT_DELTA(stored_location[0], actual_location[0], 1e-6);
            TS_ASSERT_DELTA(stored_location[1], actual_location[1], 1e-6);
        }

        // Clear memory
        delete p_mesh;
    }

    void TestInterpolateSolutionFromCellMovementReproducesFieldWhenCellsHaveNotMoved()
    {
        // InterpolateSolutionFromCellMovement() builds a mesh from cell centres, deforms a copy of
        // the FE mesh according to how those centres have moved since mOldCellLocations was last
        // recorded, then interpolates the previous solution from the deformed mesh back onto the
        // (fixed) original FE mesh nodes.
        //
        // We check a mathematically exact special case here, rather than trusting a value obtained
        // by simply running the code once: if every cell's centre is unchanged since the last
        // timestep, every interpolated displacement is exactly zero (a weighted average of zero
        // vectors), so the "deformed" mesh is geometrically identical to the original FE mesh, and
        // the method must return the input solution field completely unchanged.
        //
        // A non-zero *uniform* cell displacement was considered instead of this (it is also exactly
        // hand-computable, since piecewise-linear FE interpolation reproduces an affine function
        // exactly) but is unsafe to use: the method's second interpolation stage restricts its
        // containing-element search to each node's own originally-incident elements, and a corner
        // of a rectangular FE mesh is incident to only one or two elements spanning at most its
        // local 90-degree angle. The four corners of the box need that restricted search to succeed
        // in four different (mutually exclusive, one per quadrant) directions, so *any* non-zero
        // uniform translation of the whole mesh is guaranteed to fail this search at one corner or
        // more - hitting the method's assert(0) on a Debug build. The zero-displacement case has no
        // such failure mode: every node's own position is trivially "contained" (as a vertex, with
        // weight 1) by its own unmoved incident elements.
        MAKE_PTR_ARGS(UniformSourceParabolicPde<2>, p_pde, (0.0, 0.0, 1.0, 0.0));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(4.0, 2.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0));
        p_pde_modifier->SetDependentVariableName("variable");
        p_pde_modifier->SetMoveSolutionWithCells(true);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        // Four cells forming a box that comfortably contains the FE mesh box [0,4]x[0,2] in its
        // convex hull, so every FE mesh node has a containing element in the cell mesh.
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, -2.0, -2.0));
        nodes.push_back(new Node<2>(1, false,  6.0, -2.0));
        nodes.push_back(new Node<2>(2, false,  6.0,  4.0));
        nodes.push_back(new Node<2>(3, false, -2.0,  4.0));
        NodesOnlyMesh<2> cell_mesh;
        cell_mesh.ConstructNodesWithoutMesh(nodes, 20.0);

        std::vector<CellPtr> cells;
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, cell_mesh.GetNumNodes());

        NodeBasedCellPopulation<2> cell_population(cell_mesh, cells);

        // Cells have not moved: their recorded "old" locations equal their current locations.
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            p_pde_modifier->mOldCellLocations[*cell_iter] = cell_population.GetLocationOfCellCentre(*cell_iter);
        }

        // A real, non-constant "previous" solution field: each FE node's own x-coordinate.
        unsigned num_fe_nodes = p_pde_modifier->mpFeMesh->GetNumNodes();
        std::vector<double> solution_values(num_fe_nodes);
        for (unsigned node_index=0; node_index<num_fe_nodes; node_index++)
        {
            solution_values[node_index] = p_pde_modifier->mpFeMesh->GetNode(node_index)->rGetLocation()[0];
        }
        Vec solution = PetscTools::CreateVec(solution_values);
        p_pde_modifier->mSolution = solution;

        Vec interpolated_solution = p_pde_modifier->InterpolateSolutionFromCellMovement(cell_population);
        ReplicatableVector interpolated_solution_repl(interpolated_solution);

        for (unsigned node_index=0; node_index<num_fe_nodes; node_index++)
        {
            TS_ASSERT_DELTA(interpolated_solution_repl[node_index], solution_values[node_index], 1e-9);
        }

        PetscTools::Destroy(interpolated_solution);

        // Tidy up (ConstructNodesWithoutMesh() copies the nodes rather than taking ownership)
        delete nodes[0];
        delete nodes[1];
        delete nodes[2];
        delete nodes[3];
    }

    void TestUpdateAtEndOfTimeStepWhenMovingSolutionWithCells()
    {
        // Coverage: exercises the mMoveSolutionWithCells branch via a real UpdateAtEndOfTimeStep()
        // call (elsewhere only hit directly/via SetupSolve()). Needs an averaged source PDE, since
        // UpdateAtEndOfTimeStep() asserts HasAveragedSourcePde(). Cells don't move here, so
        // displacement is exactly zero - the same safe case as the identity test above.
        double constant_coefficient = 0.0;
        double linear_coefficient = -1.0;
        double diffusion_coefficient = 1.0;
        double rate_coefficient = 0.1;
        bool scale_by_cell_volume = false;

        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(4.0, 2.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Four cells inside the FE mesh box [0,4]x[0,2] - unlike
        // TestInterpolateSolutionFromCellMovementReproducesFieldWhenCellsHaveNotMoved() above, this
        // test also exercises UpdateCellPdeElementMap(), which requires every cell centre to have a
        // containing element in mpFeMesh itself.
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 1.0, 0.5));
        nodes.push_back(new Node<2>(1, false, 3.0, 0.5));
        nodes.push_back(new Node<2>(2, false, 1.0, 1.5));
        nodes.push_back(new Node<2>(3, false, 3.0, 1.5));
        NodesOnlyMesh<2> cell_mesh;
        cell_mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, cell_mesh.GetNumNodes());
        for (unsigned i=0; i<cells.size(); i++)
        {
            cells[i]->GetCellData()->SetItem("variable", 1.0);
        }

        NodeBasedCellPopulation<2> cell_population(cell_mesh, cells);

        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient, scale_by_cell_volume));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid, 1.0));
        p_pde_modifier->SetDependentVariableName("variable");
        p_pde_modifier->SetMoveSolutionWithCells(true);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

        p_pde_modifier->SetupSolve(cell_population, "TestUpdateAtEndOfTimeStepWhenMovingSolutionWithCells");
        TS_ASSERT_EQUALS(p_pde_modifier->mOldCellLocations.size(), cell_population.GetNumRealCells());

        // Kill one cell before UpdateAtEndOfTimeStep(), so a correctly-refreshed mOldCellLocations
        // must shrink to match the reduced population - a stale or no-op refresh would still show
        // the pre-kill count of 4. This is what actually distinguishes "the refresh ran" from "the
        // values already happened to be right" (see NOTE below).
        cell_population.Begin()->Kill();
        cell_population.RemoveDeadCells();
        cell_population.Update();
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 3u);

        SimulationTime::Instance()->IncrementTimeOneStep();
        // NOTE: with zero displacement, InterpolateSolutionFromCellMovement() (line 88) is a
        // mathematical no-op, so this only proves the call doesn't crash - not that its result is
        // actually used. Proving that would need nonzero cell displacement, which risks the
        // corner-search failure documented in the identity test above; not attempted here.
        TS_ASSERT_THROWS_NOTHING(p_pde_modifier->UpdateAtEndOfTimeStep(cell_population));

        // Unlike line 88, this IS a meaningful check: mOldCellLocations must have shrunk to 3 to
        // match the post-kill population, which only happens if the refresh actually ran.
        TS_ASSERT_EQUALS(p_pde_modifier->mOldCellLocations.size(), cell_population.GetNumRealCells());
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            c_vector<double, 2> stored_location = p_pde_modifier->mOldCellLocations[*cell_iter];
            c_vector<double, 2> actual_location = cell_population.GetLocationOfCellCentre(*cell_iter);
            TS_ASSERT_DELTA(stored_location[0], actual_location[0], 1e-6);
            TS_ASSERT_DELTA(stored_location[1], actual_location[1], 1e-6);
        }

        // Tidy up (ConstructNodesWithoutMesh() copies the nodes rather than taking ownership)
        delete nodes[0];
        delete nodes[1];
        delete nodes[2];
        delete nodes[3];
    }

    void TestMeshBasedSquareMonolayer()
    {
        HoneycombMeshGenerator generator(10,10,0);
        std::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        std::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                       cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 50u);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        double constant_coefficient = 0.0;
        double linear_coefficient = -1.0;
        double diffusion_coefficient = 1.0;
        double rate_coefficient = 0.1;
        bool scale_by_cell_volume = false;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient, scale_by_cell_volume));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedParabolicPdeWithMeshOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_55 = cell_population.GetCellUsingLocationIndex(55);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[0], 5.5, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[1], 4.3301, 1e-4);

        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable"), 0.0799, 1e-4);
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable_grad_x"), -0.0734, 1e-4);
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable_grad_y"), -0.001, 1e-4);

        // Check this was for bcs on box boundary
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(), true);
    }

    void TestMeshBasedSquareMonolayerWithBCsOnCells()
    {
        HoneycombMeshGenerator generator(10,10,0);
        std::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        std::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                       cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 50u);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        double constant_coefficient = 0.0;
        double linear_coefficient = -1.0;
        double diffusion_coefficient = 1.0;
        double rate_coefficient = 0.1;
        bool scale_by_cell_volume = false;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient, scale_by_cell_volume));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        // Change where BCS are applied here we want the bcs applied on the boundary of the cells.
        // Good for compact tissues
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(), true);
        p_pde_modifier->SetBcsOnBoxBoundary(false);
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(), false);

        // Note we also set the radius to be 0.8 so have a confluent PDE domain.
        TS_ASSERT_EQUALS(p_pde_modifier->GetTypicalCellRadius(), 0.5);
        p_pde_modifier->SetTypicalCellRadius(0.8);
        TS_ASSERT_EQUALS(p_pde_modifier->GetTypicalCellRadius(), 0.8);

        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestAveragedParabolicPdeWithMeshOnSquareBcsOnCells");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_55 = cell_population.GetCellUsingLocationIndex(55);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[0], 5.5, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[1], 4.3301, 1e-4);

        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable"), 0.1028, 1e-4);
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable_grad_x"), -0.0819, 1e-4);
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable_grad_y"), 0.0008, 1e-4);

    }

void TestMeshBasedSquareMonolayerWithBCsOnBoundingSpehre()
    {
        HoneycombMeshGenerator generator(10,10,0);
        std::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        std::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                       cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 50u);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        double constant_coefficient = 0.0;
        double linear_coefficient = -1.0;
        double diffusion_coefficient = 1.0;
        double rate_coefficient = 0.1;
        bool scale_by_cell_volume = false;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient, scale_by_cell_volume));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        // Change where BCS are applied here we want the bcs applied on the boundary of the cells.
        // Good for compact tissues
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(), true);
        p_pde_modifier->SetBcsOnBoxBoundary(false);
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoxBoundary(), false);

        // Apply BCs on the bounding sphere.
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoundingSphere(), false);
        p_pde_modifier->SetBcsOnBoundingSphere(true);
        TS_ASSERT_EQUALS(p_pde_modifier->AreBcsSetOnBoundingSphere(), true);

        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);

        p_pde_modifier->SetupSolve(cell_population,"TestAveragedParabolicPdeWithMeshOnSquareBcsOnBoundingSphere");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_55 = cell_population.GetCellUsingLocationIndex(55);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[0], 5.5, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[1], 4.3301, 1e-4);

        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable"), 0.0914, 1e-4);
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable_grad_x"), -0.0764, 1e-4);
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable_grad_y"), 0.0004, 1e-4);

    }


    // Only difference from above tests (TestMeshBasedSquareMonolayer) is the use of Neuman BCs here.
    // Note we can only do this on the box not on the cells or the bounding sphere
    void TestMeshBasedSquareMonolayerWithNeumanBcs()
    {
        HoneycombMeshGenerator generator(10,10,0);
        std::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        std::shared_ptr<AbstractCellProperty> p_apoptotic_property =
            cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }

            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 50u);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        double constant_coefficient = 0.0;
        double linear_coefficient = -1.0;
        double diffusion_coefficient = 1.0;
        double rate_coefficient = 0.1;
        bool scale_by_cell_volume = false;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient, scale_by_cell_volume));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedParabolicPdeWithNeumannWithMeshOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_55 = cell_population.GetCellUsingLocationIndex(55);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[0], 5.5, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[1], 4.3301, 1e-4);

        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable"), 0.121, 1e-4);
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable_grad_x"), -0.1028, 1e-4);
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable_grad_y"), -0.0053, 1e-4);
    }

    void TestNodeBasedSquareMonolayer()
    {
        HoneycombMeshGenerator generator(10,10,0);
        std::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        std::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 50u);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        double constant_coefficient = 0.0;
        double linear_coefficient = -1.0;
        double diffusion_coefficient = 1.0;
        double rate_coefficient = 0.1;
        bool scale_by_cell_volume = false;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient, scale_by_cell_volume));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedParabolicPdeWithNodeOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_55 = cell_population.GetCellUsingLocationIndex(55);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[0], 5.5, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[1], 4.3301, 1e-4);

        // Exactly the same as all Off Lattice models as same cell cetres
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable"), 0.0799, 1e-4);

        // Clear memory
        delete p_mesh;
    }

    void TestVertexBasedSquareMonolayer()
    {
        HoneycombVertexMeshGenerator generator(10,10);
        std::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

        p_mesh->Translate(-0.5,-sqrt(3.0)/3); // Shift so cells are on top of those in the above centre based tests

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        std::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }

            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 50u);

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        double constant_coefficient = 0.0;
        double linear_coefficient = -1.0;
        double diffusion_coefficient = 1.0;
        double rate_coefficient = 0.1;
        bool scale_by_cell_volume = false;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient, scale_by_cell_volume));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedParabolicPdeWithVertexOnSquare");

        // Run for 10 timesteps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_55 = cell_population.GetCellUsingLocationIndex(55);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[0], 5.5, 1e-4);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[1], 4.3301, 1e-4);

        // Exactly the same as all Off Lattice models as same cell cetres
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable"), 0.0799, 1e-4);
    }

    void TestPottsBasedSquareMonolayer()
    {
        PottsMeshGenerator<2> generator(100, 10, 4, 100, 10, 4);
        std::shared_ptr<PottsMesh<2> > p_mesh = generator.GetMesh();

        // Translate and scale so cells are on top of those in the above centre based tests.
        p_mesh->Translate(-31.5,-31.5);
        p_mesh->Scale(0.25,0.25 *sqrt(3.0)*0.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Make cells with x<5.0 apoptotic (so no source term)
        std::shared_ptr<AbstractCellProperty> p_apoptotic_property =
                        cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetCentroidOfElement(i);
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }

            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 50u);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        double constant_coefficient = 0.0;
        double linear_coefficient = -1.0;
        double diffusion_coefficient = 1.0;
        double rate_coefficient = 0.1;
        bool scale_by_cell_volume = false;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient, scale_by_cell_volume));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedParabolicPdeWithPottsOnSquare");

        // Run for 10 timesteps
        for (unsigned i=0; i<10; i++)
           {
            SimulationTime::Instance()->IncrementTimeOneStep();
               p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
               p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_55 = cell_population.GetCellUsingLocationIndex(55);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[0], 5, 1e-4); //5.5 in Off Lattice
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[1], 4.3301, 1e-4);

        // Similar to all Off Lattice models as different cell cetres
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable"), 0.0799, 1e-2);

        // Exactly the same as all On Lattice models as same cell cetres
        TS_ASSERT_DELTA(p_cell_55->GetCellData()->GetItem("variable"), 0.0848, 1e-4);
    }

    void TestCaBasedSquareMonolayer()
    {
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        std::shared_ptr<PottsMesh<2> > p_mesh = generator.GetMesh();

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
        std::shared_ptr<AbstractCellProperty> p_apoptotic_property =
            cells[0]->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<ApoptoticCellProperty>();
        for (unsigned i=0; i<cells.size(); i++)
        {
            c_vector<double,2> cell_location;
            cell_location = p_mesh->GetNode(i)->rGetLocation();
            if (cell_location(0) < 5.0)
            {
                cells[i]->AddCellProperty(p_apoptotic_property);
            }
            // Set initial condition for PDE
            cells[i]->GetCellData()->SetItem("variable",1.0);
        }
        TS_ASSERT_EQUALS(p_apoptotic_property->GetCellCount(), 50u);

        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set up simulation time for file output
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 10);

        // Create PDE and boundary condition objects
        double constant_coefficient = 0.0;
        double linear_coefficient = -1.0;
        double diffusion_coefficient = 1.0;
        double rate_coefficient = 0.1;
        bool scale_by_cell_volume = false;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde<2>, p_pde, (cell_population, constant_coefficient, linear_coefficient, diffusion_coefficient, rate_coefficient, scale_by_cell_volume));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));

        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<2> lower(-5.0, -5.0);
        ChastePoint<2> upper(15.0, 15.0);
        MAKE_PTR_ARGS(ChasteCuboid<2>, p_cuboid, (lower, upper));

        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false, p_cuboid));
        p_pde_modifier->SetDependentVariableName("variable");

        // For coverage, output the solution gradient
        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->SetupSolve(cell_population,"TestAveragedParabolicPdeWithCaOnSquare");

        // Run for 10 time steps
        for (unsigned i=0; i<10; i++)
        {
            SimulationTime::Instance()->IncrementTimeOneStep();
            p_pde_modifier->UpdateAtEndOfTimeStep(cell_population);
            p_pde_modifier->UpdateAtEndOfOutputTimeStep(cell_population);
        }

        // Test the solution at some fixed points to compare with other cell populations
        CellPtr p_cell_55 = cell_population.GetCellUsingLocationIndex(55);
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[0], 5, 1e-4); //5.5 in Off Lattice
        TS_ASSERT_DELTA(cell_population.GetLocationOfCellCentre(p_cell_55)[1], 4.3301, 1e-4);

        // Similar to all Off Lattice models as different cell cetres
        TS_ASSERT_DELTA( p_cell_55->GetCellData()->GetItem("variable"), 0.0799, 1e-2);

        // Exactly the same as all On Lattice models as same cell cetres
        TS_ASSERT_DELTA(p_cell_55->GetCellData()->GetItem("variable"), 0.0848, 1e-4);
    }
};

#endif /*TESTPARABOLICBOXDOMAINPDEMODIFIER_HPP_*/
