/*

Copyright (c) 2005-2025, University of Oxford.
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

#ifndef TESTIMMERSEDBOUNDARYSIMULATIONMODIFIER_HPP_
#define TESTIMMERSEDBOUNDARYSIMULATIONMODIFIER_HPP_

// Needed for the test environment
#include "AbstractCellBasedTestSuite.hpp"

// Includes from trunk
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "FluidSource.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryLinearInteractionForce.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <limits>

// This test is never run in parallel
#include "FakePetscSetup.hpp"

class TestImmersedBoundarySimulationModifier : public AbstractCellBasedTestSuite
{
public:
    void TestConstructorAndGetAndSetMethods()
    {
        ImmersedBoundarySimulationModifier<2> modifier;

        // Test GetNodeNeighbourUpdateFrequency() and SetNodeNeighbourUpdateFrequency()
        TS_ASSERT_EQUALS(modifier.GetNodeNeighbourUpdateFrequency(), 1u);
        modifier.SetNodeNeighbourUpdateFrequency(2);
        TS_ASSERT_EQUALS(modifier.GetNodeNeighbourUpdateFrequency(), 2u);

        // Test GetReynoldsNumber() and SetReynoldsNumber()
        TS_ASSERT_DELTA(modifier.GetReynoldsNumber(), 1e-4, 1e-6);
        modifier.SetReynoldsNumber(1e-5);
        TS_ASSERT_DELTA(modifier.GetReynoldsNumber(), 1e-5, 1e-6);
    }

    void TestOutputParametersWithImmersedBoundarySimulationModifier()
    {
        std::string output_directory = "TestOutputParametersWithImmersedBoundarySimulationModifier";
        OutputFileHandler output_file_handler(output_directory, false);

        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_modifier);
        TS_ASSERT_EQUALS(p_modifier->GetIdentifier(), "ImmersedBoundarySimulationModifier-2");

        out_stream modifier_parameter_file = output_file_handler.OpenOutputFile("ImmersedBoundarySimulationModifier.parameters");
        p_modifier->OutputSimulationModifierParameters(modifier_parameter_file);
        modifier_parameter_file->close();

        // Compare the generated file in test output with a reference copy in the source code
        FileFinder generated = output_file_handler.FindFile("ImmersedBoundarySimulationModifier.parameters");
        FileFinder reference("cell_based/test/data/TestSimulationModifierOutputParameters/ImmersedBoundarySimulationModifier.parameters",
                             RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated, reference);
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestSetupConstantMemberVariables()
    {
        // Set up SimulationTime - needed by SetupConstantMemberVariables()
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        ImmersedBoundarySimulationModifier<2> modifier;

        TS_ASSERT_DELTA(modifier.mGridSpacingX, 0.0, 1e-6);
        TS_ASSERT_DELTA(modifier.mGridSpacingY, 0.0, 1e-6);
        TS_ASSERT_DELTA(modifier.mFftNorm, 0.0, 1e-6);
        TS_ASSERT(modifier.mpArrays == NULL);
        TS_ASSERT(modifier.mpFftInterface == NULL);

        // Test that the correct exception is throw if trying to use the wrong cell population class
        HoneycombVertexMeshGenerator vertex_generator(2, 2);
        MutableVertexMesh<2, 2>* p_vertex_mesh = vertex_generator.GetMesh().get();
        std::vector<CellPtr> vertex_cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> vertex_cells_generator;
        vertex_cells_generator.GenerateBasicRandom(vertex_cells, p_vertex_mesh->GetNumElements(), p_diff_type);
        VertexBasedCellPopulation<2> vertex_cell_population(*p_vertex_mesh, vertex_cells);

        TS_ASSERT_THROWS_THIS(modifier.SetupConstantMemberVariables(vertex_cell_population),
                              "Cell population must be immersed boundary");

        // Now test SetupConstantMemberVariables() using an immersed boundary cell population
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        TS_ASSERT_THROWS_NOTHING(modifier.SetupConstantMemberVariables(cell_population));
        TS_ASSERT_DELTA(modifier.mGridSpacingX, 1.0 / 128.0, 1e-6);
        TS_ASSERT_DELTA(modifier.mGridSpacingY, 1.0 / 128.0, 1e-6);
        TS_ASSERT_DELTA(modifier.mFftNorm, 128.0 * 128.0, 1e-6);
        TS_ASSERT(modifier.mpArrays != NULL);
        TS_ASSERT(modifier.mpFftInterface != NULL);
    }

    void TestClearForcesAndSources()
    {
        // Set up SimulationTime - needed by SetupConstantMemberVariables()
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create an immersed boundary cell population where each node holds a non-zero applied force
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            c_vector<double, 2> force;
            force[0] = i * 0.01;
            force[1] = 2 * i * 0.01;

            p_mesh->GetNode(i)->ClearAppliedForce();
            p_mesh->GetNode(i)->AddAppliedForceContribution(force);
        }

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        // For coverage of ClearForcesAndSources(), specify that the cell population has active sources
        cell_population.SetIfPopulationHasActiveSources(true);

        // Call SetupConstantMemberVariables() first, to pass the mesh to the simulation modifier
        ImmersedBoundarySimulationModifier<2> modifier;
        modifier.SetupConstantMemberVariables(cell_population);

        // Test ClearForcesAndSources() correctly resets the applied force on each node
        modifier.ClearForcesAndSources();

        for (unsigned i = 0; i < p_mesh->GetNumNodes(); ++i)
        {
            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetAppliedForce()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetAppliedForce()[1], 0.0, 1e-6);
        }

        // Test ClearForcesAndSources() correctly resets mpArrays
        multi_array<double, 3>& r_force_grids = modifier.mpArrays->rGetModifiableForceGrids();
        multi_array<double, 3>& r_rhs_grid = modifier.mpArrays->rGetModifiableRightHandSideGrids();

        for (unsigned x = 0; x < modifier.mpMesh->GetNumGridPtsX(); ++x)
        {
            for (unsigned y = 0; y < modifier.mpMesh->GetNumGridPtsY(); ++y)
            {
                TS_ASSERT_DELTA(r_rhs_grid[2][x][y], 0.0, 1e-6);
                for (unsigned dim = 0; dim < 2; ++dim)
                {
                    TS_ASSERT_DELTA(r_force_grids[dim][x][y], 0.0, 1e-6);
                }
            }
        }
    }

    void TestPropagateForcesToFluidGrid()
    {
        {
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up and apply a force to the node
            c_vector<double, 2> force;
            force[0] = 0.01;
            force[1] = 0.01;
            nodes[0]->AddAppliedForceContribution(force);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetupConstantMemberVariables(cell_population);

            // Check that force grids are zero before force propagation
            auto& r_force_grids = modifier.mpArrays->rGetModifiableForceGrids();
            for (unsigned grid = 0; grid <= 1; ++grid)
            {
                for (unsigned x = 0; x < modifier.mpMesh->GetNumGridPtsX(); ++x)
                {
                    for (unsigned y = 0; y < modifier.mpMesh->GetNumGridPtsY(); ++y)
                    {
                        TS_ASSERT_EQUALS(r_force_grids[grid][x][y], 0.0);
                    }
                }
            }

            modifier.PropagateForcesToFluidGrid();

            // Check that there are only 16 non-zero entries per grid
            unsigned non_zero_entries_count = 0;
            for (unsigned grid = 0; grid <= 1; ++grid)
            {
                for (unsigned x = 0; x < modifier.mpMesh->GetNumGridPtsX(); ++x)
                {
                    for (unsigned y = 0; y < modifier.mpMesh->GetNumGridPtsY(); ++y)
                    {
                        if (r_force_grids[grid][x][y] != 0.0)
                        {
                            non_zero_entries_count += 1;
                        }
                    }
                }
            }
            TS_ASSERT_EQUALS(non_zero_entries_count, 32u);

            // Check the non-zero entries' values
            TS_ASSERT_DELTA(r_force_grids[0][4][4], 0.0020, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][4][5], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][4][6], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][4][7], 0.0020, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][5][4], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][5][5], 0.0707, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][5][6], 0.0707, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][5][7], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][6][4], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][6][5], 0.0707, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][6][6], 0.0707, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][6][7], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][7][4], 0.0020, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][7][5], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][7][6], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[0][7][7], 0.0020, 0.0001);

            TS_ASSERT_DELTA(r_force_grids[1][4][4], 0.0020, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][4][5], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][4][6], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][4][7], 0.0020, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][5][4], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][5][5], 0.0707, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][5][6], 0.0707, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][5][7], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][6][4], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][6][5], 0.0707, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][6][6], 0.0707, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][6][7], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][7][4], 0.0020, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][7][5], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][7][6], 0.0121, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[1][7][7], 0.0020, 0.0001);

            SimulationTime::Instance()->Destroy();
        }

        { // Coverage for zero field sums
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up and apply a force to the node
            c_vector<double, 2> force;
            force[0] = 0.01;
            force[1] = 0.01;
            nodes[0]->AddAppliedForceContribution(force);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetupConstantMemberVariables(cell_population);
            modifier.SetZeroFieldSums(true);
            modifier.PropagateForcesToFluidGrid();
        }
    }

    void TestPropagateFluidSourcesToGrid()
    {
        {
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);
            auto& r_balancing_fluid_sources = mesh.rGetBalancingFluidSources();
            r_balancing_fluid_sources.clear();

            auto p_fluid_source = std::make_shared<FluidSource<2>>(0, 0.55, 0.55);
            p_fluid_source->SetStrength(0.1);
            elems.back()->SetFluidSource(p_fluid_source);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetupConstantMemberVariables(cell_population);

            modifier.PropagateFluidSourcesToGrid();

            // Check that there are only 16 non-zero entries per grid
            auto& r_force_grids = modifier.mpArrays->rGetModifiableRightHandSideGrids();

            unsigned non_zero_entries_count = 0;
            for (unsigned x = 0; x < modifier.mpMesh->GetNumGridPtsX(); ++x)
            {
                for (unsigned y = 0; y < modifier.mpMesh->GetNumGridPtsY(); ++y)
                {
                    if (abs(r_force_grids[2][x][y]) > 0.0001)
                    {
                        non_zero_entries_count += 1;
                    }
                }
            }
            TS_ASSERT_EQUALS(non_zero_entries_count, 16u);

            // Check the non-zero entries' values
            TS_ASSERT_DELTA(r_force_grids[2][4][4], 0.0536, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][4][5], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][4][6], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][4][7], 0.0536, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][5][4], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][5][5], 1.8213, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][5][6], 1.8213, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][5][7], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][6][4], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][6][5], 1.8213, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][6][6], 1.8213, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][6][7], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][7][4], 0.0536, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][7][5], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][7][6], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][7][7], 0.0536, 0.0001);

            SimulationTime::Instance()->Destroy();

        }
        { // With zero field sums
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);
            auto& r_balancing_fluid_sources = mesh.rGetBalancingFluidSources();
            r_balancing_fluid_sources.clear();

            auto p_fluid_source = std::make_shared<FluidSource<2>>(0, 0.55, 0.55);
            p_fluid_source->SetStrength(0.1);
            elems.back()->SetFluidSource(p_fluid_source);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetZeroFieldSums(true);
            modifier.SetupConstantMemberVariables(cell_population);

            modifier.PropagateFluidSourcesToGrid();

            // Check that there are only 16 non-zero entries per grid
            auto& r_force_grids = modifier.mpArrays->rGetModifiableRightHandSideGrids();

            unsigned non_zero_entries_count = 0;
            for (unsigned x = 0; x < modifier.mpMesh->GetNumGridPtsX(); ++x)
            {
                for (unsigned y = 0; y < modifier.mpMesh->GetNumGridPtsY(); ++y)
                {
                    if (abs(r_force_grids[2][x][y]) > 0.0001)
                    {
                        non_zero_entries_count += 1;
                    }
                }
            }
            TS_ASSERT_EQUALS(non_zero_entries_count, 16u);

            // Check the non-zero entries' values
            TS_ASSERT_DELTA(r_force_grids[2][4][4], 0.0536, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][4][5], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][4][6], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][4][7], 0.0536, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][5][4], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][5][5], 1.8213, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][5][6], 1.8213, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][5][7], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][6][4], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][6][5], 1.8213, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][6][6], 1.8213, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][6][7], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][7][4], 0.0536, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][7][5], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][7][6], 0.3124, 0.0001);
            TS_ASSERT_DELTA(r_force_grids[2][7][7], 0.0536, 0.0001);

            SimulationTime::Instance()->Destroy();

        }
    }

    void TestWarnings()
    {
        // For coverage
        { // Grid point / noise skip
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetAdditiveNormalNoise(true);
            modifier.SetNoiseSkip(3);
            modifier.SetupConstantMemberVariables(cell_population);

            SimulationTime::Instance()->Destroy();

        }
        { // Too many grid points
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 150, 150);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetAdditiveNormalNoise(true);
            modifier.SetupConstantMemberVariables(cell_population);

            SimulationTime::Instance()->Destroy();

        }
    }

    void TestSolveNavierStokesSpectral()
    {
        { // Coverage for zero field sums
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.1, 0.1));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 4, 4);

            const unsigned int num_grid_pts = 4;
            multi_array<double, 3>& r_vel_grids = mesh.rGetModifiable2dVelocityGrids();

            for (unsigned int dim = 0; dim < 2; dim++) {
                for (unsigned x = 0; x < num_grid_pts; ++x) {
                    for (unsigned y = 0; y < num_grid_pts; ++y) {
                        TS_ASSERT_DELTA(r_vel_grids[dim][x][y], 0.0, 0.00001);
                    }
                }
            }

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetAdditiveNormalNoise(true);
            modifier.SetZeroFieldSums(true);
            modifier.SetupConstantMemberVariables(cell_population);

            for (unsigned int dim = 0; dim < 2; dim++) {
                for (unsigned x = 0; x < num_grid_pts; ++x) {
                    for (unsigned y = 0; y < num_grid_pts; ++y) {
                        TS_ASSERT_DELTA(r_vel_grids[dim][x][y], 0.0, 0.00001);
                    }
                }
            }

            modifier.SolveNavierStokesSpectral();


            std::vector<double> expected_vals = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            // Check values of m_
            unsigned int counter = 0;
            for (unsigned int dim = 0; dim < 2; dim++) {
                for (unsigned x = 0; x < num_grid_pts; ++x) {
                    for (unsigned y = 0; y < num_grid_pts; ++y) {
                        TS_ASSERT_DELTA(r_vel_grids[dim][x][y], expected_vals[counter], 0.00001);
                        counter++;
                    }
                }
            }
            SimulationTime::Instance()->Destroy();

        }
    }

    void TestUpdateFluidVelocityGrids()
    {
        // This calls all other untested methods

        { // Coverage for normal noise
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 150, 150);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetAdditiveNormalNoise(true);
            modifier.SetupConstantMemberVariables(cell_population);
            modifier.UpdateAtEndOfTimeStep(cell_population);


            const unsigned int num_grid_pts = 4;

            std::vector<double> expected_vals = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0 };

            // Check values of m_
            multi_array<double, 3>& r_vel_grids = mesh.rGetModifiable2dVelocityGrids();
            unsigned int counter = 0;
            for (unsigned int dim = 0; dim < 2; dim++) {
                for (unsigned x = 0; x < num_grid_pts; ++x) {
                    for (unsigned y = 0; y < num_grid_pts; ++y) {
                        TS_ASSERT_DELTA(r_vel_grids[dim][x][y], expected_vals[counter], 0.00001);
                        counter++;
                    }
                }
            }
            SimulationTime::Instance()->Destroy();

        }
    }

    void TestArchiving()
    {
        // Create a file for archiving
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "ImmersedBoundarySimulationModifier.arch";

        // Separate scope to write the archive
        {
            // Initialise a growth modifier and set a non-standard mature target area
            ImmersedBoundarySimulationModifier<2> modifier;

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Serialize
            output_arch << modifier;
        }

        // Separate scope to read the archive
        {
            ImmersedBoundarySimulationModifier<2> modifier;

            // Restore the modifier
            std::ifstream ifs(archive_filename.c_str());
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> modifier;
        }
    }

    void TestDelta1D()
    {
        {
            // Division by zero
            ImmersedBoundarySimulationModifier<2> modifier;
            TS_ASSERT(std::isnan(modifier.Delta1D(0.0, 0.0)))
        }
        {
            ImmersedBoundarySimulationModifier<2> modifier;
            TS_ASSERT_DELTA(modifier.Delta1D(0.0, 0.5), 0.5, 1e-9)
        }
        {
            ImmersedBoundarySimulationModifier<2> modifier;
            TS_ASSERT_DELTA(modifier.Delta1D(1.0, 0.5), 0.0, 1e-9)
        }
        {
            ImmersedBoundarySimulationModifier<2> modifier;
            TS_ASSERT_DELTA(modifier.Delta1D(0.5, 0.5), 0.25, 1e-9)
        }
    }

    void TestUpwind2d()
    {
        // Input field = constant 0 everywhere
        {
            // Set up SimulationTime - needed by SetupConstantMemberVariables()
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetupConstantMemberVariables(cell_population);

            multi_array<double, 3>& vel_grids   = mesh.rGetModifiable2dVelocityGrids();
            multi_array<double, 3>& rhs_grids   = modifier.mpArrays->rGetModifiableRightHandSideGrids();

            for (unsigned dim = 0; dim < 2; ++dim)
            {
                for (unsigned x = 0; x < 10; ++x)
                {
                    for (unsigned y = 0; y < 10; ++y)
                    {
                        TS_ASSERT_DELTA(vel_grids[dim][x][y], 0.0, 0.00001);
                        TS_ASSERT_DELTA(rhs_grids[dim][x][y], 0.0, 0.00001);
                    }
                }
            }

            modifier.Upwind2d(vel_grids, rhs_grids);

            for (unsigned dim = 0; dim < 2; ++dim)
            {
                for (unsigned x = 0; x < 10; ++x)
                {
                    for (unsigned y = 0; y < 10; ++y)
                    {
                        TS_ASSERT_DELTA(vel_grids[dim][x][y], 0.0, 0.00001);
                        TS_ASSERT_DELTA(rhs_grids[dim][x][y], 0.0, 0.00001);
                    }
                }
            }

            SimulationTime::Instance()->Destroy();
        }

        // Input field = constant 1.0 everywhere
        {
            // Set up SimulationTime - needed by SetupConstantMemberVariables()
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetupConstantMemberVariables(cell_population);

            multi_array<double, 3>& vel_grids   = mesh.rGetModifiable2dVelocityGrids();
            multi_array<double, 3>& rhs_grids   = modifier.mpArrays->rGetModifiableRightHandSideGrids();

            for (unsigned dim = 0; dim < 2; ++dim)
            {
                for (unsigned x = 0; x < 10; ++x)
                {
                    for (unsigned y = 0; y < 10; ++y)
                    {
                        vel_grids[dim][x][y] = 1.0;
                    }
                }
            }

            for (unsigned dim = 0; dim < 2; ++dim)
            {
                for (unsigned x = 0; x < 10; ++x)
                {
                    for (unsigned y = 0; y < 10; ++y)
                    {
                        TS_ASSERT_DELTA(rhs_grids[dim][x][y], 0.0, 0.00001);
                    }
                }
            }

            modifier.Upwind2d(vel_grids, rhs_grids);

            for (unsigned dim = 0; dim < 2; ++dim)
            {
                for (unsigned x = 0; x < 10; ++x)
                {
                    for (unsigned y = 0; y < 10; ++y)
                    {
                        TS_ASSERT_DELTA(rhs_grids[dim][x][y], 0.0, 0.00001);
                    }
                }
            }
            SimulationTime::Instance()->Destroy();
        }

        // Test propagation
        {
            // Set up SimulationTime - needed by SetupConstantMemberVariables()
            SimulationTime::Instance()->SetStartTime(0.0);
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

            // Create a single node, single element mesh
            std::vector<Node<2>*> nodes;
            nodes.push_back(new Node<2>(0, true, 0.55, 0.55));
            nodes.push_back(new Node<2>(1, true, 0.2, 0.2));
            nodes.push_back(new Node<2>(2, true, 0.1, 0.2));

            std::vector<ImmersedBoundaryElement<2, 2>*> elems;
            elems.push_back(new ImmersedBoundaryElement<2, 2>(0, nodes));

            ImmersedBoundaryMesh<2,2> mesh(nodes, elems, {}, 10, 10);

            // Set up a cell population
            std::vector<CellPtr> cells;
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            CellsGenerator<UniformCellCycleModel, 2> cells_generator;
            cells_generator.GenerateBasicRandom(cells, mesh.GetNumElements(), p_diff_type);
            ImmersedBoundaryCellPopulation<2> cell_population(mesh, cells);

            // Set up simulation modifier
            ImmersedBoundarySimulationModifier<2> modifier;
            modifier.SetupConstantMemberVariables(cell_population);

            multi_array<double, 3>& vel_grids   = mesh.rGetModifiable2dVelocityGrids();
            multi_array<double, 3>& rhs_grids   = modifier.mpArrays->rGetModifiableRightHandSideGrids();

            for (unsigned dim = 0; dim < 2; ++dim)
            {
                for (unsigned x = 0; x < 10; ++x)
                {
                    for (unsigned y = 0; y < 10; ++y)
                    {
                        TS_ASSERT_DELTA(vel_grids[dim][x][y], 0.0, 0.00001);
                        TS_ASSERT_DELTA(rhs_grids[dim][x][y], 0.0, 0.00001);
                    }
                }
            }

            vel_grids[0][1][1] = 0.1;
            vel_grids[0][1][5] = -0.1;
            vel_grids[1][5][1] = 0.1;
            vel_grids[1][5][5] = -0.1;

            modifier.Upwind2d(vel_grids, rhs_grids);

            TS_ASSERT_DELTA(rhs_grids[0][1][1], 0.1, 0.001);
            TS_ASSERT_DELTA(rhs_grids[0][1][5], -0.1, 0.001);
            TS_ASSERT_DELTA(rhs_grids[0][5][1], 0.0, 0.001);
            TS_ASSERT_DELTA(rhs_grids[0][5][5], -0.0, 0.001);

            TS_ASSERT_DELTA(rhs_grids[1][1][1], 0.0, 0.001);
            TS_ASSERT_DELTA(rhs_grids[1][1][5], -0.0, 0.001);
            TS_ASSERT_DELTA(rhs_grids[1][5][1], 0.1, 0.001);
            TS_ASSERT_DELTA(rhs_grids[1][5][5], -0.1, 0.001);

            SimulationTime::Instance()->Destroy();
        }
    }

    void TestAddImmersedBoundaryForce()
    {
        // Set up SimulationTime - needed by SetupConstantMemberVariables()
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create an immersed boundary cell population
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        // Call SetupConstantMemberVariables() to 'initialise' the simulation modifier
        ImmersedBoundarySimulationModifier<2> modifier;
        modifier.SetupConstantMemberVariables(cell_population);

        // Add two immersed boundary force objects to the simulation modifier
        MAKE_PTR(ImmersedBoundaryLinearMembraneForce<2>, p_boundary_force);
        modifier.AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetElementSpringConst(1.0 * 1e7);
        p_boundary_force->SetLaminaSpringConst(1.0 * 1e7);

        MAKE_PTR(ImmersedBoundaryLinearInteractionForce<2>, p_cell_cell_force);
        modifier.AddImmersedBoundaryForce(p_cell_cell_force);
        p_cell_cell_force->SetSpringConst(1.0 * 1e6);

        // Test AddImmersedBoundaryForceContributions() returns the correct applied force on some nodes
        // Note: testing of the force calculations themselves occurs in TestImmersedBoundaryForces
        modifier.AddImmersedBoundaryForceContributions();

        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetAppliedForce()[0], -1.28057e-09, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetAppliedForce()[1], 0.0, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetAppliedForce()[0], -2.32831e-10, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetAppliedForce()[1], 0.0, 1e-3);
    }

    void TestAdditiveNormalNoiseGettersAndSetters()
    {
        ImmersedBoundarySimulationModifier<2> modifier;

        // Noise enabled
        TS_ASSERT_EQUALS(modifier.GetAdditiveNormalNoise(), false);
        modifier.SetAdditiveNormalNoise(true);
        TS_ASSERT_EQUALS(modifier.GetAdditiveNormalNoise(), true);

        // Noise strength
        TS_ASSERT_EQUALS(modifier.GetNoiseStrength(), 0.0);
        modifier.SetNoiseStrength(1.0);
        TS_ASSERT_EQUALS(modifier.GetNoiseStrength(), 1.0);

        // Noise skip
        TS_ASSERT_EQUALS(modifier.GetNoiseSkip(), 1u);
        modifier.SetNoiseSkip(4);
        TS_ASSERT_EQUALS(modifier.GetNoiseSkip(), 4u);

        // Noise length scale
        TS_ASSERT_DELTA(modifier.GetNoiseLengthScale(), 0.1, 1e-6);
        modifier.SetNoiseLengthScale(4);
        TS_ASSERT_DELTA(modifier.GetNoiseLengthScale(), 4.0, 1e-6);

        // Noise zero field sums
        TS_ASSERT_EQUALS(modifier.GetZeroFieldSums(), false);
        modifier.SetZeroFieldSums(true);
        TS_ASSERT_EQUALS(modifier.GetZeroFieldSums(), true);
    }

    void TestAddNormalNoise()
    {
        // Set up SimulationTime - needed by SetupConstantMemberVariables()
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create an immersed boundary cell population with 64 grid points
        ImmersedBoundaryPalisadeMeshGenerator gen(2, 100, 0.2, 2.0, 0.15, true, false, false, 32);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        ImmersedBoundarySimulationModifier<2> modifier;
        modifier.SetAdditiveNormalNoise(true);
        modifier.SetupConstantMemberVariables(cell_population);
        modifier.SetNoiseStrength(0.1);

        auto& r_force_grids = modifier.mpArrays->rGetModifiableForceGrids();

        // Helper function to sum absolute values in force grids
        auto sum_force_grid_abs_values = [](multi_array<double, 3>& rForceGrid)
        {
            double total = 0.0;
            for (unsigned dim = 0; dim < rForceGrid.shape()[0]; ++dim)
            {
                for (unsigned x = 0; x < rForceGrid.shape()[1]; ++x)
                {
                    for (unsigned y = 0; y < rForceGrid.shape()[2]; ++y)
                    {
                        total += abs(rForceGrid[dim][x][y]);
                    }
                }
            }
            return total;
        };

        // Test total forces pre and post noise
        double totalForcePreNoise = sum_force_grid_abs_values(r_force_grids);
        modifier.AddNormalNoise();
        double totalForcePostNoise = sum_force_grid_abs_values(r_force_grids);

        TS_ASSERT_DIFFERS(totalForcePreNoise, totalForcePostNoise);
    }

    void TestZeroFieldSums()
    {
        // Set up SimulationTime - needed by SetupConstantMemberVariables()
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create an immersed boundary cell population with 64 grid points
        ImmersedBoundaryPalisadeMeshGenerator gen(2, 100, 0.2, 2.0, 0.15, true, false, false, 32);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        ImmersedBoundarySimulationModifier<2> modifier;
        modifier.SetupConstantMemberVariables(cell_population);

        auto& r_force_grids = modifier.mpArrays->rGetModifiableForceGrids();

        // Helper function to initialise values in force grids
        auto initialise_force_grid = [](multi_array<double, 3>& rForceGrid)
        {
            for (unsigned dim = 0; dim < rForceGrid.shape()[0]; ++dim)
            {
                for (unsigned x = 0; x < rForceGrid.shape()[1]; ++x)
                {
                    for (unsigned y = 0; y < rForceGrid.shape()[2]; ++y)
                    {
                        rForceGrid[dim][x][y] = 1.0;
                    }
                }
            }
        };

        auto checkWholeGridEquals = [](multi_array<double, 3>& grid, double value)
        {
            bool all_equal = true;
            for (unsigned dim = 0; dim < grid.shape()[0]; ++dim)
            {
                for (unsigned x = 0; x < grid.shape()[1]; ++x)
                {
                    for (unsigned y = 0; y < grid.shape()[2]; ++y)
                    {
                        if (grid[dim][x][y] != value)
                        {
                            all_equal = false;
                            break;
                        }
                    }
                }
            }
            return all_equal;
        };

        // Test total forces pre and post noise
        initialise_force_grid(r_force_grids);
        TS_ASSERT(checkWholeGridEquals(r_force_grids, 1.0));
        modifier.ZeroFieldSums(r_force_grids);
        TS_ASSERT(checkWholeGridEquals(r_force_grids, 0.0));
    }
};

#endif /*TESTIMMERSEDBOUNDARYSIMULATIONMODIFIER_HPP_*/
