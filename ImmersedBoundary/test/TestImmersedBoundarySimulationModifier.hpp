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

// Needed for the test environment
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

// Includes from trunk
#include "CellsGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FileComparison.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"

// Includes from Immersed Boundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"
#include "ImmersedBoundaryCellCellInteractionForce.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

///\todo Improve testing
class TestImmersedBoundarySimulationModifier : public AbstractCellBasedTestSuite
{
public:

    void TestConstructorAndGetAndSetMethods() throw(Exception)
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

    void TestOutputParametersWithImmersedBoundarySimulationModifier() throw(Exception)
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
        FileFinder reference("projects/ImmersedBoundary/test/data/TestOutputParametersWithImmersedBoundarySimulationModifier/ImmersedBoundarySimulationModifier.parameters",
                RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated, reference);
        TS_ASSERT(comparer.CompareFiles());
    }

    void TestSetupConstantMemberVariables() throw(Exception)
    {
        // Set up SimulationTime - needed by SetupConstantMemberVariables()
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        ImmersedBoundarySimulationModifier<2> modifier;

        TS_ASSERT_EQUALS(modifier.mNumGridPtsX, 0u);
        TS_ASSERT_EQUALS(modifier.mNumGridPtsY, 0u);
        TS_ASSERT_DELTA(modifier.mGridSpacingX, 0.0, 1e-6);
        TS_ASSERT_DELTA(modifier.mGridSpacingY, 0.0, 1e-6);
        TS_ASSERT_DELTA(modifier.mFftNorm, 0.0, 1e-6);
        TS_ASSERT(modifier.mpArrays == NULL);
        TS_ASSERT(modifier.mpFftInterface == NULL);

        // Test that the correct exception is throw if trying to use the wrong cell population class
        HoneycombVertexMeshGenerator vertex_generator(2, 2);
        MutableVertexMesh<2,2>* p_vertex_mesh = vertex_generator.GetMesh();
        std::vector<CellPtr> vertex_cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> vertex_cells_generator;
        vertex_cells_generator.GenerateBasicRandom(vertex_cells, p_vertex_mesh->GetNumElements(), p_diff_type);
        VertexBasedCellPopulation<2> vertex_cell_population(*p_vertex_mesh, vertex_cells);

        TS_ASSERT_THROWS_THIS(modifier.SetupConstantMemberVariables(vertex_cell_population),
                "Cell population must be immersed boundary");

        // Now test SetupConstantMemberVariables() using an immersed boundary cell population
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();
        std::vector<CellPtr> cells;
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        TS_ASSERT_THROWS_NOTHING(modifier.SetupConstantMemberVariables(cell_population));
        TS_ASSERT_EQUALS(modifier.mNumGridPtsX, 256u);
        TS_ASSERT_EQUALS(modifier.mNumGridPtsY, 256u);
        TS_ASSERT_DELTA(modifier.mGridSpacingX, 0.0039, 1e-4);
        TS_ASSERT_DELTA(modifier.mGridSpacingY, 0.0039, 1e-4);
        TS_ASSERT_DELTA(modifier.mFftNorm, 65536.0, 1e-6);
        TS_ASSERT(modifier.mpArrays != NULL);
        TS_ASSERT(modifier.mpFftInterface != NULL);
    }

    void TestClearForcesAndSources() throw(Exception)
    {
        // Set up SimulationTime - needed by SetupConstantMemberVariables()
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create an immersed boundary cell population where each node holds a non-zero applied force
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            c_vector<double, 2> force;
            force[0] = i*0.01;
            force[1] = 2*i*0.01;

            p_mesh->GetNode(i)->ClearAppliedForce();
            p_mesh->GetNode(i)->AddAppliedForceContribution(force);
        }

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        // For coverage of ClearForcesAndSources(), specify that the cell population has active sources
        cell_population.SetIfPopulationHasActiveSources(true);

        // Call SetupConstantMemberVariables() first, to pass the mesh to the simulation modifier
        ImmersedBoundarySimulationModifier<2> modifier;
        modifier.SetupConstantMemberVariables(cell_population);

        // Test ClearForcesAndSources() correctly resets the applied force on each node
        modifier.ClearForcesAndSources();

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetAppliedForce()[0], 0.0, 1e-6);
            TS_ASSERT_DELTA(p_mesh->GetNode(i)->rGetAppliedForce()[1], 0.0, 1e-6);
        }

        // Test ClearForcesAndSources() correctly resets mpArrays
        multi_array<double, 3>& r_force_grids = modifier.mpArrays->rGetModifiableForceGrids();
        multi_array<double, 3>& r_rhs_grid = modifier.mpArrays->rGetModifiableRightHandSideGrids();

        for (unsigned x=0; x<modifier.mNumGridPtsX; x++)
        {
            for (unsigned y=0; y<modifier.mNumGridPtsY; y++)
            {
                TS_ASSERT_DELTA(r_rhs_grid[2][x][y], 0.0, 1e-6);
                for (unsigned dim=0; dim<2; dim++)
                {
                    TS_ASSERT_DELTA(r_force_grids[dim][x][y], 0.0, 1e-6);
                }
            }
        }
    }

    void TestPropagateForcesToFluidGrid() throw(Exception)
    {
        ///\todo Test this method
    }

    void TestPropagateFluidSourcesToGrid() throw(Exception)
    {
        ///\todo Test this method
    }

    void TestSolveNavierStokesSpectral() throw(Exception)
    {
        ///\todo Test this method
    }

    void TestUpdateFluidVelocityGrids() throw(Exception)
    {
        ///\todo Test this method
    }

    void TestDelta1D() throw(Exception)
    {
        ///\todo Test this method
    }

    void TestUpwind2d() throw(Exception)
    {
        ///\todo Test this method
    }

    void TestSetMemberVariablesForTesting() throw(Exception)
    {
        ImmersedBoundarySimulationModifier<2> modifier;

        modifier.SetMemberVariablesForTesting(15, 15);

        TS_ASSERT_EQUALS(modifier.mNumGridPtsX, 15u);
        TS_ASSERT_EQUALS(modifier.mNumGridPtsY, 15u);
        TS_ASSERT_DELTA(modifier.mGridSpacingX, 0.0666, 1e-3);
        TS_ASSERT_DELTA(modifier.mGridSpacingY, 0.0666, 1e-3);
        TS_ASSERT_DELTA(modifier.mFftNorm, 15.0, 1e-4);
    }

    void TestUpdateAtEndOfTimeStep() throw(Exception)
    {
        ///\todo Test this method
    }

    void TestSetupSolve() throw(Exception)
    {
        ///\todo Test this method
    }

    void TestAddImmersedBoundaryForce() throw(Exception)
    {
        // Set up SimulationTime - needed by SetupConstantMemberVariables()
        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(2.0, 2);

        // Create an immersed boundary cell population
        ImmersedBoundaryPalisadeMeshGenerator gen(5, 100, 0.2, 2.0, 0.15, true);
        ImmersedBoundaryMesh<2,2>* p_mesh = gen.GetMesh();
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);

        // Call SetupConstantMemberVariables() to 'initialise' the simulation modifier
        ImmersedBoundarySimulationModifier<2> modifier;
        modifier.SetupConstantMemberVariables(cell_population);

        // Add two immersed boundary force objects to the simulation modifier
        MAKE_PTR(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force);
        modifier.AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetSpringConstant(1.0 * 1e7);

        MAKE_PTR(ImmersedBoundaryCellCellInteractionForce<2>, p_cell_cell_force);
        modifier.AddImmersedBoundaryForce(p_cell_cell_force);
        p_cell_cell_force->SetSpringConstant(1.0 * 1e6);

        // Test AddImmersedBoundaryForceContributions() returns the correct applied force on some nodes
        // Note: testing of the force calculations themselves occurs in TestImmersedBoundaryForces
        modifier.AddImmersedBoundaryForceContributions();

        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetAppliedForce()[0], -125.2290, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(0)->rGetAppliedForce()[1], 6352.7140, 1e-3);

        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetAppliedForce()[0], -1235.1356, 1e-3);
        TS_ASSERT_DELTA(p_mesh->GetNode(5)->rGetAppliedForce()[1], 16800.5590, 1e-3);
    }
};
