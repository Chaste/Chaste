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
#include "OffLatticeSimulation.hpp"
#include "OutputFileHandler.hpp"
#include "SmartPointers.hpp"
#include "UniformlyDistributedCellCycleModel.hpp"

// Includes from Immersed Boundary
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryHoneycombMeshGenerator.hpp"
#include "ImmersedBoundaryMembraneElasticityForce.hpp"
#include "ImmersedBoundaryCellCellInteractionForce.hpp"

#include "ForwardEulerNumericalMethod.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/make_shared.hpp>

#include "Debug.hpp"
#include "Timer.hpp"

// Simulation does not run in parallel
#include "FakePetscSetup.hpp"

class TestProfiling : public AbstractCellBasedTestSuite
{
public:

    /**
     * A helper defining a profiling test with settable variables.  Each simulation defines a hexagonal packing which
     * relaxes for a fixed number of time steps.
     *
     * @param numGridPts  the number of fluid grid points in each dimension
     * @param numNodesPerSide  the number of nodes on each side of each hexagon
     */
    void Profile(unsigned numGridPts, unsigned numNodesPerSide)
    {
        /**
         * @param numElementsX  the number of cells from left to right along the domain
         * @param numElementsY  the number of cells from top to bottom up the domain
         * @param numNodesPerCell  the number of nodes per cell (defaults to 100)
         * @param proportionalGap  the proportion of space between elements
         * @param padding  the minimum padding around the edge of the generated mesh
         */
        ImmersedBoundaryHoneycombMeshGenerator gen(5, 4, numNodesPerSide, 0.025, 0.1);
        ImmersedBoundaryMesh<2, 2>* p_mesh = gen.GetMesh();

        p_mesh->SetNumGridPtsXAndY(numGridPts);

        PRINT_VARIABLE(p_mesh->GetSpacingRatio());

        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<UniformlyDistributedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);

        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetIfPopulationHasActiveSources(false);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2,2> >());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(true);

        // Add main immersed boundary simulation modifier
        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_main_modifier);
        simulator.AddSimulationModifier(p_main_modifier);

        // Add force laws
        MAKE_PTR(ImmersedBoundaryMembraneElasticityForce<2>, p_boundary_force);
        p_main_modifier->AddImmersedBoundaryForce(p_boundary_force);
        p_boundary_force->SetSpringConstant(1e7);

        MAKE_PTR(ImmersedBoundaryCellCellInteractionForce<2>, p_cell_cell_force);
        p_main_modifier->AddImmersedBoundaryForce(p_cell_cell_force);
        p_cell_cell_force->SetSpringConstant(1e6);

        std::string output_directory = "numerics_paper/profiling_" + boost::lexical_cast<std::string>(numGridPts);
        simulator.SetOutputDirectory(output_directory);

        // Set simulation properties
        double dt = 0.01;
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(2000 * dt);
        simulator.Solve();

        // Because this will be run multiple times, we destroy the SimulationTime singleton
        SimulationTime::Destroy();
    }

    void TestProfile512() throw(Exception)
    {
        Timer timer;
        timer.Reset();

        Profile(512, 50);

        PRINT_VARIABLE(timer.GetElapsedTime());
    }

    void TestProfile1024() throw(Exception)
    {
        Timer timer;
        timer.Reset();

        Profile(1024, 100);

        PRINT_VARIABLE(timer.GetElapsedTime());
    }

    void TestProfile2048() throw(Exception)
    {
        Timer timer;
        timer.Reset();

        Profile(2048, 200);

        PRINT_VARIABLE(timer.GetElapsedTime());
    }
};
