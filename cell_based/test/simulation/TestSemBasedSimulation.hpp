/*

Copyright (c) 2005-2024, University of Oxford.
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

#ifndef TESTSEMBASEDSIMULATION_HPP_
#define TESTSEMBASEDSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <cmath>

#include "CheckpointArchiveTypes.hpp"
#include "OffLatticeSimulation.hpp"
#include "SemMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "ToroidalHoneycombMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "ChemotacticForce.hpp"
#include "RandomCellKiller.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulationWithMyStoppingEvent.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "CellIdWriter.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "SemBasedCellPopulation.hpp"
#include "SemInterCellularForce.hpp"
#include "SemIntraCellularForce.hpp"
#include "NoCellCycleModel.hpp"
#include "NodeLocationWriter.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestSemBasedSimulation : public AbstractCellBasedWithTimingsTestSuite
{
public:

    void TestSemBasedSimulationExample()
    {
        // Create a simple 2D SemBasedCellPopulation
        SemMeshGenerator generator;
        generator.GenerateSingleCell({0.0, 0.0}, {0.5, 0.5}, {8, 8});
        auto p_mesh = generator.GetMesh();

        // Assertions
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 64);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.0);

        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 64);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1);


        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1);
        std::cout << cell_population.GetNumRealCells() << std::endl;
        simulator.SetOutputDirectory("TestSemBasedSimulation");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(2.5);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<2>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

        std::cout << cell_population.GetNumRealCells() << std::endl;
        // Create some force laws and pass them to the simulation
        MAKE_PTR(SemInterCellularForce<2>, p_inter_cellular_force);
        simulator.AddForce(p_inter_cellular_force);

        std::cout << cell_population.GetNumRealCells() << std::endl;
        // Create some force laws and pass them to the simulation
        MAKE_PTR(SemIntraCellularForce<2>, p_intra_cellular_force);
        simulator.AddForce(p_intra_cellular_force);
        
        std::cout << cell_population.GetNumRealCells() << std::endl;
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1);

        // Run the simulation
        simulator.Solve();
        
        std::cout << cell_population.GetNumRealCells() << std::endl;

    }

};

#endif /*TESTSEMBASEDSIMULATION_HPP_*/
