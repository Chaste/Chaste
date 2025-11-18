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

#ifndef TESTSEMBASEDSIMULATION_HPP_
#define TESTSEMBASEDSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <cstdio>
#include <cmath>

#include "CheckpointArchiveTypes.hpp"
#include "OffLatticeSimulation.hpp"
#include "SemMeshGenerator.hpp"
#include "SemSingleElementMeshGenerator.hpp"
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
#include "SemRegionalForce.hpp"
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

    void xTestSemBasedSimulationExample2D()
    {
        SemSingleElementMeshGenerator<2> generator({5,8}, 0.5);
        auto p_mesh = generator.GetMesh();

        c_vector<double, 4> boxCollectionDomain{};
        boxCollectionDomain[0] = -1.0;
        boxCollectionDomain[1] =  2.0;
        boxCollectionDomain[2] = -1.0;
        boxCollectionDomain[3] =  2.0;

        p_mesh->SetUpBoxCollection(0.1, boxCollectionDomain);

        // Assertions
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 40);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.0);

        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 40);
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
        MAKE_PTR(SemRegionalForce<2>, p_sem_force);
        simulator.AddForce(p_sem_force);

        // Run the simulation
        simulator.Solve();
        
        std::cout << cell_population.GetNumRealCells() << std::endl;

    }

    void TestSemBasedSimulationExample3D()
    {
        SemSingleElementMeshGenerator<3> generator({ 5, 8, 5 }, 0.5);
        auto p_mesh = generator.GetMesh();

        // Set the most central nodes to region 1
        {
            c_vector<double, 3> centroid = p_mesh->GetCentroidOfElement(0u);
            auto p_elem_0 = p_mesh->GetElement(0);
            for (unsigned i = 0; i < p_elem_0->GetNumNodes(); ++i)
            {
                auto p_node = p_elem_0->GetNode(i);
                if (norm_2(p_mesh->GetVectorFromAtoB(centroid, p_node->rGetLocation())) < 0.2)
                {
                    p_node->SetRegion(1u);
                }
            }
        }

        c_vector<double, 6> boxCollectionDomain{};
        boxCollectionDomain[0] = -1.0;
        boxCollectionDomain[1] =  2.0;
        boxCollectionDomain[2] = -1.0;
        boxCollectionDomain[3] =  2.0;
        boxCollectionDomain[4] = -1.0;
        boxCollectionDomain[5] =  2.0;

        p_mesh->SetUpBoxCollection(0.1, boxCollectionDomain);

        // Assertions
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), 1);
        TS_ASSERT_EQUALS(p_mesh->GetNumNodes(), 200);

        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements());
        SemBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetDampingConstantNormal(1.0);

        TS_ASSERT_EQUALS(cell_population.GetNumElements(), 1);
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), 200);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1);


        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cell_population);
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 1);
        std::cout << cell_population.GetNumRealCells() << std::endl;
        simulator.SetOutputDirectory("TestSemBasedSimulation3D");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(2.5);
        simulator.SetNumericalMethod(boost::make_shared<ForwardEulerNumericalMethod<3>>());
        simulator.GetNumericalMethod()->SetUseUpdateNodeLocation(false);

        std::cout << cell_population.GetNumRealCells() << std::endl;
        // Create some force laws and pass them to the simulation
        MAKE_PTR(SemRegionalForce<3>, p_sem_force);
        simulator.AddForce(p_sem_force);

        // Run the simulation
        simulator.Solve();

        std::cout << cell_population.GetNumRealCells() << std::endl;
    }

};

#endif /*TESTSEMBASEDSIMULATION_HPP_*/
