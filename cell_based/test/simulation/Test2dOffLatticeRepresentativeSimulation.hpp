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

#ifndef TEST2DOFFLATTICEREPRESENTATIVESIMULATION_HPP_
#define TEST2DOFFLATTICEREPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include "CheckpointArchiveTypes.hpp"

#include "OffLatticeSimulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "RandomCellKiller.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"

/**
 * This class consists of a single test, in which a 2D model
 * of a growing monolayer of cells is simulated for a fixed
 * period of time.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class Test2DMonolayerRepresentativeSimulation : public CxxTest::TestSuite
{
public:

    void Test2DMonolayerRepresentativeSimulationForProfiling()
    {
        // Set start time
        SimulationTime::Instance()->SetStartTime(0.0);

        // Create a simple mesh
        HoneycombMeshGenerator generator(5, 5, 0);
        MutableMesh<2,2>* p_mesh = generator.GetCircularMesh(3.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel();
            p_model->SetMaxTransitGenerations(UINT_MAX);
            p_model->SetTransitCellG1Duration(1.0);
            p_model->SetStemCellG1Duration(1.0);

            double birth_time = -RandomNumberGenerator::Instance()->ranf()*
                                (  p_model->GetStemCellG1Duration()
                                 + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create a cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Test2DMonolayerRepresentativeSimulationForProfiling");
        simulator.SetEndTime(50.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Run simulation
        simulator.Solve();

        // Tidy up
        SimulationTime::Destroy();
        RandomNumberGenerator::Destroy();
    }
};

#endif /*TEST2DOFFLATTICEREPRESENTATIVESIMULATION_HPP_*/
