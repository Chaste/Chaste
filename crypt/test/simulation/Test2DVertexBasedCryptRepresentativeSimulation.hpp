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

#ifndef TEST2DVERTEXCRYPTREPRESENTATIVESIMULATION_HPP_
#define TEST2DVERTEXCRYPTREPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based or crypt headers
#include "CellBasedSimulationArchiver.hpp"

#include "CryptSimulation2d.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "VertexCryptBoundaryForce.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "SloughingCellKiller.hpp"
#include "CellBasedEventHandler.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FakePetscSetup.hpp"

/**
 * This class consists of a single test, in which a 2D vertex model
 * of a colonic crypt is simulated for a fixed period of time.
 *
 * This test is used for profiling, to establish the run time
 * variation as the code is developed.
 */
class Test2DVertexCryptRepresentativeSimulation : public CxxTest::TestSuite
{
public:

    void Test2DVertexCryptRepresentativeSimulationForProfiling()
    {
        // Set start time
        SimulationTime::Instance()->SetStartTime(0.0);

        // Create mesh
        unsigned crypt_width = 18;
        unsigned crypt_height = 25;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Make crypt shorter for sloughing
        double crypt_length = 20.0;

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            SimpleWntCellCycleModel* p_model = new SimpleWntCellCycleModel;
            p_model->SetDimension(2);

            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                             ( p_model->GetTransitCellG1Duration()
                                                + p_model->GetSG2MDuration() );

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        // Create cell population
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Set up Wnt gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(crypt);
        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Create crypt simulation from cell population and force law
        CryptSimulation2d simulator(crypt);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(30.0);
        simulator.SetOutputDirectory("Test2DVertexCryptRepresentativeSimulationForProfiling");

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // ...and with that the target area modifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add a cell killer
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Tidy up
        WntConcentration<2>::Destroy();
    }
};

#endif /*TEST2DMONOLAYERREPRESENTATIVESIMULATION_HPP_*/
