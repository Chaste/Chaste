/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#ifndef TEST2DVERTEXCRYPTREPRESENTATIVESIMULATION_HPP_
#define TEST2DVERTEXCRYPTREPRESENTATIVESIMULATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "TissueSimulationArchiver.hpp"

#include "VertexCryptSimulation2d.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "VertexCryptBoundaryForce.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "SloughingCellKiller.hpp"
#include "CellBasedEventHandler.hpp"
#include "WildTypeCellMutationState.hpp"

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

    void Test2DVertexCryptRepresentativeSimulationForProfiling() throw (Exception)
    {
        // Set start time
        SimulationTime::Instance()->SetStartTime(0.0);

        // Create mesh
        unsigned crypt_width = 18;
        unsigned crypt_height = 25;
        CylindricalHoneycombVertexMeshGenerator generator(crypt_width, crypt_height, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create cells
        std::vector<TissueCell> cells;
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        for (unsigned elem_index=0; elem_index<p_mesh->GetNumElements(); elem_index++)
        {
            double birth_time = - RandomNumberGenerator::Instance()->ranf()*
                                 ( TissueConfig::Instance()->GetTransitCellG1Duration()
                                    + TissueConfig::Instance()->GetSG2MDuration() );

            SimpleWntCellCycleModel* p_model = new SimpleWntCellCycleModel;
            p_model->SetDimension(2);
            p_model->SetCellProliferativeType(TRANSIT);

            TissueCell cell(p_state, p_model);
            cell.SetBirthTime(birth_time);
            cells.push_back(cell);
        }

        // Create tissue
        VertexBasedTissue<2> crypt(*p_mesh, cells);

        // Set up Wnt gradient
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetTissue(crypt);

        // Create force law
        NagaiHondaForce<2> force_law;
        std::vector<AbstractForce<2>*> force_collection;
        force_collection.push_back(&force_law);

        // Create crypt simulation from tissue and force law
        VertexCryptSimulation2d simulator(crypt, force_collection);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(30.0);
        simulator.SetOutputDirectory("Test2DVertexCryptRepresentativeSimulationForProfiling");

        // Make crypt shorter for sloughing
        TissueConfig::Instance()->SetCryptLength(20.0);
        SloughingCellKiller<2> sloughing_cell_killer(&crypt);
        simulator.AddCellKiller(&sloughing_cell_killer);

        // Run simulation
        simulator.Solve();

        // Tidy up
        WntConcentration<2>::Destroy();
    }

};

#endif /*TEST2DMONOLAYERREPRESENTATIVESIMULATION_HPP_*/
