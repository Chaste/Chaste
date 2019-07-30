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


#ifndef TESTCELLSRN_HPP_
#define TESTCELLSRN_HPP_

#include <cxxtest/TestSuite.h>

#include <fstream>
#include <iostream>

#include "Cell.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "Goldbeter1991SrnModel.hpp"
#include "SmartPointers.hpp"
#include "TransitCellProliferativeType.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "FakePetscSetup.hpp"

class TestCellSrn: public AbstractCellBasedTestSuite
{
public:

    void TestGoldbeter1991OdeSteadyState()
    {
        // Keep running until we reach steady state
        SimulationTime* p_simulation_time = SimulationTime::Instance();

        // run until 100, with dt=0.01
        double t1=100;
        double dt=0.01;
        unsigned num_steps=(unsigned) t1/dt;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(t1, num_steps+1);

        boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_transit_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        // Create a cell-cycle model
        FixedG1GenerationalCellCycleModel* p_cell_model = new FixedG1GenerationalCellCycleModel();
        Goldbeter1991SrnModel* p_srn_model = new Goldbeter1991SrnModel();
        CellPtr p_tn_cell(new Cell(p_healthy_state, p_cell_model, p_srn_model, false, CellPropertyCollection()));
        p_tn_cell->SetCellProliferativeType(p_transit_type);
        p_tn_cell->InitialiseCellCycleModel();
        p_tn_cell->InitialiseSrnModel();

        // Run the cell simulation until t1
        while (!p_simulation_time->IsFinished())
        {
            p_simulation_time->IncrementTimeOneStep();
            if (p_tn_cell->ReadyToDivide())
            {
                p_tn_cell->Divide();
            }
        }

        // Get final state variables
        double current_time = SimulationTime::Instance()->GetTime();
        std::cout << "Finished ODE - " << "mSimulatedToTime : " << p_srn_model->GetSimulatedToTime() << ", current_time : " << current_time << std::endl;

        // Direct access to state variables
        double C = dynamic_cast<Goldbeter1991SrnModel*>(p_tn_cell->GetSrnModel())->GetC();
        TS_ASSERT_DELTA(C, 0.5470, 1e-4);
        double M = dynamic_cast<Goldbeter1991SrnModel*>(p_tn_cell->GetSrnModel())->GetM();
        TS_ASSERT_DELTA(M, 0.2936, 1e-4);
        double X = dynamic_cast<Goldbeter1991SrnModel*>(p_tn_cell->GetSrnModel())->GetX();
        TS_ASSERT_DELTA(X, 0.0067, 1e-4);

        // Indirect access to state vector
        C = dynamic_cast<Goldbeter1991SrnModel*>(p_tn_cell->GetSrnModel())->GetProteinConcentrations()[0];
        TS_ASSERT_DELTA(C, 0.5470, 1e-4);
        M = dynamic_cast<Goldbeter1991SrnModel*>(p_tn_cell->GetSrnModel())->GetProteinConcentrations()[1];
        TS_ASSERT_DELTA(M, 0.2936, 1e-4);
        X = dynamic_cast<Goldbeter1991SrnModel*>(p_tn_cell->GetSrnModel())->GetProteinConcentrations()[2];
        TS_ASSERT_DELTA(X, 0.0067, 1e-4);

        std::cout << "Finished ODE - " << "C : " << C << ", M : " << M  << ", X : " << X << std::endl;
    }
};

#endif /* TESTCELLSRN_HPP_ */
