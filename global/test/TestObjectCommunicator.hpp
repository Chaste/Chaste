/*

Copyright (c) 2005-2012, University of Oxford.
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

#ifndef _TESTOBJECTCOMMUNICATOR_HPP_
#define _TESTOBJECTCOMMUNICATOR_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <string>


#include "Cell.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AbstractCellCycleModel.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "CellData.hpp"
#include "WildTypeCellMutationState.hpp"

#include "OutputFileHandler.hpp"

#include "SmartPointers.hpp"
#include "CommandLineArguments.hpp"
#include "Node.hpp"

#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "ObjectCommunicator.hpp"

class TestObjectCommunicator: public CxxTest::TestSuite
{

public:

    void TestSendingNode() throw(Exception)
    {
        ObjectCommunicator communicator;
        MPI_Status status;

        if (PetscTools::GetMyRank() == 0)
        {
              Node<2>* p_node = new Node<2>(1, false, 0.0, 1.0);

              // Send a copy to every other process
              for(unsigned i=1; i<PetscTools::GetNumProcs(); i++)
              {
                  communicator.SendObject<Node<2> >(p_node, i, 123);
              }

              delete p_node;
        }
        else
        {
            // Receive node and check the correct index has been sent.
            Node<2>* p_recv_node = communicator.RecvObject<Node<2> >(0, 123, status);
            TS_ASSERT_EQUALS(p_recv_node->GetIndex(), 1u);

            // Check the correct location has been set.
            c_vector<double, 2> location;
            location[0] = 0.0;
            location[1] = 1.0;
            TS_ASSERT_EQUALS(p_recv_node->rGetLocation()[0], location[0]);
            TS_ASSERT_EQUALS(p_recv_node->rGetLocation()[1], location[1]);

            // Check the correct boundary property has been recevied.
            TS_ASSERT(!p_recv_node->IsBoundaryNode());

        }
    }

    void TestSendingCell() throw(Exception)
    {
        ObjectCommunicator communicator;
        MPI_Status status;

        if (PetscTools::GetMyRank() == 0)
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // Create mutation state
            boost::shared_ptr<AbstractCellProperty> p_healthy_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

            // Create cell-cycle model
            StochasticDurationGenerationBasedCellCycleModel* p_cell_model = new StochasticDurationGenerationBasedCellCycleModel();
            p_cell_model->SetCellProliferativeType(STEM);

            // Create cell property collection
            CellPropertyCollection collection;
            MAKE_PTR(CellLabel, p_label);
            collection.AddProperty(p_label);

            // Create cell
            CellPtr p_cell(new Cell(p_healthy_state, p_cell_model, false, collection));
            p_cell->InitialiseCellCycleModel();
            p_simulation_time->IncrementTimeOneStep();

            // Send to the other process
            for(unsigned i=1; i<PetscTools::GetNumProcs(); i++)
            {
                communicator.SendObject<Cell>(&(*p_cell), i, 123);
            }

        }
        else
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            // Initialize a cell
            Cell* p_cell_2 = communicator.RecvObject<Cell>(0, 123, status);

            TS_ASSERT_EQUALS(p_cell_2->GetCellCycleModel()->GetIdentifier(), "StochasticDurationGenerationBasedCellCycleModel");
            TS_ASSERT(p_cell_2->HasCellProperty<CellLabel>());
        }

        SimulationTime::Destroy();
    }

};

#endif //_TESTOBJECTCOMMUNICATOR_HPP_
