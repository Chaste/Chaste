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

#ifndef TESTODEBASEDSRNMODELS_HPP_
#define TESTODEBASEDSRNMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "AbstractSrnModel.hpp"
#include "NullSrnModel.hpp"
#include "DeltaNotchSrnModel.hpp"
#include "Goldbeter1991SrnModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"
#include "UniformCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestOdeBasedSrnModels : public AbstractCellBasedTestSuite
{
public:

    void TestDeltaNotchSrnCorrectBehaviour()
    {
        TS_ASSERT_THROWS_NOTHING(DeltaNotchSrnModel srn_model);

        DeltaNotchSrnModel* p_srn_model = new DeltaNotchSrnModel();

        // Create a vector of initial conditions
        std::vector<double> starter_conditions;
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.5);
        p_srn_model->SetInitialConditions(starter_conditions);

        UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model, false, CellPropertyCollection()));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->GetCellData()->SetItem("mean delta", 1.0);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();

        // Now updated to initial conditions
        TS_ASSERT_DELTA(p_srn_model->GetNotch(), 0.5, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetDelta(), 0.5, 1e-4);

        // Now update the SRN
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        double end_time = 10.0;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        while (p_simulation_time->GetTime() < end_time)
        {
            p_simulation_time->IncrementTimeOneStep();
            p_srn_model->SimulateToCurrentTime();
        }

        // Test converged to steady state
        TS_ASSERT_DELTA(p_srn_model->GetNotch(), 0.9900, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetDelta(), 0.0101, 1e-4);
    }

    void TestDeltaNotchSrnCreateCopy()
    {
        // Test with DeltaNotchSrnModel
        DeltaNotchSrnModel* p_model= new DeltaNotchSrnModel;

        // Set ODE system
        std::vector<double> state_variables;
        state_variables.push_back(2.0);
        state_variables.push_back(3.0);
        p_model->SetOdeSystem(new DeltaNotchOdeSystem(state_variables));

        p_model->SetInitialConditions(state_variables);

        // Create a copy
        DeltaNotchSrnModel* p_model2 = static_cast<DeltaNotchSrnModel*> (p_model->CreateSrnModel());

        // Check correct initializations
        TS_ASSERT_EQUALS(p_model2->GetNotch(), 2.0);
        TS_ASSERT_EQUALS(p_model2->GetDelta(), 3.0);

        // Destroy models
        delete p_model;
        delete p_model2;
    }

    void TestArchiveDeltaNotchSrnModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "delta_notch_srn.arch";

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractSrnModel* p_srn_model = new DeltaNotchSrnModel;

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);

            // We must create a cell to be able to initialise the cell srn model's ODE system
            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->GetCellData()->SetItem("mean delta", 10.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();
            p_cell->SetBirthTime(0.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Read mean Delta from CellData
            static_cast<DeltaNotchSrnModel*>(p_srn_model)->UpdateDeltaNotch();
            TS_ASSERT_DELTA(static_cast<DeltaNotchSrnModel*>(p_srn_model)->GetMeanNeighbouringDelta(), 10.0, 1e-12);

            output_arch << p_srn_model;

            // Note that here, deletion of the cell-cycle model and srn is handled by the cell destructor
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractSrnModel* p_srn_model;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_srn_model;

            TS_ASSERT_DELTA(static_cast<DeltaNotchSrnModel*>(p_srn_model)->GetMeanNeighbouringDelta(), 10.0, 1e-12);

            delete p_srn_model;
        }
    }

    void TestGoldbeter1991SrnCorrectBehaviour()
    {
        TS_ASSERT_THROWS_NOTHING(Goldbeter1991SrnModel srn_model);

        Goldbeter1991SrnModel* p_srn_model = new Goldbeter1991SrnModel();

        // Create a vector of initial conditions
        std::vector<double> starter_conditions;
        starter_conditions.push_back(0.5);
        starter_conditions.push_back(0.6);
        starter_conditions.push_back(0.7);
        p_srn_model->SetInitialConditions(starter_conditions);

        UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model, false, CellPropertyCollection()));
        p_cell->SetCellProliferativeType(p_diff_type);
        p_cell->GetCellData()->SetItem("mean delta", 1.0);
        p_cell->InitialiseCellCycleModel();
        p_cell->InitialiseSrnModel();

        // Now updated to initial conditions
        TS_ASSERT_DELTA(p_srn_model->GetC(), 0.5, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetM(), 0.6, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetX(), 0.7, 1e-4);

        // Now update the SRN
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        double end_time = 10.0;
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        while (p_simulation_time->GetTime() < end_time)
        {
            p_simulation_time->IncrementTimeOneStep();
            p_srn_model->SimulateToCurrentTime();
        }

        // Test converged to steady state
        TS_ASSERT_DELTA(p_srn_model->GetC(), 5.5642, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetM(), 4.7817, 1e-4);
        TS_ASSERT_DELTA(p_srn_model->GetX(), 2.1160, 1e-4);
    }

    void TestGoldbeter1991SrnCreateCopy()
    {
        // Test with Goldbeter1991SrnModel
        Goldbeter1991SrnModel* p_model= new Goldbeter1991SrnModel;

        // Set ODE system
        std::vector<double> state_variables;
        state_variables.push_back(2.0);
        state_variables.push_back(3.0);
        state_variables.push_back(4.0);

        p_model->SetOdeSystem(new Goldbeter1991OdeSystem(state_variables));

        // Create a copy
        Goldbeter1991SrnModel* p_model2 = static_cast<Goldbeter1991SrnModel*> (p_model->CreateSrnModel());

        // Check correct initializations
        TS_ASSERT_EQUALS(p_model2->GetC(), 2.0);
        TS_ASSERT_EQUALS(p_model2->GetM(), 3.0);
        TS_ASSERT_EQUALS(p_model2->GetX(), 4.0);

        // Destroy models
        delete p_model;
        delete p_model2;
    }

    void TestArchiveGoldbeter1991SrnModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "delta_notch_srn.arch";

        // Create an output archive
        {
            double end_time = 10.0;
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, 100);

            UniformCellCycleModel* p_cc_model = new UniformCellCycleModel();

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractSrnModel* p_srn_model = new Goldbeter1991SrnModel;

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);

            // We must create a cell to be able to initialise the cell srn model's ODE system
            CellPtr p_cell(new Cell(p_healthy_state, p_cc_model, p_srn_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->InitialiseCellCycleModel();
            p_cell->InitialiseSrnModel();
            p_cell->SetBirthTime(0.0);

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Now update the SRN so the state variables are different from ICS
            while (p_simulation_time->GetTime() < end_time)
            {
                p_simulation_time->IncrementTimeOneStep();
                p_srn_model->SimulateToCurrentTime();
            }

            TS_ASSERT_DELTA(static_cast<Goldbeter1991SrnModel*>(p_srn_model)->GetC(), 3.6955, 1e-4);
            TS_ASSERT_DELTA(static_cast<Goldbeter1991SrnModel*>(p_srn_model)->GetM(), 7.4468, 1e-4);
            TS_ASSERT_DELTA(static_cast<Goldbeter1991SrnModel*>(p_srn_model)->GetX(), 12.7484, 1e-4);

            output_arch << p_srn_model;

            // Note that here, deletion of the cell-cycle model and srn is handled by the cell destructor
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractSrnModel* p_srn_model;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_srn_model;

            TS_ASSERT_DELTA(static_cast<Goldbeter1991SrnModel*>(p_srn_model)->GetC(), 3.6955, 1e-4);
            TS_ASSERT_DELTA(static_cast<Goldbeter1991SrnModel*>(p_srn_model)->GetM(), 7.4468, 1e-4);
            TS_ASSERT_DELTA(static_cast<Goldbeter1991SrnModel*>(p_srn_model)->GetX(), 12.7484, 1e-4);

            // Destroy model
            delete p_srn_model;
        }
    }

    void TestSrnModelOutputParameters()
    {
        std::string output_directory = "TestSrnModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);


        // Test with NullSrnModel
        {
            NullSrnModel srn_model;

            TS_ASSERT_EQUALS(srn_model.GetIdentifier(), "NullSrnModel");

            out_stream parameter_file = output_file_handler.OpenOutputFile("null_srn_results.parameters");
            srn_model.OutputSrnModelParameters(parameter_file);
            parameter_file->close();

            {
                // Compare the generated file in test output with a reference copy in the source code.
                FileFinder generated = output_file_handler.FindFile("null_srn_results.parameters");
                FileFinder reference("cell_based/test/data/TestSrnModels/null_srn_results.parameters",
                                     RelativeTo::ChasteSourceRoot);
                FileComparison comparer(generated, reference);
                TS_ASSERT(comparer.CompareFiles());
            }
        }

        // Test with DeltaNotchSrnModel
        {
            DeltaNotchSrnModel srn_model;

            TS_ASSERT_EQUALS(srn_model.GetIdentifier(), "DeltaNotchSrnModel");

            out_stream parameter_file = output_file_handler.OpenOutputFile("delta_notch_srn_results.parameters");
            srn_model.OutputSrnModelParameters(parameter_file);
            parameter_file->close();

            {
                // Compare the generated file in test output with a reference copy in the source code.
                FileFinder generated = output_file_handler.FindFile("delta_notch_srn_results.parameters");
                FileFinder reference("cell_based/test/data/TestSrnModels/delta_notch_srn_results.parameters",
                                     RelativeTo::ChasteSourceRoot);
                FileComparison comparer(generated, reference);
                TS_ASSERT(comparer.CompareFiles());
            }
        }

        // Test with Goldbeter1991SrnModel
        {
            Goldbeter1991SrnModel srn_model;

            TS_ASSERT_EQUALS(srn_model.GetIdentifier(), "Goldbeter1991SrnModel");

            out_stream parameter_file = output_file_handler.OpenOutputFile("gb_1991_srn_results.parameters");
            srn_model.OutputSrnModelParameters(parameter_file);
            parameter_file->close();

            {
                // Compare the generated file in test output with a reference copy in the source code.
                FileFinder generated = output_file_handler.FindFile("gb_1991_srn_results.parameters");
                FileFinder reference("cell_based/test/data/TestSrnModels/gb_1991_srn_results.parameters",
                                     RelativeTo::ChasteSourceRoot);
                FileComparison comparer(generated, reference);
                TS_ASSERT(comparer.CompareFiles());
            }
        }
    }
};

#endif /* TESTODEBASEDSRNMODELS_HPP_ */
