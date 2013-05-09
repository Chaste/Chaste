/*

Copyright (c) 2005-2013, University of Oxford.
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
#ifndef TESTDELTANOTCHCELLCYCLEMODEL_HPP_
#define TESTDELTANOTCHCELLCYCLEMODEL_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>

#include "DeltaNotchCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestDeltaNotchCellCycleModel : public AbstractCellBasedTestSuite
{
public:

    ///\todo test correct behaviour of ODE system and state variables
    void TestCorrectBehaviour() throw(Exception)
    {
        TS_ASSERT_THROWS_NOTHING(DeltaNotchCellCycleModel cell_model3);

        DeltaNotchCellCycleModel* p_stem_model = new DeltaNotchCellCycleModel;
        p_stem_model->SetDimension(2);

        // Change G1 duration for this model
        p_stem_model->SetStemCellG1Duration(1.0);

        DeltaNotchCellCycleModel* p_transit_model = new DeltaNotchCellCycleModel;
        p_transit_model->SetDimension(3);

        // Change G1 duration for this model
        p_stem_model->SetTransitCellG1Duration(1.0);  ///\todo Is this a copy and paste error?

        DeltaNotchCellCycleModel* p_diff_model = new DeltaNotchCellCycleModel;
        p_diff_model->SetDimension(1);

        MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        CellPtr p_stem_cell(new Cell(p_healthy_state, p_stem_model));
        p_stem_cell->SetCellProliferativeType(p_stem_type);
        p_stem_cell->GetCellData()->SetItem("mean delta", 0.0);
        p_stem_cell->InitialiseCellCycleModel();

        CellPtr p_transit_cell(new Cell(p_healthy_state, p_transit_model));
        p_transit_cell->SetCellProliferativeType(p_transit_type);
        p_transit_cell->GetCellData()->SetItem("mean delta", 0.0);
        p_transit_cell->InitialiseCellCycleModel();

        CellPtr p_diff_cell(new Cell(p_healthy_state, p_diff_model));
        p_diff_cell->SetCellProliferativeType(p_diff_type);
        p_diff_cell->GetCellData()->SetItem("mean delta", 0.0);
        p_diff_cell->InitialiseCellCycleModel();

        SimulationTime* p_simulation_time = SimulationTime::Instance();
        unsigned num_steps = 100;
        double end_time = 2.0*(p_stem_model->GetStemCellG1Duration() + p_stem_model->GetSG2MDuration());
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(end_time, num_steps);

        for (unsigned i=0; i<num_steps; i++)
        {
            p_simulation_time->IncrementTimeOneStep();

            // The numbers for the G1 durations below are taken from the first three random numbers generated
            CheckReadyToDivideAndPhaseIsUpdated(p_stem_model, 3.19525);
            CheckReadyToDivideAndPhaseIsUpdated(p_transit_model, 3.18569);
            CheckReadyToDivideAndPhaseIsUpdated(p_diff_model, 132);  // any old number
        }

        // Coverage of Get Methods
        TS_ASSERT_DELTA(p_diff_model->GetNotch(), 0.0, 1e-4);
        TS_ASSERT_DELTA(p_diff_model->GetDelta(), 1.0, 1e-4);
        TS_ASSERT_DELTA(p_diff_model->GetMeanNeighbouringDelta(), 0.0, 1e-4);

        // Setting mean delta via cell
        p_diff_cell->GetCellData()->SetItem("mean delta", 4.2);
        p_diff_model->UpdateDeltaNotch();
        TS_ASSERT_DELTA(p_diff_model->GetMeanNeighbouringDelta(), 4.2, 1e-4);
    }

    ///\todo test archiving of ODE system and state variables
    void TestArchiveDeltaNotchCellCycleModel()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "delta_notch_cell_cycle.arch";

        double random_number_test = 0;

        // Create an output archive
        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(2.0, 4);

            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel;
            p_model->SetDimension(2);

            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);

            CellPtr p_cell(new Cell(p_healthy_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->GetCellData()->SetItem("mean delta", 0.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->SetBirthTime(-1.1);
            p_simulation_time->IncrementTimeOneStep();
            p_simulation_time->IncrementTimeOneStep();

            p_cell->ReadyToDivide(); // updates phases

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            TS_ASSERT_DELTA(p_model->GetSDuration(), 5.0, 1e-12);

            CellPtr const p_const_cell = p_cell;
            output_arch << p_const_cell;

            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.1, 1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(), 2.1, 1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), G_ONE_PHASE);

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            random_number_test = p_gen->ranf();
            SimulationTime::Destroy();
        }

        {
            SimulationTime* p_simulation_time = SimulationTime::Instance();
            p_simulation_time->SetStartTime(0.0);
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            p_gen->Reseed(128);

            CellPtr p_cell;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_cell;

            TS_ASSERT_DELTA(RandomNumberGenerator::Instance()->ranf(), random_number_test, 1e-7);

            AbstractCellCycleModel* p_model = p_cell->GetCellCycleModel();

            // Test that the cell-cycle model was restored correctly
            TS_ASSERT_DELTA(p_model->GetBirthTime(), -1.1, 1e-12);
            TS_ASSERT_DELTA(p_model->GetAge(), 2.1, 1e-12);
            TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), G_ONE_PHASE);
            TS_ASSERT_DELTA(p_model->GetSDuration(), 5.0, 1e-12);
        }
    }

    void TestCellCycleModelOutputParameters()
    {
        std::string output_directory = "TestCellCycleModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with DeltaNotchCellCycleModel
        DeltaNotchCellCycleModel cell_cycle_model;
        TS_ASSERT_EQUALS(cell_cycle_model.GetIdentifier(), "DeltaNotchCellCycleModel");

        out_stream parameter_file = output_file_handler.OpenOutputFile("delta_notch_results.parameters");
        cell_cycle_model.OutputCellCycleModelParameters(parameter_file);
        parameter_file->close();

        {
            // Compare the generated file in test output with a reference copy in the source code.
            FileFinder generated = output_file_handler.FindFile("delta_notch_results.parameters");
            FileFinder reference("cell_based/test/data/TestCellCycleModels/delta_notch_results.parameters",
                    RelativeTo::ChasteSourceRoot);
            FileComparison comparer(generated, reference);
            TS_ASSERT(comparer.CompareFiles());
        }
    }

    void TestCreateCopyCellCycleModel()
    {
        // Test with DeltaNotchCellCycleModel
        DeltaNotchCellCycleModel* p_model= new DeltaNotchCellCycleModel;

        // Set ODE system
        std::vector<double> state_variables;
        state_variables.push_back(1.0);
        state_variables.push_back(1.0);
        p_model->SetOdeSystem(new DeltaNotchOdeSystem(state_variables));

        // Set model parameters.
        p_model->SetBirthTime(2.0);
        p_model->SetDimension(2);
        p_model->SetGeneration(2);
        p_model->SetMaxTransitGenerations(10);

        // Create a copy
        DeltaNotchCellCycleModel* p_model2 = static_cast<DeltaNotchCellCycleModel*> (p_model->CreateCellCycleModel());

        // Check correct initializations
        TS_ASSERT_EQUALS(p_model2->GetBirthTime(),2);
        TS_ASSERT_EQUALS(p_model2->GetDimension(), 2u);
        TS_ASSERT_EQUALS(p_model2->GetGeneration(), 2u);
        TS_ASSERT_EQUALS(p_model2->GetMaxTransitGenerations(), 10u);

        // Destruct models
        delete p_model;
        delete p_model2;
    }
};

#endif /* TESTDELTANOTCHCELLCYCLEMODEL_HPP_ */
