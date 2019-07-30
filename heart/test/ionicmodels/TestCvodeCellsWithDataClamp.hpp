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

#ifndef TESTCVODECELLSWITHDATACLAMP_HPP_
#define TESTCVODECELLSWITHDATACLAMP_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCvodeCellWithDataClamp.hpp"
#include "Shannon2004CvodeDataClamp.hpp"
#include "RandomNumberGenerator.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "ArchiveLocationInfo.hpp"
#include "OutputFileHandler.hpp"
#include "FixedModifier.hpp" // We are using this test to test archiving cells with modifiers too.
#include <Timer.hpp>

#include "FakePetscSetup.hpp"

class TestCvodeCellsWithDataClamp : public CxxTest::TestSuite
{
private:
#ifdef CHASTE_CVODE
    boost::shared_ptr<CellShannon2004FromCellMLCvodeDataClamp> mpModel;
#endif

public:
    void TestInterpolatorTimesAndGenerateReferenceTrace()
    {
#ifdef CHASTE_CVODE
        OutputFileHandler handler("CvodeCellsWithDataClamp");

        boost::shared_ptr<AbstractIvpOdeSolver> p_empty_solver;
        boost::shared_ptr<AbstractStimulusFunction> p_empty_stimulus;

        // N.B. Because we use the Shannon model as a lot of examples,
        // here it is actually a Shannon->WithModifiers->WithDataClamp->CvodeCell
        // (the WithModifiers doesn't need to be there to use the data clamp!)
        mpModel.reset(new CellShannon2004FromCellMLCvodeDataClamp(p_empty_solver,p_empty_stimulus));

        TS_ASSERT_EQUALS(mpModel->HasParameter("membrane_data_clamp_current_conductance"), true);

        mpModel->SetMaxSteps(5000);
        mpModel->UseCellMLDefaultStimulus();

        // Run a simulation without clamping switched on

        Timer::Reset();
        double end_time = 400.0;
        OdeSolution solution = mpModel->Compute(0, end_time, 0.2);
        Timer::Print("OdeSolution");
        std::vector<double> expt_times = solution.rGetTimes();
        std::vector<double> expt_data = solution.GetAnyVariable("membrane_voltage");
        solution.WriteToFile("CvodeCellsWithDataClamp","shannon_original_no_clamp", "ms", 1, false); // false to clean

        TS_ASSERT_THROWS_THIS(mpModel->TurnOnDataClamp(),
            "Before calling TurnOnDataClamp(), please provide experimental data via the SetExperimentalData() method.");

        // Test the interpolation methods.
        {
            mpModel->SetExperimentalData(expt_times, expt_data);

            // Note - unless the data clamp is switched on the below method just returns DOUBLE_UNSET to save time interpolating.
            double time = 100.0;
            TS_ASSERT_EQUALS(mpModel->GetExperimentalVoltageAtTimeT(time), DOUBLE_UNSET);

            // So now turn on the data clamp
            mpModel->TurnOnDataClamp();

# if CHASTE_SUNDIALS_VERSION >= 20400
            double tol = 5e-3; // mV
#else
            double tol = 0.2; // mV
#endif
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(time), -8.55863245e+01, tol);

            // So turn it off again
            mpModel->TurnOffDataClamp();
            TS_ASSERT_DELTA(mpModel->GetParameter("membrane_data_clamp_current_conductance"), 0.0, 1e-12);
            mpModel->TurnOnDataClamp(200.0);
            TS_ASSERT_DELTA(mpModel->GetParameter("membrane_data_clamp_current_conductance"), 200.0, 1e-12);
            mpModel->TurnOffDataClamp();
            TS_ASSERT_DELTA(mpModel->GetParameter("membrane_data_clamp_current_conductance"), 0.0, 1e-12);
            mpModel->TurnOnDataClamp();
            TS_ASSERT_DELTA(mpModel->GetParameter("membrane_data_clamp_current_conductance"), 100.0, 1e-12); // the default

            // Test a couple of times where no interpolation is needed (on data points).
            time = 116.0;
            double v_at_116 = 1.53670634e+01;
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(time), v_at_116, tol);

            time = 116.2;
            double v_at_116_2 = 1.50089546e+01;
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(time), v_at_116_2, tol);

            // Now test a time where interpolation is required.
            time = 116.1;
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(time), 0.5*(v_at_116 + v_at_116_2), tol);

            // Test ends
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(0.0), expt_data[0], 1e-4);
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(end_time), expt_data.back(), 1e-4);

            // Test exceptions
            TS_ASSERT_THROWS_CONTAINS(mpModel->GetExperimentalVoltageAtTimeT(-1e-12),
                                      "is outside the times stored in the data clamp");
            TS_ASSERT_THROWS_CONTAINS(mpModel->GetExperimentalVoltageAtTimeT(end_time+1e-12),
                                      "is outside the times stored in the data clamp");

            //std::cout << "membrane_data_clamp_current_conductance = " << mpModel->GetParameter("membrane_data_clamp_current_conductance") << std::endl << std::flush;
            //std::cout << "mpModel->GetExperimentalVoltageAtTimeT(time) = " << mpModel->GetExperimentalVoltageAtTimeT(time) << std::endl << std::flush;

            unsigned how_many = 10000u;
            Timer::Reset();
            for (unsigned i=0; i<how_many; i++)
            {
                mpModel->GetExperimentalVoltageAtTimeT(time);
            }
            Timer::PrintAndReset("GetExperimentalVoltageAtTimeT");
        }

        // Generate some noisier data - more like a real cell.
        out_stream experimental_voltage_results_file = handler.OpenOutputFile("Shannon_noisy_data.dat");
        double experimental_noise_sd = 0.25; // directly from Teun's AP data trace of a stationary-looking bit
        for (unsigned i=0; i<expt_data.size(); i++)
        {
            double random_number = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
            expt_data[i] +=  (10+experimental_noise_sd*random_number); // make experimental data very different to see if clamp works
            *experimental_voltage_results_file << expt_times[i] << "\t" << expt_data[i] << std::endl;
        }
        experimental_voltage_results_file->close();

        // In this half of the test, try running simulations with the clamp.
        {
            // Now solve with data clamp to the noisy data
            mpModel->SetExperimentalData(expt_times, expt_data);
            mpModel->ResetToInitialConditions();
            mpModel->TurnOnDataClamp();

            mpModel->SetParameter("membrane_data_clamp_current_conductance",0.001);

            std::vector<double> data_clamp_times;
            std::vector<double> data_clamp_voltage;

            for (int i = -4; i < 3; i++)
            {
                double clamp_conductance = pow(10,i);
                std::cout << "clamp_conductance = " << clamp_conductance << std::endl << std::flush;
                std::stringstream output_file;
                output_file << "Shannon_test_solution_with_data_clamp_conductance_exponent_" << i << ".dat";

                mpModel->ResetToInitialConditions();
                mpModel->SetParameter("membrane_data_clamp_current_conductance",clamp_conductance);
                solution = mpModel->Compute(0, 400, 0.2);
                data_clamp_times = solution.rGetTimes();
                data_clamp_voltage = solution.GetAnyVariable("membrane_voltage");

                out_stream data_clamp_voltage_results_file = handler.OpenOutputFile(output_file.str());
                for (unsigned j=0; j<data_clamp_voltage.size(); j++)
                {
                    *data_clamp_voltage_results_file << data_clamp_times[j] << "\t" << data_clamp_voltage[j] << "\n";
                }
                data_clamp_voltage_results_file->close();
            }
        }
#else
        std::cout << "Cvode is not enabled.\n";
#endif
    }

    void TestArchivingCvodeCellsWithDataClamp()
    {
        // We also hijack this test to test the archiving and restoration of modifiers.
#ifdef CHASTE_CVODE
        //Archive
        OutputFileHandler handler("archive", false);
        handler.SetArchiveDirectory();
        std::string archive_filename = ArchiveLocationInfo::GetProcessUniqueFilePath("shannon_with_data_clamp.arch");

        bool data_clamp_state;
        bool data_available;
        double fixed_modifier_value = 56.0;
        std::vector<double> times;
        std::vector<double> voltages;

        // Save
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Using friend status to directly look at member variables.
            data_clamp_state = mpModel->mDataClampIsOn;
            data_available = mpModel->mDataAvailable;
            times = mpModel->mExperimentalTimes;
            voltages = mpModel->mExperimentalVoltages;

            // Check we are actually checking something!
            TS_ASSERT_EQUALS(data_clamp_state, true);
            TS_ASSERT_EQUALS(data_available, true);
            TS_ASSERT(times.size()>10u);
            TS_ASSERT_EQUALS(times.size(), voltages.size());
            TS_ASSERT_EQUALS(mpModel->HasModifier("membrane_slow_delayed_rectifier_potassium_current"), true);
            TS_ASSERT_DELTA(mpModel->GetModifier("membrane_slow_delayed_rectifier_potassium_current")->Calc(0,1), 0.0, 1e-12);

            boost::shared_ptr<AbstractModifier> p_fixed(new FixedModifier(fixed_modifier_value));
            mpModel->SetModifier("membrane_slow_delayed_rectifier_potassium_current", p_fixed);

            TS_ASSERT_DELTA(mpModel->GetModifier("membrane_slow_delayed_rectifier_potassium_current")->Calc(0,1), fixed_modifier_value, 1e-12);

            // Archive as an AbstractCvodeCell pointer (not via boost shared pointer)
            AbstractCvodeCell* const p_cell = mpModel.get();
            output_arch <<  p_cell;
        }

        // This should free up the memory and delete cell model.
        mpModel.reset();

        // Load
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            AbstractCvodeCell* p_cell;
            input_arch >> p_cell;

            TS_ASSERT_EQUALS(p_cell->GetNumberOfStateVariables(), 39u);

            // Check modifiers were archived correctly
            if (dynamic_cast<AbstractCardiacCellWithModifiers<AbstractCvodeCellWithDataClamp>* >(p_cell) == NULL)
            {
                // Pointer could not be cast as the right kind, so throw error.
                TS_ASSERT(false);
            }

            AbstractCardiacCellWithModifiers<AbstractCvodeCellWithDataClamp>* p_modifiers_cell = static_cast<AbstractCardiacCellWithModifiers<AbstractCvodeCellWithDataClamp>*>(p_cell);
            TS_ASSERT_EQUALS(p_modifiers_cell->HasModifier("membrane_slow_delayed_rectifier_potassium_current"),true);
            boost::shared_ptr<AbstractModifier> p_modifier = p_modifiers_cell->GetModifier("membrane_slow_delayed_rectifier_potassium_current");
            TS_ASSERT_DELTA(p_modifier->Calc(2 /*param*/,3 /*time*/), fixed_modifier_value, 1e-12); // Fixed modifier returns 56 from above..

            // Check data clamp was archived correctly.
            if (dynamic_cast<AbstractCvodeCellWithDataClamp*>(p_cell) == NULL)
            {
                // Pointer could not be cast as the right kind, so throw error.
                TS_ASSERT(false);
            }

            AbstractCvodeCellWithDataClamp* p_data_clamp_cell = static_cast<AbstractCvodeCellWithDataClamp*>(p_cell);

            TS_ASSERT_EQUALS(p_data_clamp_cell->mDataClampIsOn, data_clamp_state);
            TS_ASSERT_EQUALS(p_data_clamp_cell->mDataAvailable, data_available);
            TS_ASSERT_EQUALS(p_data_clamp_cell->mExperimentalTimes.size(), times.size());
            TS_ASSERT_EQUALS(p_data_clamp_cell->mExperimentalVoltages.size(), voltages.size());

            for (unsigned i=0; i<times.size(); i++)
            {
                TS_ASSERT_DELTA(p_data_clamp_cell->mExperimentalTimes[i], times[i], 1e-12);
                TS_ASSERT_DELTA(p_data_clamp_cell->mExperimentalVoltages[i], voltages[i], 1e-12);
            }

            // Check we have a functioning unarchived cell (this will check that all internal member variables that pointed to modifiers still work).
            p_cell->Compute(0, 100, 1);

            delete p_cell;
        }
#else
       std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
   }
};

#endif /* TESTABSTRACTCVODECELLWITHDATACLAMP_HPP_ */
