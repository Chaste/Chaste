/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "CellMLLoader.hpp"
#include "AbstractCvodeCellWithDataClamp.hpp"
#include "RandomNumberGenerator.hpp"
#include <Timer.hpp>

#include "FakePetscSetup.hpp"

class TestCvodeCellsWithDataClamp : public CxxTest::TestSuite
{
private:
#ifdef CHASTE_CVODE
    boost::shared_ptr<AbstractCvodeCellWithDataClamp> mpModel;
#endif

public:
    void TestInterpolatorTimesAndGenerateReferenceTrace() throw(Exception)
    {
#ifdef CHASTE_CVODE
        // Test generation of a clamp model.

        FileFinder shannon_cellml("heart/src/odes/cellml/Shannon2004.cellml", RelativeTo::ChasteSourceRoot);
        OutputFileHandler handler("CvodeCellsWithDataClamp");
        std::vector<std::string> options;
        options.push_back("--cvode-data-clamp");
        options.push_back("--expose-annotated-variables");

        CellMLLoader loader(shannon_cellml, handler, options);
        mpModel =  boost::static_pointer_cast<AbstractCvodeCellWithDataClamp>(loader.LoadCvodeCell());

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

            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(time), -8.55863245e+01, 2e-3);

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
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(time), v_at_116, 2e-3);
            time = 116.2;
            double v_at_116_2 = 1.50089546e+01;
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(time), v_at_116_2, 2e-3);

            // Now test a time where interpolation is required.
            time = 116.1;
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(time), 0.5*(v_at_116 + v_at_116_2), 2e-3);

            // Test ends
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(0.0), expt_data[0], 1e-4);
            TS_ASSERT_DELTA(mpModel->GetExperimentalVoltageAtTimeT(end_time), expt_data.back(), 1e-4);

            // Test exceptions
            TS_ASSERT_THROWS_CONTAINS(mpModel->GetExperimentalVoltageAtTimeT(-1e-12),
                                      "is outside the times stored in the data clamp");
            TS_ASSERT_THROWS_CONTAINS(mpModel->GetExperimentalVoltageAtTimeT(end_time+1e-12),
                                      "is outside the times stored in the data clamp");

            std::cout << "membrane_data_clamp_current_conductance = " << mpModel->GetParameter("membrane_data_clamp_current_conductance") << std::endl << std::flush;
            std::cout << "mpModel->GetExperimentalVoltageAtTimeT(time) = " << mpModel->GetExperimentalVoltageAtTimeT(time) << std::endl << std::flush;

            unsigned how_many = 10000;
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

        std::cout << "Just created an AbstractCvodeCellWithDataClamp!" << std::endl << std::flush;

        // In this half of the test, try running simulations with the clamp.
        {
            // Now solve with data clamp to the noisy data
            mpModel->SetExperimentalData(expt_times, expt_data);
            mpModel->ResetToInitialConditions();
            mpModel->TurnOnDataClamp();

            std::cout << "membrane_data_clamp_current_conductance = " << mpModel->GetParameter("membrane_data_clamp_current_conductance") << std::endl << std::flush;
            mpModel->SetParameter("membrane_data_clamp_current_conductance",0.001);
            std::cout << "membrane_data_clamp_current_conductance = " << mpModel->GetParameter("membrane_data_clamp_current_conductance") << std::endl << std::flush;

            double clamp_conductance;
            std::string output_file;
            std::vector<double> data_clamp_times;
            std::vector<double> data_clamp_voltage;
            out_stream data_clamp_voltage_results_file;

            for (int i = -4; i < 3; i++)
            {
                clamp_conductance = pow(10,i);
                std::cout << "clamp_conductance = " << clamp_conductance << std::endl << std::flush;
                output_file = "Shannon_test_solution_with_data_clamp_conductance_exponent_" + boost::lexical_cast<std::string>(i) + ".dat";

                mpModel->ResetToInitialConditions();
                mpModel->SetParameter("membrane_data_clamp_current_conductance",clamp_conductance);
                solution = mpModel->Compute(0, 400, 0.2);
                data_clamp_times = solution.rGetTimes();
                data_clamp_voltage = solution.GetAnyVariable("membrane_voltage");
                data_clamp_voltage_results_file = handler.OpenOutputFile(output_file);
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
};

#endif /* TESTABSTRACTCVODECELLWITHDATACLAMP_HPP_ */
