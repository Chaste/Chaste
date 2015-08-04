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
#include "CheckpointArchiveTypes.hpp"
#include "ArchiveLocationInfo.hpp"
#include <Timer.hpp>

#include "FakePetscSetup.hpp"

class TestCvodeCellsWithDataClamp : public CxxTest::TestSuite
{
private:
#ifdef CHASTE_CVODE
    boost::shared_ptr<AbstractCvodeCell> mpModel;
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
        mpModel = loader.LoadCvodeCell();

        boost::shared_ptr<AbstractCvodeCellWithDataClamp> p_model = boost::static_pointer_cast<AbstractCvodeCellWithDataClamp>(mpModel);

        TS_ASSERT_EQUALS(p_model->HasParameter("membrane_data_clamp_current_conductance"), true);

        p_model->SetMaxSteps(5000);
        p_model->UseCellMLDefaultStimulus();

        // Run a simulation without clamping switched on

        Timer::Reset();
        double end_time = 400.0;
        OdeSolution solution = p_model->Compute(0, end_time, 0.2);
        Timer::Print("OdeSolution");
        std::vector<double> expt_times = solution.rGetTimes();
        std::vector<double> expt_data = solution.GetAnyVariable("membrane_voltage");
        solution.WriteToFile("CvodeCellsWithDataClamp","shannon_original_no_clamp", "ms", 1, false); // false to clean

        TS_ASSERT_THROWS_THIS(p_model->TurnOnDataClamp(),
            "Before calling TurnOnDataClamp(), please provide experimental data via the SetExperimentalData() method.");

        // Test the interpolation methods.
        {
            p_model->SetExperimentalData(expt_times, expt_data);

            // Note - unless the data clamp is switched on the below method just returns DOUBLE_UNSET to save time interpolating.
            double time = 100.0;
            TS_ASSERT_EQUALS(p_model->GetExperimentalVoltageAtTimeT(time), DOUBLE_UNSET);

            // So now turn on the data clamp
            p_model->TurnOnDataClamp();

            TS_ASSERT_DELTA(p_model->GetExperimentalVoltageAtTimeT(time), -8.55863245e+01, 2e-3);

            // So turn it off again
            p_model->TurnOffDataClamp();
            TS_ASSERT_DELTA(p_model->GetParameter("membrane_data_clamp_current_conductance"), 0.0, 1e-12);
            p_model->TurnOnDataClamp(200.0);
            TS_ASSERT_DELTA(p_model->GetParameter("membrane_data_clamp_current_conductance"), 200.0, 1e-12);
            p_model->TurnOffDataClamp();
            TS_ASSERT_DELTA(p_model->GetParameter("membrane_data_clamp_current_conductance"), 0.0, 1e-12);
            p_model->TurnOnDataClamp();
            TS_ASSERT_DELTA(p_model->GetParameter("membrane_data_clamp_current_conductance"), 100.0, 1e-12); // the default

            // Test a couple of times where no interpolation is needed (on data points).
            time = 116.0;
            double v_at_116 = 1.53670634e+01;
            TS_ASSERT_DELTA(p_model->GetExperimentalVoltageAtTimeT(time), v_at_116, 2e-3);
            time = 116.2;
            double v_at_116_2 = 1.50089546e+01;
            TS_ASSERT_DELTA(p_model->GetExperimentalVoltageAtTimeT(time), v_at_116_2, 2e-3);

            // Now test a time where interpolation is required.
            time = 116.1;
            TS_ASSERT_DELTA(p_model->GetExperimentalVoltageAtTimeT(time), 0.5*(v_at_116 + v_at_116_2), 2e-3);

            // Test ends
            TS_ASSERT_DELTA(p_model->GetExperimentalVoltageAtTimeT(0.0), expt_data[0], 1e-4);
            TS_ASSERT_DELTA(p_model->GetExperimentalVoltageAtTimeT(end_time), expt_data.back(), 1e-4);

            // Test exceptions
            TS_ASSERT_THROWS_CONTAINS(p_model->GetExperimentalVoltageAtTimeT(-1e-12),
                                      "is outside the times stored in the data clamp");
            TS_ASSERT_THROWS_CONTAINS(p_model->GetExperimentalVoltageAtTimeT(end_time+1e-12),
                                      "is outside the times stored in the data clamp");

            //std::cout << "membrane_data_clamp_current_conductance = " << p_model->GetParameter("membrane_data_clamp_current_conductance") << std::endl << std::flush;
            //std::cout << "mpModel->GetExperimentalVoltageAtTimeT(time) = " << p_model->GetExperimentalVoltageAtTimeT(time) << std::endl << std::flush;

            unsigned how_many = 10000;
            Timer::Reset();
            for (unsigned i=0; i<how_many; i++)
            {
                p_model->GetExperimentalVoltageAtTimeT(time);
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
            p_model->SetExperimentalData(expt_times, expt_data);
            p_model->ResetToInitialConditions();
            p_model->TurnOnDataClamp();

            p_model->SetParameter("membrane_data_clamp_current_conductance",0.001);

            double clamp_conductance;
            std::vector<double> data_clamp_times;
            std::vector<double> data_clamp_voltage;
            out_stream data_clamp_voltage_results_file;

            for (int i = -4; i < 3; i++)
            {
                clamp_conductance = pow(10,i);
                std::cout << "clamp_conductance = " << clamp_conductance << std::endl << std::flush;
                std::stringstream output_file;
                output_file << "Shannon_test_solution_with_data_clamp_conductance_exponent_" << i << ".dat";

                p_model->ResetToInitialConditions();
                p_model->SetParameter("membrane_data_clamp_current_conductance",clamp_conductance);
                solution = p_model->Compute(0, 400, 0.2);
                data_clamp_times = solution.rGetTimes();
                data_clamp_voltage = solution.GetAnyVariable("membrane_voltage");
                data_clamp_voltage_results_file = handler.OpenOutputFile(output_file.str());
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

    void TestArchivingCvodeCellsWithDataClamp() throw(Exception)
    {
#ifdef CHASTE_CVODE
       //Archive
       OutputFileHandler handler("archive", false);
       handler.SetArchiveDirectory();
       std::string archive_filename =  ArchiveLocationInfo::GetProcessUniqueFilePath("shannon_with_data_clamp.arch");

       bool data_clamp_state;
       bool data_available;
       std::vector<double> times;
       std::vector<double> voltages;

       // Save
       {
           std::ofstream ofs(archive_filename.c_str());
           boost::archive::text_oarchive output_arch(ofs);

           boost::shared_ptr<AbstractCvodeCellWithDataClamp> p_data_clamp_cell =
                   boost::static_pointer_cast<AbstractCvodeCellWithDataClamp>(mpModel);

           // Archive as an AbstractCvodeCell pointer (not via boost shared pointer)
           AbstractCvodeCell* const p_cell = mpModel.get();
           output_arch <<  p_cell;

           // Using friend status to directly look at member variables.
           data_clamp_state = p_data_clamp_cell->mDataClampIsOn;
           data_available = p_data_clamp_cell->mDataAvailable;
           times = p_data_clamp_cell->mExperimentalTimes;
           voltages = p_data_clamp_cell->mExperimentalVoltages;

           // Check we are actually checking something!
           TS_ASSERT_EQUALS(data_clamp_state, true);
           TS_ASSERT_EQUALS(data_available, true);
           TS_ASSERT(times.size()>10u);
           TS_ASSERT_EQUALS(times.size(), voltages.size());
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

           for(unsigned i=0; i<times.size(); i++)
           {
               TS_ASSERT_DELTA(p_data_clamp_cell->mExperimentalTimes[i], times[i], 1e-12);
               TS_ASSERT_DELTA(p_data_clamp_cell->mExperimentalVoltages[i], voltages[i], 1e-12);
           }
           delete p_cell;
       }
#else
       std::cout << "Cvode is not enabled.\n";
#endif // CHASTE_CVODE
   }
};

#endif /* TESTABSTRACTCVODECELLWITHDATACLAMP_HPP_ */
