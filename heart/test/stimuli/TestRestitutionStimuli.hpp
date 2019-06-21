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

#ifndef _TESTRESTITUTIONSTIMULI_HPP_
#define _TESTRESTITUTIONSTIMULI_HPP_

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "Exception.hpp"
#include "S1S2Stimulus.hpp"
#include "SteadyStateRestitutionStimulus.hpp"
#include "OutputFileHandler.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class TestRestitutionStimuli : public CxxTest::TestSuite
{
private:
    void DoesSteadyStateStimulusPerformCorrectly(SteadyStateRestitutionStimulus* pStimulus, double magnitude, double duration_of_stimulus,
                                                 double startTime, const std::vector<double>& pacing_cycle_lengths, unsigned number_of_pulses)
    {
        double time = startTime;
        for (unsigned pacing_length_index = 0; pacing_length_index<pacing_cycle_lengths.size(); pacing_length_index++)
        {
            for (unsigned pulse_index=0; pulse_index<number_of_pulses; pulse_index++)
            {
                double pulse_time = time;
                TS_ASSERT_DELTA(pStimulus->GetStimulus(pulse_time-0.01),0.0, 1e-9);
                TS_ASSERT_DELTA(pStimulus->GetStimulus(pulse_time+0.01),magnitude, 1e-9);
                TS_ASSERT_DELTA(pStimulus->GetStimulus(pulse_time+duration_of_stimulus-0.01),magnitude, 1e-9);
                TS_ASSERT_DELTA(pStimulus->GetStimulus(pulse_time+duration_of_stimulus+0.01),0.0, 1e-9);

                time = time + pacing_cycle_lengths[pacing_length_index];
            }
        }
    }

    void DoesStimulusPerformCorrectly(S1S2Stimulus* pStimulus, double magnitude, double duration_of_stimulus,
                                      double duration_of_s1, double s1_period, double startTime, const std::vector<double>& s2_periods)
    {
        for (unsigned experiment_index=0 ; experiment_index<s2_periods.size() ; experiment_index++)
        {
            std::cout << "Experiment = " << experiment_index << "\n" << std::flush;
            pStimulus->SetS2ExperimentPeriodIndex(experiment_index);

            // First S1 Pulse
            TS_ASSERT_DELTA(pStimulus->GetStimulus(0.0), 0.0, 1e-9);
            TS_ASSERT_DELTA(pStimulus->GetStimulus(startTime), magnitude, 1e-9);
            TS_ASSERT_DELTA(pStimulus->GetStimulus(startTime + 0.1), magnitude, 1e-9);
            TS_ASSERT_DELTA(pStimulus->GetStimulus(startTime + duration_of_stimulus), magnitude, 1e-9);
            TS_ASSERT_DELTA(pStimulus->GetStimulus(startTime + duration_of_stimulus + 0.1), 0.0, 1e-9);

            // First S2 Pulse (or final S1 pulse for physiologists)
            TS_ASSERT_DELTA(pStimulus->GetStimulus(duration_of_s1 + startTime - 0.1), 0.0, 1e-9);
            TS_ASSERT_DELTA(pStimulus->GetStimulus(duration_of_s1 + startTime), magnitude, 1e-9);
            TS_ASSERT_DELTA(pStimulus->GetStimulus(duration_of_s1 + startTime + 0.1), magnitude, 1e-9);
            TS_ASSERT_DELTA(pStimulus->GetStimulus(duration_of_s1 + startTime + duration_of_stimulus), magnitude, 1e-9);
            TS_ASSERT_DELTA(pStimulus->GetStimulus(duration_of_s1 + startTime + duration_of_stimulus + 0.1), 0.0, 1e-9);

            // Second S2 Pulse (or S2 pulse for physiologists!)
            double time;
            time = duration_of_s1 + s2_periods[experiment_index] + startTime - 0.1;
            TS_ASSERT_DELTA(pStimulus->GetStimulus(time), 0.0, 1e-9);
            time = duration_of_s1 + s2_periods[experiment_index] + startTime;
            TS_ASSERT_DELTA(pStimulus->GetStimulus(time), magnitude, 1e-9);
            time = duration_of_s1 + s2_periods[experiment_index] + startTime + 0.1;
            TS_ASSERT_DELTA(pStimulus->GetStimulus(time), magnitude, 1e-9);
            time = duration_of_s1 + s2_periods[experiment_index] + duration_of_stimulus + startTime - 0.1;
            TS_ASSERT_DELTA(pStimulus->GetStimulus(time), magnitude, 1e-9);
            time = duration_of_s1 + s2_periods[experiment_index] + duration_of_stimulus + startTime;
            TS_ASSERT_DELTA(pStimulus->GetStimulus(time), magnitude, 1e-9);
            time = duration_of_s1 + 2*s2_periods[experiment_index]+ duration_of_stimulus + startTime + 0.1;
            TS_ASSERT_DELTA(pStimulus->GetStimulus(time), 0.0, 1e-9);
        }
    }

public:

    void TestS1S2StimulusSetup(void)
    {
        double magnitude = 50;
        double duration_of_stimulus = 5;
        double duration_of_s1 = 10000;
        double s1_period = 1000;
        double startTime = 50;
        std::vector<double> s2_periods;
        s2_periods.push_back(1000);
        s2_periods.push_back(900);
        s2_periods.push_back(800);
        s2_periods.push_back(700);

        S1S2Stimulus stimulus(magnitude, duration_of_stimulus, duration_of_s1, s1_period, startTime, s2_periods);

        DoesStimulusPerformCorrectly(&stimulus, magnitude, duration_of_stimulus,
                                     duration_of_s1, s1_period, startTime, s2_periods);

        TS_ASSERT_THROWS_THIS(stimulus.SetS2ExperimentPeriodIndex(s2_periods.size()),
                              "There are fewer S2 frequency values than the one you have requested.");

        TS_ASSERT_EQUALS(stimulus.GetNumS2FrequencyValues(),s2_periods.size());
    }

    void TestSteadyStateRestitutionCellStimulusSetup(void)
    {
        double magnitude = 50;
        double duration_of_stimulus = 5;
        double startTime = 50;

        std::vector<double> pacing_cycle_lengths;
        pacing_cycle_lengths.push_back(1000);
        pacing_cycle_lengths.push_back(900);
        pacing_cycle_lengths.push_back(800);
        pacing_cycle_lengths.push_back(700);

        unsigned number_of_pulses = 10;

        SteadyStateRestitutionStimulus stimulus(magnitude, duration_of_stimulus, startTime, pacing_cycle_lengths, number_of_pulses);

        DoesSteadyStateStimulusPerformCorrectly(&stimulus, magnitude, duration_of_stimulus,
                                                startTime, pacing_cycle_lengths, number_of_pulses);

    }


    void TestArchivingS1S2Stimulus(void)
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "s1s2_stimulus.arch";

        double magnitude = 50;
        double duration_of_stimulus = 5;
        double duration_of_s1 = 10000;
        double s1_period = 1000;
        double startTime = 50;
        std::vector<double> s2_periods;
        s2_periods.push_back(1000);
        s2_periods.push_back(900);
        s2_periods.push_back(800);
        s2_periods.push_back(700);

        // Create and archive stimulus
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            AbstractStimulusFunction* const p_stimulus = new S1S2Stimulus(magnitude, duration_of_stimulus, duration_of_s1, s1_period, startTime, s2_periods);

            // Should always archive a pointer
            output_arch << p_stimulus;

            delete p_stimulus;
        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create a pointer
            AbstractStimulusFunction* p_loaded_stim;
            input_arch >> p_loaded_stim;

            DoesStimulusPerformCorrectly(static_cast<S1S2Stimulus*>(p_loaded_stim), magnitude, duration_of_stimulus,
                                                  duration_of_s1, s1_period, startTime, s2_periods);

            delete p_loaded_stim;
        }
    }

    void TestArchivingSteadyStateStimulus(void)
    {
        OutputFileHandler handler("archive",false);
        std::string archive_filename;
        archive_filename = handler.GetOutputDirectoryFullPath() + "dynamic_stimulus.arch";

        double magnitude = 50;
        double duration_of_stimulus = 5;
        double startTime = 50;

        std::vector<double> pacing_cycle_lengths;
        pacing_cycle_lengths.push_back(1000);
        pacing_cycle_lengths.push_back(900);
        pacing_cycle_lengths.push_back(800);
        pacing_cycle_lengths.push_back(700);

        unsigned number_of_pulses = 10;


        // Create and archive stimulus
        {
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            AbstractStimulusFunction* const p_stimulus = new SteadyStateRestitutionStimulus(magnitude, duration_of_stimulus, startTime, pacing_cycle_lengths, number_of_pulses);

            // Should always archive a pointer
            output_arch << p_stimulus;

            delete p_stimulus;
        }

        // Restore
        {
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Create a pointer
            AbstractStimulusFunction* p_loaded_stim;
            input_arch >> p_loaded_stim;

            DoesSteadyStateStimulusPerformCorrectly(static_cast<SteadyStateRestitutionStimulus*>(p_loaded_stim),
                                                    magnitude, duration_of_stimulus, startTime,
                                                    pacing_cycle_lengths, number_of_pulses);

            delete p_loaded_stim;
        }
    }
};

#endif //_TESTRESTITUTIONCELLSTIMULI_HPP_
