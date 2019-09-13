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

#include "SteadyStateRunner.hpp"
#include "ZeroStimulus.hpp"

#ifdef CHASTE_CVODE

void SteadyStateRunner::RunToSteadyStateImplementation()
{
    // Get necessary things from stimulus current
    boost::shared_ptr<RegularStimulus> p_reg_stim = boost::static_pointer_cast<RegularStimulus>(mpModel->GetStimulusFunction());
    boost::shared_ptr<ZeroStimulus> p_zero_stim(new ZeroStimulus);
    const double pacing_cycle_length = p_reg_stim->GetPeriod(); //ms
    double stimulus_duration = p_reg_stim->GetDuration(); // ms
    double stimulus_start_time = p_reg_stim->GetStartTime(); // ms
    double stimulus_end_time = stimulus_start_time + stimulus_duration;
    double maximum_time_step = pacing_cycle_length; // ms

    bool force_reset_setting = mpModel->GetForceReset();
    bool minimal_reset_setting = mpModel->GetMinimalReset();

    mpModel->SetMaxSteps(1e5); // Per pace.
    mpModel->SetForceReset(true); // Best way to deal with discontinuities in the RHS according to CVODE

    // Set up vectors to monitor progress
    std::vector<double> old_state_vars;
    std::vector<double> new_state_vars;
    CopyToStdVector(mpModel->rGetStateVariables(), old_state_vars);

    const unsigned num_paces_to_analyse = (mTwoPaceScan) ? 2u : 1u;
    for (unsigned i = 0; i < mMaxNumPaces; i = i + num_paces_to_analyse)
    {
        // Look at two paces in the hope of detecting alternans and still saying we're steady-ish.
        for (unsigned j = 0; j < num_paces_to_analyse; j++)
        {
            // Pre-stimulus (can skip if stimulus is applied at t=0)
            if (stimulus_start_time > 0)
            {
                mpModel->SetMinimalReset(true); // We are just carrying on from last stop (not on stimulus) so don't need a solver reset.
                mpModel->SetStimulusFunction(p_zero_stim);
                mpModel->Solve((pacing_cycle_length) * (double)(i + j), (pacing_cycle_length) * (double)(i + j) + stimulus_start_time, maximum_time_step);
                mpModel->SetForceReset(true); // Reset every solve call after this (next two) to deal with discontinuities in stimulus.
            }

            mpModel->SetStimulusFunction(p_reg_stim); // RegularStimulus applies at stimulus_start_time<=t<=stimulus_end_time (includes bounds).
            mpModel->Solve((pacing_cycle_length) * (double)(i + j) + stimulus_start_time, (pacing_cycle_length) * (double)(i + j) + stimulus_end_time, maximum_time_step);

            // Post-stimulus
            mpModel->SetStimulusFunction(p_zero_stim); // added this because the RegularStimulus behaviour is to apply the stimulus at t=stimulus_end_time, which might confuse matters.
            mpModel->Solve((pacing_cycle_length) * (double)(i + j) + stimulus_end_time, (pacing_cycle_length) * (double)(i + j + 1), maximum_time_step);
        }

        this->mNumEvaluations += num_paces_to_analyse;
        CopyToStdVector(mpModel->rGetStateVariables(), new_state_vars);

        // Calculate the change in the norm of the state variables
        double temp = 0;
        for (unsigned j = 0; j < old_state_vars.size(); j++)
        {
            temp += fabs(new_state_vars[j] - old_state_vars[j]);
        }

        if (temp < 1e-6)
        {
            break; // say we are converged enough to steady state to stop here.
        }

        old_state_vars = new_state_vars;
    }

    // Reset stimulus to normal
    mpModel->SetStimulusFunction(p_reg_stim);
    mpModel->SetForceReset(force_reset_setting);
    mpModel->SetMinimalReset(minimal_reset_setting);
}

#endif // CHASTE_CVODE
