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

#include "SteadyStateRunner.hpp"

#ifdef CHASTE_CVODE

void SteadyStateRunner::RunToSteadyStateImplementation()
{
    // Get necessary things from stimulus current
    boost::shared_ptr<RegularStimulus> p_reg_stim =
                         boost::static_pointer_cast<RegularStimulus>(mpModel->GetStimulusFunction());
    const double pacing_cycle_length = p_reg_stim->GetPeriod(); //ms
    double maximum_time_step = p_reg_stim->GetDuration();  // ms

    // Set up vectors to monitor progress
    std::vector<double> old_state_vars;
    std::vector<double> new_state_vars;
    CopyToStdVector(mpModel->rGetStateVariables(),old_state_vars);

    mpModel->SetMaxSteps(1e5); // Per pace.

    /*
     * This was an idea about solving roughly to start with and then refining later,
     * but it overrides anything set by the user and is therefore a bit dangerous.
     * It also didn't seem to make very much difference...
     */
    //bool rough_solution = true;
    //mpModel->SetTolerances(1e-4,1e-6);

    // Use optimisations from #1912 - since we don't interfere with the model between Solve() calls don't reset.
    mpModel->SetMinimalReset(true);

    const unsigned num_paces_to_analyse = (mTwoPaceScan) ? 2u : 1u;

    for (unsigned i=0; i<mMaxNumPaces; i=i+num_paces_to_analyse)
    {
        // Look at two paces in the hope of detecting alternans and still saying we're steady-ish.
        mpModel->Solve((pacing_cycle_length)*(double)(i), (pacing_cycle_length)*(double)(i+num_paces_to_analyse), maximum_time_step);
        this->mNumEvaluations += num_paces_to_analyse;
        CopyToStdVector(mpModel->rGetStateVariables(),new_state_vars);

        // Calculate the change in the norm of the state variables
        double temp = 0;
        for (unsigned j=0; j<old_state_vars.size(); j++)
        {
            temp += fabs(new_state_vars[j] - old_state_vars[j]);
        }

//        if (rough_solution && temp < 1e-4)
//        {
//            mpModel->SetTolerances(1e-6,1e-8);
//            rough_solution = false;
//        }

        if (temp < 1e-6)
        {
            break;
        }

        old_state_vars = new_state_vars;
    }

    // Since we don't know what we are going to do with it next re-set it on every solve call.
    mpModel->SetMinimalReset(false);
}

#endif // CHASTE_CVODE

