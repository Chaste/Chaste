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

#include "S1S2Stimulus.hpp"
#include "Exception.hpp"

S1S2Stimulus::S1S2Stimulus(double magnitude, double stimulusDuration, double s1Duration, double s1Period, double startTime, std::vector<double> s2Periods)
    : MultiStimulus()
{
    mNumS2FrequencyValues = s2Periods.size();
    for (unsigned i=0; i<mNumS2FrequencyValues; i++)
    {
        boost::shared_ptr<MultiStimulus> p_experiment_stimulus(new MultiStimulus());

        // The 0.5*(s1Period-stimulusDuration) taken off the stop time is to avoid getting a double stimulus at the instant s1Duration+startTime (#2442).
        boost::shared_ptr<RegularStimulus> p_s1(new RegularStimulus(magnitude, stimulusDuration, s1Period, startTime, s1Duration+startTime-0.5*(s1Period-stimulusDuration)));
        boost::shared_ptr<RegularStimulus> p_s2(new RegularStimulus(magnitude, stimulusDuration, s2Periods[i], s1Duration+startTime, s1Duration + 2*s2Periods[i]+startTime));

        p_experiment_stimulus->AddStimulus(p_s1);
        p_experiment_stimulus->AddStimulus(p_s2);

        // We hijack the mStimulus of the MultiStimulus class to
        // record different stimuli for each S2 frequency.
        this->AddStimulus(p_experiment_stimulus);
    }
    mS2Index=0u;
}

double S1S2Stimulus::GetStimulus(double time)
{
    return this->mStimuli[mS2Index]->GetStimulus(time);
}

void S1S2Stimulus::SetS2ExperimentPeriodIndex(unsigned index)
{
    if (index < mNumS2FrequencyValues)
    {
        mS2Index = index;
    }
    else
    {
        EXCEPTION("There are fewer S2 frequency values than the one you have requested.");
    }
}

unsigned S1S2Stimulus::GetNumS2FrequencyValues()
{
    return mNumS2FrequencyValues;
}


#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(S1S2Stimulus)


