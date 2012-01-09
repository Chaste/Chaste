/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/
#include "S1S2Stimulus.hpp"

S1S2Stimulus::S1S2Stimulus(double magnitude, double stimulusDuration, double s1Duration, double s1Period, double startTime, std::vector<double> s2Periods)
    : MultiStimulus()
{
    mNumS2FrequencyValues = s2Periods.size();
    for (unsigned i=0; i<mNumS2FrequencyValues; i++)
    {
        boost::shared_ptr<MultiStimulus> p_experiment_stimulus(new MultiStimulus());

        boost::shared_ptr<RegularStimulus> p_s1(new RegularStimulus(magnitude, stimulusDuration, s1Period, startTime, s1Duration+startTime));
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
CHASTE_CLASS_EXPORT(S1S2Stimulus);


