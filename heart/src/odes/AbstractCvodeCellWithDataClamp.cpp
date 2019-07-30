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

#ifdef CHASTE_CVODE

#include <sstream>
#include <iostream>
#include <cmath>

#include "AbstractCvodeCellWithDataClamp.hpp"
#include "Exception.hpp"
#include "HeartConfig.hpp"
#include "VectorHelperFunctions.hpp"
#include <algorithm>


AbstractCvodeCellWithDataClamp::AbstractCvodeCellWithDataClamp(boost::shared_ptr<AbstractIvpOdeSolver> pSolver/* unused */,
                                                               unsigned numberOfStateVariables,
                                                               unsigned voltageIndex,
                                                               boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCvodeCell(pSolver, numberOfStateVariables, voltageIndex, pIntracellularStimulus),
      mDataAvailable(false),
      mDataClampIsOn(false)
{
}

AbstractCvodeCellWithDataClamp::~AbstractCvodeCellWithDataClamp()
{
}

void AbstractCvodeCellWithDataClamp::SetExperimentalData(std::vector<double> experimentalTimes,
                                                         std::vector<double> experimentalVoltages)
{
    assert(experimentalVoltages.size()==experimentalTimes.size());
    assert(experimentalVoltages.size()>1u);
    mExperimentalTimes = experimentalTimes;
    mExperimentalVoltages = experimentalVoltages;
    mDataAvailable = true;
}

double AbstractCvodeCellWithDataClamp::GetExperimentalVoltageAtTimeT(const double& rTime) // Linearly interpolate between data points
{
    if (!mDataClampIsOn)
    {
        // For speed we want to be able to skip this method.
        return DOUBLE_UNSET;
    }

    if (rTime > mExperimentalTimes.back() || rTime < mExperimentalTimes[0])
    {
        EXCEPTION("The requested time " << rTime << " is outside the times stored in the data clamp.");
    }

    std::vector<double>::iterator upper_vec =
            std::upper_bound(mExperimentalTimes.begin(), mExperimentalTimes.end(), rTime);
    // upper_vec points to first element that is (strictly) > rTime (so it can't be first one if rTime is in range)

    // Special case - here the lower bound is equal to upper bound and both are pointing to end!
    if (upper_vec == mExperimentalTimes.end())
    {
        return mExperimentalVoltages.back();
    }

    unsigned upper_index = upper_vec - mExperimentalTimes.begin();

    assert(upper_index>0u); // Sanity check to avoid seg faults below.
    // if this assert ever fails we may need the code below, but I can't see how it could!
//    if (upper_index == 0u)
//    {
//        // Special case - here the lower bound is equal to upper bound and no interpolation needed
//        return mExperimentalVoltages[0];
//    }

    double lower_time = mExperimentalTimes[upper_index-1];
    double increment = rTime - lower_time;
    double time_step = mExperimentalTimes[upper_index] - mExperimentalTimes[upper_index-1];
    double voltage_step = mExperimentalVoltages[upper_index] - mExperimentalVoltages[upper_index-1];
    return mExperimentalVoltages[upper_index-1] + increment * voltage_step / time_step;
}

void AbstractCvodeCellWithDataClamp::TurnOffDataClamp()
{
    this->SetParameter("membrane_data_clamp_current_conductance",0);
    mDataClampIsOn = false;
}

void AbstractCvodeCellWithDataClamp::TurnOnDataClamp(double conductance)
{
    if (!mDataAvailable)
    {
        EXCEPTION("Before calling TurnOnDataClamp(), please provide experimental data via the SetExperimentalData() method.");
    }
    assert(mExperimentalTimes.size() > 1u);
    this->SetParameter("membrane_data_clamp_current_conductance", conductance);
    mDataClampIsOn = true;
}

#endif // CHASTE_CVODEWITHDATACLAMP
