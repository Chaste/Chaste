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

#include "AbstractSteadyStateRunner.hpp"

#ifdef CHASTE_CVODE

AbstractSteadyStateRunner::AbstractSteadyStateRunner(boost::shared_ptr<AbstractCvodeCell> pModel)
        : mpModel(pModel),
          mNumEvaluations(0u),
          mMaxNumPaces(10000u),
          mSuppressOutput(false)
{
}

bool AbstractSteadyStateRunner::RunToSteadyState()
{
    // Get timing information from the cell stimulus (only works for regular stimuli)
    if (!boost::dynamic_pointer_cast<RegularStimulus>(mpModel->GetStimulusFunction()))
    {
        EXCEPTION("Steady State approximations only work for models with RegularStimulus objects.");
    }

    RunToSteadyStateImplementation();

    if (mNumEvaluations < mMaxNumPaces)
    {
        if (!mSuppressOutput)
        {
            std::cout << "Steady state detected after " << GetNumEvaluations() << " paces." << std::endl;
        }
        return true;
    }
    else
    {
        if (!mSuppressOutput)
        {
            std::cout << "Steady state not reached after " << GetNumEvaluations() << " paces." << std::endl;
        }
        WARNING("Model " << mpModel->GetSystemName() << " did not reach steady state within " << mMaxNumPaces << " paces.");
        return false;
    }
}

void AbstractSteadyStateRunner::SuppressOutput(bool suppress)
{
    mSuppressOutput = suppress;
}

unsigned AbstractSteadyStateRunner::GetNumEvaluations()
{
    return mNumEvaluations;
}

void AbstractSteadyStateRunner::SetMaxNumPaces(unsigned numPaces)
{
    if (numPaces == 0u)
    {
        EXCEPTION("Please set a maximum number of paces that is positive");
    }
    mMaxNumPaces = numPaces;
}

#endif // CHASTE_CVODE
