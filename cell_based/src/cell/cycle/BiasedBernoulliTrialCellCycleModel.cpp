/*

Copyright (c) 2005-2023, University of Oxford.
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

#include "BiasedBernoulliTrialCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

BiasedBernoulliTrialCellCycleModel::BiasedBernoulliTrialCellCycleModel()
    : AbstractCellCycleModel(),
      mMaxDivisionProbability(0.1),
      mMinimumDivisionAge(1.0)
{
}

BiasedBernoulliTrialCellCycleModel::BiasedBernoulliTrialCellCycleModel(const BiasedBernoulliTrialCellCycleModel& rModel)
   : AbstractCellCycleModel(rModel),
     mMaxDivisionProbability(rModel.mMaxDivisionProbability),
     mMinimumDivisionAge(rModel.mMinimumDivisionAge)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

bool BiasedBernoulliTrialCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        if (GetAge() > mMinimumDivisionAge)
        {
            double dt = SimulationTime::Instance()->GetTimeStep();
            double p = mpCell->GetCellData()->GetItem("bias");
            assert(p >= 0.0 && p <= 1.0);
            if (!(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
            {
                RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                if (p_gen->ranf() < p*mMaxDivisionProbability*dt)
                {
                    mReadyToDivide = true;
                }
            }
        }
    }
    return mReadyToDivide;
}

AbstractCellCycleModel* BiasedBernoulliTrialCellCycleModel::CreateCellCycleModel()
{
    return new BiasedBernoulliTrialCellCycleModel(*this);
}

void BiasedBernoulliTrialCellCycleModel::SetMaxDivisionProbability(double maxDivisionProbability)
{
    assert(maxDivisionProbability >= 0.0 && maxDivisionProbability <= 1.0);
    mMaxDivisionProbability = maxDivisionProbability;
}

double BiasedBernoulliTrialCellCycleModel::GetMaxDivisionProbability()
{
    return mMaxDivisionProbability;
}

void BiasedBernoulliTrialCellCycleModel::SetMinimumDivisionAge(double minimumDivisionAge)
{
    assert(minimumDivisionAge >= 0.0);
    mMinimumDivisionAge = minimumDivisionAge;
}

double BiasedBernoulliTrialCellCycleModel::GetMinimumDivisionAge()
{
    return mMinimumDivisionAge;
}

double BiasedBernoulliTrialCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 1.0/mMaxDivisionProbability;
}

double BiasedBernoulliTrialCellCycleModel::GetAverageStemCellCycleTime()
{
    return 1.0/mMaxDivisionProbability;
}

void BiasedBernoulliTrialCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MaxDivisionProbability>" << mMaxDivisionProbability << "</MaxDivisionProbability>\n";
    *rParamsFile << "\t\t\t<MinimumDivisionAge>" << mMinimumDivisionAge << "</MinimumDivisionAge>\n";

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(BiasedBernoulliTrialCellCycleModel)
