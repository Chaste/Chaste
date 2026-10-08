/*

Copyright (c) 2005-2026, University of Oxford.
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

#include "SemSpatiallyCorrelatedRandomForce.hpp"

#include "Warnings.hpp"

template<unsigned DIM>
SemSpatiallyCorrelatedRandomForce<DIM>::SemSpatiallyCorrelatedRandomForce()
    : AbstractSemRandomForce<DIM>(),
      mCorrelationLength(1.0),
      mSmallCorrelationLengthWarningThreshold(1e-12),
      mRandomSeed(0u),
      mHasRandomSeed(false),
      mGeneratorsNeedRefresh(true)
{
    for (unsigned dim = 0; dim < DIM; ++dim)
    {
        mLowerCorner[dim] = 0.0;
        mUpperCorner[dim] = 1.0;
        mPeriodicity[dim] = false;
    }
}

template<unsigned DIM>
void SemSpatiallyCorrelatedRandomForce<DIM>::RefreshRandomFieldGenerators()
{
    assert(mCorrelationLength > 0.0);

    const double generator_length_scale = 1.0 / mCorrelationLength;

    for (unsigned dim = 0; dim < DIM; ++dim)
    {
        mpRandomFieldGenerators[dim] = std::make_unique<OffLatticeRandomFieldGenerator<DIM> >(
            mLowerCorner,
            mUpperCorner,
            mPeriodicity,
            generator_length_scale);

        if (mHasRandomSeed)
        {
            mpRandomFieldGenerators[dim]->SetRandomSeed(mRandomSeed + dim);
        }
    }

    mGeneratorsNeedRefresh = false;
}

template<unsigned DIM>
std::vector<c_vector<double, DIM> > SemSpatiallyCorrelatedRandomForce<DIM>::GenerateUnitNoise(const std::vector<Node<DIM>*>& rNodes)
{
    if (mCorrelationLength < mSmallCorrelationLengthWarningThreshold)
    {
        WARN_ONCE_ONLY("SemSpatiallyCorrelatedRandomForce correlation length is very small; the random field may be poorly resolved.");
    }

    if (mGeneratorsNeedRefresh)
    {
        RefreshRandomFieldGenerators();
    }

    std::vector<c_vector<double, DIM> > noise(rNodes.size());
    for (unsigned dim = 0; dim < DIM; ++dim)
    {
        const std::vector<double> field_samples = mpRandomFieldGenerators[dim]->SampleRandomField(rNodes);
        assert(field_samples.size() == rNodes.size());

        for (unsigned node_index = 0; node_index < rNodes.size(); ++node_index)
        {
            noise[node_index][dim] = field_samples[node_index];
        }
    }

    return noise;
}

template<unsigned DIM>
double SemSpatiallyCorrelatedRandomForce<DIM>::GetCorrelationLength() const
{
    return mCorrelationLength;
}

template<unsigned DIM>
void SemSpatiallyCorrelatedRandomForce<DIM>::SetCorrelationLength(double correlationLength)
{
    assert(correlationLength > 0.0);
    mCorrelationLength = correlationLength;
    mGeneratorsNeedRefresh = true;
}

template<unsigned DIM>
double SemSpatiallyCorrelatedRandomForce<DIM>::GetSmallCorrelationLengthWarningThreshold() const
{
    return mSmallCorrelationLengthWarningThreshold;
}

template<unsigned DIM>
void SemSpatiallyCorrelatedRandomForce<DIM>::SetSmallCorrelationLengthWarningThreshold(double threshold)
{
    assert(threshold >= 0.0);
    mSmallCorrelationLengthWarningThreshold = threshold;
}

template<unsigned DIM>
std::array<double, DIM> SemSpatiallyCorrelatedRandomForce<DIM>::GetLowerCorner() const
{
    return mLowerCorner;
}

template<unsigned DIM>
void SemSpatiallyCorrelatedRandomForce<DIM>::SetLowerCorner(std::array<double, DIM> lowerCorner)
{
    mLowerCorner = lowerCorner;
    mGeneratorsNeedRefresh = true;
}

template<unsigned DIM>
std::array<double, DIM> SemSpatiallyCorrelatedRandomForce<DIM>::GetUpperCorner() const
{
    return mUpperCorner;
}

template<unsigned DIM>
void SemSpatiallyCorrelatedRandomForce<DIM>::SetUpperCorner(std::array<double, DIM> upperCorner)
{
    mUpperCorner = upperCorner;
    mGeneratorsNeedRefresh = true;
}

template<unsigned DIM>
std::array<bool, DIM> SemSpatiallyCorrelatedRandomForce<DIM>::GetPeriodicity() const
{
    return mPeriodicity;
}

template<unsigned DIM>
void SemSpatiallyCorrelatedRandomForce<DIM>::SetPeriodicity(std::array<bool, DIM> periodicity)
{
    mPeriodicity = periodicity;
    mGeneratorsNeedRefresh = true;
}

template<unsigned DIM>
void SemSpatiallyCorrelatedRandomForce<DIM>::SetRandomSeed(unsigned randomSeed)
{
    mRandomSeed = randomSeed;
    mHasRandomSeed = true;
    mGeneratorsNeedRefresh = true;
}

template<unsigned DIM>
void SemSpatiallyCorrelatedRandomForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CorrelationLength>" << mCorrelationLength << "</CorrelationLength>\n";
    *rParamsFile << "\t\t\t<SmallCorrelationLengthWarningThreshold>" << mSmallCorrelationLengthWarningThreshold << "</SmallCorrelationLengthWarningThreshold>\n";
    *rParamsFile << "\t\t\t<HasRandomSeed>" << mHasRandomSeed << "</HasRandomSeed>\n";
    if (mHasRandomSeed)
    {
        *rParamsFile << "\t\t\t<RandomSeed>" << mRandomSeed << "</RandomSeed>\n";
    }

    AbstractSemRandomForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class SemSpatiallyCorrelatedRandomForce<1>;
template class SemSpatiallyCorrelatedRandomForce<2>;
template class SemSpatiallyCorrelatedRandomForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemSpatiallyCorrelatedRandomForce)
