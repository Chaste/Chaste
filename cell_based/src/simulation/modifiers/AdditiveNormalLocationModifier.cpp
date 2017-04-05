/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "AdditiveNormalLocationModifier.hpp"
#include "RandomNumberGenerator.hpp"

template <unsigned DIM>
AdditiveNormalLocationModifier<DIM>::AdditiveNormalLocationModifier()
        : AbstractCellBasedSimulationModifier<DIM>(),
          mMean(0.0),
          mStdDev(0.0)
{
}

template <unsigned DIM>
AdditiveNormalLocationModifier<DIM>::~AdditiveNormalLocationModifier()
{
}

template <unsigned DIM>
void AdditiveNormalLocationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    for (unsigned node_idx = 0; node_idx < rCellPopulation.GetNumNodes(); ++node_idx)
    {
        c_vector<double, DIM>& r_location = rCellPopulation.GetNode(node_idx)->rGetModifiableLocation();

        r_location[0] = fmod(r_location[0] + p_gen->NormalRandomDeviate(mMean, mStdDev) + 1.0, 1.0);
        r_location[1] = fmod(r_location[1] + p_gen->NormalRandomDeviate(mMean, mStdDev) + 1.0, 1.0);
    }
}

template <unsigned DIM>
void AdditiveNormalLocationModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory)
{
}

template <unsigned DIM>
void AdditiveNormalLocationModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<Mean>" << mMean << "</Mean>\n";
    *rParamsFile << "\t\t\t<StdDev>" << mStdDev << "</StdDev>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

template <unsigned DIM>
double AdditiveNormalLocationModifier<DIM>::GetMean() const
{
    return mMean;
}

template <unsigned DIM>
void AdditiveNormalLocationModifier<DIM>::SetMean(double mean)
{
    mMean = mean;
}

template <unsigned DIM>
double AdditiveNormalLocationModifier<DIM>::GetStdDev() const
{
    return mStdDev;
}

template <unsigned DIM>
void AdditiveNormalLocationModifier<DIM>::SetStdDev(double stdDev)
{
    mStdDev = stdDev;
}

// Explicit instantiation
template class AdditiveNormalLocationModifier<1>;
template class AdditiveNormalLocationModifier<2>;
template class AdditiveNormalLocationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AdditiveNormalLocationModifier)
