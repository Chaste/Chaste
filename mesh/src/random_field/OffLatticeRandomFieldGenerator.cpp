/*

Copyright (c) 2005-2025, University of Oxford.
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

#include "OffLatticeRandomFieldGenerator.hpp"

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <memory>
#include <numeric>

#include "Exception.hpp"
#include "RandomNumberGenerator.hpp"
#include "RandomFieldUtilityFunctions.hpp"

template <unsigned SPACE_DIM>
OffLatticeRandomFieldGenerator<SPACE_DIM>::OffLatticeRandomFieldGenerator(
    std::array<double, SPACE_DIM> lowerCorner,
    std::array<double, SPACE_DIM> upperCorner,
    std::array<bool, SPACE_DIM> periodicity,
    double lengthScale,
    double boxWidth)
    : mLowerCorner(lowerCorner),
      mUpperCorner(upperCorner),
      mPeriodicity(periodicity),
      mLengthScale(lengthScale),
      mpBoxCollection(nullptr)
{
    // Reset the box width if the default value is being used
    if (boxWidth == DOUBLE_UNSET)
    {
        const double max_box_width = 4.0 * lengthScale;
        boxWidth = std::min(max_box_width, upperCorner[0] - lowerCorner[0]);
    }

    // Set up the box collection
    c_vector<double, 2 * SPACE_DIM> domain_size;
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        domain_size[2 * dim] = lowerCorner[dim];
        domain_size[2 * dim + 1] = upperCorner[dim];
    }

    ///\todo: Obsolete box collection should take a std::array<bool, DIM>
    bool periodic_x = periodicity[0];
    bool periodic_y = SPACE_DIM > 1 ? periodicity[1] : false;
    bool periodic_z = SPACE_DIM > 2 ? periodicity[2] : false;

    mpBoxCollection = std::make_unique<ObsoleteBoxCollection<SPACE_DIM>>(
        boxWidth,
        domain_size,
        periodic_x,
        periodic_y,
        periodic_z
    );

    mpBoxCollection->SetupLocalBoxesHalfOnly();
    mOpenSimplex = OpenSimplex2S();
}

template <unsigned SPACE_DIM>
void OffLatticeRandomFieldGenerator<SPACE_DIM>::SetRandomSeed(
    const unsigned seed)
{
    mOpenSimplex = OpenSimplex2S(seed);
}

template <unsigned SPACE_DIM>
std::vector<double> OffLatticeRandomFieldGenerator<SPACE_DIM>::SampleRandomField(
    const std::vector<Node<SPACE_DIM>*>& rNodes)
{
    return this->SampleRandomFieldAtTime(rNodes, 100.0 * mLengthScale * RandomNumberGenerator::Instance()->ranf());
}

template <unsigned SPACE_DIM>
std::vector<double> OffLatticeRandomFieldGenerator<SPACE_DIM>::SampleRandomFieldAtTime(
    const std::vector<Node<SPACE_DIM>*>& rNodes,
    const double time)
{
    std::vector<double> samples(rNodes.size());
    for (unsigned i = 0; i < samples.size(); ++i)
    {
        const c_vector<double, SPACE_DIM> node_location = rNodes[i]->rGetLocation();
        switch (SPACE_DIM)
        {
            case 1:
            {
                samples[i] = random_field::Reshape(mOpenSimplex.noise2_XBeforeY(node_location[0] * mLengthScale, time + 0.5));
                break;
            }
            case 2:
            {
                samples[i] = random_field::Reshape(mOpenSimplex.noise3_XYBeforeZ(node_location[0] * mLengthScale, node_location[1] * mLengthScale, time));
                break;
            }
            case 3:
            {
                samples[i] = random_field::Reshape(mOpenSimplex.noise4_XYBeforeZW(node_location[0] * mLengthScale, node_location[1] * mLengthScale, node_location[2] * mLengthScale, time));
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }
    }

    return samples;
}
// Explicit instantiation
template class OffLatticeRandomFieldGenerator<1>;
template class OffLatticeRandomFieldGenerator<2>;
template class OffLatticeRandomFieldGenerator<3>;
