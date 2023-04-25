/*

Copyright (c) 2005-2018, University of Oxford.
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

#include "UniformGridRandomFieldGenerator.hpp"

#include <iomanip>
#include <fstream>
#include <numeric>

#include "Exception.hpp"
#include "FileFinder.hpp"
#include "Node.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

#include "Debug.hpp"

template <unsigned SPACE_DIM>
UniformGridRandomFieldGenerator<SPACE_DIM>::UniformGridRandomFieldGenerator(std::array<double, SPACE_DIM> lowerCorner,
                                                                            std::array<double, SPACE_DIM> upperCorner,
                                                                            std::array<unsigned, SPACE_DIM> numGridPts,
                                                                            std::array<bool, SPACE_DIM> periodicity,
                                                                            double lengthScale)
        : mLowerCorner(lowerCorner),
          mUpperCorner(upperCorner),
          mNumGridPts(numGridPts),
          mPeriodicity(periodicity),
          mLengthScale(lengthScale),
          mCacheDir("CachedRandomFields/")
{
    // Calculate the total number of grid points
    mNumTotalGridPts = std::accumulate(mNumGridPts.begin(), mNumGridPts.end(), 1u, std::multiplies<unsigned>());

    // Check parameters are sensible
    {
        for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
        {
            assert(mLowerCorner[dim] < mUpperCorner[dim]);
            assert(mNumGridPts[dim] > 0);
        }

        assert(mLengthScale > 0.0);
    }

    // Calculate the grid spacings
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        mGridSpacing[dim] = (mUpperCorner[dim] - mLowerCorner[dim]) / mNumGridPts[dim];
        mOneOverGridSpacing[dim] = mNumGridPts[dim] / (mUpperCorner[dim] - mLowerCorner[dim]);
    }
}

template <unsigned SPACE_DIM>
double UniformGridRandomFieldGenerator<SPACE_DIM>::GetSquaredDistAtoB(
    const c_vector<double, SPACE_DIM>& rLocation1,
    const c_vector<double, SPACE_DIM>& rLocation2) const noexcept
{
    double dist_squared = 0.0;

    // Loop over dimensions
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        // Get the (non-periodic) absolute difference
        double delta_this_dim = std::fabs(rLocation2[dim] - rLocation1[dim]);

        // If this dimension is periodic, recalculate the distance if necessary
        if (mPeriodicity[dim])
        {
            const double domain_width_this_dim = mUpperCorner[dim] - mLowerCorner[dim];

            delta_this_dim = std::min(delta_this_dim, domain_width_this_dim - delta_this_dim);
        }

        dist_squared += delta_this_dim * delta_this_dim;
    }

    return dist_squared;
}

template <unsigned SPACE_DIM>
std::string UniformGridRandomFieldGenerator<SPACE_DIM>::GetFilenameFromParams() const noexcept
{
    std::stringstream name;
    name << std::fixed << std::setprecision(3);

    switch(SPACE_DIM)
    {
        case 1:
        {
            name << "x_"
                 << mLowerCorner[0] << "_"
                 << mUpperCorner[0] << "_"
                 << mNumGridPts[0] << "_"
                 << mPeriodicity[0] << "_"
                 << mLengthScale;
            break;
        }
        case 2:
        {
            name << "xy_"
                 << mLowerCorner[0] << "_" << mLowerCorner[1] << "_"
                 << mUpperCorner[0] << "_" << mUpperCorner[1] << "_"
                 << mNumGridPts[0] << "_" << mNumGridPts[1] << "_"
                 << mPeriodicity[0] << "_" << mPeriodicity[1] << "_"
                 << mLengthScale;
            break;
        }
        case 3:
        {
            name << "xyz_"
                 << mLowerCorner[0] << "_" << mLowerCorner[1] << "_" << mLowerCorner[2] << "_"
                 << mUpperCorner[0] << "_" << mUpperCorner[1] << "_" << mUpperCorner[2] << "_"
                 << mNumGridPts[0] << "_" << mNumGridPts[1] << "_" << mNumGridPts[2] << "_"
                 << mPeriodicity[0] << "_" << mPeriodicity[1] << "_" << mPeriodicity[2] << "_"
                 << mLengthScale;
            break;
        }
        default:
            NEVER_REACHED;
    }

    // Add file extension
    name << ".rfg";
    return name.str();
}

template <unsigned SPACE_DIM>
std::vector<double> UniformGridRandomFieldGenerator<SPACE_DIM>::SampleRandomField() const noexcept
{ 
    return this->SampleRandomFieldAtTime(rand());
}

template <unsigned SPACE_DIM>
std::vector<double> UniformGridRandomFieldGenerator<SPACE_DIM>::SampleRandomFieldAtTime(double time) const noexcept
{ 
    // TODO: randomise seed
    OpenSimplex2S os(123);
    auto reshape = [](const double val) {
        double distFromHalf = 2.0 * (0.5 - std::abs(0.5 - std::abs(val)));
        double strength = (1.0 - std::abs(val)) + distFromHalf * 0.2;
        return val * (1.0 - 0.18 * strength);
    };

    std::vector<double> samples(mNumTotalGridPts);
    
    switch(SPACE_DIM) {
        case 1:
            for (unsigned int x = 0; x < mNumGridPts[0]; x++) {
                samples[x] = reshape(reshape(os.noise2_XBeforeY(x * mLengthScale, time)));
            }
            break;
        
        case 2:
            for (unsigned int x = 0; x < mNumGridPts[0]; x++) {
                for (unsigned int y = 0; y < mNumGridPts[1]; y++) {
                    samples[mNumGridPts[1] * y + x] = reshape(reshape(os.noise3_XYBeforeZ(x * mLengthScale, y * mLengthScale, time)));
                }
            }
            break;

        case 3:
            for (unsigned int x = 0; x < mNumGridPts[0]; x++) {
                for (unsigned int y = 0; y < mNumGridPts[1]; y++) {
                    for (unsigned int z = 0; z < mNumGridPts[2]; z++) {
                        samples[mNumGridPts[2] * z * y + mNumGridPts[1] * y + x] = reshape(reshape(os.noise4_XYBeforeZW(x * mLengthScale, y * mLengthScale, z * mLengthScale, time)));
                    }
                }
            }
            break;
            
        default:
            NEVER_REACHED;
            break;

    }
    
    return samples;
}
template <unsigned SPACE_DIM>
double UniformGridRandomFieldGenerator<SPACE_DIM>::Interpolate(const std::vector<double>& rRandomField,
                                                               const c_vector<double, SPACE_DIM>& rLocation) const noexcept
{
    assert(mNumTotalGridPts == rRandomField.size());

    // Find the nearest node
    std::array<long, SPACE_DIM> lower_left;

    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        lower_left[dim] = std::floor((rLocation[dim] - mLowerCorner[dim]) / mGridSpacing[dim]);

        if (lower_left[dim] < 0)
        {
            lower_left[dim] = 0;
            WARN_ONCE_ONLY("Interpolating outside random field grid: does the random field need to be larger?");
        }
        else if (lower_left[dim] >= mNumGridPts[dim])
        {
            lower_left[dim] = mNumGridPts[dim] - 1;
            WARN_ONCE_ONLY("Interpolating outside random field grid: does the random field need to be larger?");
        }
    }

    double interpolated_value{};

    // Perform the interpolation
    switch(SPACE_DIM)
    {
        case 1:
        {
            // The value of the field at points either side of rLocation
            std::array<long, SPACE_DIM> upper_idx = lower_left;
            upper_idx[0] =  (upper_idx[0] + 1) % mNumGridPts[0];

            const double field_at_lower = rRandomField[GetLinearIndex(lower_left)];
            const double field_at_upper = rRandomField[GetLinearIndex(upper_idx)];

            // Perform a simple linear interpolation
            const std::array<double, SPACE_DIM> lower_location = GetPositionUsingGridIndex(lower_left);
            const double interpolant = (rLocation[0] - lower_location[0]) * mOneOverGridSpacing[0];

            interpolated_value = field_at_lower * (1.0 - interpolant) + field_at_upper * interpolant;
            break;
        }
        case 2:
        {
            // The value of the field at the four points of the square containing rLocation
            std::array<long, SPACE_DIM> x_upper = lower_left;
            x_upper[0] = (x_upper[0] + 1) % mNumGridPts[0];

            std::array<long, SPACE_DIM> y_upper = lower_left;
            y_upper[1] = (y_upper[1] + 1) % mNumGridPts[1];

            std::array<long, SPACE_DIM> xy_upper = x_upper;
            xy_upper[1] = y_upper[1];

            const double field_xy_lower = rRandomField[GetLinearIndex(lower_left)];
            const double field_x_upper = rRandomField[GetLinearIndex(x_upper)];
            const double field_y_upper = rRandomField[GetLinearIndex(y_upper)];
            const double field_xy_upper = rRandomField[GetLinearIndex(xy_upper)];

            // Perform a simple bilinear interpolation
            const std::array<double, SPACE_DIM> lower_location = GetPositionUsingGridIndex(lower_left);
            const double dist_x_lower = rLocation[0] - lower_location[0];
            const double dist_x_upper = mGridSpacing[0] - dist_x_lower;
            const double dist_y_lower = rLocation[1] - lower_location[1];
            const double dist_y_upper = mGridSpacing[1] - dist_y_lower;

            interpolated_value = mOneOverGridSpacing[0] * mOneOverGridSpacing[1] * (field_xy_lower * dist_x_upper * dist_y_upper +
                                                                                    field_x_upper * dist_x_lower * dist_y_upper +
                                                                                    field_y_upper * dist_x_upper * dist_y_lower +
                                                                                    field_xy_upper * dist_x_lower * dist_y_lower);

            break;
        }
        case 3:
        {
            // \todo: perform a suitable interpolation rather than just near neighbour
            interpolated_value = rRandomField[GetLinearIndex(lower_left)];

            break;
        }
        default:
            NEVER_REACHED;
    }

    return interpolated_value;
}

template <unsigned SPACE_DIM>
long UniformGridRandomFieldGenerator<SPACE_DIM>::GetLinearIndex(std::array<long, SPACE_DIM> gridIndex) const noexcept
{
    long linear_index;

    // \todo: double check I'm not calculating the transpose of each point
    switch(SPACE_DIM)
    {
        case 1:
        {
            linear_index = gridIndex[0];
            break;
        }
        case 2:
        {
            linear_index = gridIndex[0] +
                           gridIndex[1] * mNumGridPts[0];
            break;
        }
        case 3:
        {
            linear_index = gridIndex[0] +
                           gridIndex[1] * mNumGridPts[0] +
                           gridIndex[2] * mNumGridPts[0] * mNumGridPts[1];
            break;
        }
        default:
            NEVER_REACHED;
    }

    return linear_index;
}

template <unsigned SPACE_DIM>
std::array<double, SPACE_DIM> UniformGridRandomFieldGenerator<SPACE_DIM>::GetPositionUsingGridIndex(std::array<long, SPACE_DIM> gridIndex) const noexcept
{
    std::array<double, SPACE_DIM> position = {{}};

    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        position[dim] = mLowerCorner[dim] + mGridSpacing[dim] * gridIndex[dim];
    }

    return position;
}

// Explicit instantiation
template class UniformGridRandomFieldGenerator<1>;
template class UniformGridRandomFieldGenerator<2>;
template class UniformGridRandomFieldGenerator<3>;
