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
#include "RandomFieldUtilityFunctions.hpp"

template <unsigned SPACE_DIM>
UniformGridRandomFieldGenerator<SPACE_DIM>::UniformGridRandomFieldGenerator(
    std::array<double, SPACE_DIM> lowerCorner,
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
    mNumTotalGridPts = std::accumulate(mNumGridPts.begin(),
                                       mNumGridPts.end(),
                                       1u,
                                       std::multiplies<unsigned>());

    // Check parameters are sensible
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        assert(mLowerCorner[dim] < mUpperCorner[dim]);
        assert(mNumGridPts[dim] > 0);
    }
    assert(mLengthScale > 0.0);

    // Calculate the grid spacings
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        mGridSpacing[dim] = (mUpperCorner[dim] - mLowerCorner[dim]) / mNumGridPts[dim];
        mOneOverGridSpacing[dim] = mNumGridPts[dim] / (mUpperCorner[dim] - mLowerCorner[dim]);
    }

    mOpenSimplex = OpenSimplex2S(0);
}

template <unsigned SPACE_DIM>
void UniformGridRandomFieldGenerator<SPACE_DIM>::SetRandomSeed(const unsigned seed)
{
    mOpenSimplex = OpenSimplex2S(seed);
}

template <unsigned SPACE_DIM>
std::vector<double> UniformGridRandomFieldGenerator<SPACE_DIM>::SampleRandomField()
{
    return this->SampleRandomFieldAtTime(rand());
}

template <unsigned SPACE_DIM>
std::vector<double> UniformGridRandomFieldGenerator<SPACE_DIM>::SampleRandomFieldAtTime(double time)
{
    std::vector<double> samples(mNumTotalGridPts);

    switch(SPACE_DIM)
    {
        case 1:
        {
            for (unsigned x = 0; x < mNumGridPts[0]; x++)
            {
                samples[x] = random_field::Reshape(mOpenSimplex.noise2_XBeforeY(x * mLengthScale, time));
            }
            break;
        }
        case 2:
        {
            for (unsigned x = 0; x < mNumGridPts[0]; x++)
            {
                for (unsigned y = 0; y < mNumGridPts[1]; y++)
                {
                    samples[mNumGridPts[1] * y + x] = random_field::Reshape(mOpenSimplex.noise3_XYBeforeZ(x * mLengthScale, y * mLengthScale, time));
                }
            }
            break;
        }
        case 3:
        {
            for (unsigned x = 0; x < mNumGridPts[0]; x++)
            {
                for (unsigned y = 0; y < mNumGridPts[1]; y++)
                {
                    for (unsigned z = 0; z < mNumGridPts[2]; z++)
                    {
                        samples[mNumGridPts[2] * z * y + mNumGridPts[1] * y + x] = random_field::Reshape(mOpenSimplex.noise4_XYBeforeZW(x * mLengthScale, y * mLengthScale, z * mLengthScale, time));
                    }
                }
            }
            break;
        }
        default:
            NEVER_REACHED;
            break;
    }

    return samples;
}

template <unsigned SPACE_DIM>
double UniformGridRandomFieldGenerator<SPACE_DIM>::Interpolate(const std::vector<double>& rRandomField,
                                                               const c_vector<double, SPACE_DIM>& rLocation) const
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
            // The value of the field at the eight points of the cube containing rLocation
            std::array<long, SPACE_DIM> x_upper = lower_left;
            x_upper[0] = (x_upper[0] + 1) % mNumGridPts[0];

            std::array<long, SPACE_DIM> y_upper = lower_left;
            y_upper[1] = (y_upper[1] + 1) % mNumGridPts[1];

            std::array<long, SPACE_DIM> z_upper = lower_left;
            z_upper[2] = (z_upper[2] + 1) % mNumGridPts[2];

            std::array<long, SPACE_DIM> xy_upper = x_upper;
            xy_upper[1] = y_upper[1];

            std::array<long, SPACE_DIM> xz_upper = x_upper;
            xz_upper[2] = z_upper[2];

            std::array<long, SPACE_DIM> yz_upper = y_upper;
            yz_upper[2] = z_upper[2];

            std::array<long, SPACE_DIM> xyz_upper = xy_upper;
            xyz_upper[2] = z_upper[2];

            const double field_xyz_lower = rRandomField[GetLinearIndex(lower_left)];
            const double field_x_upper = rRandomField[GetLinearIndex(x_upper)];
            const double field_y_upper = rRandomField[GetLinearIndex(y_upper)];
            const double field_z_upper = rRandomField[GetLinearIndex(z_upper)];
            const double field_xy_upper = rRandomField[GetLinearIndex(xy_upper)];
            const double field_xz_upper = rRandomField[GetLinearIndex(xz_upper)];
            const double field_yz_upper = rRandomField[GetLinearIndex(yz_upper)];
            const double field_xyz_upper = rRandomField[GetLinearIndex(xyz_upper)];

            // Perform a trilinear interpolation
            const std::array<double, SPACE_DIM> lower_location = GetPositionUsingGridIndex(lower_left);
            const double dist_x_lower = rLocation[0] - lower_location[0];
            const double dist_y_lower = rLocation[1] - lower_location[1];
            const double dist_z_lower = rLocation[2] - lower_location[2];

            const double xd = dist_x_lower / mGridSpacing[0];
            const double yd = dist_y_lower / mGridSpacing[1];
            const double zd = dist_z_lower / mGridSpacing[2];

            const double c00 = field_xyz_lower * (1.0 - xd) + (field_x_upper * xd);
            const double c01 = field_z_upper * (1.0 - xd) + (field_xz_upper * xd);
            const double c10 = field_y_upper * (1.0 - xd) + (field_xy_upper * xd);
            const double c11 = field_yz_upper * (1.0 - xd) + (field_xyz_upper * xd);

            const double c0 = c00 * (1.0 - yd) + (c10 * yd);
            const double c1 = c01 * (1.0 - yd) + (c11 * yd);

            const double c = c0 * (1.0 - zd) + (c1 * zd);
            interpolated_value = c;

            break;
        }
        default:
            NEVER_REACHED;
    }

    return interpolated_value;
}

template <unsigned SPACE_DIM>
long UniformGridRandomFieldGenerator<SPACE_DIM>::GetLinearIndex(std::array<long, SPACE_DIM> gridIndex) const
{
    long linear_index;

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
std::array<double, SPACE_DIM> UniformGridRandomFieldGenerator<SPACE_DIM>::GetPositionUsingGridIndex(std::array<long, SPACE_DIM> gridIndex) const
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
