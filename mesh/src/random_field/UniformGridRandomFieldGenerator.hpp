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

#ifndef UNIFORM_GRID_RANDOM_FIELD_GENERATOR_HPP_
#define UNIFORM_GRID_RANDOM_FIELD_GENERATOR_HPP_

#include <array>
#include <vector>

#include "ChasteSerialization.hpp"
#include "UblasVectorInclude.hpp"

#include "OpenSimplex2S.hpp"

/**
 * A UniformGridRandomFieldGenerator for adding spatially-correlated noise to a mesh.
 */
template<unsigned SPACE_DIM>
class UniformGridRandomFieldGenerator
{
private:

    friend class boost::serialization::access;
    friend class TestSamplesFromCachedRandomField;
    friend class TestUniformGridRandomFieldGenerator;

    /** Coordinates of the lower corner of the grid. */
    std::array<double, SPACE_DIM> mLowerCorner;

    /** Coordinates of the upper corner of the grid. */
    std::array<double, SPACE_DIM> mUpperCorner;

    /** Number of grid points required in each dimension. */
    std::array<unsigned, SPACE_DIM> mNumGridPts;

    /** Whether the grid is periodic in each dimension. */
    std::array<bool, SPACE_DIM> mPeriodicity;

    /** Length scale over which the noise is to be correlated. */
    double mLengthScale;

    /**
     * The product of mNumGridPts; the total number of grid points in the
     * cuboid.
     */
    unsigned mNumTotalGridPts;

    /**
     * Store the calculated grid spacings to avoid recalculation during
     * interpolation.
     */
    std::array<double, SPACE_DIM> mGridSpacing;

    /**
     * Store the calculated one over grid spacings to avoid recalculation during
     * interpolation.
     */
    std::array<double, SPACE_DIM> mOneOverGridSpacing;

    /**
     * The directory name, relative to $CHASTE_TEST_OUTPUT, in which cached
     * random fields will be saved.
     */
    const std::string mCacheDir;

    /**
     * Number of eigenvalues calculated, such that their sum minimally exceeds
     * mTraceProportion * mNumTotalGridPts.
     */
    unsigned mNumEigenvals;

    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mLowerCorner;
        archive & mUpperCorner;
        archive & mNumGridPts;
        archive & mPeriodicity;
        archive & mLengthScale;
    }

    /**
     * Helper method for Interpolate(). Get the linear index (along the flat
     * vector representing a random field instance) given the grid index
     * (x,y,z).
     *
     * @param gridIndex the (x,y,z) coordinate to translate into a linear index
     *
     * @return the linear index along the flat vector representing the (x,y,z)
     *         grid coordinate
     */
    long GetLinearIndex(std::array<long, SPACE_DIM> gridIndex) const;

    /**
     * Helper method for Interpolate(). Get the position in space of the grid
     * point with (x,y,z) grid index gridIndex.
     *
     * @param gridIndex the (x,y,z) coordinate to calculate the position of
     *
     * @return the position in space of the (x,y,z) coordinate
     */
    std::array<double, SPACE_DIM> GetPositionUsingGridIndex(std::array<long, SPACE_DIM> gridIndex) const;

    /**
     * Noise generator
     */
    OpenSimplex2S mOpenSimplex;

public:

    /**
     * Constructor that takes all parameters as arguments.
     *
     * @param lowerCorner the lower corner of the rectangular grid
     * @param upperCorner the upper corner of the rectangular grid
     * @param numGridPts the number of grid points in each dimension
     * @param periodicity whether the grid is periodic in each dimension
     * @param lengthScale the length scale of the correlation when calculating the covariance matrix
     */
    UniformGridRandomFieldGenerator(std::array<double, SPACE_DIM> lowerCorner,
                                    std::array<double, SPACE_DIM> upperCorner,
                                    std::array<unsigned, SPACE_DIM> numGridPts,
                                    std::array<bool, SPACE_DIM> periodicity,
                                    double lengthScale);

    /**
     * Sample an instance of the random field.
     *
     * Calls SampleRandomFieldAtTime() at a time chosen uniformly at random from
     * 0 to 1.
     *
     * @return A vector representing an instance of the random field.
     */
    std::vector<double> SampleRandomField();

    /**
     * Sample an instance of the random field at a given time.
     *
     * @param time time at which to sample the random field
     *
     * @return A vector representing an instance of the random field.
     */
    std::vector<double> SampleRandomFieldAtTime(double time);

    /**
     * Interpolate from the random field by returning the value at the node of
     * the random field closest to the given location.
     *
     * @param rRandomField the random field, assumed obtained from
     *                     SampleRandomField()
     * @param rLocation the location to which we identify the value at the
     *                  closest node in the random field
     * @return the value of the random field closest to rLocation
     */
    double Interpolate(const std::vector<double>& rRandomField,
                       const c_vector<double, SPACE_DIM>& rLocation) const;

    /**
     * Set the random seed used for the generator.
     *
     * @param seed the new seed
     */
    void SetRandomSeed(const unsigned seed);
};

#endif /*UNIFORM_GRID_RANDOM_FIELD_GENERATOR_HPP_*/
