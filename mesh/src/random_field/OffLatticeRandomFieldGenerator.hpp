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

#ifndef OFF_LATTICE_RANDOM_FIELD_GENERATOR_HPP_
#define OFF_LATTICE_RANDOM_FIELD_GENERATOR_HPP_

#include <array>
#include <memory>
#include <vector>

#include "ChasteSerialization.hpp"
#include "Node.hpp"
#include "ObsoleteBoxCollection.hpp"
#include "UblasVectorInclude.hpp"

#include "OpenSimplex2S.hpp"

/**
 * An OffLatticeRandomFieldGenerator for adding spatially-correlated noise to a mesh with off-lattice nodes.
 */
template<unsigned SPACE_DIM>
class OffLatticeRandomFieldGenerator
{
private:

    friend class boost::serialization::access;
    friend class TestOffLatticeRandomFieldGenerator;

    /** Coordinates of the lower corner of the grid. */
    std::array<double, SPACE_DIM> mLowerCorner;

    /** Coordinates of the upper corner of the grid. */
    std::array<double, SPACE_DIM> mUpperCorner;

    /** Whether the grid is periodic in each dimension. */
    std::array<bool, SPACE_DIM> mPeriodicity;

    /** Length scale over which the noise is to be correlated. */
    double mLengthScale;

    /**
     * An owning pointer to a box collection for efficiently keeping track of
     * nearby nodes.
     */
    std::unique_ptr<ObsoleteBoxCollection<SPACE_DIM>> mpBoxCollection;

    /** The number of nodes passed in at the most recent Update(). */
    unsigned mNumNodesAtLastUpdate;

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
        archive & mPeriodicity;
        archive & mLengthScale;
    }

    /**
     * Noise generator
     */
    OpenSimplex2S mOpenSimplex;

public:

    /**
     * Constructor that takes all parameters as arguments, and sets up an
     * appropriate box collection.
     *
     * @param lowerCorner the lower corner of the rectangular grid
     * @param upperCorner the upper corner of the rectangular grid
     * @param periodicity whether the grid is periodic in each dimension
     * @param lengthScale the length scale of the correlation when calculating
     *                    the covariance matrix
     * @param boxWidth the box size for the box collection. Defaults to
     *                 DOUBLE_UNSET in which case it is reset appropriately.
     */
    OffLatticeRandomFieldGenerator(std::array<double, SPACE_DIM> lowerCorner,
                                   std::array<double, SPACE_DIM> upperCorner,
                                   std::array<bool, SPACE_DIM> periodicity,
                                   double lengthScale,
                                   double boxWidth=DOUBLE_UNSET);

    /**
     * Default constructor.
     */
    virtual ~OffLatticeRandomFieldGenerator() = default;

    /**
     * Sample an instance of the random field.
     *
     * Calls SampleRandomFieldAtTime() at a time chosen uniformly at random from
     * 0 to 100 * lengthScale.
     *
     * @param rNodes a vector of nodes at which to sample the random field
     *
     * @return A vector representing an instance of the random field.
     */
    std::vector<double> SampleRandomField(const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Sample an instance of the random field at a given time. First, draw
     * mNumTotalGridPts random numbers from N(0,1), and then create an
     * appropriate linear combination of the eigenvectors.
     *
     * @param rNodes a vector nodes at which to sample the random field
     * @param time time at which to sample the random field
     *
     * @return A vector representing an instance of the random field.
     */
    std::vector<double> SampleRandomFieldAtTime(const std::vector<Node<SPACE_DIM>*>& rNodes, double time);

    /**
     * Set the random seed used for the generator.
     *
     * @param seed the new seed
     */
    void SetRandomSeed(const unsigned seed);
};

#endif /*OFF_LATTICE_RANDOM_FIELD_GENERATOR_HPP_*/
