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

#ifndef SEMSPATIALLYCORRELATEDRANDOMFORCE_HPP_
#define SEMSPATIALLYCORRELATEDRANDOMFORCE_HPP_

#include <array>
#include <memory>

#include "AbstractSemRandomForce.hpp"
#include "OffLatticeRandomFieldGenerator.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * SEM random force using the existing off-lattice spatial random field
 * generator to correlate the unit noise between nearby nodes.
 */
template<unsigned DIM>
class SemSpatiallyCorrelatedRandomForce : public AbstractSemRandomForce<DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSemRandomForce<DIM> >(*this);
        archive & mCorrelationLength;
        archive & mLowerCorner;
        archive & mUpperCorner;
        archive & mPeriodicity;
        archive & mSmallCorrelationLengthWarningThreshold;
        archive & mRandomSeed;
        archive & mHasRandomSeed;
    }

    /** Physical correlation length exposed to SEM users. */
    double mCorrelationLength;

    /** Lower corner of the random-field domain. */
    std::array<double, DIM> mLowerCorner;

    /** Upper corner of the random-field domain. */
    std::array<double, DIM> mUpperCorner;

    /** Periodicity of the random-field domain. */
    std::array<bool, DIM> mPeriodicity;

    /** Warn when the requested correlation length is below this value. */
    double mSmallCorrelationLengthWarningThreshold;

    /** Optional random seed for deterministic field generation. */
    unsigned mRandomSeed;

    /** Whether mRandomSeed has been explicitly set. */
    bool mHasRandomSeed;

    /** Lazily-created random-field generators, one for each force component. */
    std::array<std::unique_ptr<OffLatticeRandomFieldGenerator<DIM> >, DIM> mpRandomFieldGenerators;

    /** Whether the generators need rebuilding before the next use. */
    bool mGeneratorsNeedRefresh;

    /**
     * Build or rebuild the random-field generators.
     */
    void RefreshRandomFieldGenerators();

protected:

    /**
     * Generate spatially correlated unit noise.
     *
     * @param rNodes the SEM nodes to sample noise for
     *
     * @return a vector of DIM-dimensional unit noise vectors
     */
    std::vector<c_vector<double, DIM> > GenerateUnitNoise(const std::vector<Node<DIM>*>& rNodes) override;

public:

    /**
     * Constructor.
     */
    SemSpatiallyCorrelatedRandomForce();

    /**
     * Destructor.
     */
    virtual ~SemSpatiallyCorrelatedRandomForce() = default;

    /**
     * @return the correlation length
     */
    double GetCorrelationLength() const;

    /**
     * @param correlationLength the new correlation length
     */
    void SetCorrelationLength(double correlationLength);

    /**
     * @return the small-correlation-length warning threshold
     */
    double GetSmallCorrelationLengthWarningThreshold() const;

    /**
     * @param threshold the new small-correlation-length warning threshold
     */
    void SetSmallCorrelationLengthWarningThreshold(double threshold);

    /**
     * @return the lower corner of the random-field domain
     */
    std::array<double, DIM> GetLowerCorner() const;

    /**
     * @param lowerCorner the lower corner of the random-field domain
     */
    void SetLowerCorner(std::array<double, DIM> lowerCorner);

    /**
     * @return the upper corner of the random-field domain
     */
    std::array<double, DIM> GetUpperCorner() const;

    /**
     * @param upperCorner the upper corner of the random-field domain
     */
    void SetUpperCorner(std::array<double, DIM> upperCorner);

    /**
     * @return the random-field periodicity
     */
    std::array<bool, DIM> GetPeriodicity() const;

    /**
     * @param periodicity the random-field periodicity
     */
    void SetPeriodicity(std::array<bool, DIM> periodicity);

    /**
     * Set the random seed used by the component random fields.
     *
     * @param randomSeed the random seed
     */
    void SetRandomSeed(unsigned randomSeed);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile) override;
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SemSpatiallyCorrelatedRandomForce)

#endif /*SEMSPATIALLYCORRELATEDRANDOMFORCE_HPP_*/
