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

#ifndef ABSTRACTIMMERSEDBOUNDARYFORCE_HPP_
#define ABSTRACTIMMERSEDBOUNDARYFORCE_HPP_

#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * An abstract immersed boundary force class, for use in  immersed boundary
 * cell-based simulations.
 *
 * NOTE: while this class shares similarities with AbstractForce, it is not
 * passed to a simulation object, but instead to an
 * ImmersedBoundarySimulationModifier as all actual immersed boundary logic
 * is handled by the modifer, rather than the typical simulation
 * infrastructure.
 *
 */
template<unsigned DIM>
class AbstractImmersedBoundaryForce : public Identifiable
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive& mAdditiveNormalNoise;
        archive& mNormalNoiseMean;
        archive& mNormalNoiseStdDev;
    }

protected:

    /**
     * Add the random noise to nodes, if mAdditiveNormalNoise is set. This
     * method handles randomization of the forces in such a way as to ensure no
     * net change to forces across the whole domain.
     *
     * @param rCellPopulation an immersed boundary cell population
     */
    void AddNormalNoiseToNodes(
        ImmersedBoundaryCellPopulation<DIM>& rCellPopulation);

    /** Whether to apply multiplicative normal noise to the calculated force. */
    bool mAdditiveNormalNoise;

    /**
     * The mean of the Normal distribution from which random noise variations
     * are drawn. Default value is 1.0
     */
    double mNormalNoiseMean;

    /**
     * The standard deviation of the Normal distribution from which random noise
     * variations are drawn. Default value is 0.0
     */
    double mNormalNoiseStdDev;

public:

    /**
     * Default constructor.
     */
    AbstractImmersedBoundaryForce();

    /**
     * Destructor.
     */
    virtual ~AbstractImmersedBoundaryForce();

    /**
     * Calculates the force on each immersed boundary node. As this method is
     * pure virtual, it must be overridden in subclasses.
     *
     * @param rNodePairs reference to a vector set of node pairs between which
     *     to contribute the force
     * @param rCellPopulation an immersed boundary cell population
     */
    virtual void AddImmersedBoundaryForceContribution(
        std::vector<std::pair<Node<DIM>*, Node<DIM>*> >& rNodePairs,
        ImmersedBoundaryCellPopulation<DIM>& rCellPopulation)=0;

    /**
     * Outputs force parameters to file. As this method is pure virtual, it must
     * be overridden in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputImmersedBoundaryForceParameters(out_stream& rParamsFile)=0;

    /** @return mAdditiveNormalNoise */
    bool GetAdditiveNormalNoise() const;

    /**
     * Enable or disable additive normal noise
     *
     * @param additiveNormalNoise whether to include multiplicative normal noise
     */
    void SetAdditiveNormalNoise(bool additiveNormalNoise);

    /**
     *
     * @return mNormalNoiseMean
     * */
    double GetNormalNoiseMean() const;

    /**
     * Set the mean value of the normal noise
     *
     * @param normalNoiseMean the new value of mNormalNoiseMean
     */
    void SetNormalNoiseMean(double normalNoiseMean);

    /** @return mNormalNoiseStdDev */
    double GetNormalNoiseStdDev() const;

    /**
     * Set the standard deviation of the normal noise
     *
     * @param normalNoiseStdDev the new value of mNormalNoiseStdDev
     */
    void SetNormalNoiseStdDev(double normalNoiseStdDev);
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractImmersedBoundaryForce)

#endif /*ABSTRACTIMMERSEDBOUNDARYFORCE_HPP_*/
