/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef LINEARSPRINGWITHVARIABLESPRINGCONSTANTSFORCE_HPP
#define LINEARSPRINGWITHVARIABLESPRINGCONSTANTSFORCE_HPP

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "GeneralisedLinearSpringForce.hpp"

/**
 * A subclass of GeneralisedLinearSpringForce with variable spring constants.
 */
template<unsigned DIM>
class LinearSpringWithVariableSpringConstantsForce : public GeneralisedLinearSpringForce<DIM>
{
    friend class TestLinearSpringWithVariableSpringConstantsForce;

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
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<DIM> >(*this);
        archive & mUseEdgeBasedSpringConstant;
        archive & mUseMutantSprings;
        archive & mMutantMutantMultiplier;
        archive & mNormalMutantMultiplier;
        archive & mUseBCatSprings;
        archive & mUseApoptoticSprings;
        archive & mBetaCatSpringScaler;
        archive & mApoptoticSpringTensionStiffness;
        archive & mApoptoticSpringCompressionStiffness;
    }

protected:

    /** Whether to use spring constant proportional to cell-cell contact length/area (defaults to false). */
    bool mUseEdgeBasedSpringConstant;

    /** Whether to use different stiffnesses depending on whether either cell is a mutant. */
    bool mUseMutantSprings;

    /** Multiplier for spring stiffness if mutant. */
    double mMutantMutantMultiplier;

    /** Multiplier for spring stiffness if mutant. */
    double mNormalMutantMultiplier;

    /** Use springs which are dependent on beta-catenin levels. */
    bool mUseBCatSprings;

    /** Use springs which are dependent on whether cells are apoptotic. */
    bool mUseApoptoticSprings;

    /** Scaling factor for beta catenin to spring strength. */
    double mBetaCatSpringScaler;

    /** Non-dimensionalized 'stiffness' of a apoptotic cell under tension. */
    double mApoptoticSpringTensionStiffness;

    /** Non-dimensionalized 'stiffness' of a apoptotic cell under compression. */
    double mApoptoticSpringCompressionStiffness;

public:

    /**
     * Constructor.
     */
    LinearSpringWithVariableSpringConstantsForce();

    /**
     * Destructor.
     */
    virtual ~LinearSpringWithVariableSpringConstantsForce();

    /**
     * Set whether to use an edge-based spring constant.
     *
     * @param useEdgeBasedSpringConstant whether to use an edge-based spring constant
     */
    void SetEdgeBasedSpringConstant(bool useEdgeBasedSpringConstant);

    /**
     * Use different spring strengths depending on two cells:
     * Normal-normal, Normal-mutant, mutant-mutant
     *
     * @param useMutantSprings whether to use mutant springs
     * @param mutantMutantMultiplier the multiplier for springs connecting two mutant cells
     * @param normalMutantMultiplier the multiplier for springs connecting a mutant cell with a normal cell
     */
    void SetMutantSprings(bool useMutantSprings, double mutantMutantMultiplier=2, double normalMutantMultiplier=1.5);

    /**
     * Use the amount of beta-catenin on an edge to find spring constant.
     *
     * @param useBCatSprings whether to use beta-catenin-dependent spring stiffness
     */
    void SetBetaCateninSprings(bool useBCatSprings);

    /**
     * Set spring stiffness to be dependent on whether cells are apoptotic
     *
     * @param useApoptoticSprings whether to have apoptosis-dependent spring stiffness
     */
    void SetApoptoticSprings(bool useApoptoticSprings);

    /**
     * Return a multiplication factor for the spring constant, which
     * may depend on whether the given pair of neighbouring cells are
     * e.g. undergoing apoptosis, have mutations, or experience variable
     * levels of beta catenin.
     *
     * @param nodeAGlobalIndex index of one neighbouring node
     * @param nodeBGlobalIndex index of the other neighbouring node
     * @param rCellPopulation the cell population
     * @param isCloserThanRestLength whether the neighbouring nodes lie closer than the rest length of their connecting spring
     *
     * @return the multiplication factor.
     */
    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                      unsigned nodeBGlobalIndex,
                                                      AbstractCellPopulation<DIM>& rCellPopulation,
                                                      bool isCloserThanRestLength);

    /**
     * Overridden AddForceContribution method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * @return mBetaCatSpringScaler
     */
    double GetBetaCatSpringScaler();

    /**
     * Set mBetaCatSpringScaler.
     *
     * @param betaCatSpringScaler the new value of mBetaCatSpringScaler
     */
    void SetBetaCatSpringScaler(double betaCatSpringScaler);

    /**
     * @return mApoptoticSpringTensionStiffness
     */
    double GetApoptoticSpringTensionStiffness();

    /**
     * Set mApoptoticSpringTensionStiffness.
     *
     * @param apoptoticSpringTensionStiffness the new value of mApoptoticSpringTensionStiffness
     */
    void SetApoptoticSpringTensionStiffness(double apoptoticSpringTensionStiffness);

    /**
     * @return mApoptoticSpringCompressionStiffness
     */
    double GetApoptoticSpringCompressionStiffness();

    /**
     * Set mApoptoticSpringCompressionStiffness.
     *
     * @param apoptoticSpringCompressionStiffness the new value of mApoptoticSpringCompressionStiffness
     */
    void SetApoptoticSpringCompressionStiffness(double apoptoticSpringCompressionStiffness);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LinearSpringWithVariableSpringConstantsForce)

#endif /*LINEARSPRINGWITHVARIABLESPRINGCONSTANTSFORCE_HPP*/
